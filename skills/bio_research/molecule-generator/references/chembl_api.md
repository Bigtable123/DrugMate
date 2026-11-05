# ChEMBL REST API Reference

## Base URL

```
https://www.ebi.ac.uk/chembl/api/data/
```

Target/activity API requests use an explicit `.json` endpoint and include `format=json`
as a query parameter. The HTTP helper normalizes the suffix centrally for
activity queries, target lookups, and target searches; an existing `.json`
suffix is preserved. Caller-owned query dictionaries are not modified.

Target/activity requests use a 60-second socket timeout; the lightweight
status probe uses 10 seconds. These are per-operation waits, not a deadline
for retrieving all pages. Raising the pipeline's outer step timeout does not
change these HTTP timeouts. API reads make at most three total attempts
(the initial request plus two retries), with 1s and 2s backoff, for HTTP
408/429/5xx, URL/connection errors, socket timeouts including response reads,
and incomplete HTTP response bodies. Other 4xx and malformed JSON are not
retried. Each failed attempt is logged, and the final exception propagates;
partial pages never become a complete library. The status probe remains a
single 10-second check, and its failure alone does not prove that activity
retrieval is unavailable. On retrieval failure, preserve the raw capture
and follow [Pocket2Mol recovery](../SKILL.md#chembl-unavailable-pocket2mol-recovery)
instead of repeatedly restarting the same unavailable source.

---

## Key Endpoints

### Status Check

```
GET /status.json
```

Returns API status. This is an advisory availability check: its failure must
not block a target lookup by name/UniProt or activity retrieval with a known
target ID. Those endpoints may still work while status is unavailable.

---

### Target Search (by name)

```
GET /target/search.json?q={query}&limit=10&format=json
```

**Response shape:**
```json
{
  "page_meta": {"limit": 10, "offset": 0, "total_count": 5, "next": null},
  "targets": [
    {
      "target_chembl_id": "CHEMBL301",
      "pref_name": "Cyclin-dependent kinase 2",
      "organism": "Homo sapiens",
      "target_type": "SINGLE PROTEIN",
      "target_components": [{"accession": "P24941"}]
    }
  ]
}
```

**Selection priority:** Prefer `target_type == "SINGLE PROTEIN"` + `organism == "Homo sapiens"`.

---

### Target Lookup (by UniProt)

```
GET /target.json?target_components__accession={uniprot_id}&limit=5&format=json
```

Same response shape as search.

---

### Target Validation (by ChEMBL ID)

```
GET /target/{chembl_target_id}.json?format=json
```

Returns single target object (not wrapped in array).

---

### Activity Query

```
GET /activity.json?target_chembl_id={id}&standard_type={type}&standard_value__lte={threshold}&standard_units=nM&limit=100&offset=0&format=json
```

**Key parameters:**

| Parameter | Description |
|-----------|-------------|
| `target_chembl_id` | ChEMBL target identifier |
| `standard_type` | Activity type: IC50, Ki, Kd, EC50 |
| `standard_value__lte` | Upper bound for activity value (nM) |
| `standard_units` | "nM" |
| `limit` | Max results per page (max 1000) |
| `offset` | Pagination offset |

The client accepts activity rows whose normalized `standard_relation` is one of
`=`, `<`, or `<=`, and rejects unsupported relations such as `>`/`>=`. The
relation filter is applied locally so the raw capture can audit every row
returned by the numeric ChEMBL query.

**Response shape:**
```json
{
  "page_meta": {"limit": 100, "offset": 0, "total_count": 500, "next": "/activity.json?offset=100&..."},
  "activities": [
    {
      "molecule_chembl_id": "CHEMBL12345",
      "canonical_smiles": "CC(=O)Nc1ccc(O)cc1",
      "standard_value": 150.0,
      "standard_type": "IC50",
      "standard_units": "nM",
      "standard_relation": "=",
      "activity_id": 1234567,
      "assay_chembl_id": "CHEMBL67890",
      "target_chembl_id": "CHEMBL301"
    }
  ]
}
```

---

## Pagination

ChEMBL exposes `limit`/`offset` pagination through `page_meta.next`:
- Default limit: 20
- Max limit: 1000
- `page_meta.next` is null when no more pages

**Required implementation contract:**
```python
params = {"offset": 0, "limit": 100, **activity_filters}
retrieved = 0
expected = None
while True:
    data = api_get("/activity.json", params)
    activities = data["activities"]
    retrieved += len(activities)
    if expected is None:
        expected = data["page_meta"]["total_count"]
    next_page = data["page_meta"].get("next")
    if not next_page:
        break
    params.update(query_parameters_from(next_page))
    time.sleep(0.5)  # rate limiting

if retrieved != expected:
    raise RuntimeError("incomplete ChEMBL result set")
```

- Follow the exact `page_meta.next` continuation rather than synthesizing the
  next offset.
- A safety page limit may prevent runaway/cyclic responses, but if `next` is
  still present at that limit the retrieval is `complete=false`,
  `truncated=true`, and must fail rather than supply a partial active set.
- Persist `expected_records`, `retrieved_records`, page count, `complete`,
  `truncated`, `next`, and any error beside the raw pages.
- A missing/malformed `total_count`, malformed/cyclic `next`, or a final count
  mismatch is an explicit incomplete-retrieval failure.
- Treat target, activity type, units, threshold, and endpoint as immutable while
  following `next`; a continuation link or returned row that contradicts the
  original query is quarantined/fails integrity checks rather than silently
  changing the cohort.

---

## Rate Limiting

- No authentication required
- Informal rate limit: ~10 requests/second
- Returns HTTP 429 on excessive requests
- Recommended: 0.5s sleep between paginated requests
- Retry strategy: at most three total attempts, with 1s then 2s backoff, for
  HTTP 408/429/5xx and transient transport/read errors as described above

---

## Common Target ChEMBL IDs

| Target | ChEMBL ID | UniProt |
|--------|-----------|---------|
| CDK2 | CHEMBL301 | P24941 |
| EGFR | CHEMBL203 | P00533 |
| BRAF | CHEMBL5145 | P15056 |
| ABL1 | CHEMBL1862 | P00519 |
| PI3Kα | CHEMBL3267 | P42336 |
| JAK2 | CHEMBL2971 | O60674 |

---

## Notes

- Preserve `standard_relation` for `=`, `<`, and `<=`; never reinterpret a
  censored value as an exact point.
- Deduplicate by `canonical_smiles` (same compound may have multiple assay entries)
- Some entries have `pchembl_value` (negative log molar) which is pre-computed
- For repeated measurements, compute an ordinary median from exact (`=`)
  observations only. Keep every censored observation in the audit trail. If no
  exact observation exists, retain the largest reported upper bound (`<x` or
  `<=x`) as a conservative censored summary; do not manufacture an exact value.
- Carry relation, explicit lower/upper bound fields, inclusivity, aggregation
  method, measurement counts, activity/assay IDs, and the original measurement
  list through standardized-SMILES deduplication, the library manifest, and SDF
  tags.
- Preserve rows marked by `data_validity_comment` or `potential_duplicate` in
  the raw and per-compound audit trail, but exclude them from analytical
  aggregation. Count the two flags independently because one row may carry
  both. A SMILES group containing only quarantined rows is omitted and counted
  as a quarantined compound group rather than aborting other valid groups.
- Preserve assay type, BAO format, document IDs, and assay-stratified summaries.
  Assay heterogeneity is an audit warning; it does not silently change the
  primary exact-median/censored-upper-bound selection policy.
- Target search may return non-human orthologues; filter by organism
