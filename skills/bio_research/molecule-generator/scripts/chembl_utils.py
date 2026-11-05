"""ChEMBL API utilities for target resolution and active compound retrieval.

Uses only stdlib urllib (no requests dependency), following the same HTTP patterns
as structure-prep/scripts/pdb_utils.py.
"""
from __future__ import annotations

import json
import logging
import math
import time
import urllib.error
import urllib.request
from http.client import IncompleteRead
from urllib.parse import parse_qs, urlencode, urlparse
from datetime import datetime, timezone
from pathlib import Path


CHEMBL_API_BASE = "https://www.ebi.ac.uk/chembl/api/data"
CHEMBL_STATUS_URL = f"{CHEMBL_API_BASE}/status.json"
logger = logging.getLogger(__name__)

_REQUEST_TIMEOUT = 60
_PAGE_SLEEP = 0.5
_MAX_ATTEMPTS = 3
_MAX_ACTIVITY_PAGES = 1000
ACTIVE_STANDARD_RELATIONS = {"=", "<", "<="}


def check_chembl_api_available() -> bool:
    """Check if ChEMBL REST API is reachable."""
    try:
        req = urllib.request.Request(CHEMBL_STATUS_URL, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=10) as resp:
            return resp.status == 200
    except urllib.error.HTTPError as exc:
        exc.close()
        return False
    except Exception:
        return False


def _api_get(endpoint: str, params: dict | None = None) -> dict:
    """Read a JSON endpoint with at most three attempts for transient faults.

    Use the explicit JSON route for every resource, including target/search
    and paginated activity requests. The socket timeout applies to each
    network operation, not to the complete paginated retrieval.
    """
    parsed = urlparse(endpoint)
    path = parsed.path.strip("/")
    if not path.endswith(".json"):
        path += ".json"
    query = parse_qs(parsed.query, keep_blank_values=True)
    query.update(params or {})
    query["format"] = "json"
    url = f"{CHEMBL_API_BASE}/{path}?{urlencode(query, doseq=True)}"

    for attempt in range(1, _MAX_ATTEMPTS + 1):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=_REQUEST_TIMEOUT) as resp:
                payload = resp.read()
        except urllib.error.HTTPError as exc:
            # HTTPError owns the failed response stream too.
            exc.close()
            retryable = exc.code in {408, 429} or 500 <= exc.code <= 599
            if not retryable or attempt == _MAX_ATTEMPTS:
                logger.warning(
                    "ChEMBL request failed on attempt %d/%d: %s: %s",
                    attempt, _MAX_ATTEMPTS, url, exc,
                )
                raise
            error = exc
        except (urllib.error.URLError, TimeoutError, ConnectionError, IncompleteRead) as exc:
            if attempt == _MAX_ATTEMPTS:
                logger.warning(
                    "ChEMBL request failed on attempt %d/%d: %s: %s: %s",
                    attempt, _MAX_ATTEMPTS, url, type(exc).__name__, exc,
                )
                raise
            error = exc
        else:
            # Invalid JSON is a response/data error, not a transport retry.
            return json.loads(payload.decode("utf-8"))

        delay = 2 ** (attempt - 1)
        logger.warning(
            "ChEMBL request attempt %d/%d failed; retrying in %ss: %s: %s: %s",
            attempt, _MAX_ATTEMPTS, delay, url, type(error).__name__, error,
        )
        time.sleep(delay)

    raise AssertionError("ChEMBL retry loop exited unexpectedly")


def resolve_target_to_chembl_id(
    target_name: str | None = None,
    chembl_target_id: str | None = None,
    uniprot_id: str | None = None,
) -> dict:
    """Resolve any identifier to a ChEMBL target ID.

    Resolution order:
      1. chembl_target_id provided → validate via direct lookup
      2. uniprot_id → query by accession
      3. target_name → free-text search, pick best human match

    Returns: {"chembl_id": str, "pref_name": str, "organism": str, "target_type": str}
    Raises ValueError if resolution fails.
    """
    if chembl_target_id:
        try:
            data = _api_get(f"target/{chembl_target_id}")
            return {
                "chembl_id": data.get("target_chembl_id", chembl_target_id),
                "pref_name": data.get("pref_name", ""),
                "organism": data.get("organism", ""),
                "target_type": data.get("target_type", ""),
            }
        except Exception as e:
            raise ValueError(f"Failed to validate ChEMBL target {chembl_target_id}: {e}")

    if uniprot_id:
        try:
            data = _api_get("target", {"target_components__accession": uniprot_id, "limit": 5})
            targets = data.get("targets", [])
            if not targets:
                raise ValueError(f"No ChEMBL target found for UniProt {uniprot_id}")
            for t in targets:
                if "homo sapiens" in (t.get("organism", "") or "").lower():
                    return {
                        "chembl_id": t["target_chembl_id"],
                        "pref_name": t.get("pref_name", ""),
                        "organism": t.get("organism", ""),
                        "target_type": t.get("target_type", ""),
                    }
            t = targets[0]
            return {
                "chembl_id": t["target_chembl_id"],
                "pref_name": t.get("pref_name", ""),
                "organism": t.get("organism", ""),
                "target_type": t.get("target_type", ""),
            }
        except ValueError:
            raise
        except Exception as e:
            raise ValueError(f"Failed to resolve UniProt {uniprot_id}: {e}")

    if target_name:
        try:
            # Strategy 1: Try gene_name filter for short names (e.g., "CDK2", "EGFR")
            targets = []
            if len(target_name) <= 10 and target_name.isalnum():
                try:
                    data = _api_get("target", {
                        "target_components__target_component_synonyms__component_synonym__iexact": target_name,
                        "target_type": "SINGLE PROTEIN",
                        "organism": "Homo sapiens",
                        "limit": 5,
                    })
                    targets = data.get("targets", [])
                except Exception:
                    pass

            # Strategy 2: Try pref_name contains (for full names)
            if not targets:
                try:
                    data = _api_get("target", {
                        "pref_name__icontains": target_name,
                        "target_type": "SINGLE PROTEIN",
                        "organism": "Homo sapiens",
                        "limit": 10,
                    })
                    targets = data.get("targets", [])
                except Exception:
                    pass

            # Strategy 3: Free-text search as last resort
            if not targets:
                data = _api_get("target/search", {"q": target_name, "limit": 20})
                targets = data.get("targets", [])

            if not targets:
                raise ValueError(f"No ChEMBL target found for name '{target_name}'")

            # Priority: SINGLE PROTEIN + Homo sapiens > other types
            type_priority = {"SINGLE PROTEIN": 0, "PROTEIN COMPLEX": 1, "PROTEIN FAMILY": 2, "SELECTIVITY GROUP": 3}
            def _score(t):
                is_human = 1 if "homo sapiens" in (t.get("organism", "") or "").lower() else 0
                type_rank = type_priority.get(t.get("target_type", ""), 9)
                return (-is_human, type_rank)

            targets.sort(key=_score)
            t = targets[0]
            return {
                "chembl_id": t["target_chembl_id"],
                "pref_name": t.get("pref_name", ""),
                "organism": t.get("organism", ""),
                "target_type": t.get("target_type", ""),
            }
        except ValueError:
            raise
        except Exception as e:
            raise ValueError(f"Failed to search for target '{target_name}': {e}")

    raise ValueError("Must provide target_name, chembl_target_id, or uniprot_id")


def fetch_active_compounds(
    chembl_target_id: str,
    activity_type: str = "IC50",
    threshold_nm: float = 1000.0,
    max_compounds: int = 50,
    raw_output_file: str | Path | None = None,
    *,
    max_pages: int = _MAX_ACTIVITY_PAGES,
) -> list[dict]:
    """Fetch active compounds from ChEMBL with activity below threshold.

    Follow ChEMBL's ``page_meta.next`` links until the result set is complete.
    A safety limit is retained to prevent an unbounded or cyclic API response,
    but reaching it with another page pending is an explicit retrieval failure.

    Repeated exact measurements are summarized by their median.  Censored
    measurements (``<``/``<=``) are never inserted into that exact median; if
    no exact measurement exists, the least-potent reported upper bound is kept
    as a conservative, explicitly censored summary.
    """
    if max_pages < 1:
        raise ValueError("max_pages must be at least 1")
    threshold_value = float(threshold_nm)
    if not math.isfinite(threshold_value) or threshold_value <= 0:
        raise ValueError("threshold_nm must be a finite positive number")
    threshold_query_value: int | float = (
        int(threshold_value) if threshold_value.is_integer() else threshold_value
    )

    params = {
        "target_chembl_id": chembl_target_id,
        "standard_type": activity_type,
        # Preserve a caller-supplied fractional-nM cutoff; integer coercion
        # would silently make the retrieval cohort stricter.
        "standard_value__lte": threshold_query_value,
        "standard_units": "nM",
        "limit": 100,
        "offset": 0,
    }

    all_activities = []
    raw_pages = []
    request_params = dict(params)
    seen_next_links: set[str] = set()
    expected_records: int | None = None
    retrieved_records = 0
    next_page: str | None = None
    quality_audit = {
        "flagged_data_validity_records": 0,
        "flagged_potential_duplicate_records": 0,
        "quarantined_quality_records": 0,
        "eligible_compound_groups": 0,
        "quarantined_compound_groups": 0,
        "excluded_data_validity_records": 0,
        "excluded_potential_duplicate_records": 0,
        "excluded_query_mismatch_records": 0,
    }
    quality_eligible_records = 0

    while True:
        try:
            data = _api_get("activity", dict(request_params))
        except Exception as exc:
            _write_activity_capture(
                raw_output_file,
                query=params,
                pages=raw_pages,
                expected_records=expected_records,
                retrieved_records=retrieved_records,
                accepted_records=quality_eligible_records,
                complete=False,
                truncated=False,
                max_pages=max_pages,
                next_page=next_page,
                error=str(exc),
                failed_offset=request_params.get("offset"),
                quality_audit=quality_audit,
            )
            raise RuntimeError(
                "ChEMBL activity request failed at offset "
                f"{request_params.get('offset')}: {exc}"
            ) from exc

        raw_pages.append(data)

        activities = data.get("activities", [])
        retrieved_records += len(activities)
        page_meta = data.get("page_meta", {}) or {}
        raw_total_count = page_meta.get("total_count")
        page_expected = _as_nonnegative_int(raw_total_count)
        if page_expected is None:
            error = (
                "ChEMBL page_meta.total_count is missing or invalid: "
                f"{raw_total_count!r}"
            )
            _write_activity_capture(
                raw_output_file,
                query=params,
                pages=raw_pages,
                expected_records=expected_records,
                retrieved_records=retrieved_records,
                accepted_records=quality_eligible_records,
                complete=False,
                truncated=bool(page_meta.get("next")),
                max_pages=max_pages,
                next_page=page_meta.get("next"),
                error=error,
                quality_audit=quality_audit,
            )
            raise RuntimeError(error)
        if expected_records is None:
            expected_records = page_expected
        elif page_expected != expected_records:
            error = (
                "ChEMBL page_meta.total_count changed during pagination: "
                f"{expected_records} -> {page_expected}"
            )
            _write_activity_capture(
                raw_output_file,
                query=params,
                pages=raw_pages,
                expected_records=expected_records,
                retrieved_records=retrieved_records,
                accepted_records=quality_eligible_records,
                complete=False,
                truncated=bool(page_meta.get("next")),
                max_pages=max_pages,
                next_page=page_meta.get("next"),
                error=error,
                quality_audit=quality_audit,
            )
            raise RuntimeError(error)

        for act in activities:
            validity_comment = str(act.get("data_validity_comment") or "").strip()
            potential_duplicate = _truthy_flag(act.get("potential_duplicate"))
            # Count ChEMBL's two quality flags independently.  A record may
            # carry both, so the two flag counts are not mutually exclusive;
            # ``quarantined_quality_records`` is their union.
            if validity_comment:
                quality_audit["flagged_data_validity_records"] += 1
                quality_audit["excluded_data_validity_records"] += 1
            if potential_duplicate:
                quality_audit["flagged_potential_duplicate_records"] += 1
                # ChEMBL marks these records because counting them alongside
                # their source-equivalent row can overweight one experiment.
                quality_audit["excluded_potential_duplicate_records"] += 1
            quality_quarantined = bool(validity_comment or potential_duplicate)
            if quality_quarantined:
                quality_audit["quarantined_quality_records"] += 1
            if not _activity_row_matches_query(
                act,
                target_chembl_id=chembl_target_id,
                activity_type=activity_type,
                units="nM",
            ):
                quality_audit["excluded_query_mismatch_records"] += 1
                continue
            smiles = act.get("canonical_smiles")
            value = act.get("standard_value")
            relation = _normalize_standard_relation(act.get("standard_relation"))
            if relation not in ACTIVE_STANDARD_RELATIONS:
                continue
            if not smiles or value is None:
                continue
            try:
                value_float = float(value)
            except (TypeError, ValueError):
                continue
            if not math.isfinite(value_float) or value_float < 0:
                continue
            if value_float > threshold_value:
                continue

            activity = {
                "chembl_id": act.get("molecule_chembl_id", ""),
                "smiles": smiles,
                "activity_nm": value_float,
                "activity_type": act.get("standard_type", activity_type),
                "standard_relation": relation,
                "assay_chembl_id": act.get("assay_chembl_id", ""),
                "activity_chembl_id": act.get("activity_id", ""),
                "target_chembl_id": act.get("target_chembl_id", ""),
                "standard_units": act.get("standard_units", "nM"),
                "assay_type": act.get("assay_type"),
                "assay_description": act.get("assay_description"),
                "bao_format": act.get("bao_format"),
                "bao_label": act.get("bao_label"),
                "document_chembl_id": act.get("document_chembl_id"),
                "document_journal": act.get("document_journal"),
                "document_year": act.get("document_year"),
                "data_validity_comment": act.get("data_validity_comment"),
                "data_validity_description": act.get(
                    "data_validity_description"
                ),
                "potential_duplicate": potential_duplicate,
                "standard_flag": act.get("standard_flag"),
                "pchembl_value": act.get("pchembl_value"),
                "source_id": act.get("src_id"),
                "aggregation_included": not quality_quarantined,
                "aggregation_exclusion_reasons": [
                    reason
                    for reason, present in (
                        ("data_validity_comment", bool(validity_comment)),
                        ("potential_duplicate", potential_duplicate),
                    )
                    if present
                ],
            }
            activity.update(_activity_bound_fields(value_float, relation))
            all_activities.append(activity)
            if not quality_quarantined:
                quality_eligible_records += 1

        next_value = page_meta.get("next")
        next_page = str(next_value) if next_value else None
        if next_page is None:
            break

        if len(raw_pages) >= max_pages:
            error = (
                f"ChEMBL activity retrieval reached safety limit of {max_pages} "
                "pages while page_meta.next still advertised more records"
            )
            _write_activity_capture(
                raw_output_file,
                query=params,
                pages=raw_pages,
                expected_records=expected_records,
                retrieved_records=retrieved_records,
                accepted_records=quality_eligible_records,
                complete=False,
                truncated=True,
                max_pages=max_pages,
                next_page=next_page,
                error=error,
                quality_audit=quality_audit,
            )
            raise RuntimeError(error)

        if next_page in seen_next_links:
            error = f"ChEMBL pagination cycle detected at {next_page}"
            _write_activity_capture(
                raw_output_file,
                query=params,
                pages=raw_pages,
                expected_records=expected_records,
                retrieved_records=retrieved_records,
                accepted_records=quality_eligible_records,
                complete=False,
                truncated=True,
                max_pages=max_pages,
                next_page=next_page,
                error=error,
                quality_audit=quality_audit,
            )
            raise RuntimeError(error)
        seen_next_links.add(next_page)
        try:
            request_params = _params_from_next_page(next_page, request_params)
        except Exception as exc:
            error = f"Cannot follow ChEMBL page_meta.next {next_page!r}: {exc}"
            _write_activity_capture(
                raw_output_file,
                query=params,
                pages=raw_pages,
                expected_records=expected_records,
                retrieved_records=retrieved_records,
                accepted_records=quality_eligible_records,
                complete=False,
                truncated=True,
                max_pages=max_pages,
                next_page=next_page,
                error=error,
                quality_audit=quality_audit,
            )
            raise RuntimeError(error) from exc
        time.sleep(_PAGE_SLEEP)

    complete = retrieved_records == expected_records
    if not complete:
        error = (
            "ChEMBL pagination ended without page_meta.next, but the retrieved "
            f"record count does not match total_count ({retrieved_records} != "
            f"{expected_records})"
        )
        _write_activity_capture(
            raw_output_file,
            query=params,
            pages=raw_pages,
            expected_records=expected_records,
            retrieved_records=retrieved_records,
            accepted_records=quality_eligible_records,
            complete=False,
            truncated=False,
            max_pages=max_pages,
            next_page=None,
            error=error,
            quality_audit=quality_audit,
        )
        raise RuntimeError(error)

    group_audit = _quality_group_audit(all_activities)
    quality_audit.update(group_audit)
    deduped = _deduplicate_by_smiles(all_activities)
    _write_activity_capture(
        raw_output_file,
        query=params,
        pages=raw_pages,
        expected_records=expected_records,
        retrieved_records=retrieved_records,
        accepted_records=quality_eligible_records,
        complete=True,
        truncated=False,
        max_pages=max_pages,
        next_page=None,
        quality_audit=quality_audit,
    )

    deduped.sort(key=lambda x: x["activity_nm"])
    return deduped[:max_compounds]


def _write_activity_capture(
    output_file: str | Path | None,
    *,
    query: dict,
    pages: list[dict],
    expected_records: int | None,
    retrieved_records: int,
    accepted_records: int,
    complete: bool,
    truncated: bool,
    max_pages: int,
    next_page: str | None,
    error: str | None = None,
    failed_offset: object = None,
    quality_audit: dict[str, int] | None = None,
) -> None:
    """Persist raw pages and an explicit completeness contract when requested."""
    if not output_file:
        return
    clean_query = {
        **query,
        "standard_relation_allowed": sorted(ACTIVE_STANDARD_RELATIONS),
    }
    retrieval = {
        "expected_records": expected_records,
        "retrieved_records": retrieved_records,
        "accepted_records": accepted_records,
        "pages_retrieved": len(pages),
        "complete": bool(complete),
        "truncated": bool(truncated),
        "max_pages": max_pages,
        "next": next_page,
        **dict(quality_audit or {}),
    }
    payload = {
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
        "api_base": CHEMBL_API_BASE,
        "query": clean_query,
        # Keep these top-level fields easy for orchestrators to inspect while
        # retaining the structured retrieval block for forward extension.
        "expected_records": expected_records,
        "retrieved_records": retrieved_records,
        "complete": bool(complete),
        "truncated": bool(truncated),
        "retrieval": retrieval,
        "pages": pages,
    }
    if error:
        payload["error"] = error
    if failed_offset is not None:
        payload["failed_offset"] = failed_offset
    Path(output_file).write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )


def _params_from_next_page(next_page: str, current: dict) -> dict:
    """Build the next request from ChEMBL's advertised pagination link."""
    parsed = urlparse(next_page)
    query = parse_qs(parsed.query, keep_blank_values=True)
    if not query:
        raise RuntimeError(
            f"ChEMBL page_meta.next has no query parameters: {next_page}"
        )
    normalized_path = parsed.path.rstrip("/")
    if not (
        normalized_path.endswith("/activity")
        or normalized_path.endswith("/activity.json")
    ):
        raise RuntimeError(
            "ChEMBL page_meta.next changed the activity endpoint: "
            f"{parsed.path!r}"
        )
    params = dict(current)
    immutable = {
        key: value for key, value in current.items()
        if key not in {"limit", "offset"}
    }
    for key, values in query.items():
        if not values:
            continue
        if len(values) != 1:
            raise RuntimeError(
                "ChEMBL page_meta.next repeated query parameter "
                f"{key!r}: {values!r}"
            )
        value: object = values[-1]
        if key == "format":
            if str(value).strip().casefold() != "json":
                raise RuntimeError(
                    "ChEMBL page_meta.next changed response format: "
                    f"{value!r}"
                )
            continue
        if key in {"limit", "offset"}:
            try:
                value = int(str(value))
            except ValueError as exc:
                raise RuntimeError(
                    f"ChEMBL page_meta.next has invalid {key}: {value}"
                ) from exc
            if key == "offset" and int(value) <= int(current.get("offset", -1)):
                raise RuntimeError(
                    "ChEMBL page_meta.next did not advance the offset: "
                    f"{current.get('offset')} -> {value}"
                )
            params[key] = value
            continue
        if key not in immutable:
            raise RuntimeError(
                f"ChEMBL page_meta.next introduced unexpected query key {key!r}"
            )
        if str(value) != str(immutable[key]):
            raise RuntimeError(
                "ChEMBL page_meta.next changed immutable query parameter "
                f"{key!r}: {immutable[key]!r} -> {value!r}"
            )
    return params


def _activity_row_matches_query(
    activity: dict,
    *,
    target_chembl_id: str,
    activity_type: str,
    units: str,
) -> bool:
    """Reject a response row that drifts outside the immutable activity query."""
    expected = {
        "target_chembl_id": str(target_chembl_id).strip().casefold(),
        "standard_type": str(activity_type).strip().casefold(),
        "standard_units": str(units).strip().casefold(),
    }
    for field, wanted in expected.items():
        observed = activity.get(field)
        # Some test doubles/older ChEMBL responses omit redundant filter fields;
        # a present contradictory value is the integrity failure.
        if observed not in (None, "") and str(observed).strip().casefold() != wanted:
            return False
    return True


def _truthy_flag(value: object) -> bool:
    if value in (None, "", 0, False):
        return False
    if isinstance(value, str):
        normalized = value.strip().casefold()
        if normalized in {"0", "false", "no", "n"}:
            return False
        if normalized in {"1", "true", "yes", "y"}:
            return True
    return bool(value)


def _as_nonnegative_int(value: object) -> int | None:
    if value in (None, ""):
        return None
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return None
    return parsed if parsed >= 0 else None


def _normalize_standard_relation(value: object) -> str:
    relation = str(value or "").strip()
    relation = relation.replace("&lt;", "<").replace("&gt;", ">")
    return relation


def _activity_bound_fields(value_nm: float, relation: str) -> dict:
    """Return an explicit interval without inventing a lower censoring bound."""
    value = float(value_nm)
    if relation == "=":
        return {
            "activity_lower_bound_nm": value,
            "activity_upper_bound_nm": value,
            "activity_lower_bound_inclusive": True,
            "activity_upper_bound_inclusive": True,
            "activity_is_censored": False,
        }
    if relation in {"<", "<="}:
        return {
            "activity_lower_bound_nm": None,
            "activity_upper_bound_nm": value,
            "activity_lower_bound_inclusive": None,
            "activity_upper_bound_inclusive": relation == "<=",
            "activity_is_censored": True,
        }
    raise ValueError(f"unsupported ChEMBL activity relation: {relation!r}")


def _measurement_record(compound: dict) -> dict:
    value = float(compound["activity_nm"])
    relation = _normalize_standard_relation(compound.get("standard_relation")) or "="
    record = {
        "value_nm": value,
        "standard_relation": relation,
        "chembl_id": compound.get("chembl_id", ""),
        "assay_chembl_id": compound.get("assay_chembl_id", ""),
        "activity_chembl_id": compound.get("activity_chembl_id", ""),
        "activity_type": compound.get("activity_type", ""),
        "target_chembl_id": compound.get("target_chembl_id", ""),
        "standard_units": compound.get("standard_units", "nM"),
        "assay_type": compound.get("assay_type"),
        "assay_description": compound.get("assay_description"),
        "bao_format": compound.get("bao_format"),
        "bao_label": compound.get("bao_label"),
        "document_chembl_id": compound.get("document_chembl_id"),
        "document_journal": compound.get("document_journal"),
        "document_year": compound.get("document_year"),
        "data_validity_comment": compound.get("data_validity_comment"),
        "data_validity_description": compound.get(
            "data_validity_description"
        ),
        "potential_duplicate": _truthy_flag(
            compound.get("potential_duplicate")
        ),
        "standard_flag": compound.get("standard_flag"),
        "pchembl_value": compound.get("pchembl_value"),
        "source_id": compound.get("source_id"),
    }
    bounds = _activity_bound_fields(value, relation)
    record.update({
        "lower_bound_nm": bounds["activity_lower_bound_nm"],
        "upper_bound_nm": bounds["activity_upper_bound_nm"],
        "lower_bound_inclusive": bounds["activity_lower_bound_inclusive"],
        "upper_bound_inclusive": bounds["activity_upper_bound_inclusive"],
        "is_censored": bounds["activity_is_censored"],
    })
    return record


def _flatten_measurements(compounds: list[dict]) -> list[dict]:
    measurements = []
    for compound in compounds:
        nested = compound.get("activity_measurements")
        if isinstance(nested, list) and nested:
            for item in nested:
                if not isinstance(item, dict) or item.get("value_nm") is None:
                    continue
                relation = _normalize_standard_relation(
                    item.get("standard_relation")
                ) or "="
                normalized = {
                    **item,
                    "value_nm": float(item["value_nm"]),
                    "standard_relation": relation,
                }
                bounds = _activity_bound_fields(normalized["value_nm"], relation)
                normalized.update({
                    "lower_bound_nm": bounds["activity_lower_bound_nm"],
                    "upper_bound_nm": bounds["activity_upper_bound_nm"],
                    "lower_bound_inclusive": bounds[
                        "activity_lower_bound_inclusive"
                    ],
                    "upper_bound_inclusive": bounds[
                        "activity_upper_bound_inclusive"
                    ],
                    "is_censored": bounds["activity_is_censored"],
                    "potential_duplicate": _truthy_flag(
                        normalized.get("potential_duplicate")
                    ),
                })
                measurements.append(normalized)
        elif compound.get("activity_nm") is not None:
            measurements.append(_measurement_record(compound))
    return sorted(
        measurements,
        key=lambda item: (
            float(item["value_nm"]),
            str(item.get("standard_relation") or ""),
            str(item.get("chembl_id") or ""),
            str(item.get("assay_chembl_id") or ""),
            str(item.get("activity_chembl_id") or ""),
        ),
    )


def aggregate_activity_measurements(compounds: list[dict]) -> dict:
    """Aggregate repeated measurements without treating censoring as equality.

    Exact measurements, when present, define a conventional median point
    estimate.  Censored measurements remain attached for audit but do not enter
    that median.  If every observation is censored, the largest reported upper
    bound is retained so the summary does not claim greater potency than any
    input record establishes.
    """
    if not compounds:
        raise ValueError("cannot aggregate an empty activity group")
    all_measurements = _flatten_measurements(compounds)
    if not all_measurements:
        raise ValueError("activity group contains no numeric measurements")

    measurements: list[dict] = []
    for measurement in all_measurements:
        reasons = []
        if str(measurement.get("data_validity_comment") or "").strip():
            reasons.append("data_validity_comment")
        if _truthy_flag(measurement.get("potential_duplicate")):
            reasons.append("potential_duplicate")
        measurement["aggregation_included"] = not reasons
        measurement["aggregation_exclusion_reasons"] = reasons
        if not reasons:
            measurements.append(measurement)
    if not measurements:
        raise ValueError(
            "activity group contains no quality-eligible measurements"
        )

    exact = [m for m in measurements if m["standard_relation"] == "="]
    if exact:
        exact_values = sorted(float(item["value_nm"]) for item in exact)
        n_exact = len(exact_values)
        aggregate_value = (
            exact_values[n_exact // 2]
            if n_exact % 2
            else (exact_values[n_exact // 2 - 1] + exact_values[n_exact // 2]) / 2
        )
        selected_measurement = min(
            exact,
            key=lambda item: abs(float(item["value_nm"]) - aggregate_value),
        )
        aggregate_relation = "="
        aggregation_method = "median_exact_measurements"
    else:
        # At the same upper bound, <= is the weaker/more conservative claim.
        selected_measurement = max(
            measurements,
            key=lambda item: (
                float(item["upper_bound_nm"]),
                item["standard_relation"] == "<=",
            ),
        )
        aggregate_value = float(selected_measurement["upper_bound_nm"])
        aggregate_relation = str(selected_measurement["standard_relation"])
        aggregation_method = "conservative_largest_censored_upper_bound"

    eligible_compounds = [
        item for item in compounds
        if not str(item.get("data_validity_comment") or "").strip()
        and not _truthy_flag(item.get("potential_duplicate"))
    ]
    if not eligible_compounds:
        # Nested, already-aggregated records carry their quality decision on
        # each measurement rather than necessarily on the wrapper itself.
        eligible_compounds = [
            item for item in compounds
            if any(
                isinstance(measurement, dict)
                and measurement.get("aggregation_included") is not False
                for measurement in (item.get("activity_measurements") or [])
            )
        ]
    if not eligible_compounds:
        raise ValueError("activity group has no eligible representative")
    selected_activity_id = str(selected_measurement.get("activity_chembl_id") or "")
    matching_compounds = [
        item for item in eligible_compounds
        if selected_activity_id
        and str(item.get("activity_chembl_id") or "") == selected_activity_id
    ]
    representative = min(
        matching_compounds or eligible_compounds,
        key=lambda item: abs(
            float(item.get("activity_nm", aggregate_value)) - aggregate_value
        ),
    ).copy()
    aggregate_bounds = _activity_bound_fields(aggregate_value, aggregate_relation)
    values = [float(item["value_nm"]) for item in measurements]
    observed_values = [float(item["value_nm"]) for item in all_measurements]
    assay_ids = sorted({
        str(item.get("assay_chembl_id"))
        for item in all_measurements
        if item.get("assay_chembl_id")
    })
    chembl_ids = sorted({
        str(item.get("chembl_id"))
        for item in all_measurements
        if item.get("chembl_id")
    })
    eligible_assay_ids = sorted({
        str(item.get("assay_chembl_id"))
        for item in measurements
        if item.get("assay_chembl_id")
    })
    quarantined_assay_ids = sorted({
        str(item.get("assay_chembl_id"))
        for item in all_measurements
        if not item.get("aggregation_included") and item.get("assay_chembl_id")
    })
    assay_types = sorted({
        str(item.get("assay_type")) for item in all_measurements
        if item.get("assay_type") not in (None, "")
    })
    bao_formats = sorted({
        str(item.get("bao_format")) for item in all_measurements
        if item.get("bao_format") not in (None, "")
    })
    data_validity_count = sum(
        bool(str(item.get("data_validity_comment") or "").strip())
        for item in all_measurements
    )
    potential_duplicate_count = sum(
        _truthy_flag(item.get("potential_duplicate"))
        for item in all_measurements
    )
    selected_provenance_fields = (
        "assay_type", "assay_description", "bao_format", "bao_label",
        "document_chembl_id", "document_journal", "document_year",
        "data_validity_comment", "data_validity_description",
        "potential_duplicate", "standard_flag", "pchembl_value", "source_id",
        "target_chembl_id", "standard_units", "activity_type",
    )
    representative.update({
        "activity_nm": aggregate_value,
        "standard_relation": aggregate_relation,
        **aggregate_bounds,
        "activity_aggregation_method": aggregation_method,
        "activity_measurement_count": len(all_measurements),
        "activity_aggregation_measurement_count": len(measurements),
        "activity_quarantined_measurement_count": (
            len(all_measurements) - len(measurements)
        ),
        "activity_data_validity_quarantined_count": data_validity_count,
        "activity_potential_duplicate_quarantined_count": (
            potential_duplicate_count
        ),
        "activity_exact_measurement_count": len(exact),
        "activity_censored_measurement_count": len(measurements) - len(exact),
        "activity_range_nm": [min(values), max(values)],
        "activity_observed_range_nm": [
            min(observed_values), max(observed_values),
        ],
        "activity_measurements": all_measurements,
        "activity_assay_types": assay_types,
        "activity_bao_formats": bao_formats,
        "activity_assay_heterogeneity": (
            len(assay_ids) > 1 or len(assay_types) > 1 or len(bao_formats) > 1
        ),
        "activity_assay_strata": _activity_assay_strata(all_measurements),
        "activity_document_chembl_ids": sorted({
            str(item.get("document_chembl_id")) for item in all_measurements
            if item.get("document_chembl_id") not in (None, "")
        }),
        "assay_chembl_ids": assay_ids,
        "activity_eligible_assay_chembl_ids": eligible_assay_ids,
        "activity_quarantined_assay_chembl_ids": quarantined_assay_ids,
        "chembl_ids": chembl_ids,
        "assay_chembl_id": selected_measurement.get("assay_chembl_id", ""),
        "chembl_id": selected_measurement.get("chembl_id", ""),
        "activity_chembl_id": selected_measurement.get("activity_chembl_id", ""),
        **{
            field: selected_measurement.get(field)
            for field in selected_provenance_fields
        },
    })
    return representative


def _activity_assay_strata(measurements: list[dict]) -> list[dict]:
    """Describe assay-type/BAO strata without changing the primary summary."""
    groups: dict[tuple[str, str], list[dict]] = {}
    for measurement in measurements:
        key = (
            str(measurement.get("assay_type") or "unknown"),
            str(measurement.get("bao_format") or "unknown"),
        )
        groups.setdefault(key, []).append(measurement)

    strata = []
    for (assay_type, bao_format), group in sorted(groups.items()):
        eligible = [item for item in group if item.get("aggregation_included")]
        exact = [item for item in eligible if item.get("standard_relation") == "="]
        censored = [item for item in eligible if item.get("standard_relation") != "="]
        summary_value = None
        summary_relation = None
        summary_method = None
        if exact:
            values = sorted(float(item["value_nm"]) for item in exact)
            middle = len(values) // 2
            summary_value = (
                values[middle]
                if len(values) % 2
                else (values[middle - 1] + values[middle]) / 2
            )
            summary_relation = "="
            summary_method = "median_exact_measurements"
        elif censored:
            selected = max(
                censored,
                key=lambda item: (
                    float(item["upper_bound_nm"]),
                    item.get("standard_relation") == "<=",
                ),
            )
            summary_value = float(selected["upper_bound_nm"])
            summary_relation = str(selected.get("standard_relation"))
            summary_method = "conservative_largest_censored_upper_bound"
        strata.append({
            "assay_type": assay_type,
            "bao_format": bao_format,
            "assay_chembl_ids": sorted({
                str(item.get("assay_chembl_id")) for item in group
                if item.get("assay_chembl_id") not in (None, "")
            }),
            "document_chembl_ids": sorted({
                str(item.get("document_chembl_id")) for item in group
                if item.get("document_chembl_id") not in (None, "")
            }),
            "measurement_count": len(group),
            "aggregation_measurement_count": len(eligible),
            "quarantined_measurement_count": len(group) - len(eligible),
            "exact_measurement_count": len(exact),
            "censored_measurement_count": len(censored),
            "summary_nm": summary_value,
            "summary_relation": summary_relation,
            "summary_method": summary_method,
        })
    return strata


def _deduplicate_by_smiles(compounds: list[dict]) -> list[dict]:
    """Aggregate repeated measurements by SMILES with censoring preserved."""
    grouped = {}
    for comp in compounds:
        smi = comp["smiles"]
        grouped.setdefault(smi, []).append(comp)

    aggregated = []
    for measurements in grouped.values():
        try:
            aggregated.append(aggregate_activity_measurements(measurements))
        except ValueError as exc:
            if "quality-eligible" not in str(exc) and "eligible representative" not in str(exc):
                raise
    return aggregated


def _quality_group_audit(compounds: list[dict]) -> dict[str, int]:
    """Count SMILES groups retained or wholly quarantined by quality flags."""
    grouped: dict[str, list[dict]] = {}
    for compound in compounds:
        grouped.setdefault(str(compound.get("smiles") or ""), []).append(compound)
    quarantined = sum(
        not any(
            not str(item.get("data_validity_comment") or "").strip()
            and not _truthy_flag(item.get("potential_duplicate"))
            for item in group
        )
        for group in grouped.values()
    )
    return {
        "eligible_compound_groups": len(grouped) - quarantined,
        "quarantined_compound_groups": quarantined,
    }
