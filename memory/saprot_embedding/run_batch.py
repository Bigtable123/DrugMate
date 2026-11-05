#!/usr/bin/env python3
"""Encode actual task inputs with local SaProt weights. See README.md.

No input files are modified. The model sees one polymer fragment at a time;
protein and pocket means are residue-weighted across fragments and chains.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import fcntl
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import traceback

import numpy as np
from Bio.Data.PDBData import protein_letters_3to1_extended

def configured_path(name):
    """Read an optional deployment path without choosing a host-specific default."""
    value = os.environ.get(name, '').strip()
    return Path(value).expanduser() if value else None


DEFAULT_MANIFEST = configured_path('SAPROT_MANIFEST')
DEFAULT_OUTPUT = configured_path('SAPROT_OUTPUT')
DEFAULT_MODEL = configured_path('SAPROT_MODEL')
DEFAULT_FOLDSEEK = configured_path('SAPROT_FOLDSEEK')
STANDARD = dict(zip(
    'ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),
    'ARNDCQEGHILKMFPSTWYV'))
PARENT_NAMES = {v: k for k, v in STANDARD.items()}
AA_ALPHABET = set(PARENT_NAMES)
TRUE = {'true', '1', 'yes'}


def sha_file(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def atomic_bytes(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix='.' + path.name, dir=path.parent)
    try:
        with os.fdopen(fd, 'wb') as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def atomic_json(path, value):
    atomic_bytes(path, (json.dumps(value, indent=2, allow_nan=False) + '\n').encode())


def read_json(path):
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return None


def case_directory(output, row):
    parts = [row['task_type'], row['dataset'], row['target']]
    if any(not p or p in {'.', '..'} or '/' in p or '\\' in p for p in parts):
        raise ValueError('Unsafe task/dataset/target identifier')
    return Path(output).joinpath(*parts)


def parse_pdb(path):
    """Keep original order/IDs, choose one altloc, and preserve TER fragments.

    ATOM residues with CA are retained, including unknown modifications (# AA).
    HETATM amino acids are included only for recognized parent amino acids on
    chains containing ATOM protein residues. Cofactors/ions are excluded.
    """
    groups = {}
    segment = 0
    previous_chain = None
    atom_chains = set()
    for line in Path(path).read_text().splitlines():
        if line.startswith('ENDMDL'):
            break
        if line.startswith('TER'):
            segment += 1
            previous_chain = None
            continue
        if line[:6] not in {'ATOM  ', 'HETATM'}:
            continue
        if len(line) < 54:
            raise ValueError('Truncated PDB atom record')
        chain = line[21]
        if previous_chain is not None and chain != previous_chain:
            segment += 1
        previous_chain = chain
        name, resname = line[12:16].strip(), line[17:20].strip()
        element = line[76:78].strip().upper() or name.lstrip('0123456789')[:1].upper()
        if element in {'H', 'D'}:
            continue
        xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        if not all(np.isfinite(xyz)):
            raise ValueError('Non-finite atom coordinates')
        occupancy = float(line[54:60].strip() or 0)
        bfactor = float(line[60:66].strip() or 0)
        key = (segment, chain, line[22:26], line[26:27])
        groups.setdefault(key, []).append(dict(
            line=line, name=name, resname=resname, altloc=line[16],
            xyz=xyz, occupancy=occupancy, bfactor=bfactor, element=element,
            is_atom=line.startswith('ATOM  ')))
        if line.startswith('ATOM  ') and name == 'CA':
            atom_chains.add(chain)
    residues = []
    no_ca = []
    for (segment, chain, number, icode), atoms in groups.items():
        candidates = [a for a in atoms if a['name'] == 'CA' and (
            a['is_atom'] or (chain in atom_chains and a['resname'] in protein_letters_3to1_extended))]
        if not candidates:
            if any(a['is_atom'] and a['resname'] in protein_letters_3to1_extended for a in atoms):
                no_ca.append(dict(chain=chain, resseq=number.strip(), icode=icode.strip()))
            continue
        # Prefer greatest CA occupancy, then blank, A, and lexical altloc.
        ca = min(candidates, key=lambda a: (-a['occupancy'], a['altloc'] != ' ', a['altloc'] != 'A', a['altloc']))
        selected = {}
        for atom in atoms:
            if atom['resname'] != ca['resname'] or atom['altloc'] not in {' ', ca['altloc']}:
                continue
            old = selected.get(atom['name'])
            if old is None or (atom['altloc'] == ca['altloc'], atom['occupancy']) > (old['altloc'] == ca['altloc'], old['occupancy']):
                selected[atom['name']] = atom
        aa = protein_letters_3to1_extended.get(ca['resname'], '#')
        if aa not in AA_ALPHABET:
            aa = '#'
        residues.append(dict(
            index=len(residues), segment=segment, chain=chain,
            resseq=number.strip(), insertion_code=icode.strip(), resname=ca['resname'],
            aa=aa, altloc=ca['altloc'].strip(), record='ATOM' if ca['is_atom'] else 'HETATM',
            ca=ca['xyz'], bfactor=float(np.mean([a['bfactor'] for a in selected.values()])),
            missing_backbone=sorted({'N', 'CA', 'C'} - selected.keys()),
            atoms=list(selected.values())))
    if not residues:
        raise ValueError('No encodable protein residues with CA')
    return residues, no_ca


def normalized_segment_pdb(residues):
    """Canonical parent names/numbering make Foldseek output order checkable."""
    lines = []
    serial = 0
    for position, residue in enumerate(residues, 1):
        parent = PARENT_NAMES.get(residue['aa'], 'GLY')
        for atom in residue['atoms']:
            serial += 1
            original = atom['line'].ljust(80)
            line = ('ATOM  ' + f'{serial:5d}' + original[11:16] + ' ' + parent +
                    ' A' + f'{position:4d}' + ' ' + original[27:])
            # Foldseek treats inconsistent legacy SEGID fields as chain breaks,
            # even when PDB chain IDs agree (common in repaired terminal atoms).
            # This file already represents exactly one explicit input fragment.
            line = line[:72] + '    ' + line[76:]
            lines.append(line.rstrip() + '\n')
    return ''.join(lines) + 'TER\nEND\n'


def box_mask(residues, box):
    center = np.asarray(box['center'], dtype=float)
    size = np.asarray(box['size'], dtype=float)
    if center.shape != (3,) or size.shape != (3,) or not np.isfinite(center).all() or not np.isfinite(size).all() or (size <= 0).any():
        raise ValueError('Invalid docking box center/size')
    return [any(np.all(np.abs(np.asarray(a['xyz']) - center) <= size / 2) for a in r['atoms']) for r in residues]


def preparation_recipe(row, args, foldseek_identity):
    path = Path(row['protein_path'])
    current_hash = sha_file(path)
    if row.get('protein_sha256') and current_hash != row['protein_sha256']:
        raise ValueError('Protein changed since selection manifest: ' + str(path))
    if row.get('input_validation_passed', 'true').lower() not in TRUE:
        raise ValueError('Selection manifest marks this input as invalid')
    box = row.get('known_pocket_box_path', '') if row['task_type'] == 'docking' else ''
    if row['task_type'] == 'docking' and not box:
        raise ValueError('Docking input is missing known_pocket_box_path')
    effective = row.get('effective_receptor_pdbqt_path', '')
    if effective and row.get('effective_receptor_pdbqt_sha256') and sha_file(effective) != row['effective_receptor_pdbqt_sha256']:
        raise ValueError('Effective docking receptor changed since selection')
    mask = row.get('plddt_mask', '').lower() in TRUE
    threshold = float(row.get('plddt_threshold') or 70) if mask else None
    return dict(
        pipeline_code_sha256=sha_file(__file__), protein_sha256=current_hash, source_record=row,
        task_type=row['task_type'], box_sha256=sha_file(box) if box else None,
        plddt_mask=mask, plddt_threshold=threshold,
        foldseek=foldseek_identity, parser='parent_amino_acids_TER_fragments_v1',
        unknown_aa_policy='AA #; GLY geometry surrogate only for structure conversion',
        incomplete_backbone_policy='retain residue with masked 3Di #',
        foldseek_timeout_seconds=args.foldseek_timeout)


def prepare_case(row, args, foldseek_identity):
    directory = case_directory(args.output, row)
    directory.mkdir(parents=True, exist_ok=True)
    recipe = preparation_recipe(row, args, foldseek_identity)
    fingerprint = stable_hash(recipe)
    cached = read_json(directory / 'prepared.json')
    if not args.overwrite and cached and cached.get('preparation_fingerprint') == fingerprint:
        return cached, True
    residues, no_ca = parse_pdb(row['protein_path'])
    fragments = []
    for residue in residues:
        if not fragments or fragments[-1]['original_segment'] != residue['segment']:
            fragments.append(dict(id=f'fragment_{len(fragments):04d}', original_segment=residue['segment'],
                                  chain=residue['chain'], residue_indices=[]))
        fragments[-1]['residue_indices'].append(residue['index'])
    geometry = ['#'] * len(residues)
    warnings = []
    with tempfile.TemporaryDirectory(prefix='.foldseek_', dir=directory) as temporary:
        tmp = Path(temporary)
        expected = {}
        for fragment in fragments:
            valid = [residues[i] for i in fragment['residue_indices'] if not residues[i]['missing_backbone']]
            if not valid:
                warnings.append(fragment['id'] + ': no complete N/CA/C backbone; all 3Di masked')
                continue
            name = fragment['id'] + '.pdb'
            (tmp / name).write_text(normalized_segment_pdb(valid))
            expected[name + '_A'] = valid
        if expected:
            output = tmp / 'descriptors.tsv'
            command = [str(args.foldseek), 'structureto3didescriptor', '-v', '0', '--threads', '1',
                       '--chain-name-mode', '1', '--file-include', r'.*\.pdb$', str(tmp), str(output)]
            proc = subprocess.run(command, capture_output=True, text=True, timeout=args.foldseek_timeout)
            if proc.returncode:
                raise RuntimeError('Foldseek failed: ' + (proc.stderr + proc.stdout)[-2000:])
            found = set()
            for line in output.read_text().splitlines():
                fields = line.split('\t')
                if len(fields) < 3:
                    raise ValueError('Malformed Foldseek descriptor')
                descriptor, aa, di = fields[:3]
                descriptor = Path(descriptor.split(' ')[0]).name
                if descriptor not in expected or descriptor in found:
                    raise ValueError('Unexpected or duplicate Foldseek fragment: ' + descriptor)
                found.add(descriptor)
                valid = expected[descriptor]
                sequence = ''.join(r['aa'] if r['aa'] != '#' else 'G' for r in valid)
                if aa != sequence or len(di) != len(valid):
                    raise ValueError(f'Foldseek residue mismatch in {descriptor}: expected {len(valid)}, got {len(aa)}; refusing ambiguous residue mapping')
                for residue, token in zip(valid, di):
                    geometry[residue['index']] = token.lower() if token.lower() in 'pynwrqhgdlvtmfsaeikc' else '#'
            for descriptor in expected.keys() - found:
                if len(expected[descriptor]) > 2:
                    raise ValueError('Foldseek omitted fragment: ' + descriptor)
                warnings.append(descriptor + ': tiny fragment omitted by Foldseek; 3Di masked')
    for residue, di in zip(residues, geometry):
        residue['geometry_3di'] = di
        residue['confidence_masked'] = bool(recipe['plddt_mask'] and residue['bfactor'] < recipe['plddt_threshold'])
        residue['model_3di'] = '#' if residue['confidence_masked'] else di
        residue['sa_token'] = residue['aa'] + residue['model_3di']
    box = read_json(row['known_pocket_box_path']) if row['task_type'] == 'docking' else None
    pocket = box_mask(residues, box) if box is not None else [False] * len(residues)
    if row['task_type'] == 'docking' and not any(pocket):
        raise ValueError('No protein residues in the actual docking box')
    for residue, inside in zip(residues, pocket):
        residue['in_docking_box'] = inside
        del residue['atoms']
    for fragment in fragments:
        fragment['sa_sequence'] = ''.join(residues[i]['sa_token'] for i in fragment['residue_indices'])
    prepared = dict(
        preparation_fingerprint=fingerprint, preparation_recipe=recipe,
        source_record=row, residue_count=len(residues), fragments=fragments, residues=residues,
        pocket_box=box, pocket_residue_count=sum(pocket),
        omitted_protein_residues_without_CA=no_ca, warnings=warnings,
        modified_amino_acid_count=sum(r['resname'] not in STANDARD for r in residues),
        included_HETATM_amino_acid_count=sum(r['record'] == 'HETATM' for r in residues),
        incomplete_backbone_residue_count=sum(bool(r['missing_backbone']) for r in residues),
        unknown_amino_acid_count=sum(r['aa'] == '#' for r in residues),
        confidence_masked_residue_count=sum(r['confidence_masked'] for r in residues),
        geometry_masked_residue_count=sum(r['geometry_3di'] == '#' for r in residues))
    atomic_json(directory / 'prepared.json', prepared)
    return prepared, False


def window_ranges(length, size, overlap):
    if not 0 <= overlap < size:
        raise ValueError('Window overlap must satisfy 0 <= overlap < window size')
    if length <= size:
        return [(0, length)]
    starts = list(range(0, length - size + 1, size - overlap))
    if starts[-1] != length - size:
        starts.append(length - size)
    return [(start, start + size) for start in starts]


def normalize(array):
    norm = np.linalg.norm(array, axis=-1, keepdims=True)
    return array / np.maximum(norm, 1e-12)


class Encoder:
    def __init__(self, args):
        import torch
        import transformers
        from transformers import EsmModel, EsmTokenizer
        self.torch = torch
        torch.set_num_threads(args.threads)
        torch.set_num_interop_threads(1)
        device = ('cuda' if torch.cuda.is_available() else 'cpu') if args.device == 'auto' else args.device
        self.device = torch.device(device)
        if self.device.type == 'cuda' and not torch.cuda.is_available():
            raise RuntimeError('CUDA requested, but no usable GPU is visible; use --device cpu')
        dtype_name = ('float16' if self.device.type == 'cuda' else 'float32') if args.dtype == 'auto' else args.dtype
        if self.device.type == 'cpu' and dtype_name == 'float16':
            raise ValueError('CPU float16 is unsupported; use --dtype float32')
        self.dtype = getattr(torch, dtype_name)
        self.tokenizer = EsmTokenizer.from_pretrained(str(args.model), local_files_only=True)
        weight_files = sorted(args.model.glob('*.safetensors')) or sorted(args.model.glob('pytorch_model*.bin'))
        if not weight_files:
            raise ValueError('No Hugging Face weights found in model directory')
        identity_files = weight_files + [args.model / name for name in [
            'config.json', 'vocab.txt', 'tokenizer_config.json', 'special_tokens_map.json'] if (args.model / name).is_file()]
        self.recipe = dict(model_path=str(args.model.resolve()),
                           model_files_sha256={p.name: sha_file(p) for p in identity_files},
                           torch_version=torch.__version__, transformers_version=transformers.__version__,
                           device_type=self.device.type, inference_dtype=dtype_name,
                           storage_dtype=args.storage_dtype, window_size=args.window_size,
                           window_overlap=args.window_overlap, window_merge='uniform mean in overlaps',
                           pooling='mean over residues; exclude CLS/EOS/PAD',
                           context='independent TER/chain fragments; all fragments retained')
        print(f'Loading model on {self.device} ({dtype_name})', flush=True)
        self.model = EsmModel.from_pretrained(str(args.model), local_files_only=True,
                                            add_pooling_layer=False, torch_dtype=self.dtype).eval().to(self.device)
        self.dimension = self.model.config.hidden_size
        if args.window_size > 1022:
            raise ValueError('Use --window-size <= 1022 for the validated encoding recipe')

    def encode(self, prepared, args):
        torch = self.torch
        result = np.zeros((prepared['residue_count'], self.dimension), dtype=np.float32)
        windows = []
        with torch.inference_mode():
            for fragment in prepared['fragments']:
                tokens = [prepared['residues'][i]['sa_token'] for i in fragment['residue_indices']]
                total = np.zeros((len(tokens), self.dimension), dtype=np.float32)
                counts = np.zeros(len(tokens), dtype=np.int32)
                for start, end in window_ranges(len(tokens), args.window_size, args.window_overlap):
                    # Supply known tokens directly: one input token per residue.
                    ids = self.tokenizer.convert_tokens_to_ids(tokens[start:end])
                    if any(i is None or i == self.tokenizer.unk_token_id for i in ids):
                        raise ValueError('SA token not present in the model vocabulary')
                    ids = [self.tokenizer.cls_token_id] + ids + [self.tokenizer.eos_token_id]
                    input_ids = torch.tensor([ids], dtype=torch.long, device=self.device)
                    hidden = self.model(input_ids=input_ids, attention_mask=torch.ones_like(input_ids),
                                        output_hidden_states=False, output_attentions=False).last_hidden_state
                    values = hidden[0, 1:-1].float().cpu().numpy()
                    if values.shape != (end - start, self.dimension) or not np.isfinite(values).all():
                        raise ValueError('Invalid model output shape or non-finite embedding')
                    total[start:end] += values
                    counts[start:end] += 1
                    windows.append(dict(fragment=fragment['id'], start=start, end=end))
                if (counts == 0).any():
                    raise ValueError('Windowing left unencoded residues')
                result[fragment['residue_indices']] = total / counts[:, None]
        return result, windows


def output_valid(directory, fingerprint):
    completion = read_json(directory / 'complete.json')
    if not completion or completion.get('encoding_fingerprint') != fingerprint:
        return False
    return all((directory / name).is_file() and sha_file(directory / name) == digest
               for name, digest in completion.get('output_sha256', {}).items()) and set(completion.get('output_sha256', {})) == {'embeddings.npz', 'metadata.json'}


def encode_case(row, prepared, encoder, args):
    directory = case_directory(args.output, row)
    fingerprint = stable_hash(dict(preparation_fingerprint=prepared['preparation_fingerprint'], encoder=encoder.recipe))
    if not args.overwrite and output_valid(directory, fingerprint):
        return 'skipped'
    embeddings, windows = encoder.encode(prepared, args)
    residues = prepared['residues']
    chain_ids = list(dict.fromkeys(r['chain'] for r in residues))
    chain_means = np.stack([embeddings[[r['chain'] == chain for r in residues]].mean(axis=0) for chain in chain_ids])
    protein_mean = embeddings.mean(axis=0)
    mask = np.asarray([r['in_docking_box'] for r in residues], dtype=bool)
    pocket_means = embeddings[mask].mean(axis=0, keepdims=True) if mask.any() else np.empty((0, encoder.dimension), np.float32)
    arrays = dict(
        residue_embeddings=embeddings.astype(args.storage_dtype),
        protein_embedding=protein_mean, protein_embedding_normalized=normalize(protein_mean),
        chain_embeddings=chain_means, chain_embeddings_normalized=normalize(chain_means),
        chain_ids=np.asarray(chain_ids), pocket_embeddings=pocket_means,
        pocket_embeddings_normalized=normalize(pocket_means), pocket_mask=mask,
        residue_chain=np.asarray([r['chain'] for r in residues]),
        residue_number=np.asarray([r['resseq'] for r in residues]),
        residue_insertion_code=np.asarray([r['insertion_code'] for r in residues]),
        residue_segment=np.asarray([r['segment'] for r in residues], dtype=np.int32),
        residue_name=np.asarray([r['resname'] for r in residues]),
        residue_aa=np.asarray([r['aa'] for r in residues]),
        ca_coordinates=np.asarray([r['ca'] for r in residues], dtype=np.float32),
        bfactor_or_plddt=np.asarray([r['bfactor'] for r in residues], dtype=np.float32),
        structural_token_mask=np.asarray([r['model_3di'] == '#' for r in residues]))
    if not np.isfinite(arrays['residue_embeddings']).all():
        raise ValueError('Storage dtype caused non-finite embeddings')
    buffer = io.BytesIO()
    np.savez_compressed(buffer, **arrays)
    atomic_bytes(directory / 'embeddings.npz', buffer.getvalue())
    metadata = dict(
        schema_version=1, encoding_fingerprint=fingerprint,
        preparation_fingerprint=prepared['preparation_fingerprint'], source_record=row,
        encoder=encoder.recipe, residue_count=len(residues), embedding_dimension=encoder.dimension,
        chain_ids=chain_ids, pocket_residue_count=int(mask.sum()),
        residue_mapping='prepared.json residues; index equals residue_embeddings row',
        window_ranges=windows, long_fragment_windowing_used=any(
            len(f['residue_indices']) > args.window_size for f in prepared['fragments']),
        confidence_masked_residue_count=prepared['confidence_masked_residue_count'],
        incomplete_backbone_residue_count=prepared['incomplete_backbone_residue_count'],
        geometry_masked_residue_count=prepared['geometry_masked_residue_count'],
        included_HETATM_amino_acid_count=prepared['included_HETATM_amino_acid_count'],
        omitted_protein_residues_without_CA=prepared['omitted_protein_residues_without_CA'],
        warnings=prepared['warnings'])
    atomic_json(directory / 'metadata.json', metadata)
    atomic_json(directory / 'complete.json', dict(
        encoding_fingerprint=fingerprint,
        output_sha256={name: sha_file(directory / name) for name in ['embeddings.npz', 'metadata.json']},
        completed_at_unix=time.time()))
    return 'success'


def select_rows(args):
    with args.manifest.open(newline='') as handle:
        rows = list(csv.DictReader(handle))
    required = {'task_type', 'dataset', 'target', 'protein_path'}
    if not rows or not required.issubset(rows[0]):
        raise ValueError('Manifest is empty or missing required columns')
    rows = [r for r in rows if (args.task_type == 'all' or r['task_type'] == args.task_type)
            and (not args.dataset or r['dataset'] in args.dataset)
            and (not args.target or r['target'] in args.target)]
    if args.limit:
        rows = rows[:args.limit]
    identities = [(r['task_type'], r['dataset'], r['target']) for r in rows]
    if len(set(identities)) != len(identities):
        raise ValueError('Duplicate task/dataset/target records in manifest')
    if not rows:
        raise ValueError('No records matched the requested filters')
    for row in rows:
        case_directory(args.output, row)
    return rows


def arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, default=DEFAULT_MANIFEST, help='Input CSV; or set SAPROT_MANIFEST')
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT, help='Output root; or set SAPROT_OUTPUT')
    parser.add_argument('--model', type=Path, default=DEFAULT_MODEL, help='Local model directory; or set SAPROT_MODEL')
    parser.add_argument('--foldseek', type=Path, default=DEFAULT_FOLDSEEK, help='Foldseek executable; or set SAPROT_FOLDSEEK')
    parser.add_argument('--device', default='auto', help='auto, cpu, cuda, or cuda:N')
    parser.add_argument('--dtype', choices=['auto', 'float32', 'float16', 'bfloat16'], default='auto')
    parser.add_argument('--storage-dtype', choices=['float16', 'float32'], default='float16')
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--prepare-workers', type=int, default=4)
    parser.add_argument('--foldseek-timeout', type=int, default=120)
    parser.add_argument('--window-size', type=int, default=1022)
    parser.add_argument('--window-overlap', type=int, default=128)
    parser.add_argument('--task-type', choices=['all', 'docking', 'pocket_prediction'], default='all')
    parser.add_argument('--dataset', action='append', help='Repeat to select multiple datasets')
    parser.add_argument('--target', action='append', help='Repeat to select multiple target IDs')
    parser.add_argument('--limit', type=int)
    parser.add_argument('--prepare-only', action='store_true', help='Validate/map/convert all selected structures without loading SaProt')
    parser.add_argument('--overwrite', action='store_true', help='Recompute matching cached preparation/embeddings')
    args = parser.parse_args(argv)
    required_paths = ['manifest', 'output', 'foldseek']
    if not args.prepare_only:
        required_paths.append('model')
    for name in required_paths:
        if getattr(args, name) is None:
            parser.error(f'Provide --{name} or set SAPROT_{name.upper()}')
    for name in ['manifest', 'output', 'model', 'foldseek']:
        value = getattr(args, name)
        if value is not None:
            setattr(args, name, value.expanduser())
    if min(args.threads, args.prepare_workers, args.foldseek_timeout, args.window_size) < 1:
        parser.error('Thread/worker/time/window values must be positive')
    if args.limit is not None and args.limit < 1:
        parser.error('--limit must be positive')
    if args.window_size > 1022 or not 0 <= args.window_overlap < args.window_size:
        parser.error('Require window-size <= 1022 and 0 <= overlap < window-size')
    return args


def main(argv=None):
    args = arguments(argv)
    args.output.mkdir(parents=True, exist_ok=True)
    started = time.time()
    rows = select_rows(args)
    results = {str(case_directory(args.output, r)): dict(
        task_type=r['task_type'], dataset=r['dataset'], target=r['target'], status='pending',
        output_directory=str(case_directory(args.output, r))) for r in rows}
    prepared = {}

    def save_summary():
        counts = {}
        for result in results.values():
            counts[result['status']] = counts.get(result['status'], 0) + 1
        atomic_json(args.output / 'run_summary.json', dict(
            manifest=str(args.manifest.resolve()), manifest_sha256=sha_file(args.manifest),
            prepare_only=args.prepare_only, selected=len(rows), status_counts=counts,
            elapsed_seconds=round(time.time() - started, 2), records=list(results.values())))
        buffer = io.StringIO()
        writer = csv.DictWriter(buffer, fieldnames=['task_type', 'dataset', 'target', 'status', 'output_directory', 'error'], extrasaction='ignore')
        writer.writeheader()
        writer.writerows(results.values())
        atomic_bytes(args.output / 'embedding_index.csv', buffer.getvalue().encode())

    with (args.output / '.run.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('Another writer is using this output directory')
        foldseek_identity = dict(binary_sha256=sha_file(args.foldseek), version=subprocess.check_output(
            [str(args.foldseek), 'version'], text=True, timeout=30).strip())
        print(f'Selected {len(rows)} inputs; output={args.output}', flush=True)
        try:
            with concurrent.futures.ThreadPoolExecutor(max_workers=args.prepare_workers) as pool:
                futures = {pool.submit(prepare_case, r, args, foldseek_identity): r for r in rows}
                for completed, future in enumerate(concurrent.futures.as_completed(futures), 1):
                    row = futures[future]
                    key = str(case_directory(args.output, row))
                    try:
                        plan, cached = future.result()
                        prepared[key] = plan
                        results[key]['status'] = 'prepared_cached' if cached else 'prepared'
                        results[key]['residue_count'] = plan['residue_count']
                    except Exception as exc:
                        results[key].update(status='failed_preparation', error=f'{type(exc).__name__}: {exc}')
                        print(f'PREP FAILED {row["target"]}: {exc}', flush=True)
                    if completed % 25 == 0 or completed == len(rows):
                        print(f'Prepared {completed}/{len(rows)}; usable={len(prepared)}', flush=True)
                        save_summary()
            if not args.prepare_only and prepared:
                encoder = Encoder(args)
                for number, row in enumerate(rows, 1):
                    key = str(case_directory(args.output, row))
                    if key not in prepared:
                        continue
                    began = time.time()
                    try:
                        status = encode_case(row, prepared[key], encoder, args)
                        results[key].update(status=status, seconds=round(time.time() - began, 3))
                    except Exception as exc:
                        results[key].update(status='failed_encoding', error=f'{type(exc).__name__}: {exc}')
                        print(traceback.format_exc(), file=sys.stderr, flush=True)
                        if encoder.device.type == 'cuda':
                            encoder.torch.cuda.empty_cache()
                    print(f'[{number}/{len(rows)}] {row["task_type"]}/{row["dataset"]}/{row["target"]}: {results[key]["status"]} ({time.time() - began:.2f}s)', flush=True)
                    if number % 10 == 0 or number == len(rows) or results[key]['status'].startswith('failed'):
                        save_summary()
        finally:
            save_summary()
    failures = sum(r['status'].startswith('failed') for r in results.values())
    print(f'Finished: selected={len(rows)}, failed={failures}; see {args.output / "run_summary.json"}', flush=True)
    return 1 if failures else 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print('Interrupted. Rerun the same command to resume.', file=sys.stderr)
        sys.exit(130)
