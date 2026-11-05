#!/usr/bin/env python3
"""Append target-assets protein/cavity embeddings; never update frozen v1 files.

Reuses run_batch's exact PDB/Foldseek/model recipe. A cavity vector is the
mean of full-protein residue embeddings selected by distance to cavity points,
not an independent learned embedding of the empty-space grid.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import fcntl
import io
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import time
import traceback

import numpy as np
import scipy
from scipy.spatial import cKDTree
import run_batch as core

ASSETS = core.configured_path('SAPROT_TARGET_ASSETS')
RBCAVITY = core.configured_path('SAPROT_RBCAVITY')
RBT_ROOT = core.configured_path('SAPROT_RBT_ROOT')
COHORTS = {'ADB': 'aligndockbench', 'PB': 'posebusters', 'PL': 'plinder'}


def read_cavity(path):
    """Read legacy little-endian rDock serialization, with bounded reads.

    RbtDockingSite: string title, min/max (3 doubles), border double, int nCav.
    Each RbtCavity: 3-double grid step, int nCoords, nCoords * 3 doubles.
    The remaining distance grid is not used for residue selection. Every real
    file is independently read/exported by native rbcavity before encoding.
    """
    raw = Path(path).read_bytes()
    position = 0

    def take(fmt):
        nonlocal position
        size = struct.calcsize('<' + fmt)
        if position + size > len(raw):
            raise ValueError('Truncated cavity.as')
        result = struct.unpack_from('<' + fmt, raw, position)
        position += size
        return result

    def title(expected):
        nonlocal position
        length, = take('i')
        if length != len(expected) or raw[position:position+length] != expected:
            raise ValueError('Unsupported cavity.as title/endian/layout')
        position += length

    title(b'RbtDockingSite')
    lower = np.array(take('3d'))
    upper = np.array(take('3d'))
    border, = take('d')
    count, = take('i')
    if not np.isfinite([*lower, *upper, border]).all() or (lower > upper).any() or border < 0 or not 1 <= count <= 10000:
        raise ValueError('Invalid cavity header')
    cavities = []
    for index in range(count):
        step = np.array(take('3d'))
        n, = take('i')
        if not np.isfinite(step).all() or (step <= 0).any() or n <= 0 or n > (len(raw)-position)//24:
            raise ValueError('Invalid/truncated cavity coordinate list')
        points = np.frombuffer(raw, dtype='<f8', count=n*3, offset=position).reshape(n, 3).copy()
        position += n*24
        if not np.isfinite(points).all() or len(np.unique(points, axis=0)) != n:
            raise ValueError('Invalid/duplicate cavity grid coordinates')
        cavities.append(dict(id=index+1, step=step, points=points))
    grid_offset = position
    title(b'RbtBaseGrid')
    all_points = np.concatenate([c['points'] for c in cavities])
    if not np.allclose(all_points.min(0), lower, atol=1e-8, rtol=0) or not np.allclose(all_points.max(0), upper, atol=1e-8, rtol=0):
        raise ValueError('Cavity bounds disagree with coordinates')
    return cavities, dict(format='rDock legacy binary little-endian int32/float64',
                         lower=lower.tolist(), upper=upper.tolist(), border_angstrom=border,
                         distance_grid_offset=grid_offset, distance_grid_bytes=len(raw)-grid_offset)


def native_validate(row, cavities, directory, args):
    """Compare every decoded coordinate with native rbcavity's InsightII dump."""
    with tempfile.TemporaryDirectory(prefix='.native_cavity_', dir=directory) as temporary:
        temp = Path(temporary)
        (temp/'cavity.as').symlink_to(Path(row['cavity_path']).resolve())
        (temp/'cavity.prm').write_text('RBT_PARAMETER_FILE_V1.00\nTITLE native cavity validation\nRECEPTOR_FILE '+row['protein_mol2_path']+'\n')
        command = [str(args.rbcavity), '-rcavity.prm', '-ras', '-d']
        proc = subprocess.run(command, cwd=temp, env=dict(os.environ, RBT_ROOT=str(args.rbt_root), RBT_HOME=str(temp)),
                              capture_output=True, text=True, timeout=args.foldseek_timeout)
        if proc.returncode or 'DOCKING SITE' not in proc.stdout:
            raise ValueError('Native rbcavity could not read cavity.as: '+(proc.stdout+proc.stderr)[-1000:])
        grids = sorted(temp.glob('cavity_cav*.grd'))
        if len(grids) != len(cavities):
            raise ValueError('Native cavity count differs')
        for cavity in cavities:
            lines = (temp/f'cavity_cav{cavity["id"]}.grd').read_text().splitlines()
            if lines[:2] != ['RBT FFT GRID', '(1F15.10)']:
                raise ValueError('Unsupported native grid export')
            intervals = np.array(list(map(int, lines[3].split())))
            limits = np.array(list(map(int, lines[4].split())))
            if len(limits) != 7 or limits[0] != 1:
                raise ValueError('Unsupported InsightII axis order')
            lower, upper = limits[1::2], limits[2::2]
            shape = upper-lower+1
            if (shape <= 0).any() or not np.array_equal(intervals, shape-1):
                raise ValueError('Invalid InsightII grid dimensions')
            values = np.fromstring(' '.join(lines[5:]), sep=' ')
            if len(values) != int(np.prod(shape)) or not np.isin(values, [0., 1.]).all():
                raise ValueError('Invalid InsightII occupancy data')
            zyx = np.array(np.unravel_index(np.flatnonzero(values), tuple(shape[::-1]))).T
            points = (zyx[:, ::-1]+lower)*cavity['step']
            expected = cavity['points']
            if len(points) != len(expected) or cKDTree(expected).query(points)[0].max() > 1e-6:
                raise ValueError('Decoded cavity coordinates differ from native rDock')
        return dict(binary=str(args.rbcavity), binary_sha256=core.sha_file(args.rbcavity),
                    command=command, exit_code=proc.returncode,
                    validation='all cavity coordinate sets equal native InsightII occupied points',
                    checked_cavity_count=len(cavities),
                    native_warning_count=sum('WARNING' in line for line in proc.stdout.splitlines()))


def residue_distances(residues, points):
    coordinates = np.asarray([a['xyz'] for r in residues for a in r['atoms']])
    sizes = np.array([len(r['atoms']) for r in residues])
    atom_distances = cKDTree(points).query(coordinates)[0]
    return np.minimum.reduceat(atom_distances, np.r_[0, np.cumsum(sizes)[:-1]])


def save_npz(path, **arrays):
    buffer = io.BytesIO()
    np.savez_compressed(buffer, **arrays)
    core.atomic_bytes(path, buffer.getvalue())


def prepare(row, args, foldseek_identity):
    if row['cavity_path'] and (args.rbcavity is None or args.rbt_root is None):
        raise ValueError('Cavity assets require --rbcavity and --rbt-root '
                         '(or SAPROT_RBCAVITY and SAPROT_RBT_ROOT)')
    plan, cached = core.prepare_case(row, args, foldseek_identity)
    directory = core.case_directory(args.output, row)
    recipe = dict(core_preparation_fingerprint=plan['preparation_fingerprint'],
                  extension_code_sha256=core.sha_file(__file__), cavity_sha256=row['cavity_sha256'],
                  protein_mol2_sha256=row['protein_mol2_sha256'], shell_angstrom=args.cavity_distance,
                  selection='any selected protein heavy atom within shell of any original cavity point; whole residue',
                  missing_policy='missing cavity stays missing; no ligand-derived substitute',
                  native_rbcavity_sha256=core.sha_file(args.rbcavity) if row['cavity_path'] else None,
                  scipy_version=scipy.__version__)
    fingerprint = core.stable_hash(recipe)
    pocket = core.read_json(directory/'pocket_prepared.json')
    if not args.overwrite and pocket and pocket.get('preparation_fingerprint') == fingerprint:
        geometry = directory/'cavity_points.npz'
        if not row['cavity_path'] or (geometry.is_file() and core.sha_file(geometry) == pocket.get('geometry_sha256')):
            return plan, pocket, cached
    n = plan['residue_count']
    if not row['cavity_path']:
        pocket = dict(preparation_fingerprint=fingerprint, recipe=recipe, status='missing_cavity_asset',
                      pocket_mask=[False]*n, component_masks=[], component_ids=[],
                      min_heavy_atom_distance_angstrom=None, native_validation=None,
                      cavity_count=0, point_count=0, pocket_residue_count=0)
    else:
        if core.sha_file(row['cavity_path']) != row['cavity_sha256'] or core.sha_file(row['protein_mol2_path']) != row['protein_mol2_sha256']:
            raise ValueError('Cavity/MOL2 changed since input manifest creation')
        cavities, header = read_cavity(row['cavity_path'])
        native = native_validate(row, cavities, directory, args)
        residues, _ = core.parse_pdb(row['protein_path'])
        if [(r['segment'],r['chain'],r['resseq'],r['insertion_code'],r['resname']) for r in residues] != [(r['segment'],r['chain'],r['resseq'],r['insertion_code'],r['resname']) for r in plan['residues']]:
            raise ValueError('Pocket mapping and model residue mapping disagree')
        distances = np.stack([residue_distances(residues, c['points']) for c in cavities])
        masks = distances <= args.cavity_distance
        union = masks.any(0)
        if not union.any():
            raise ValueError('No protein residues near cavity; possible coordinate-frame mismatch')
        save_npz(directory/'cavity_points.npz', points=np.concatenate([c['points'] for c in cavities]),
                 component_ids=np.concatenate([np.full(len(c['points']),c['id'],dtype=np.int32) for c in cavities]),
                 grid_steps=np.stack([c['step'] for c in cavities]))
        pocket = dict(preparation_fingerprint=fingerprint, recipe=recipe, status='available', header=header,
                      pocket_mask=union.tolist(), component_masks=masks.tolist(),
                      component_ids=[c['id'] for c in cavities], cavity_count=len(cavities),
                      point_count=sum(len(c['points']) for c in cavities), pocket_residue_count=int(union.sum()),
                      component_residue_counts=masks.sum(1).tolist(),
                      min_heavy_atom_distance_angstrom=distances.min(0).tolist(),
                      component_min_heavy_atom_distance_angstrom=distances.tolist(),
                      native_validation=native, geometry_sha256=core.sha_file(directory/'cavity_points.npz'))
    core.atomic_json(directory/'pocket_prepared.json', pocket)
    return plan, pocket, False


def encode(row, plan, pocket, encoder, args):
    directory = core.case_directory(args.output, row)
    fingerprint = core.stable_hash(dict(protein=plan['preparation_fingerprint'], pocket=pocket['preparation_fingerprint'], encoder=encoder.recipe))
    names = {'embeddings.npz','metadata.json','prepared.json','pocket_prepared.json'}
    if pocket['status'] == 'available':
        names |= {'cavity_points.npz','pocket_embeddings.npz'}
    complete = core.read_json(directory/'complete.json')
    if not args.overwrite and complete and complete.get('encoding_fingerprint') == fingerprint and set(complete.get('output_sha256',{})) == names:
        if all((directory/name).is_file() and core.sha_file(directory/name) == complete['output_sha256'][name] for name in names):
            return 'skipped'
    embeddings, windows = encoder.encode(plan,args)
    residues = plan['residues']
    chains = list(dict.fromkeys(r['chain'] for r in residues))
    chain_means = np.stack([embeddings[[r['chain']==c for r in residues]].mean(0) for c in chains])
    protein_mean = embeddings.mean(0)
    mask = np.array(pocket['pocket_mask'],dtype=bool)
    means = embeddings[mask].mean(0,keepdims=True) if mask.any() else np.empty((0,encoder.dimension),np.float32)
    arrays = dict(residue_embeddings=embeddings.astype(args.storage_dtype),
                  protein_embedding=protein_mean,protein_embedding_normalized=core.normalize(protein_mean),
                  chain_embeddings=chain_means,chain_embeddings_normalized=core.normalize(chain_means),chain_ids=np.array(chains),
                  pocket_embeddings=means,pocket_embeddings_normalized=core.normalize(means),pocket_mask=mask,
                  residue_chain=np.array([r['chain'] for r in residues]),residue_number=np.array([r['resseq'] for r in residues]),
                  residue_insertion_code=np.array([r['insertion_code'] for r in residues]),
                  residue_segment=np.array([r['segment'] for r in residues],dtype=np.int32),
                  residue_name=np.array([r['resname'] for r in residues]),residue_aa=np.array([r['aa'] for r in residues]),
                  ca_coordinates=np.array([r['ca'] for r in residues],dtype=np.float32),
                  bfactor_or_plddt=np.array([r['bfactor'] for r in residues],dtype=np.float32),
                  structural_token_mask=np.array([r['model_3di']=='#' for r in residues]))
    if not np.isfinite(arrays['residue_embeddings']).all():
        raise ValueError('Nonfinite stored embedding')
    save_npz(directory/'embeddings.npz',**arrays)
    if mask.any():
        component_masks = np.array(pocket['component_masks'],dtype=bool)
        valid = component_masks.any(1)
        component_means = np.stack([embeddings[m].mean(0) for m in component_masks[valid]])
        save_npz(directory/'pocket_embeddings.npz',pocket_embeddings=means,pocket_embeddings_normalized=core.normalize(means),
                 residue_embeddings=embeddings[mask].astype(args.storage_dtype),residue_indices=np.flatnonzero(mask),
                 component_embeddings=component_means,component_embeddings_normalized=core.normalize(component_means),
                 component_ids=np.array(pocket['component_ids'],dtype=np.int32)[valid],component_masks=component_masks[valid])
    metadata = dict(schema_version=1,encoding_fingerprint=fingerprint,source_record=row,encoder=encoder.recipe,
                    residue_count=len(residues),embedding_dimension=encoder.dimension,chain_ids=chains,
                    window_ranges=windows,long_fragment_windowing_used=any(len(f['residue_indices'])>args.window_size for f in plan['fragments']),
                    protein_preparation_fingerprint=plan['preparation_fingerprint'],pocket_preparation_fingerprint=pocket['preparation_fingerprint'],
                    pocket_status=pocket['status'],pocket_residue_count=int(mask.sum()),
                    pocket_definition=pocket['recipe']['selection'],cavity_shell_angstrom=args.cavity_distance,
                    pocket_embedding_semantics='mean of full-protein SaProt residue embeddings near original cavity points; not a learned grid-shape encoder',
                    cavity_count=pocket['cavity_count'],cavity_point_count=pocket['point_count'],
                    missing_backbone_residue_count=plan['incomplete_backbone_residue_count'],
                    omitted_protein_residues_without_CA=plan['omitted_protein_residues_without_CA'],warnings=plan['warnings'])
    core.atomic_json(directory/'metadata.json',metadata)
    core.atomic_json(directory/'complete.json',dict(encoding_fingerprint=fingerprint,
                     output_sha256={name:core.sha_file(directory/name) for name in sorted(names)},completed_at_unix=time.time()))
    return 'success'


def inventory(args):
    rows = []
    for directory in sorted(args.assets.iterdir()):
        if not directory.is_dir() or directory.name.split('_')[0] not in COHORTS:
            continue
        if args.target and directory.name not in args.target:
            continue
        metadata = json.loads((directory/'metadata.json').read_text())
        def asset(name):
            path = directory/name
            if not path.is_file():
                return '', ''
            actual = core.sha_file(path)
            recorded = metadata.get('files',{}).get(name,{}).get('sha256')
            if recorded and actual != recorded:
                raise ValueError('Asset differs from source metadata: '+str(path))
            return str(path.resolve()),actual
        protein,protein_sha = asset('protein.pdb')
        cavity,cavity_sha = asset('cavity.as')
        mol2,mol2_sha = asset('protein.mol2')
        if not protein or (cavity and not mol2):
            raise ValueError('Missing required protein or native-validation MOL2: '+str(directory))
        rows.append(dict(task_type='molecule_generation',dataset=COHORTS[directory.name.split('_')[0]],target=directory.name,
                         protein_path=protein,protein_sha256=protein_sha,cavity_path=cavity,cavity_sha256=cavity_sha,
                         protein_mol2_path=mol2,protein_mol2_sha256=mol2_sha,
                         source_metadata_path=str(directory/'metadata.json'),source_metadata_sha256=core.sha_file(directory/'metadata.json'),
                         source_cohort=metadata.get('cohort'),plddt_mask='False',plddt_threshold='',
                         structure_variant='generation_target_asset',input_validation_passed='True'))
    if not rows:
        raise ValueError('No targets selected')
    return rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--assets',type=Path,default=ASSETS,help='Target assets root; or set SAPROT_TARGET_ASSETS')
    parser.add_argument('--output',type=Path,default=core.DEFAULT_OUTPUT,help='Output root; or set SAPROT_OUTPUT')
    parser.add_argument('--model',type=Path,default=core.DEFAULT_MODEL,help='Local model directory; or set SAPROT_MODEL')
    parser.add_argument('--foldseek',type=Path,default=core.DEFAULT_FOLDSEEK,help='Foldseek executable; or set SAPROT_FOLDSEEK')
    parser.add_argument('--rbcavity',type=Path,default=RBCAVITY,help='rbcavity executable; or set SAPROT_RBCAVITY')
    parser.add_argument('--rbt-root',type=Path,default=RBT_ROOT,help='rDock data directory; or set SAPROT_RBT_ROOT')
    parser.add_argument('--device',default='auto')
    parser.add_argument('--dtype',choices=['auto','float32','float16','bfloat16'],default='auto')
    parser.add_argument('--storage-dtype',choices=['float16','float32'],default='float16')
    parser.add_argument('--threads',type=int,default=8)
    parser.add_argument('--prepare-workers',type=int,default=4)
    parser.add_argument('--foldseek-timeout',type=int,default=120)
    parser.add_argument('--window-size',type=int,default=1022)
    parser.add_argument('--window-overlap',type=int,default=128)
    parser.add_argument('--cavity-distance',type=float,default=4.0)
    parser.add_argument('--target',action='append')
    parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--overwrite',action='store_true')
    args=parser.parse_args()
    required_paths = {'assets': 'SAPROT_TARGET_ASSETS', 'output': 'SAPROT_OUTPUT',
                      'foldseek': 'SAPROT_FOLDSEEK'}
    if not args.prepare_only:
        required_paths['model'] = 'SAPROT_MODEL'
    for name, variable in required_paths.items():
        if getattr(args, name) is None:
            parser.error(f'Provide --{name} or set {variable}')
    for name in ['assets', 'output', 'model', 'foldseek', 'rbcavity', 'rbt_root']:
        value = getattr(args, name)
        if value is not None:
            setattr(args, name, value.expanduser())
    if min(args.threads,args.prepare_workers,args.foldseek_timeout,args.window_size)<1 or args.window_size>1022 or not 0<=args.window_overlap<args.window_size or not np.isfinite(args.cavity_distance) or args.cavity_distance<=0:
        parser.error('Invalid threads/window/cavity distance')
    namespace=args.output/'molecule_generation'
    namespace.mkdir(parents=True,exist_ok=True)
    started=time.time()
    with (namespace/'.run.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        rows=inventory(args)
        buffer=io.StringIO();writer=csv.DictWriter(buffer,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        core.atomic_bytes(namespace/'input_manifest.csv',buffer.getvalue().encode())
        states={r['target']:dict(task_type=r['task_type'],dataset=r['dataset'],target=r['target'],status='pending',
                                cavity_status='available' if r['cavity_path'] else 'missing_cavity_asset',
                                output_directory=str(core.case_directory(args.output,r))) for r in rows}
        def summary():
            counts={}
            for r in states.values(): counts[r['status']]=counts.get(r['status'],0)+1
            core.atomic_json(namespace/'run_summary.json',dict(selected=len(rows),prepare_only=args.prepare_only,status_counts=counts,
                             cavity_assets=sum(bool(r['cavity_path']) for r in rows),missing_cavity_assets=sum(not r['cavity_path'] for r in rows),
                             elapsed_seconds=round(time.time()-started,2),records=list(states.values())))
            buffer=io.StringIO();writer=csv.DictWriter(buffer,fieldnames=['task_type','dataset','target','status','cavity_status','output_directory','error'],extrasaction='ignore')
            writer.writeheader();writer.writerows(states.values());core.atomic_bytes(namespace/'embedding_index.csv',buffer.getvalue().encode())
        try:
            print(f'Selected {len(rows)} generation targets; cavities={sum(bool(r["cavity_path"]) for r in rows)}',flush=True)
            foldseek_identity=dict(binary_sha256=core.sha_file(args.foldseek),version=subprocess.check_output([str(args.foldseek),'version'],text=True,timeout=30).strip())
            plans={}
            with concurrent.futures.ThreadPoolExecutor(max_workers=args.prepare_workers) as pool:
                futures={pool.submit(prepare,r,args,foldseek_identity):r for r in rows}
                for number,future in enumerate(concurrent.futures.as_completed(futures),1):
                    row=futures[future];state=states[row['target']]
                    try:
                        plan,pocket,cached=future.result();plans[row['target']]=(plan,pocket)
                        state.update(status='prepared_cached' if cached else 'prepared',residue_count=plan['residue_count'],pocket_residue_count=pocket['pocket_residue_count'])
                    except Exception as exc:
                        state.update(status='failed_preparation',error=f'{type(exc).__name__}: {exc}')
                        print('PREPARATION FAILED '+row['target']+': '+str(exc),flush=True)
                    if number%10==0 or number==len(rows):
                        print(f'Prepared {number}/{len(rows)}; usable={len(plans)}',flush=True);summary()
            if not args.prepare_only and plans:
                encoder=core.Encoder(args)
                for number,row in enumerate(rows,1):
                    if row['target'] not in plans: continue
                    state=states[row['target']];began=time.time()
                    try:
                        state.update(status=encode(row,*plans[row['target']],encoder,args),seconds=round(time.time()-began,3))
                    except Exception as exc:
                        state.update(status='failed_encoding',error=f'{type(exc).__name__}: {exc}')
                        print(traceback.format_exc(),flush=True)
                        if encoder.device.type=='cuda': encoder.torch.cuda.empty_cache()
                    print(f'[{number}/{len(rows)}] {row["target"]}: {state["status"]} ({time.time()-began:.2f}s)',flush=True)
                    if number%5==0 or number==len(rows): summary()
        finally:
            summary()
        failures=sum(s['status'].startswith('failed') for s in states.values())
        print(f'Finished: selected={len(rows)}, failed={failures}; summary={namespace/"run_summary.json"}',flush=True)
        return 1 if failures else 0


if __name__=='__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print('Interrupted; rerun same command to resume.',file=sys.stderr)
        sys.exit(130)
