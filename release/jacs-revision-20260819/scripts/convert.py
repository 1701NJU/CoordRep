#!/usr/bin/env python3
"""
Convert coordination complex structures to canonical CoordRep strings.

Supports:
    - Single XYZ file (+optional BO file)
    - Batch conversion of a tmQM dataset directory
    - Batch conversion of CIF files from COD

Usage:
    # Single structure
    python -m libcoordrep.scripts.convert --xyz complex.xyz --bo complex.BO

    # Batch tmQM directory
    python -m libcoordrep.scripts.convert --tmqm-dir tmQM-master/tmQM --output results.jsonl

    # Batch CIF directory
    python -m libcoordrep.scripts.convert --cif-dir data/cod_cif/ --output results.jsonl --limit 1000

    # Pipe-friendly: output only the CoordRep string
    python -m libcoordrep.scripts.convert --xyz complex.xyz --quiet
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from libcoordrep.coordrep import encode, batch_encode, CoordRepConfig
from libcoordrep.coordrep.io.tmqm_reader import TMQMReader


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert structures to canonical CoordRep strings"
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--xyz", type=str, help="Path to a single XYZ file")
    group.add_argument("--tmqm-dir", type=str, help="Path to tmQM data directory")
    group.add_argument("--cif-dir", type=str, help="Path to CIF file directory")

    parser.add_argument("--bo", type=str, default=None,
                        help="Path to Wiberg bond-order file (for --xyz mode)")
    parser.add_argument("--output", "-o", type=str, default=None,
                        help="Output file path (.jsonl for batch, stdout if omitted)")
    parser.add_argument("--limit", type=int, default=None,
                        help="Max number of structures to process (batch mode)")
    parser.add_argument("--bo-threshold", type=float, default=0.3,
                        help="Bond-order threshold for coordination detection (default: 0.3)")
    parser.add_argument("--quiet", "-q", action="store_true",
                        help="Only print the CoordRep string (single-file mode)")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Print detailed information")

    return parser.parse_args()


def convert_single(args):
    """Convert a single XYZ file."""
    config = CoordRepConfig(bo_threshold_fixed=args.bo_threshold)

    cc = encode(xyz_path=args.xyz, bo_path=args.bo, config=config)
    cc_can = cc.canonicalize()
    coordrep_str = cc_can.to_string()

    if args.quiet:
        print(coordrep_str)
        return

    print(f"Source:    {cc.source_id}")
    print(f"Metal:    {cc.metal.element} (ox={cc.metal.oxidation}, d={cc.metal.dcount})")

    if cc.graph:
        print(f"CN:       {len(cc.graph.donor_indices)}")

    if cc.shape:
        print(f"Shape:    {cc.shape}")

    print(f"Ligands:  {len(cc.ligands)}")
    for lig in cc.ligands:
        print(f"  {lig.lig_id}: {lig.smiles} (dent={lig.dent}, donors={lig.donor_elements})")

    if cc.constraints.trans_pairs:
        print(f"Trans:    {len(cc.constraints.trans_pairs)} pairs")
    if cc.constraints.cis_pairs:
        print(f"Cis:      {len(cc.constraints.cis_pairs)} pairs")

    if cc.issues:
        print(f"\nIssues ({len(cc.issues)}):")
        for issue in cc.issues:
            print(f"  {issue}")

    print(f"\nCoordRep: {coordrep_str}")

    if args.output:
        record = cc_can.to_dict()
        record['coordrep'] = coordrep_str
        with open(args.output, 'w') as f:
            json.dump(record, f, indent=2)
        print(f"\nSaved to {args.output}")


def convert_tmqm_batch(args):
    """Batch-convert a tmQM dataset directory."""
    config = CoordRepConfig(bo_threshold_fixed=args.bo_threshold)
    reader = TMQMReader(args.tmqm_dir)

    molecules = list(reader.iter_molecules(limit=args.limit))
    print(f"Loaded {len(molecules)} molecules from tmQM")

    complexes = batch_encode(molecules, config=config, show_progress=True)

    success = 0
    output_records = []
    for cc in complexes:
        cc_can = cc.canonicalize()
        coordrep_str = cc_can.to_string()
        is_ok = cc.metal.element != '?'
        if is_ok:
            success += 1

        record = {
            'source_id': cc.source_id,
            'metal': cc.metal.element,
            'cn': len(cc.graph.donor_indices) if cc.graph else 0,
            'coordrep': coordrep_str,
            'valid': is_ok,
        }
        output_records.append(record)

        if args.verbose and is_ok:
            print(f"  {cc.source_id}: {coordrep_str[:80]}...")

    print(f"\nResults: {success}/{len(complexes)} succeeded "
          f"({100*success/len(complexes):.1f}%)")

    if args.output:
        with open(args.output, 'w') as f:
            for rec in output_records:
                f.write(json.dumps(rec) + '\n')
        print(f"Saved to {args.output}")
    else:
        for rec in output_records:
            if rec['valid']:
                print(rec['coordrep'])


def convert_cif_batch(args):
    """Batch-convert CIF files."""
    from libcoordrep.coordrep.io.cif_reader import CIFReader

    config = CoordRepConfig(bo_threshold_fixed=args.bo_threshold)
    reader = CIFReader(args.cif_dir)
    molecules = list(reader.iter_molecules(limit=args.limit))
    print(f"Loaded {len(molecules)} molecules from CIF directory")

    complexes = batch_encode(molecules, config=config, show_progress=True)

    success = 0
    output_records = []
    for cc in complexes:
        cc_can = cc.canonicalize()
        coordrep_str = cc_can.to_string()
        is_ok = cc.metal.element != '?'
        if is_ok:
            success += 1

        record = {
            'source_id': cc.source_id,
            'metal': cc.metal.element,
            'cn': len(cc.graph.donor_indices) if cc.graph else 0,
            'coordrep': coordrep_str,
            'valid': is_ok,
        }
        output_records.append(record)

    print(f"\nResults: {success}/{len(complexes)} succeeded "
          f"({100*success/len(complexes):.1f}%)")

    if args.output:
        with open(args.output, 'w') as f:
            for rec in output_records:
                f.write(json.dumps(rec) + '\n')
        print(f"Saved to {args.output}")
    else:
        for rec in output_records:
            if rec['valid']:
                print(rec['coordrep'])


def main():
    args = parse_args()

    if args.xyz:
        convert_single(args)
    elif args.tmqm_dir:
        convert_tmqm_batch(args)
    elif args.cif_dir:
        convert_cif_batch(args)


if __name__ == '__main__':
    main()
