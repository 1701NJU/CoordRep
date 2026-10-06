"""
Invariance Tests

 canonicalize()
- bit-wise
- bit-wise
- rounded values
"""

import sys
import numpy as np
from pathlib import Path
from scipy.spatial.transform import Rotation
import copy

sys.path.insert(0, str(Path(__file__).parent.parent))

from coordrep import encode_molecule, CoordRepConfig
from coordrep.io.tmqm_reader import TMQMReader, RawMolecule, Atom


def random_rotation(coords: np.ndarray, seed: int = None) -> np.ndarray:
    """"""
    if seed is not None:
        np.random.seed(seed)
    
    R = Rotation.random().as_matrix()
    return coords @ R.T


def random_permutation(atoms: list, seed: int = None) -> list:
    """"""
    if seed is not None:
        np.random.seed(seed)
    
    indices = np.random.permutation(len(atoms))
    new_atoms = []
    
    for new_idx, old_idx in enumerate(indices):
        old_atom = atoms[old_idx]
        new_atoms.append(Atom(
            index=new_idx,
            element=old_atom.element,
            x=old_atom.x,
            y=old_atom.y,
            z=old_atom.z
        ))
    
    return new_atoms, indices


def add_noise(atoms: list, sigma: float = 0.01, seed: int = None) -> list:
    """"""
    if seed is not None:
        np.random.seed(seed)
    
    new_atoms = []
    for atom in atoms:
        noise = np.random.normal(0, sigma, 3)
        new_atoms.append(Atom(
            index=atom.index,
            element=atom.element,
            x=atom.x + noise[0],
            y=atom.y + noise[1],
            z=atom.z + noise[2]
        ))
    
    return new_atoms


def test_rotation_invariance(mol: RawMolecule, n_rotations: int = 100) -> dict:
    """"""
    config = CoordRepConfig.default()
    
    original_cc = encode_molecule(mol, config)
    original_string = original_cc.canonicalize().to_string()
    
    matches = 0
    mismatches = []
    
    for i in range(n_rotations):
        coords = mol.get_coords()
        rotated_coords = random_rotation(coords, seed=i)
        
        new_atoms = []
        for j, atom in enumerate(mol.atoms):
            new_atoms.append(Atom(
                index=j,
                element=atom.element,
                x=rotated_coords[j, 0],
                y=rotated_coords[j, 1],
                z=rotated_coords[j, 2]
            ))
        
        rotated_mol = RawMolecule(
            mol_id=mol.mol_id,
            atoms=new_atoms,
            bond_orders=mol.bond_orders
        )
        
        rotated_cc = encode_molecule(rotated_mol, config)
        rotated_string = rotated_cc.canonicalize().to_string()
        
        if rotated_string == original_string:
            matches += 1
        else:
            mismatches.append({
                'rotation': i,
                'original': original_string[:100],
                'rotated': rotated_string[:100]
            })
    
    return {
        'total': n_rotations,
        'matches': matches,
        'match_rate': matches / n_rotations,
 'mismatches': mismatches[:5] # 5
    }


def test_permutation_invariance(mol: RawMolecule, n_permutations: int = 100) -> dict:
    """"""
    config = CoordRepConfig.default()
    
    original_cc = encode_molecule(mol, config)
    original_string = original_cc.canonicalize().to_string()
    
    matches = 0
    mismatches = []
    
    for i in range(n_permutations):
        new_atoms, perm = random_permutation(mol.atoms, seed=i)
        
        new_bo = None
        if mol.bond_orders is not None:
            new_bo = mol.bond_orders[perm][:, perm]
        
        permuted_mol = RawMolecule(
            mol_id=mol.mol_id,
            atoms=new_atoms,
            bond_orders=new_bo
        )
        
        permuted_cc = encode_molecule(permuted_mol, config)
        permuted_string = permuted_cc.canonicalize().to_string()
        
        if permuted_string == original_string:
            matches += 1
        else:
            mismatches.append({
                'permutation': i,
                'original': original_string[:100],
                'permuted': permuted_string[:100]
            })
    
    return {
        'total': n_permutations,
        'matches': matches,
        'match_rate': matches / n_permutations,
        'mismatches': mismatches[:5]
    }


def test_noise_stability(mol: RawMolecule, sigmas: list = [0.01, 0.02], n_trials: int = 50) -> dict:
    """"""
    config = CoordRepConfig.default()
    
    results = {}
    
    for sigma in sigmas:
        original_cc = encode_molecule(mol, config)
        original_shape = original_cc.shape.values_rounded if original_cc.shape else None
        
        stable = 0
        jumps = 0
        
        for i in range(n_trials):
            noisy_atoms = add_noise(mol.atoms, sigma=sigma, seed=i)
            noisy_mol = RawMolecule(
                mol_id=mol.mol_id,
                atoms=noisy_atoms,
                bond_orders=mol.bond_orders
            )
            
            noisy_cc = encode_molecule(noisy_mol, config)
            noisy_shape = noisy_cc.shape.values_rounded if noisy_cc.shape else None
            
            if original_shape is not None and noisy_shape is not None:
                if np.allclose(original_shape, noisy_shape, atol=0.1):
                    stable += 1
                else:
                    jumps += 1
            else:
                stable += 1
        
        results[f'sigma={sigma}'] = {
            'stable': stable,
            'jumps': jumps,
            'stability_rate': stable / n_trials
        }
    
    return results


def run_invariance_tests(data_dir: str, n_samples: int = 50, n_trials: int = 100):
    """"""
    print("=" * 60)
    print("CoordRep Invariance Tests")
    print("=" * 60)
    
    reader = TMQMReader(data_dir)
    reader.load()
    
    rotation_results = []
    permutation_results = []
    noise_results = []
    
    count = 0
    for mol in reader.iter_molecules(limit=n_samples):
        if count >= n_samples:
            break
        
        print(f"\n[{count + 1}/{n_samples}] Testing {mol.mol_id}...")
        
        rot_result = test_rotation_invariance(mol, n_trials)
        rotation_results.append(rot_result['match_rate'])
        print(f"  Rotation: {rot_result['match_rate'] * 100:.1f}% match")
        
        perm_result = test_permutation_invariance(mol, n_trials)
        permutation_results.append(perm_result['match_rate'])
        print(f"  Permutation: {perm_result['match_rate'] * 100:.1f}% match")
        
        noise_result = test_noise_stability(mol, sigmas=[0.01, 0.02])
        noise_results.append(noise_result)
        for key, val in noise_result.items():
            print(f"  Noise {key}: {val['stability_rate'] * 100:.1f}% stable")
        
        count += 1
    
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    
    avg_rotation = np.mean(rotation_results) if rotation_results else 0
    avg_permutation = np.mean(permutation_results) if permutation_results else 0
    
    print(f"\nRotation Invariance: {avg_rotation * 100:.1f}% average")
    print(f"Permutation Invariance: {avg_permutation * 100:.1f}% average")
    
    if noise_results:
        for sigma in ['sigma=0.01', 'sigma=0.02']:
            rates = [r[sigma]['stability_rate'] for r in noise_results if sigma in r]
            avg = np.mean(rates) if rates else 0
            print(f"Noise Stability ({sigma}): {avg * 100:.1f}% average")
    
    print("\n" + "=" * 60)
    print("Acceptance Criteria")
    print("=" * 60)
    
    rotation_pass = avg_rotation >= 0.99
    permutation_pass = avg_permutation >= 0.99
    
    print(f"✓ Rotation Invariance ≥ 99%: {'PASS' if rotation_pass else 'FAIL'}")
    print(f"✓ Permutation Invariance ≥ 99%: {'PASS' if permutation_pass else 'FAIL'}")
    
    return {
        'rotation': avg_rotation,
        'permutation': avg_permutation,
        'passed': rotation_pass and permutation_pass
    }


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-dir', type=str, required=True)
    parser.add_argument('--n-samples', type=int, default=50)
    parser.add_argument('--n-trials', type=int, default=100)
    args = parser.parse_args()
    
    result = run_invariance_tests(args.data_dir, args.n_samples, args.n_trials)
    
    sys.exit(0 if result['passed'] else 1)





