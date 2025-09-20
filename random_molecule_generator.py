"""Random molecule generator using RDKit and SELFIES.

The script builds molecules by combining a basis structure that contains
placeholder dummy atoms (e.g. ``[*:1]``) with fragment libraries for each
R-group.  Each generated molecule is returned both as SMILES and SELFIES to
simplify downstream processing in environments such as Google Colab.

Example (Colab):

.. code-block:: python

    !pip install rdkit selfies

    from random_molecule_generator import generate_random_molecules

    basis = "c1cc([*:1])ccc1[*:2]"  # para-disubstituted benzene
    r_groups = {
        1: ["[*:1]C", "[*:1]O", "[*:1]N"],
        2: ["[*:2]O", "[*:2]C(=O)O", "[*:2]N"],
    }

    molecules = generate_random_molecules(basis, r_groups, n_molecules=5, seed=0)
    for mol in molecules:
        print(mol.smiles, mol.selfies, mol.attachments)

"""

from __future__ import annotations

import random
from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from rdkit import Chem
import selfies as sf


@dataclass
class GeneratedMolecule:
    """Container storing the outputs of a randomised R-group combination."""

    smiles: str
    selfies: str
    attachments: Dict[int, str]
    mol: Chem.Mol = field(repr=False)


def _parse_smiles(smiles: str, label: str) -> Chem.Mol:
    """Parse a SMILES string and raise a helpful error if sanitisation fails."""

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse {label} SMILES: {smiles}")
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:  # RDKit can raise ValueError or SanitException
        raise ValueError(f"SMILES for {label} could not be sanitised: {smiles}") from exc
    return mol


def _find_placeholder_map_nums(mol: Chem.Mol) -> List[int]:
    """Return the sorted atom-map numbers that mark R-group positions."""

    map_nums: List[int] = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "*":
            map_num = atom.GetAtomMapNum()
            if map_num == 0:
                raise ValueError(
                    "Dummy atoms in the basis structure must use atom-map numbers "
                    "(e.g. [*:1]) to mark R-groups."
                )
            map_nums.append(map_num)
    if not map_nums:
        raise ValueError("The basis structure does not contain any dummy atoms for R-groups.")
    if len(map_nums) != len(set(map_nums)):
        raise ValueError("Each R-group in the basis structure must use a unique atom-map number.")
    return sorted(map_nums)


def _prepare_r_groups(r_group_smiles: Mapping[int, Sequence[str]]) -> Dict[int, List[Tuple[str, Chem.Mol]]]:
    """Parse and cache every fragment in the R-group library."""

    library: Dict[int, List[Tuple[str, Chem.Mol]]] = {}
    for map_num, smiles_list in r_group_smiles.items():
        if not smiles_list:
            raise ValueError(f"No substituents defined for R{map_num}.")
        parsed: List[Tuple[str, Chem.Mol]] = []
        for smi in smiles_list:
            parsed.append((smi, _parse_smiles(smi, f"R{map_num} fragment")))
        library[map_num] = parsed
    return library


def _get_single_neighbor_index(mol: Chem.RWMol, atom_idx: int) -> Tuple[int, Chem.rdchem.BondType]:
    """Return the sole neighbour of a dummy atom and the connecting bond type."""

    atom = mol.GetAtomWithIdx(atom_idx)
    neighbors = list(atom.GetNeighbors())
    if len(neighbors) != 1:
        raise ValueError(
            "Dummy atoms must have exactly one neighbour so the attachment point is well defined."
        )
    neighbor_idx = neighbors[0].GetIdx()
    bond = mol.GetBondBetweenAtoms(atom_idx, neighbor_idx)
    bond_type = bond.GetBondType() if bond is not None else Chem.rdchem.BondType.SINGLE
    return neighbor_idx, bond_type


def attach_fragment(base_mol: Chem.Mol, fragment_mol: Chem.Mol, map_num: int) -> Chem.Mol:
    """Attach a fragment (with ``[*:map_num]``) onto the corresponding basis site."""

    base = Chem.RWMol(base_mol)
    frag = Chem.RWMol(fragment_mol)

    base_dummy_idx = next(
        (atom.GetIdx() for atom in base.GetAtoms() if atom.GetSymbol() == "*" and atom.GetAtomMapNum() == map_num),
        None,
    )
    if base_dummy_idx is None:
        raise ValueError(f"Basis molecule does not contain an R{map_num} attachment point.")
    frag_dummy_idx = next(
        (atom.GetIdx() for atom in frag.GetAtoms() if atom.GetSymbol() == "*" and atom.GetAtomMapNum() == map_num),
        None,
    )
    if frag_dummy_idx is None:
        raise ValueError(f"Fragment for R{map_num} does not define an attachment point [*:{map_num}].")

    base_neighbor_idx, base_bond_type = _get_single_neighbor_index(base, base_dummy_idx)
    frag_neighbor_idx, frag_bond_type = _get_single_neighbor_index(frag, frag_dummy_idx)

    bond_type = base_bond_type
    if bond_type == Chem.rdchem.BondType.SINGLE and frag_bond_type != Chem.rdchem.BondType.SINGLE:
        bond_type = frag_bond_type
    if bond_type == Chem.rdchem.BondType.ZERO:
        bond_type = Chem.rdchem.BondType.SINGLE

    combined = Chem.RWMol(Chem.CombineMols(base, frag))
    frag_offset = base.GetNumAtoms()
    combined_base_neighbor = base_neighbor_idx
    combined_frag_neighbor = frag_offset + frag_neighbor_idx
    combined.AddBond(combined_base_neighbor, combined_frag_neighbor, bond_type)

    for idx in sorted([base_dummy_idx, frag_offset + frag_dummy_idx], reverse=True):
        combined.RemoveAtom(idx)

    mol = combined.GetMol()
    Chem.SanitizeMol(mol)
    return mol


def _assemble_molecule(
    basis: Chem.Mol,
    placeholder_map_nums: Sequence[int],
    r_group_library: Mapping[int, Sequence[Tuple[str, Chem.Mol]]],
    rng: random.Random,
) -> Tuple[Chem.Mol, Dict[int, str]]:
    """Combine the basis with randomly selected fragments for each R-group."""

    working = Chem.Mol(basis)
    attachments: Dict[int, str] = {}
    for map_num in placeholder_map_nums:
        options = r_group_library[map_num]
        smiles, fragment = rng.choice(options)
        attachments[map_num] = smiles
        working = attach_fragment(working, fragment, map_num)
    working = Chem.RemoveHs(working)
    Chem.SanitizeMol(working)
    return working, attachments


def generate_random_molecules(
    basis_smiles: str,
    r_group_smiles: Mapping[int, Sequence[str]],
    n_molecules: int,
    seed: Optional[int] = None,
) -> List[GeneratedMolecule]:
    """Generate ``n_molecules`` random combinations of the basis and R-groups."""

    if n_molecules <= 0:
        raise ValueError("n_molecules must be a positive integer.")

    rng = random.Random(seed)
    basis = _parse_smiles(basis_smiles, "basis structure")
    placeholder_map_nums = _find_placeholder_map_nums(basis)
    library = _prepare_r_groups(r_group_smiles)

    missing = sorted(set(placeholder_map_nums) - set(library))
    if missing:
        raise ValueError(
            "Missing R-group definitions for the following attachment points: "
            + ", ".join(f"R{m}" for m in missing)
        )

    molecules: List[GeneratedMolecule] = []
    for _ in range(n_molecules):
        mol, attachments = _assemble_molecule(basis, placeholder_map_nums, library, rng)
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        selfies_string = sf.encoder(smiles)
        molecules.append(GeneratedMolecule(smiles=smiles, selfies=selfies_string, attachments=attachments, mol=mol))
    return molecules


def generate_random_selfies(
    basis_smiles: str,
    r_group_smiles: Mapping[int, Sequence[str]],
    n_molecules: int,
    seed: Optional[int] = None,
) -> List[str]:
    """Convenience wrapper returning just the SELFIES strings."""

    return [entry.selfies for entry in generate_random_molecules(basis_smiles, r_group_smiles, n_molecules, seed=seed)]


def main() -> None:
    """Demonstrate the generator on a simple para-disubstituted benzene scaffold."""

    basis_smiles = "c1cc([*:1])ccc1[*:2]"
    r_groups = {
        1: ["[*:1]C", "[*:1]O", "[*:1]N"],
        2: ["[*:2]O", "[*:2]C(=O)O", "[*:2]N"],
    }

    molecules = generate_random_molecules(basis_smiles, r_groups, n_molecules=5, seed=0)
    for idx, entry in enumerate(molecules, start=1):
        attachment_summary = ", ".join(f"R{key}={value}" for key, value in sorted(entry.attachments.items()))
        print(f"{idx}: {entry.smiles} -> {entry.selfies} ({attachment_summary})")


if __name__ == "__main__":
    main()
