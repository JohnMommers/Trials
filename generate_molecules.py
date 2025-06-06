import itertools
from rdkit import Chem
from rdkit.Chem import Draw
import selfies as sf

# utility to attach a substituent at a given mapping number

def attach_substituent(base_smiles: str, sub_smiles: str, map_num: int) -> Chem.Mol:
    base = Chem.MolFromSmiles(base_smiles)
    sub = Chem.MolFromSmiles(sub_smiles)
    base_rw = Chem.RWMol(base)
    sub_rw = Chem.RWMol(sub)

    base_star_idx = None
    for atom in base_rw.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetAtomMapNum() == map_num:
            base_star_idx = atom.GetIdx()
            break
    if base_star_idx is None:
        raise ValueError(f"Placeholder {map_num} not found")

    sub_star_idx = None
    for atom in sub_rw.GetAtoms():
        if atom.GetAtomicNum() == 0:
            sub_star_idx = atom.GetIdx()
            break
    if sub_star_idx is None:
        raise ValueError("Substituent missing star atom")

    combo = Chem.CombineMols(base_rw, sub_rw)
    combo_rw = Chem.RWMol(combo)

    sub_star_idx2 = sub_star_idx + base_rw.GetNumAtoms()

    base_nbrs = [n.GetIdx() for n in combo_rw.GetAtomWithIdx(base_star_idx).GetNeighbors()]
    sub_nbrs = [n.GetIdx() for n in combo_rw.GetAtomWithIdx(sub_star_idx2).GetNeighbors()]
    if len(base_nbrs) != 1 or len(sub_nbrs) != 1:
        raise ValueError("Star atoms must have single neighbors")

    combo_rw.AddBond(base_nbrs[0], sub_nbrs[0], Chem.BondType.SINGLE)

    for idx in sorted([base_star_idx, sub_star_idx2], reverse=True):
        combo_rw.RemoveAtom(idx)
    return combo_rw.GetMol()


def build_molecule(skeleton_smiles: str, substituents: dict) -> Chem.Mol:
    current_smiles = skeleton_smiles
    for map_num in sorted(substituents.keys()):
        mol = attach_substituent(current_smiles, substituents[map_num], map_num)
        current_smiles = Chem.MolToSmiles(mol)
    return Chem.MolFromSmiles(current_smiles)


def enumerate_molecules(skeleton_smiles: str, rgroup_smiles_lists: list) -> dict:
    """Return dict of canonical SMILES to (Mol, SELFIES) for all combinations."""

    results = {}
    for combo in itertools.product(*rgroup_smiles_lists):
        substituents = {i + 1: combo[i] for i in range(len(combo))}
        mol = build_molecule(skeleton_smiles, substituents)
        smi = Chem.MolToSmiles(mol, canonical=True)
        selfie = sf.encoder(smi)
        results[smi] = (mol, selfie)
    return results


def plot_molecules(mol_dict: dict, filename: str = "molecules.png") -> None:
    mols = [pair[0] for pair in mol_dict.values()]
    legends = list(mol_dict.keys())
    img = Draw.MolsToGridImage(mols, legends=legends, molsPerRow=4, subImgSize=(200,200))
    img.save(filename)


if __name__ == "__main__":
    # Example main structure with three attachment points labeled 1-3
    skeleton_smiles = "[*:1]c1cc([*:2])ccc1[*:3]"

    # R groups are represented with SELFIES strings
    r1_selfies = ["[F]", "[Cl]"]
    r2_selfies = ["[O]", "[N]"]
    r3_selfies = ["[C]", "[S]"]

    # Convert to SMILES and prepend attachment point
    r1 = [f"[*:1]{sf.decoder(s)}" for s in r1_selfies]
    r2 = [f"[*:1]{sf.decoder(s)}" for s in r2_selfies]
    r3 = [f"[*:1]{sf.decoder(s)}" for s in r3_selfies]

    molecules = enumerate_molecules(skeleton_smiles, [r1, r2, r3])
    plot_molecules(molecules)
    print(f"Generated {len(molecules)} unique molecules. Image saved to molecules.png")
    for smi, (_, selfie) in molecules.items():
        print(smi, selfie)
