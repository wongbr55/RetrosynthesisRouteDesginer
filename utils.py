"""
File that contains common utility functions
"""

from rdkit import Chem
from rdkit.Chem.rdchem import Mol, Atom, Bond, BondType
from rdkit.Chem import Draw
from rdkit.Chem.rdChemReactions import ChemicalReaction
from typing import List, Set, Tuple

def find_atom(atom_map_num: int, mols: List[Mol]) -> Atom:
    """
    Finds Atom object that has the same atom map number as atom_map_num
    """
    for mol in mols:
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() == atom_map_num:
                return atom
    return None


def find_bond(index1: int, index2: int, mols: List[Mol]) -> Bond:
    """
    Finds Bond object that has same bond map numbers as index1 and index2
    """
    
    for mol in mols:
        for bond in mol.GetBonds():
            possible1, possible2 = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
            if (index1 == possible1 and index2 == possible2) or (index1 == possible2 and index2 == possible1):
                return bond
    return None


def highlight_reaction_core(mol, changing_atoms, changing_bonds):
    """
    Draws highlted reaction core on reactant and product
    Only used for debugging
    :param mol:
    :param changing_atoms:
    :param changing_bonds:
    :return:
    """
    atom_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomMapNum() in changing_atoms]
    bond_indices = [bond.GetIdx() for bond in mol.GetBonds() if
                    (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) in changing_bonds or (
                        bond.GetEndAtomIdx(), bond.GetBeginAtomIdx()) in changing_bonds]
    return Draw.MolToFile(mol, "testing.png", highlightAtoms=atom_indices, highlightBonds=bond_indices)


def print_atom_and_bond(reaction_core: Tuple[Set[Atom], Set[Bond]]) -> None:
    """
    Prints the map numbers of Atom and Bond objects in a reaction core
    Used for debugging purposes
    :param reaction_core:
    :return:
    """
    print("## ATOMS ##")
    for atom in reaction_core[0]:
        print("Atom :" + str(atom.GetAtomMapNum()))
    print("## BONDS ##")
    for bond in reaction_core[1]:
        print("From " + str(bond.GetBeginAtom().GetAtomMapNum()) + " to " + str(bond.GetEndAtom().GetAtomMapNum()))
