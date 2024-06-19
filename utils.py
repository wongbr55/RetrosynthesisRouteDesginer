"""
File containing class information for all defined classes (currently only ReactionCore)
"""

from rdkit.Chem.rdchem import Atom, Bond
from typing import Set, Tuple

class ReactionCore():
    
    """
    Class that represents the reaction core of a reaction
    Contains wrapper functinaility for most set operations
    """
    
    atoms: Set[Atoms]
    bonds: Set[Bonds]
    atom_map_nums: Set[int]
    bond_map_nums: Set[Tuple[int]]
    
    def __init__(self) -> None:
        self.atoms = set()
        self.bonds = set()
        self.atom_map_nums = set()
        self.bond_map_nums = set()
    
    def __str__(self) -> None:
        print("Atoms: " + str(self.atom_map_nums))
        print("Bonds: " + str(self.bond_map_nums))
    
    def add_atom_to_core(self, atom: Atom) -> None:
        """
        Adds atom to core
        Assumes that atom is not already in the current core
        """
        self.atoms.add(atom)
        self.atom_map_nums.add(atom.GetAtomMapNum())
    
    def add_bond_to_core(self, bond: Bond) -> None:
        """
        Adds bond to core
        Assumes that bond is not already in current core
        """
        self.bonds.add(Bond)
        self.bond_map_nums.add(bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum())
    
    def check_atom_in_core(self, atom: Atom) -> bool:
        """
        Checks whether an atom is in the current core
        """
        return atom.GetAtomMapNum() in self.atom_map_nums

    def check_bond_in_core(self, bond: Bond) -> bool:
        """
        Checks whether a bond is in the current core
        """
        index1, index2 = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
        return (index1, index2) in self.bond_map_nums or (index2, index1) in self.bond_map_nums

    def add_atoms(self, atom_set: Set[Atom]) -> None:
        """
        Adds atom_set to atoms and mutates atom_map_nums accoridingly
        """
        for atom in atom_set:
            if atom.GetAtomMapNum() not in self.atom_map_nums:
                self.atoms.add(atom)
                self.atom_map_nums.add(atom.GetAtomMapNum())
    
    def add_bonds(self, bond_set: Set[Bond]) -> None:
        """
        Adds bond_set to bonds and mutates bond_map_nums accoridingly
        """
        self.bonds = bond_set
        for bond in bond_set:
            index1, index2 = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
            if not (index1, index2) in self.bond_map_nums and not (index2, index1) in self.bond_map_nums:
                self.bonds.add(bond)
                self.bond_map_nums.add((index1, index2))
    