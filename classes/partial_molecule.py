"""
File containing class information for all defined classes (currently only ReactionCore)
"""
from rdkit import Chem
from rdkit.Chem.rdchem import Atom, Bond, Mol, EditableMol
from typing import Set, Tuple, List, Dict
from utils import find_bond_set

class ReactionCore():
    
    """
    A class that represents the reaction core of a reaction
    """
    atoms: Set[Atom]
    bonds: Set[Bond]
    atom_map_nums: Set[int]
    bond_map_nums: Set[Tuple[int]]
    

    def __init__(self) -> None:
        self.atoms = set()
        self.bonds = set()
        self.atom_map_nums = set()
        self.bond_map_nums = set()
        
    def __str__(self) -> str:
        return "Atoms: " + str(self.atom_map_nums) + '\n' + "Bonds: " + str(self.bond_map_nums)
        
    def add_atom(self, atom: Atom) -> None:
        """
        Adds atom to core
        Assumes that atom is not already in the current core
        """
        self.atoms.add(atom)
        self.atom_map_nums.add(atom.GetAtomMapNum())
    
    def add_bond(self, bond: Bond) -> None:
        """
        Adds bond to core
        Assumes that bond is not already in current core
        """
        self.bonds.add(bond)
        self.bond_map_nums.add((bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()))
    
    def check_atom(self, atom: Atom) -> bool:
        """
        Checks whether an atom is in the current core with an Atom object
        """
        return self.check_atom_map_num(atom.GetAtomMapNum())

    def check_atom_map_num(self, atom_map_num: int):
        """
        Checks whether an atom is in the current core with athe atom map number
        """
        return atom_map_num in self.atom_map_nums

    def check_bond(self, bond: Bond) -> bool:
        """
        Checks whether a bond is in the current core with a Bond object
        """
        index1, index2 = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
        return self.check_bond_map_num(index1, index2)
    
    def check_bond_map_num(self, index1: int, index2: int):
        """
        Checks whether a bond is in the current core with the atom map numbers
        """
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
        for bond in bond_set:
            index1, index2 = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
            if not (index1, index2) in self.bond_map_nums and not (index2, index1) in self.bond_map_nums:
                self.bonds.add(bond)
                self.bond_map_nums.add((index1, index2))
    
    def get_bond_from_atom_map_num(self, endpoint1: int, endpoint2: int) -> Bond:
        """
        Returns Bond object with endpoints endpoint1 and endpoint2
        """
        for bond in self.bonds:
            possible1, possible2 = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
            if (possible1 == endpoint1 and possible2 == endpoint2) or (possible1 == endpoint2 and possible2 == endpoint1):
                return bond
        return None
         
    def get_mol(self) -> Mol:
        """
        Returns a new Mol object for current Reaction core with proper bonds and such
        Mol object may contain multiple molecules
        """
        new_mol = EditableMol(Chem.Mol())
        map_num_to_idx = {}
        for atom in self.atoms:
            # transfers[atom.GetAtomMapNum()] = counter
            new_atom = Atom(atom.GetAtomicNum())
            new_atom.SetAtomMapNum(atom.GetAtomMapNum())
            # new_atom.SetAtomMapNum(counter)
            map_num_to_idx[new_atom.GetAtomMapNum()] = new_mol.AddAtom(new_atom)
            # counter += 1
        for bond_num in self.bond_map_nums:
            atom_num1, atom_num2 = bond_num[0], bond_num[1]
            bond = self.get_bond_from_atom_map_num(atom_num1, atom_num2)
            idx1, idx2 = map_num_to_idx[atom_num1], map_num_to_idx[atom_num2]
            # new_num1, new_num2 = transfers[atom_num1], transfers[atom_num2]
            # new_mol.AddBond(new_num1, new_num2, bond.GetBondType())
            new_mol.AddBond(idx1, idx2, bond.GetBondType())
        
        return new_mol.GetMol()
    
    def get_smarts(self) -> str:
        """
        Returns the SMARTS for the current reaction core
        """
        return Chem.MolToSmarts(self.get_mol())


class Fragment(ReactionCore):
    """
    Class that represents fragments of a product template being broken apart
    """
    
    def __init__(self, atoms: Set[Atom], bonds: Set[Bond]) -> None:
        super().__init__()
        for atom in atoms:
            self.atoms.add(atom)
            self.atom_map_nums.add(atom.GetAtomMapNum())
        for bond in bonds:
            self.bonds.add(bond)
            self.bond_map_nums.add((bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()))
    
    def get_atom_from_atom_map_num(self, atom_map_num: int) -> Atom:
        """
        Returns Atom object with atom_map_num
        """
        for atom in self.atoms:
            if atom.GetAtomMapNum() == atom_map_num:
                return atom
        return None
    
    def remove_atom(self, atom: Atom) -> None:
        """
        Removes atom from self.atoms and self.atom_map_nums
        """
        self.atoms.remove(atom)
        self.atom_map_nums.remove(atom.GetAtomMapNum())
    
    def remove_bond(self, bond: Bond) -> None:
        """
        Removes bond from self.bonds and self.bond_map_nums
        """
        endpoint1, endpoint2 = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
        self.bonds.remove(bond)
        if (endpoint1, endpoint2) in self.bond_map_nums:
            self.bond_map_nums.remove((endpoint1, endpoint2))
        elif (endpoint2, endpoint1) in self.bond_map_nums:
            self.bond_map_nums.remove((endpoint2, endpoint1))
        # otherwise only one endpoint is in self.bonds and we ignore
        
        # check to see if endpoints are still connected to anything, otherwise start to remove atoms
        all_endpoints = {end[0] for end in self.bond_map_nums}.union({end[1] for end in self.bond_map_nums})
        if endpoint1 not in all_endpoints and endpoint1 in self.atom_map_nums:
            self.remove_atom(self.get_atom_from_atom_map_num(endpoint1))
        if endpoint2 not in all_endpoints and endpoint2 in self.atom_map_nums:
            self.remove_atom(self.get_atom_from_atom_map_num(endpoint2))
    
    def fragment_with_multiple_edges(self, edges: Set[Tuple[int]]):
        """
        Fragments a molecule by cutting off all of the edges in edges
        """
        
        # remove edges
        for edge in edges:
            self.remove_bond(self.get_bond_from_atom_map_num(edge[0], edge[1]))
            
        # get the remaining fragments
        fragments = set()
        atoms_seen_so_far_map_num = set()
        for atom in self.atoms:
            if atom.GetAtomMapNum() not in atoms_seen_so_far_map_num:
                atoms_seen_so_far_map_num = self._get_fragment_for_fixed_atom(atom, atoms_seen_so_far_map_num, fragments)
        
        return fragments
    
    def _get_fragment_for_fixed_atom(self, atom: Atom, atoms_seen_so_far_map_num: Set[int], fragments: Set) -> Set[int]:
        atoms_seen_so_far_map_num.add(atom.GetAtomMapNum())
        new_fragment_atoms = set()
        new_fragment_atoms_map_num = set()
        new_fragment_bonds = set()
        
        self._add_atoms_from_curr_atom(atom, new_fragment_atoms, new_fragment_bonds, set(), new_fragment_atoms_map_num)
        new_fragment = Fragment(new_fragment_atoms, new_fragment_bonds)
        fragments.add(new_fragment)
        return atoms_seen_so_far_map_num.union(new_fragment_atoms_map_num)
    
    # def fragment(self, edge: Tuple[int]):
    #     """
    #     may not be needed
    #     From the given edge, create a new fragment object such that the current fragment object is broken
    #     At "edge" and the new and current fragment are two "new" fragments
        
    #     Returns a new Fragment object
    #     """
    #     if not self.check_bond_map_num(edge[0], edge[1]):
    #         return None

    #     endpoint1, endpoint2 = edge[0], edge[1]
    #     atom1 = self.get_atom_from_atom_map_num(endpoint1)
    #     fragment1_atom_set, fragment1_bond_set, fragment1_atom_map_num, atoms_to_exclude = set(), set(), set(), set()
    #     atoms_to_exclude.add(endpoint2)
        
    #     self._add_atoms_from_curr_atom(atom1, fragment1_atom_set, fragment1_bond_set, atoms_to_exclude, fragment1_atom_map_num)
    #     fragment2_atom_map_num = self.atom_map_nums.difference(fragment1_atom_map_num).union(atoms_to_exclude)
    #     # if there are no new fragments formed, then we just remove edge from the self and return None
    #     if len(fragment2_atom_map_num) == 1:
    #         self.bonds.remove(self.get_bond_from_atom_map_num(endpoint1, endpoint2))
    #         if (endpoint1, endpoint2) in self.bond_map_nums:
    #             self.bond_map_nums.remove((endpoint1, endpoint2))
    #         else:
    #             self.bond_map_nums.remove((endpoint2, endpoint1))
    #         return None

    #     # otherwise there are 2 distinct and valid fragments, and we seperate as we need
    #     fragment2_bond_set = self._get_bonds_atom_set(fragment2_atom_map_num)
    #     fragment2_atom_set = {self.get_atom_from_atom_map_num(atom) for atom in fragment2_atom_map_num}
        
    #     self.atoms = fragment1_atom_set
    #     self.bonds = fragment1_bond_set
    #     self.atom_map_nums = fragment1_atom_map_num
    #     self.bond_map_nums = {(bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()) for bond in fragment1_bond_set}
    #     return Fragment(fragment2_atom_set, fragment2_bond_set)
    
    def transform(self, reactant_core: ReactionCore, conversion: Dict[int, int], edges_to_not_add: Set[Tuple[int]]) -> None:
        """
        Transform the fragment into its proper reactant given reactant core and conversions
        edges_to_not_add are edges that acted as breaking edges that broke substrate into fragments
        Assume conversion of form conversion[core_map_num] = substrate_map_num
        Assume edges_to_not_add has atom map numbers of core
        """
        
        # relevant_reactant_core is part of reaction core that is relevant to current fragment
        relevant_reactant_core = ReactionCore()
        for atom in reactant_core.atoms:
            if conversion[atom.GetAtomMapNum()] in self.atom_map_nums:
                relevant_reactant_core.add_atom(atom)
        
        for bond in reactant_core.bonds:
            endpoint1, endpoint2 = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
            if relevant_reactant_core.check_atom_map_num(endpoint1) and relevant_reactant_core.check_atom_map_num(endpoint2):
                relevant_reactant_core.add_bond(bond)
        
        # all changes that need to be made to self are in relevant_reactant_core
        # we use this to see if we need to add, remove, or modify a bond
        atom_bond_index_to_add = max({conversion[key] : key for key in conversion})
        for bond in relevant_reactant_core.bonds:
            # if there is a bond in the reactants that is not in the fragment, we add it
            endpoint1, endpoint2 = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
            if not endpoint1 in conversion or not endpoint2 in conversion and \
                (endpoint1, endpoint2) not in edges_to_not_add and (endpoint2, endpoint1) not in edges_to_not_add:
            # if (endpoint1, endpoint2) not in self.bond_map_nums and (endpoint2, endpoint1) not in self.bond_map_nums:
                atom_bond_index_to_add = self._add_bond_from_core(bond, conversion, atom_bond_index_to_add)
            # otherwise if it is in the fragment but is not matching, we try to modify it
            elif (endpoint1, endpoint2) not in edges_to_not_add and (endpoint2, endpoint1) not in edges_to_not_add:
                self_endpoint1, self_endpoint2 = conversion[endpoint1], conversion[endpoint2]
                curr_bond = self.get_bond_from_atom_map_num(self_endpoint1, self_endpoint2)
                if curr_bond is None:
                    atom_bond_index_to_add = self._add_bond_from_core(bond, conversion, atom_bond_index_to_add)
                elif curr_bond.GetBondType() != bond.GetBondType():
                    curr_bond.SetBondType(bond.GetBondType())
        
        # now we check for bonds that are in fragment that are not in reactant, if so remove bond
        bonds_to_remove = set()
        sub_to_core = {conversion[key] : key for key in conversion}
        for bond in self.bonds:
            endpoint1, endpoint2 = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
            if endpoint1 in sub_to_core and endpoint2 in sub_to_core and \
                not relevant_reactant_core.check_bond_map_num(sub_to_core[endpoint1], sub_to_core[endpoint2]):
                bonds_to_remove.add(bond)
        for bond in bonds_to_remove:
            self.remove_bond(bond)

    def _add_bond_from_core(self, core_bond: Bond, core_to_sub: Dict[int, int], next_index_to_add) -> int:
        """
        Adds a bond from core to self, uses core_to_sub map to help with map number assignment
        Returns next_index_to_add + however many atoms added
        """
        atom1, atom2 = core_bond.GetBeginAtom(), core_bond.GetEndAtom()
        
        # check if either atom is not in self:
        num_added = 0
        if atom1.GetAtomMapNum() not in core_to_sub:
            core_to_sub[atom1.GetAtomMapNum()] = next_index_to_add + num_added + 1
            atom1.SetAtomMapNum(next_index_to_add + num_added + 1)
            self.add_atom(atom1)
            num_added += 1
        else:
            atom1.SetAtomMapNum(core_to_sub[atom1.GetAtomMapNum()])
        if atom2.GetAtomMapNum() not in core_to_sub:
            core_to_sub[atom1.GetAtomMapNum()] = next_index_to_add + num_added + 1
            atom2.SetAtomMapNum(next_index_to_add + num_added + 1)
            self.add_atom(atom2)
            num_added += 1
        else:
            atom2.SetAtomMapNum(core_to_sub[atom2.GetAtomMapNum()])
        self.add_bond(core_bond)
        return next_index_to_add + num_added
    
    def _find_changed_props(self, core_bond: Bond, fragment_bond: Bond) -> dict:
        prop_dict = {}
        for prop in fragment_bond.GetPropNames():
            other_prop = core_bond.GetProp(prop)
            if fragment_bond.GetProp(prop) != other_prop:
                prop_dict[prop] = other_prop
        
        return prop_dict
    
    def _find_reactant_that_matches(self, reactant_template: Set[Mol]) -> Mol:
        """
        finds a reactant from the set that matches with the current fragment
        """
        
        for reactant in reactant_template:
            reactant_atom_map_nums = {atom.GetAtomMapNum() for atom in reactant.GetAtoms()}
            if all(num in reactant_atom_map_nums for num in self.atom_map_nums):
                return reactant
    
    def _add_atoms_from_curr_atom(self, curr_atom: Atom, atoms_so_far: Set[Atom], bonds_so_far: Set[Bond], \
                                    atoms_to_exclude: Set[int], atoms_so_far_map_num: Set[int]):
        """
        Extends atoms_so_far and bonds_so_far to include everything at curr_atom and beyond excludinsg anything from atoms_to_exclude
        """
        for bond in curr_atom.GetBonds():
            if self.check_bond(bond):
                other_atom = bond.GetOtherAtom(curr_atom)
                other_atom_map_num = other_atom.GetAtomMapNum()
                if other_atom_map_num not in atoms_so_far_map_num and other_atom_map_num not in atoms_to_exclude \
                    and self.check_bond_map_num(other_atom_map_num, curr_atom.GetAtomMapNum()):
                    atoms_so_far.add(other_atom)
                    atoms_so_far_map_num.add(other_atom_map_num)
                    bonds_so_far.add(bond)
                    self._add_atoms_from_curr_atom(other_atom, atoms_so_far, bonds_so_far, atoms_to_exclude, atoms_so_far_map_num)
                elif other_atom_map_num in atoms_so_far_map_num and other_atom_map_num not in atoms_to_exclude \
                    and self.check_bond_map_num(other_atom_map_num, curr_atom.GetAtomMapNum()) and \
                        find_bond_set(other_atom_map_num, curr_atom.GetAtomMapNum(), bonds_so_far) is None:
                             bonds_so_far.add(bond)
        
    # def _get_bonds_atom_set(self, atoms: Set[int]) -> Set[Bond]:
    #     """
    #     Returns a set of Bond objects that only contain the atoms in atom
    #     Note atoms is a set of atom map numbers
    #     """
        
    #     bonds_so_far = set()
    #     for bond in self.bonds:
    #         endpoint1, endpoint2 = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
    #         if endpoint1 in atoms and endpoint2 in atoms:
    #             bonds_so_far.add(bond)
        
    #     return bonds_so_far