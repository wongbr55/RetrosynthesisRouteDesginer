"""
Classes that are involved with developing rules for fragments of substrate
"""

from partial_molecule import Fragment
from rdkit.Chem.rdchem import Atom, Bond, Mol, EditableMol


class Rule:
    """
    Abstract class that represents a "rule" for a Fragment that transforms a Fragment back to the proper reactant
    """
    
    def __init__(self) -> None:
        return None
    
    def apply_rule(self, fragment: Fragment) -> None:
        """
        Generic apply rule function that applies rule to self.fragment
        """
        raise NotImplementedError
    
    def get_bond_from_fragment(self, bond: Bond, fragment: Fragment) -> Bond:
        """
        Gets the corrosponding bond from fragment with different mapped atom numbers
        """
        raise NotImplementedError


class ModifyBond(Rule):
    """
    Rule for fragment that will modify a fixed bond's properties to match that of another
    
    fragment_bond is the bond from the fragment that needs to be modified
    reactant_bond is the bond from the reactant template that needs to be observed
    """
    
    fragment_bond: Bond
    reactant_bond: Bond
    
    def __init__(self, fragment_bond: Bond, reactant_bond: Bond) -> None:
        super().__init__()
        self.fragment_bond = fragment_bond
        self.reactant_bond = reactant_bond
    
    def apply_rule(self, fragment: Fragment) -> None:
        """
        Applies the ModifyBond rule to fragment
        """
    
        # figure out what properties changed
        changed_props = self._find_changed_props()
        bond = self.get_bond_from_fragment(self.fragment_bond, fragment)
        
        # apply new properties to the matching bond in fragment
        for prop in changed_props:
            other_prop = changed_props[prop]
            if type(other_prop) is bool:
                bond.SetBoolProp(prop, other_prop)
            elif type(other_prop) is float:
                bond.SetDoubleProp(prop, other_prop)
            elif type(other_prop) is int:
                bond.SetIntProp(prop, other_prop)
            else:
                bond.SetProp(prop, other_prop)
    
    def _find_changed_props(self) -> dict:
        
        prop_dict = {}
        for prop in self.fragment_bond.GetPropNames():
            other_prop = self.reactant_bond.GetProp(prop)
            if self.fragment_bond.GetProp(prop) != other_prop:
                prop_dict[prop] = other_prop
        
        return prop_dict


class AddBond(Rule):
    """
    Rule for fragment that will add a bond to a given fragment
    
    bond_to_add is the bond to add
    """
    
    bond_to_add: Bond
    
    def __init__(self, bond_to_add: Bond) -> None:
        super().__init__()
        self.bond_to_add = bond_to_add
    
    def apply_rule(self, fragment: Fragment) -> None:
        
        # TODO replace this with code that will actually get the proper atoms that relate to this bond
        begin_atom = ...
        end_atom = ...
        
        if not fragment.check_atom(begin_atom):
            fragment.add_atom(begin_atom)
        if not fragment.check_atom(end_atom):
            fragment.add_atom(end_atom)
        
        fragment.add_bond(self.bond_to_add)
        

class RemoveBond(Rule):
    """
    Rule for fragment that will remove a bond to a given fragment
    
    bond_to_remove is the bond to remove
    """
    
    bond_to_remove: Bond
    
    def __init__(self, bond_to_remove: Bond) -> None:
        super().__init__()
        self.bond_to_remove = bond_to_remove
        
    def apply_rule(self, fragment: Fragment) -> None:
        bond_in_fragment = self.get_bond_from_fragment(self.bond_to_remove, fragment)
        fragment.remove_bond(bond_in_fragment)