"""
Classes that are involved with developing rules for fragments of substrate
"""

from partial_molecule import Fragment
from rdkit.Chem.rdchem import Atom, Bond, Mol, EditableMol


# modify property, add back bond, remove bond, 
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
        endpoint1, endpoint2 = self.fragment_bond.GetBeginAtom().GetAtomMapNum(),self.fragment_bond.GetEndAtom().GetAtomMapNum()
        bond = fragment.get_bond_from_atom_map_num(endpoint1, endpoint2)
        
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
