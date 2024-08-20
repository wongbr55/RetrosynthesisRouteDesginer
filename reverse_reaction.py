"""
File that contains fuctionaility for reversing a reaction using a reaction template using rdkit
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdChemReactions, rdmolops, rdqueries
from rdkit.Chem.AllChem import ChemicalReaction
from rdkit.Chem.rdchem import Atom, Bond, Mol, EditableMol, BondType
import regex as re
from typing import List, Set, Tuple, Dict
from utils import highlight_reaction_core

from rxnmapper import RXNMapper

R_GROUP_PLACEHOLDER = 30
R_GROUP_PLACEHOLDER_SYMBOL = "Zn"
AR_GROUP_PLACEHOLDER = 87
AR_GROUP_PLACEHOLDER_SYMBOL = "Fr"
X_GROUP_PLACEHOLDER_SYMBOL = "Cu"

def get_reactants(substrate: str, reactant_smiles: List[str], product_smiles: List[str], debug: bool=False) -> Set[Mol]:
    """
    Uses rdkit to get the possible reactants for a given reaction using the template reactants and products strings
    """
    reaction_smarts = _get_smarts(reactant_smiles, product_smiles)
    # get the foward reaction
    foward_reaction_template = AllChem.ReactionFromSmarts(reaction_smarts)   
    backwards_reaction_template = _create_backwards_reaction(foward_reaction_template, debug)
    # run the products as the reactants for the backwards reaction
    products_without_r_groups = [Chem.AddHs(Chem.MolFromSmiles(Chem.CanonSmiles(substrate)))]
    
    if debug:
        counter = 1
        for product in products_without_r_groups:
            highlight_reaction_core(product, set(), set(), "substrate_product" + str(counter) + ".png")
            counter += 1
    # get the reactants for the current substrate    
    reactants_for_substrate = backwards_reaction_template.RunReactants(products_without_r_groups)[0]
    if debug:
        counter = 1
        for mol in reactants_for_substrate:
            highlight_reaction_core(mol, set(), set(), "new_reactant" + str(counter) + ".png")
            counter += 1

    return [Chem.MolToSmiles(mol) for mol in reactants_for_substrate]


def _get_smarts(reactant_smiles: List[str], product_smiles: List[str]) -> str:
    """
    Gets the SMARTS reaction string for reactant_smiles and product_smiles
    """
    # replace Ar groups first, weo do now so that we can get appropriate atom map
    
    # replace R groups
    replaced_r_group_reactant = [re.sub(r'R\d+', R_GROUP_PLACEHOLDER_SYMBOL, smile) for smile in reactant_smiles]
    replaced_r_group_product = [re.sub(r'R\d+', R_GROUP_PLACEHOLDER_SYMBOL, smile) for smile in product_smiles]
    if any("R" in smiles for smiles in replaced_r_group_reactant) or any ("R" in smiles for smiles in replaced_r_group_product):
        replaced_r_group_reactant = [re.sub(r'R+', R_GROUP_PLACEHOLDER_SYMBOL, smile) for smile in replaced_r_group_reactant]
        replaced_r_group_product = [re.sub(r'R+', R_GROUP_PLACEHOLDER_SYMBOL, smile) for smile in replaced_r_group_product]
    # replace X groups
    replaced_x_group_reactant = [re.sub(r'X+', X_GROUP_PLACEHOLDER_SYMBOL, smile) for smile in replaced_r_group_reactant]
    replaced_x_group_product = [re.sub(r'X+', X_GROUP_PLACEHOLDER_SYMBOL, smile) for smile in replaced_r_group_product]
    
    reactants_smarts = '.'.join([mol for mol in replaced_x_group_reactant])
    products_smarts = '.'.join([mol for mol in replaced_x_group_product])
    reaction_smarts_without_atom_map = reactants_smarts + ">>" + products_smarts
    
    # add atom map nums to reaction smarts
    rxn_mapper = RXNMapper()
    results = rxn_mapper.get_attention_guided_atom_maps([reaction_smarts_without_atom_map])
    reaction_smarts = ""
    highest_confidence = 0
    for result in results:
        if result["confidence"] > highest_confidence:
            reaction_smarts = result["mapped_rxn"]
            highest_confidence = result["confidence"]
    
    return reaction_smarts.replace(X_GROUP_PLACEHOLDER_SYMBOL, "*")
    

def _create_backwards_reaction(foward_reaction_template: ChemicalReaction, debug: bool) -> ChemicalReaction:
    """
    Takes the foward reaction template and creates the backwards one while preprocessing the reactants and products
    """
    # change query parameters to allow for better structure match
    qp = Chem.AdjustQueryParameters()
    qp.makeDummiesQueries = True
    qp.adjustDegree = False
    
    backwards_reaction_template = ChemicalReaction()
    for i in range(0, foward_reaction_template.GetNumReactantTemplates()):
        reactant = foward_reaction_template.GetReactantTemplate(i)
        new_product = _get_mol_without_r_group(Chem.AddHs(reactant, explicitOnly=True))
        modified_product = Chem.AdjustQueryProperties(new_product, qp)
        backwards_reaction_template.AddProductTemplate(modified_product)
        if debug:
            highlight_reaction_core(modified_product, set(), set(), "reactant_template" + str(i) + ".png")
    for i in range(0, foward_reaction_template.GetNumProductTemplates()):
        product = foward_reaction_template.GetProductTemplate(i)
        new_reactant = _get_mol_without_r_group(Chem.AddHs(product, explicitOnly=True))
        modified_reactant = Chem.AdjustQueryProperties(new_reactant, qp)
        backwards_reaction_template.AddReactantTemplate(modified_reactant)
        if debug:
            highlight_reaction_core(modified_reactant, set(), set(), "product_template" + str(i) + ".png")
     
    return backwards_reaction_template


def _get_mol_without_r_group(mol: Mol) -> Mol:
    """
    Creates a new Mol object from mol without the R groups, or in this case without zinc
    """
    # assemble new_mol from mol and do not add any atoms that are zinc (R group place holder)
    new_mol = EditableMol(Mol())
    mol_to_new_mol = {}
    for atom in mol.GetAtoms():
        # R group place holder
        if atom.GetAtomicNum() == R_GROUP_PLACEHOLDER:
            new_wildcard = Atom(0)
            new_wildcard.SetAtomMapNum(atom.GetAtomMapNum())
            new_wildcard.SetNoImplicit(True)
            mol_to_new_mol[atom.GetIdx()] = new_mol.AddAtom(new_wildcard)
        # otherwise normal atom
        else:    
            new_atom = Atom(atom.GetAtomicNum())
            new_atom.SetAtomMapNum(atom.GetAtomMapNum())
            new_atom.SetFormalCharge(atom.GetFormalCharge())
            new_atom.SetHybridization(atom.GetHybridization())
            new_atom.SetNoImplicit(True)
            mol_to_new_mol[atom.GetIdx()] = new_mol.AddAtom(new_atom)
    # add bonds
    for bond in mol.GetBonds():
        atom_num1, atom_num2 = bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx()
        idx1, idx2 = mol_to_new_mol[atom_num1], mol_to_new_mol[atom_num2]
        new_mol.AddBond(idx1, idx2, bond.GetBondType())
    
    new_mol = new_mol.GetMol()
    return new_mol
