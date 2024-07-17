"""
File that contains fuctionaility for reversing a reaction using a reaction template using rdkit
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import ChemicalReaction
from rdkit.Chem.rdchem import Atom, Bond, Mol, EditableMol
import regex as re
from typing import List, Set, Tuple, Dict
from utils import highlight_reaction_core

from rxnmapper import RXNMapper


def get_reactants(substrate: str, reactant_smiles: List[str], product_smiles: List[str]) -> Set[Mol]:
    """
    Uses rdkit to get the possible reactants for a given reaction using the template reactants and products strings
    """

    # process the string to remove the R group
    replaced_r_group_reactant = [re.sub(r'R\d+', 'Zn', smile) for smile in reactant_smiles]
    replaced_r_group_product = [re.sub(r'R\d+', 'Zn', smile) for smile in product_smiles]
    if any("R" in smiles for smiles in replaced_r_group_reactant) or any ("R" in smiles for smiles in replaced_r_group_product):
        replaced_r_group_reactant = [re.sub(r'R+', 'Zn', smile) for smile in replaced_r_group_reactant]
        replaced_r_group_product = [re.sub(r'R+', 'Zn', smile) for smile in replaced_r_group_product]
    
    removed_r_group_reactant_mols = [get_mol_without_r_group(Chem.MolFromSmiles(mol)) for mol in replaced_r_group_reactant]
    removed_r_group_product_mols = [get_mol_without_r_group(Chem.MolFromSmiles(mol)) for mol in replaced_r_group_product]

    # reactants_smarts = '.'.join([Chem.MolToSmiles(mol) for mol in removed_r_group_reactant_mols])
    # products_smarts = '.'.join([Chem.MolToSmiles(mol) for mol in removed_r_group_product_mols])
    reactants_smarts = '.'.join([mol for mol in replaced_r_group_reactant])
    products_smarts = '.'.join([mol for mol in replaced_r_group_product])
    reaction_smarts_without_atom_map = reactants_smarts + ">>" + products_smarts
    
    # add atom map nums to reaction smarts
    rxn_mapper = RXNMapper()
    reaction_smarts = rxn_mapper.get_attention_guided_atom_maps([reaction_smarts_without_atom_map])[0]["mapped_rxn"]
    # get the fowards and backwards reactions
    # also need to remove zinc as it was the placeholder for the R groups
    foward_reaction_template = AllChem.ReactionFromSmarts(reaction_smarts)
    # foward_reaction_template.RunReactants([Chem.MolFromSmiles("O=C(C1=CC=CC=C1/C=C/C2=CC=CC=C2)NC")])[0]
    backwards_reaction_template = ChemicalReaction()
    for i in range(0, foward_reaction_template.GetNumReactantTemplates()):
        reactant = foward_reaction_template.GetReactantTemplate(i)
        backwards_reaction_template.AddProductTemplate(get_mol_without_r_group(reactant))
        # highlight_reaction_core(reactant, set(), set(), "reactant_template" + str(i) + ".png")
    for i in range(0, foward_reaction_template.GetNumProductTemplates()):
        product = foward_reaction_template.GetProductTemplate(i)
        backwards_reaction_template.AddReactantTemplate(get_mol_without_r_group(product))
        # highlight_reaction_core(product, set(), set(), "product_template" + str(i) + ".png")
    backwards_reaction_template.Initialize()
    
    # run the products as the reactants for the backwards reaction
    products_without_r_groups = [Chem.MolFromSmiles(smile) for smile in product_smiles if "R" not in smile and\
                                re.search(r'R\d+', smile) is None] + [Chem.MolFromSmiles(substrate)]
    counter = 1
    for product in products_without_r_groups:
        # highlight_reaction_core(product, set(), set(), "substrate_product" + str(counter) + ".png")
        counter += 1
    
    # new_products = foward_reaction_template.RunReactants(Chem.MolFromSmiles("O=C(C1=CC=CC=C1/C=C/C2=CC=CC=C2)NC"))[0]
    reactants_for_substrate = backwards_reaction_template.RunReactants(products_without_r_groups)[0]
    
    # FOR DEBUGGING
    # counter = 1
    # for reactant in reactants_for_substrate:
    #     Chem.SanitizeMol(reactant)
    #     highlight_reaction_core(reactant, set(), set(), "new_reactant" + str(counter) + ".png")
    #     counter += 1
    
    return [Chem.MolToSmiles(Chem.SanitizeMol(mol)) for mol in reactants_for_substrate]


def get_mol_without_r_group(mol: Mol) -> Mol:
    """
    Creates a new Mol object from mol without the R groups, or in this case without zinc
    """
    # assemble new_mol from mol and do not add any atoms that are zinc (R group place holder)
    new_mol = EditableMol(Mol())
    mol_to_new_mol = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 30:
            new_atom = Atom(atom.GetAtomicNum())
            new_atom.SetAtomMapNum(atom.GetAtomMapNum())
            mol_to_new_mol[atom.GetIdx()] = new_mol.AddAtom(new_atom)
    for bond in mol.GetBonds():
        if bond.GetBeginAtom().GetAtomicNum() != 30 and bond.GetEndAtom().GetAtomicNum() != 30:
            atom_num1, atom_num2 = bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx()
            idx1, idx2 = mol_to_new_mol[atom_num1], mol_to_new_mol[atom_num2]
            new_mol.AddBond(idx1, idx2, bond.GetBondType())
    
    return new_mol.GetMol()
    

# if __name__ == "__main__":
    
    # wiley2_scheme_3
    # print(get_reactants("CN1C(C2=CC=CC=C2C1CC3=CC=CC=C3)=O", ["O=C(NC)C1=CC=CC=C1/C=C([R1])/[R2]"], ["O=C(N1C)C2=CC=CC=C2C1C([R1])[R2]"]))
    # print(get_reactants("CN1C(C2=CC=CC=C2C1CC3=CC=CC=C3)=O", ["O=C(C1=CC=CC=C1/C=C/C2=CC=CC=C2)NC"], ["CN1C(C2=CC=CC=C2C1CC3=CC=CC=C3)=O"]))

    # wiley38_table2
    # print(get_reactants("O=C(C1CCOCC1)C2=CC=C(C)C=C2", ["IC1CCOCC1", "OC(C1=CC([R])=CC=C1)=O"], ["O=C(C1CCOCC1)C2=CC=CC([R])=C2"]))