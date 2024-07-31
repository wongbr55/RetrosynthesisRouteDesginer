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
    # products_without_r_groups = [Chem.MolFromSmiles(Chem.CanonSmiles(smile)) for smile in product_smiles if "R" not in smile and "X" not in smile] \
    #                             + [Chem.AddHs(Chem.MolFromSmiles(Chem.CanonSmiles(substrate)))]
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
    
    # replaced_ar_group_reactant, replaced_ar_group_product = _replace_ar_groups(replaced_x_group_reactant, replaced_x_group_product)
    
    # replace Ar groups:
    # replaced_ar_group_reactant = [re.sub(r'Ar+', AR_GROUP_PLACEHOLDER_SYMBOL, smile) for smile in replaced_x_group_reactant]
    # replaced_ar_group_product = [re.sub(r'Ar+', AR_GROUP_PLACEHOLDER_SYMBOL, smile) for smile in replaced_x_group_product]
        
    reactants_smarts = '.'.join([mol for mol in replaced_x_group_reactant])
    products_smarts = '.'.join([mol for mol in replaced_x_group_product])
    reaction_smarts_without_atom_map = reactants_smarts + ">>" + products_smarts
    
    # add atom map nums to reaction smarts
    rxn_mapper = RXNMapper()
    reaction_smarts = rxn_mapper.get_attention_guided_atom_maps([reaction_smarts_without_atom_map])[0]["mapped_rxn"]
    return reaction_smarts.replace(X_GROUP_PLACEHOLDER_SYMBOL, "*")


# def _replace_ar_groups(reactant_smiles: List[str], product_smiles: List[str]) -> Tuple[List[str]]:
#     """
#     Replaces Ar groups with a 6 ring generic structure
#     Returns (reactant_smiles, product_smiles)
#     """
#     replaced_ar_group_reactant = [re.sub(r'Ar+', AR_GROUP_PLACEHOLDER_SYMBOL, smile) for smile in reactant_smiles]
#     replaced_ar_group_product = [re.sub(r'Ar+', AR_GROUP_PLACEHOLDER_SYMBOL, smile) for smile in product_smiles]
    
#     reactants = [Chem.MolFromSmiles(smiles) for smiles in replaced_ar_group_reactant]
#     products = [Chem.MolFromSmiles(smiles) for smiles in replaced_ar_group_product]
    
#     final_reactants = [_replace_ar_for_mol(mol) for mol in reactants]
#     final_products = [_replace_ar_for_mol(mol) for mol in products]
    
#     return [Chem.MolToSmiles(mol) for mol in final_reactants], [Chem.MolToSmiles(mol) for mol in final_products]
    
# def _replace_ar_for_mol(mol: Mol) -> Mol:
#     """
#     reaplaces Ar group for a single mol
#     """
#     new_mol = EditableMol(Mol())
#     for atom in mol.GetAtoms():
#         # if we have an Ar group, generate a generic 6 edge ring
#         if atom.GetAtomicNum() == AR_GROUP_PLACEHOLDER:
#             atom_ring = [Atom(R_GROUP_PLACEHOLDER) for __ in range(0, 6)]
#             atom_ring_idx = {i : new_mol.AddAtom(atom_ring[i]) for i in range(0, len(atom_ring))}
#             for i in range(0, len(atom_ring)):
#                 new_mol.AddBond(atom_ring_idx[i], atom_ring_idx[(i + 1) % len(atom_ring)], BondType.SINGLE)
        
#             neighbor = atom.GetNeighbors()[0]
#             new_mol.RemoveAtom(atom.GetIdx())
#             new_mol.AddBond(neighbor.GetIdx(), atom_ring_idx[0], BondType.SINGLE)
            
#     return new_mol.GetMol()
    

def _create_backwards_reaction(foward_reaction_template: ChemicalReaction, debug: bool) -> ChemicalReaction:
    """
    Takes the foward reaction template and creates the backwards one while preprocessing the reactants and products
    """
    # change query parameters to allow for better structure match
    qp = Chem.AdjustQueryParameters()
    qp.makeDummiesQueries = True
    qp.adjustDegree = False
    # qp.adjustSingleBondsToDegreeOneNeighbors = True
    # qp.adjustExplicitValence = True
    # qp.adjustImplicitValence = False
    
    # create the backwards reaction and label atoms with map numbers that are left unmapped
    next_index_to_add = 0
    for reactant in foward_reaction_template.GetReactants():
        next_index_to_add = max(next_index_to_add, max({atom.GetAtomMapNum() for atom in reactant.GetAtoms()}))
    for product in foward_reaction_template.GetProducts():
        next_index_to_add = max(next_index_to_add, max({atom.GetAtomMapNum() for atom in product.GetAtoms()}))
    
    backwards_reaction_template = ChemicalReaction()
    for i in range(0, foward_reaction_template.GetNumReactantTemplates()):
        reactant = foward_reaction_template.GetReactantTemplate(i)
        for atom in reactant.GetAtoms():
            if atom.GetAtomMapNum() == 0:
                atom.SetAtomMapNum(next_index_to_add + 1)
                next_index_to_add += 1
        new_product = _get_mol_without_r_group(Chem.AddHs(reactant, explicitOnly=True))
        # highlight_reaction_core(new_product, set(), set(), "og_reactant" + str(i) + ".png")
        modified_product = Chem.AdjustQueryProperties(new_product, qp)
        backwards_reaction_template.AddProductTemplate(modified_product)
        if debug:
            highlight_reaction_core(modified_product, set(), set(), "reactant_template" + str(i) + ".png")
    for i in range(0, foward_reaction_template.GetNumProductTemplates()):
        product = foward_reaction_template.GetProductTemplate(i)
        for atom in product.GetAtoms():
            if atom.GetAtomMapNum() == 0:
                atom.SetAtomMapNum(next_index_to_add + 1)
                next_index_to_add += 1
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
            mol_to_new_mol[atom.GetIdx()] = new_mol.AddAtom(new_wildcard)
        # otherwise normal atom
        else:    
            new_atom = Atom(atom.GetAtomicNum())
            new_atom.SetAtomMapNum(atom.GetAtomMapNum())
            new_atom.SetFormalCharge(atom.GetFormalCharge())
            new_atom.SetHybridization(atom.GetHybridization())
            mol_to_new_mol[atom.GetIdx()] = new_mol.AddAtom(new_atom)
    # add bonds
    # try using BondType.UNSPECIFIED???
    for bond in mol.GetBonds():
        atom_num1, atom_num2 = bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx()
        # if bond.GetBeginAtom().GetAtomicNum() != 30 and bond.GetEndAtom().GetAtomicNum() != 30:
        idx1, idx2 = mol_to_new_mol[atom_num1], mol_to_new_mol[atom_num2]
        # if bond.GetBeginAtom().GetAtomicNum() == 30 or bond.GetEndAtom().GetAtomicNum() == 30:
            # new_mol.AddBond(idx1, idx2, BondType.UNSPECIFIED)
        # else:
        new_mol.AddBond(idx1, idx2, bond.GetBondType())
    
    new_mol = new_mol.GetMol()
    Chem.SanitizeMol(new_mol, Chem.rdmolops.SANITIZE_ADJUSTHS)
    return new_mol

# if __name__ == "__main__":
    
    # wiley2_scheme_3
    # PASSES
    # print(get_reactants("CN1C(C2=CC=CC=C2C1CC3=CC=CC=C3)=O", ["O=C(NC)C1=CC=CC=C1/C=C([R1])/[R2]"], ["O=C(N1C)C2=CC=CC=C2C1C([R1])[R2]"], True))
    # print(get_reactants("CN1C(C2=CC=CC=C2C1CC3=CC=CC=C3)=O", ["O=C(C1=CC=CC=C1/C=C/C2=CC=CC=C2)NC"], ["CN1C(C2=CC=CC=C2C1CC3=CC=CC=C3)=O"]))

    # wiley38_table2
    # PASSES
    # print(get_reactants("O=C(C1CCOCC1)C2=CC=C(C)C=C2", ["IC1CCOCC1", "OC(C1=C([R])C([R])=C([R])C([R])=C1[R])=O"], ["O=C(C1=C([R])C([R])=C([R])C([R])=C1[R])C2CCOCC2"], True))
    
    # wiley6_scheme_2
    # PASSES
    # print(get_reactants("OCCC(C1=CC=C(Br)C=C1)=O", ["[R2]C1(C(C1)C2=CC=C(C=C2)[R])[R1]", "O"], ["O=C(CC([R2])(O)[R1])C1=CC=C(C=C1)[R]"]))
    
    # RSC18_scheme_2
    # DOES NOT PRODUCE VALID SMILES AND CANNOT SANITIZE, WORKS SOMETIMES
    # print(get_reactants("O=C(C1CC)N([H])C2=C1C=CC=C2", ["[R1]N(C=C1[R2])C2=C1C=C([R3])C=C2"], ["[R1]N(C(C1[R2])=O)C2=C1C=C([R3])C=C2"], True))
    
    # wiley16_scheme_2
    # PASSES
    # FAILS BECAUSE * IS NOT WORKING PROPERLY, HAVE TO CALL AdjustQueryProperties ON ALL MOLS
    # core = get_reaction_core(["[H]C1=NC2=C([X]1)C=CC=C2"], ["C12=C(C=CC=C2)N=C(N3CCOCC3)[X]1"])
    # product_core = extend_reaction_core(core[1], core[2], core[0])
    # reactant_core_smiles = Chem.MolToSmiles(core[0].get_mol(True)).split(",")
    # product_core_smiles = Chem.MolToSmiles(product_core.get_mol(True)).split(",")
    # reactant_core_smiles = [smile.replace("Cu", "X") for smile in reactant_core_smiles]
    # product_core_smiles = [smile.replace("Cu", "X") for smile in product_core_smiles]
    # print(get_reactants("CC1=CC2=C(OC(N3CCOCC3)=N2)C=C1", reactant_core_smiles, product_core_smiles))

    # print(get_reactants("CC1=CC2=C(OC(N3CCOCC3)=N2)C=C1", ["[H]C1=NC2=C([X]1)C=CC=C2", "O1CCNCC1"], ["C12=C(C=CC=C2)N=C(N3CCOCC3)[X]1"], True))
    
    
    # rxn = AllChem.ReactionFromSmarts('[c:1][#0].[#0][*:2]>>[c:1]-[*:2]')
    # counter = 1
    # for reactant in rxn.GetReactants():
    #     highlight_reaction_core(reactant, set(), set(), "reactant_template" + str(counter) + ".png")
    #     counter += 1
    # counter = 1
    # for product in rxn.GetProducts():
    #     highlight_reaction_core(product, set(), set(), "product_template" + str(counter) + ".png")
    #     counter += 1
        
    # reacts = (Chem.MolFromSmiles('*c1c(C)cccc1(O)'),Chem.MolFromSmiles('CN*'))
    # counter = 1
    # for reactant in reacts:
    #     highlight_reaction_core(reactant, set(), set(), "reactant" + str(counter) + ".png")
    #     counter += 1
    # counter = 1
    # products = rxn.RunReactants(reacts) # tuple
    # print(products)
    
    # RSC22_table_2
    # UNKOWN, PERHAPS IT IS WITH N BEING BONDED TO 2 HYDROGENS?
    # print(get_reactants("N12C=CC=CC1=NC(C3=CC=CC=C3)=C2N4C=NC=C4", ["[H]C1=C(N=C2N1C=CC=C2)C3=CC=C(C=C3)[R]", "[R2]N[R1]"], ["[R]C1=CC=C(C2=C(N3C=CC=CC3=N2)N([R2])[R1])C=C1"], True))
    
    # wiley30_table_3
    # DOES NOT ALWAYS FOLLOW REACTION TEMPLATE
    # O=S(C1=CC=C(C=C1)Cl)(C2=CC(N(C)C)=CC=C2)=O
    # print(get_reactants("O=S(C1=CC=C(C=C1)Cl)(C2=C(N(C)C)C=CC=C2)=O", ["[R2]N([R1])C1=CC=C([R3])C=C1", "O=[S-](C1=CC=C(Cl)C=C1)=O", "[Na+]"], ["O=S(C1=CC=C(Cl)C=C1)(C2=C(N([R1])[R2])C=CC([R3])=C2)=O"], True))

    # print(get_reactants("O=S(C1=CC=C(Cl)C=C1)(C2=CC(N(CC)C)=CC=C2)=O", ["[R2]N([R1])C1=CC=C([R3])C=C1", "O=[S-](C1=CC=C(Cl)C=C1)=O", "[Na+]"], ["O=S(C1=CC=C(Cl)C=C1)(C2=C(N([R1])[R2])C=CC([R3])=C2)=O"], True))


    # ja3c04864_0006
    # print(get_reactants("OC(C(C)(C)C)C/C=C/C1=CC=CC=C1", ["[R]C1=CC=C(C=C1)S(=O)(NN)=O", "O=C[R]"], ["OC([R])C/C=C/C1=CC=C([R])C=C1"], True))
    
    # ja3c06794_0004.png
    # reactants = get_reactants("ClCCCC(CC(C)(C)C)C1=CC(OCO2)=C2C=C1", ["[R]C([R])=C", "Br[C@@]([R1])([R3])[R2]", "[R4]CBr"], ["[R2]C([R1])([R3])CC([R])([R])C[R4]"], True)
    # print(reactants)
    # print([Chem.CanonSmiles(smiles) for smiles in reactants])

    # jo1c02275_0005.png
    # print(get_reactants("O=C(C1=NN(C2=CC=CC=C2)S/C1=N\C3=CC=CC=C3)C4=CC=CC=C4", ["[R1]C(/C(C(NC1=CC=CC=C1)=S)=N/N[R2])=O"], ["[R1]C(C(/C(S1)=N/C2=CC=CC=C2)=NN1[R2])=O"], True))

    # jo1c02557_0010_0011.png
    # print(get_reactants("O=C1CC(CS(=O)(C2=CC=C(C=CC=C3)C3=C2)=O)(C4=CC=CC=C4)C5=C1C=CC=C5", ["C=C(C1=CC=CC=C1C(O)=O)C2=CC=CC=C2", "[R2]S(=O)(NN)=O"], ["[R2]S(CC1(C2=CC=CC=C2)OC(C3=CC=CC=C31)=O)(=O)=O"], True))
    
    # jo1c00719_0008.png
    # print(get_reactants("N/C(C1=CC=CC=C1)=C(SC2=CC=C(Cl)C=C2)\SC3=CC=C(Cl)C=C3", ["[R]C1=CC=C(S)C=C1", "C=C(N=[N+]=[N-])C1=CC=CC=C1"], ["N/C(C1=CC=CC=C1)=C(SC2=CC=C([R])C=C2)/SC3=CC=C([R])C=C3"], True))
    
    # jo3c02090_0008.png
    # print(get_reactants("FC1=CC=C(/C=N/SC2=CC=C(C)C=C2)C=C1", ["[R]S", "FC1=CC=C(CN)C=C1"], ["[R]S/N=C/C1=CC=CC=C1"], True))
    
    # cs8b03302_0005_0007.png
    # print(get_reactants("CC1=CC=C(S(=O)(CC(OC)C2=CC=C(C)C=C2)=O)C=C1", ["CC1=CC=C(S(=O)(NN)=O)C=C1", "CC1=CC=C(C=C)C=C1", "CO"], ["CC1=CC=C(S(=O)(CC(OC)C2=CC=C(C)C=C2)=O)C=C1"]))