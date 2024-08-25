"""
File containing code for template matching a core to substrate and fragment and rule retrieval
"""

from rdkit import Chem
from rdkit.Chem.rdchem import Mol, Atom, Bond, BondType, EditableMol
from rdkit.Chem.rdChemReactions import ChemicalReaction
from typing import Set, Tuple, List, Dict

import networkx as nx
import matplotlib.pyplot as plt
from utils import node_match, edge_match, find_atom, create_atom, highlight_reaction_core, generate_graph, generate_graph_with_indicies_as_labels, ADD_BOND_FLAG, REMOVE_BOND_FLAG, MODIFY_BOND_FLAG
import core_extraction as ce

from classes.partial_molecule import Fragment, Rule

##################################################################################################
# Template matching/Fragment retrieval
##################################################################################################

def get_reactants_for_substrate(substrate: str, reactant_core:Mol, product_core: Mol, rules: Rule) -> Set[Mol]:
    """
    Matches the core template to a fixed substrate and returns proper reactants
    """
    substrate_mol = Chem.MolFromSmiles(substrate)
    # Chem.rdmolops.Kekulize(substrate_mol, clearAromaticFlags=True)
    for atom in substrate_mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    
    highlight_reaction_core(substrate_mol, set(), set(), "substrate.png")

    for rule in rules:
        print(rule)
    substrate_to_core = get_subgraph(substrate_mol, product_core)
    core_to_substrate = {substrate_to_core[key] : key for key in substrate_to_core}
    print(core_to_substrate)
    modify_rules = []
    substrate_mol = EditableMol(substrate_mol)
    for rule in rules:
        # we need to make sure rule.start_atom and rule.end_atom are in core_to_substrate
        # if they are not then these are leaving group atoms we need to add back
        if rule.start_atom in core_to_substrate:
            index1 = core_to_substrate[rule.start_atom]
        else:
            index1 = substrate_mol.AddAtom(create_atom(find_atom(rule.start_atom, [reactant_core])))
            core_to_substrate[rule.start_atom] = index1
        if rule.end_atom in core_to_substrate:
            index2 = core_to_substrate[rule.end_atom]
        else:
            index2 = substrate_mol.AddAtom(create_atom(find_atom(rule.end_atom, [reactant_core])))
            core_to_substrate[rule.end_atom] = index2
        
        # check the rule types and perform needed action
        if rule.rule_type == ADD_BOND_FLAG:
            substrate_mol.AddBond(index1, index2, rule.end_bond_type)
        elif rule.rule_type == REMOVE_BOND_FLAG:
            substrate_mol.RemoveBond(index1, index2)
        else: # rule.rule_type == MODIFY_BOND_FLAG
            modify_rules.append(rule)
    
    substrate_mol = substrate_mol.GetMol()
    # deal with modify bonds seperatly
    for rule in modify_rules:
        index1, index2 = core_to_substrate[rule.start_atom], core_to_substrate[rule.end_atom]
        substrate_mol.GetBondBetweenAtoms(index1, index2).SetBondType(rule.end_bond_type)
    
    smiles = Chem.MolToSmiles(substrate_mol).split(".")
    highlight_reaction_core(substrate_mol, set(), set(), "new_reactants.png")
    return [smile for smile in smiles]
    # return smiles
            

def get_subgraph(substrate_mol: Mol, product_template_mol: Mol) -> Dict[int, int]:
    """
    Checks to see if core_product_fragment is a subgraph of substrate_graph
    Returns {} if not, otherwise returns a dict mapping substrate to core atoms
    """
    substrate_graph, core_graph = generate_graph_with_indicies_as_labels(substrate_mol), generate_graph(product_template_mol)

    matcher = nx.algorithms.isomorphism.GraphMatcher(substrate_graph, core_graph, node_match=node_match, edge_match=edge_match)
    subgraphs = [subgraph for subgraph in matcher.subgraph_isomorphisms_iter()]
    if len(subgraphs) == 0:
        return {}
    else:
        return subgraphs[0]

if __name__ == "__main__":
    # core = ce.get_reaction_core(["O=C1C(C(OCC)=O)CCC1", "C=CC(C[R1])=C"], ["O=C1C(C(OCC)=O)(CCC(C[R1])=C)CCC1"])
    
    # wiley2_table_3
    # core = ce.get_reaction_core(["O=C(NC)C1=CC=CC=C1/C=C([R1])/[R2]"], ["O=C(N1C)C2=CC=CC=C2C1C([R1])[R2]"])
    # reactant_core, product_core = ce.extend_reaction_core(core[1], core[2], core[0])
    # highlight_reaction_core(core[1][0], set(), set(), "og_reactant.png")
    # highlight_reaction_core(core[2][0], set(), set(), "og_product.png")
    # highlight_reaction_core(reactant_core, set(), set(), "reactant_core.png")
    # highlight_reaction_core(product_core, set(), set(), "product_core.png")
    # mols = get_reactants_for_substrate("O=C(N1C)C2=CC=CC=C2C1CC3=CC=C(C(C)=O)C=C3", core[0], product_core, core[3])
    # # highlight_reaction_core(Chem.MolFromSmiles(mols), set(), set(), "new_reactant.png")
    # print(mols)
    
    # cs8bo3302
    core = ce.get_reaction_core(["CC1=CC=C(S(=O)(NN)=O)C=C1", "CC1=CC=C(C=C)C=C1", "CO"], ["CC1=CC=C(S(=O)(CC(OC)C2=CC=C(C)C=C2)=O)C=C1"])
    reactant_core, product_core = ce.extend_reaction_core(core[1], core[2], core[0])

    highlight_reaction_core(reactant_core, set(), set(), "reactant_core.png")
    highlight_reaction_core(product_core, set(), set(), "product_core.png")
    mols = get_reactants_for_substrate("CC1=CC=C(S(CC(OC)C2=CC=C(C)C=C2)(=O)=O)C=C1", core[0], product_core, core[3])
    print(mols)
    
    # wiley26_scheme_2
    # core = ce.get_reaction_core(["[R]O"], ["[R]OC1CCCO1", "[H][H]"])
    # product_core = ce.extend_reaction_core(core[1], core[2], core[0])
    # mols = get_reactants_for_substrate("COC1=CC=C(CCOC2CCCO2)C=C1", core[0], product_core)
    
    # wiley38_table_2
    # core = ce.get_reaction_core(["IC1CCOCC1", "OC(C1=CC([R])=CC=C1)=O"], ["O=C(C1CCOCC1)C2=CC=CC([R])=C2"])
    # reactant_core, product_core = ce.extend_reaction_core(core[1], core[2], core[0])
    
    # highlight_reaction_core(reactant_core[0], set(), set(), "reactant_core.png")
    # highlight_reaction_core(product_core[0], set(), set(), "product_core.png")
    # mols = get_reactants_for_substrate("O=C(C1CCOCC1)C2=CC=C(C(C)(C)C)C=C2", core[0], product_core, core[3])
    # print(mols)
    
    # rsc4_scheme_2
    # core = ce.get_reaction_core(["[R2]C(C1=CC([R1])=CC=C1)=C", "[R4]C(C1=CC([R3])=CC=C1)=C", "[R]O"], ["[R2]C(O[R])(C1=CC([R1])=CC=C1)CCC(O[R])(C2=CC=CC([R3])=C2)[R4]"])
    # product_core = ce.extend_reaction_core(core[1], core[2], core[0])
    # # print(product_core)
    # highlight_reaction_core(core[0].get_mol(), set(), set(), "reactant_core.png")
    # highlight_reaction_core(product_core.get_mol(), set(), set(), "product_core.png")
    # mols = get_reactants_for_substrate("COC(C1=CC=CC=C1)CCC(OC)C2=CC=CC=C2", core[0], product_core)
    
    # wiley27_scheme_2
    # core = ce.get_reaction_core(["[R1]N(C1=C(/C2=N\[H])C=CC([R2])=C1)C(N2C3=CC=C([R3])C=C3)=O"], ["O=C(N([R1])C(C=C([R2])C=C1)=C1C2=N3)N2C4=C3C=C([R3])C=C4", "[H][H]"])
    # reaction_core, product_core = ce.extend_reaction_core(core[1], core[2], core[0])
    # # # print(product_core)
    # highlight_reaction_core(reaction_core[0], set(), set(), "reactant_core.png")
    # highlight_reaction_core(product_core[0], set(), set(), "product_core.png")
    # mols = get_reactants_for_substrate("O=C(N(CC1=CC=CC=C1)C(C=CC=C2)=C2C3=N4)N3C5=C4C=C(CCO)C=C5", core[0], product_core, core[3])   
    # print(mols)


    # RSC22
    # core = ce.get_reaction_core(["[H]C1=C([R])N=C2N1C=CC=C2", "[R2]N[R1]"], ["[R]C1=C(N([R1])[R2])N2C=CC=CC2=N1"])
    # reactant_core, product_core = ce.extend_reaction_core(core[1], core[2], core[0])
    # highlight_reaction_core(core[1][0], set(), set(), "og_reactant.png")
    # highlight_reaction_core(core[2][0], set(), set(), "og_product.png")
    # highlight_reaction_core(reactant_core, set(), set(), "reactant_core.png")
    # highlight_reaction_core(product_core, set(), set(), "product_core.png")
    # mols = get_reactants_for_substrate("N12C=CC=CC1=NC(C3=CC=CC=C3)=C2N4C=NC=C4", core[0], product_core, core[3])
    # print(mols)
    # highlight_reaction_core(Chem.MolFromSmiles(mols), set(), set(), "new_reactant.png")

    # print(mols)   

    # product_core = ce.extend_reaction_core(core[1], core[2], core[0])
    # highlight_reaction_core(core[0].get_mol(), set(), set(), "extended_reactant_core.png")
    # highlight_reaction_core(product_core.get_mol(), set(), set(), "extended_product_core.png")
    # mols = get_reactants_for_substrate("O=C(C1CCOCC1)C2=CC=C(C)C=C2", core[0], product_core)
    # mols = get_reactants_for_substrate("C1=CC=CC=C1", core[0], product_core)
    # index = 1
    # for mol in mols:
        # highlight_reaction_core(mol, set(), set(), "new_reactant" + str(index) + ".png")
        # print(Chem.MolToSmiles(mol))
        # index += 1
    