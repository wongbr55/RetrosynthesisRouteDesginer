"""
File that contains code to evalute algorithim
Contains functionaility for both exact string mathcing using canonical SMILES and Tanimoto similiarity
"""

from rdkit import Chem, DataStructs
from typing import Set, Tuple, List, Dict
from core_extraction import get_reaction_core, extend_reaction_core
from template_rule_retrieval import get_reactants_for_substrate
import csv

def evaluate_csv(ground_truth_csv: str, reaction_template_csv: str) -> Dict[str, Tuple[List[List[int]], List[float]] | str]:
    """
    Evaluates a csv file using exact string matching and tanimoto accuracy heurstics from evaluate function
    Assume csv format is of form: Figure, Product 1 SMILES, Product 2 SMILES, Reactant 1 SMILES, Reactant 2 SMILES, Reactant 3 SMILES
    Returns a dict of evaluation metrics for each substrate
    """
    
    # reaction cores maps a file to its reactant and product cores
    # to be used when evalauting files in the ground truth
    reaction_cores = {}
    with open(reaction_template_csv, "r") as csvfile:
        reaction_templates = csv.reader(csvfile)
        next(reaction_templates)
        for row in reaction_templates:
            filename = row[0]
            products = []
            reactants = []
            if row[1] != "": products.append(row[1])
            if row[2] != "": products.append(row[2])
            if row[3] != "": reactants.append(row[3])
            if row[4] != "": reactants.append(row[4])
            if row[5] != "": reactants.append(row[5])
            try:
                core = get_reaction_core(reactants, products)
                product_core = extend_reaction_core(core[1], core[2], core[0])
                reaction_cores[filename] = (core[0], product_core)
            except:
                reaction_cores[filename] = "ERROR UPON CORE EXTRACTION"
    
    evaluations = {}
    with open(ground_truth_csv, "r") as csvfile:
        ground_truths = csv.reader(csvfile)
        next(ground_truths)
        for row in ground_truths:
            # get the reaction cores:
            cores = []
            filename = row[0]
            for file in reaction_cores:
                if file in filename and reaction_cores[file] is not str:
                    cores.append(reaction_cores[file][0])
                    cores.append(reaction_cores[file][1])
                    break
                elif file in filename and reaction_cores[file] is str:
                    break
            if cores == []:
                evaluations[filename] = "ERROR UPON CORE EXTRACTION"
            else:
            # get generated reactants and evaluate againsy ground truth
                substrate = row[1]
                ground_truth_reactants = []
                if row[3] != "": ground_truth_reactants.append(row[3])
                if row[4] != "": ground_truth_reactants.append(row[4])
                if row[5] != "": ground_truth_reactants.append(row[5])
                
                try:
                    generated_reactants = get_reactants_for_substrate(substrate, cores[0], cores[1])
                    generated_reactant_smiles = [Chem.MolToSmiles(mol) for mol in generated_reactants]
                    evaluation = evaluate(ground_truth_reactants, generated_reactant_smiles)
                    evaluations[filename] = evaluation
                except:
                    evaluations[filename] = "ERROR UPON EVALUATION"
    
    return evaluations

def evaluate(ground_truth_reactants: List[str], generated_reactants: List[str]) -> Tuple[List[List[int]], List[float]]:
    """
    Uses exact canonical string match and tanimoto similarity to evalute algorithim for a single reaction
    Returns a tuple of similarity matricies for exact string matching and tanimoto similarity (in that order)
    """
    # get canonical strings for both ground truth and generated
    canonical_ground_truth = []
    canonical_generated = []
    for smiles in ground_truth_reactants:
        canonical_ground_truth.append(canonicalize_smiles(smiles))
    for smiles in generated_reactants:
        canonical_generated.append(canonicalize_smiles(smiles))
    
    # get smilarity matricies
    exact_match_matrix = []
    tanimoto_matrix = []
    for generated in canonical_generated:
        curr_exact_match = []
        curr_tanimoto = 0
        for ground_truth in canonical_ground_truth:
            curr_exact_match.append(int(generated == ground_truth))
            similarity = tanimoto_similarity(generated, ground_truth)
            if similarity == 1:
                curr_tanimoto += 1
        exact_match_matrix.append(curr_exact_match)
        tanimoto_matrix.append(curr_tanimoto / len(canonical_ground_truth))
    
    return exact_match_matrix, tanimoto_matrix
            

def canonicalize_smiles(smiles: str) -> str:
    """
    Returns the canonicalized form of smiles
    Returns "" upon error
    """
    try:
        canon_smiles = Chem.CanonSmiles(smiles)
    except:
        canon_smiles = ""
    return canon_smiles
    

def tanimoto_similarity(smiles1: str, smiles2: str) -> float:
    """
    Returns the tanimoto similarity of smiles1 and smiles2
    Returns -1 upon error
    """
    try:
        mol1, mol2 = Chem.MolFromSmiles(smiles1), Chem.MolFromSmiles(smiles2)
        fp1, fp2 = Chem.RDKFingerprint(mol1), Chem.RDKFingerprint(mol2)
        return DataStructs.FingerprintSimilarity(fp1, fp2)
    except:
        return -1