"""
File that contains code to evalute algorithim
Contains functionaility for both exact string mathcing using canonical SMILES and Tanimoto similiarity
"""

from rdkit import Chem, DataStructs
from typing import Set, Tuple, List, Dict


def evaluate(ground_truth_reactants: List[str], generated_reactants: List[str]) -> Tuple[List[List[int]], List[List[float]]]:
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
        curr_tanimoto = []
        for ground_truth in canonical_ground_truth:
            curr_exact_match.append(int(generated == ground_truth))
            curr_tanimoto.append(tanimoto_similarity(generated, ground_truth))
        exact_match_matrix.append(curr_exact_match)
        tanimoto_matrix.append(curr_tanimoto)
    
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