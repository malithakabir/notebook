"""Ligand Extraction from mol object (RDKit PDB file reader should be used in mol object generation)
Originally written by Greg Landrum and updated by Malitha Humayun Kabir as a part of GSoC 2017 
Project : RDKit - 3Dmol.js integration
Mentors: Paul Czodrowski and Greg Landrum
Date: 28th July 2017
Email# malitha12345@gmail.com
"""

from collections import namedtuple
from rdkit import Chem

## Ligand Extract

ExtractResult = namedtuple('ExtractResult',('match','rest'))
def ExtractMolAtomsMatchingQuery(mol,func, confId,sanitize, includeAttachedHs):
    """ func should take an atom index """
    
    match = [x for x in range(mol.GetNumAtoms()) if func(x)]
    
    if includeAttachedHs:
        # bring over H atoms attached to the atoms matching the query (if necessary)
        for aid in match:
            for nbr in mol.GetAtomWithIdx(aid).GetNeighbors():
                if nbr.GetAtomicNum()==1 and nbr.GetIdx() not in match:
                    match.append(nbr.GetIdx())
                    
    ##### This chunk clones the input molecule and deletes atoms###
    res2 = Chem.RWMol(mol)
    for i in sorted(match, reverse=True):
        res2.RemoveAtom(i)
    if sanitize:
        Chem.SanitizeMol(res2)
    #### Creating new molecule...
    res = Chem.RWMol()
    
    # start with all atoms and their coordinates:
    # Should probably also handle multiple conformers
    oconf = mol.GetConformer(confId)
    nconf = Chem.Conformer(len(match))
    nconf.SetId(oconf.GetId())
    old_new_map={}
    for i,aid in enumerate(match):
        res.AddAtom(mol.GetAtomWithIdx(aid))
        nconf.SetAtomPosition(i,oconf.GetAtomPosition(aid))
        old_new_map[aid] = i
    res.AddConformer(nconf)
    
    # bonds:
    for i,aid in enumerate(match):
        for nbr in mol.GetAtomWithIdx(aid).GetNeighbors():
            if nbr.GetIdx() not in old_new_map:
                continue
            bnd = mol.GetBondBetweenAtoms(nbr.GetIdx(),aid)
            if aid != bnd.GetBeginAtomIdx():
                continue
            res.AddBond(old_new_map[aid],old_new_map[nbr.GetIdx()],bnd.GetBondType())
    if sanitize:
        Chem.SanitizeMol(res)
    return ExtractResult(res.GetMol(), res2.GetMol())
    
def ExtractMolFragment(mol, ResName, confId=-1,sanitize=True, includeAttachedHs=True):
    " see ExtractMolAtomsMatchingQuery() for kwargs "
    ids=list()
    for Idx in range(mol.GetNumAtoms()-1):
        atom=mol.GetAtomWithIdx(Idx)
        if atom.GetPDBResidueInfo().GetResidueName() == ResName:
            ids.append(atom.GetIdx())
    return ExtractMolAtomsMatchingQuery(mol, lambda x,y = ids : x in y, confId,sanitize, includeAttachedHs)
    
    
    
    
    