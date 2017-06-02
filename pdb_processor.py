
import openbabel
import pybel
import urllib2
import time

#retrievePDB will download pdb file from RCSB PDB

def retrievePDB(pdbID):
    url_string="https://files.rcsb.org/download/PDB_ID.pdb"
    url_formatted=url_string.replace("PDB_ID.pdb", pdbID)
    
    print("retrieving from "+url_formatted)
    
    response = urllib2.urlopen(url_formatted)
    html = response.read()
    
    with open(pdbID, "wb") as handle:
        handle.write(html)
    
    handle.close()
    return pdbID+' downloaded from RCSB PBD'




# getNonStandardResName only reads from loacl directory
# output non standard residue index and non standard residue name
# Water molecules automatically removed

def getNonStandardResName(pdbID):
    
    stdResList=["ala","arg","asn","asp","asx","cys",
                "glu", "gln","glx","gly","his","ile",
                "leu","lys","met","phe","pro","ser",
                "thr","trp","tyr","val"]
    stdResListUpperCase=[s.upper() for s in stdResList]
    
    mol=pybel.readfile('pdb', '4i22.pdb').next()
    res=mol.residues
    resNameAll=[x.name for x in res]
    
    stdResNameDecision=[]
    for i in resNameAll:
        stdResNameDecision.append(any([stdResName==i for stdResName in stdResListUpperCase]))
    
    nonStdResIndex=[i for i, x in enumerate(stdResNameDecision) if x==False]
    
    nonStdResName=[resNameAll[x] for x in nonStdResIndex]
    
    nonStdResNameFinal=[]
    nonStdResIndexFinal=[]
    
    for i in range(len(nonStdResName)):
        resNameTmp=nonStdResName[i]
        if resNameTmp!='HOH':
            nonStdResNameFinal.append(nonStdResName[i])
            nonStdResIndexFinal.append(nonStdResIndex[i])
    
    return [nonStdResNameFinal, nonStdResIndexFinal]





def extractLigandOpenBabel(pdbID,resNameIndex):
    
    mol=pybel.readfile('pdb', pdbID).next()
    resList=mol.residues
    resSelected=resList[resNameIndex]
    atomsList=resSelected.atoms
    
    molCreate = openbabel.OBMol()

    for atomSelected in atomsList:

        atomicNumberAssign=atomSelected.atomicnum
        coordsSel=atomSelected.coords

        a = molCreate.NewAtom()
        a.SetAtomicNum(atomicNumberAssign)   # atom
        a.SetVector(float(coordsSel[0]), float(coordsSel[1]), float(coordsSel[2])) # coordinates



    pybelmol = pybel.Molecule(molCreate)

    pybelmol.write("sdf", "outputfile_obabel.sdf", overwrite=True)
    
    return 'Ligand extraction done using OpenBabel...filename: outputfile_obabel.sdf'

    
# I am not in favor of the following approach 
# Please see my attempt in extractLigandOpenBabel
# The following exports wrong bond information

import __main__
__main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
from time import sleep
import pymol

def extractLigandPymol(pdbID,resNameSelected):
    
    pymol.finish_launching()
    
    pymol.cmd.load(pdbID)
    
    resSelect='resn RESIDUE'.replace('RESIDUE', resNameSelected)
    #pymol.cmd.png("outputfile.png", width=300, height=400)
    
    
    time.sleep(3) # delays for 3 seconds
    pymol.cmd.save('outputfile_pymol.sdf', resSelect)
    
    time.sleep(5) # delays for 5 seconds
    pymol.cmd.quit()
    
    return 'Ligand extraction done using Pymol...filename: outputfile_pymol.sdf'





