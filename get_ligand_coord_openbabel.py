


# It reads loacl file exports sdf for selected ligand
# But the exported sdf does not have bond information
# Need more work on it
# It should be the preferred approach

pdbID='4i22.pdb'
resNameIndex=297

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

pybelmol.write("sdf", "outputfile.sdf", overwrite=True)


#from rdkit.Chem import rdmolfiles
#mol=rdmolfiles.MolFromPDBFile('outputfile.sdf')