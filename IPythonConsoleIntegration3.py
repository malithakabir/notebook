
"""RDKit Conformer Browser
Derived from Greg Landrum's code by Malitha Humayun Kabir
As a part of GSoC 2017
Project : RDKit - 3Dmol.js integration
Mentors: Paul Czodrowski and Greg Landrum
Date: 11th July 2017
"""

import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from ipywidgets import interact, interactive, fixed
import ipywidgets as widgets
from IPython.display import display
import time


def processMolContainingObj(ms):
    
    try:
        moldict = dict(ms)
    except TypeError:
        if type(ms) is tuple:
            moldict=list()
            moldict.append(ms)
            moldict = dict(moldict)
        elif hasattr(ms, '__iter__') is False:
            moldict=list()
            moldict.append(('m0', ms))
            moldict = dict(moldict)
        elif type(ms) is list:
            ms_name=['m'+str(x) for x in range(len(ms))]
            ms2=[(ms_name[i],ms[i]) for i in range(len(ms))]
            moldict = dict(ms2)
    return moldict
    
    
    
def addMolToViewForScrolling(moldict, view, uid,
                             # property viewer
                             show_precalc_prop,
                             show_calc_prop,
                             # molecule and conformer selector
                             molDictKeys,
                             confId,
                             precalculated_prop,
                             calculate_prop,
                             useDrawAs,
                             drawAs,
                             color,
                             labelConf,
                             labelAtom
                            ):
    
    
    #### Molecule and conformer selection chunk
    
    # Get mol from supplied list object
    mol = moldict[molDictKeys]
    if mol.GetNumConformers()>0:
        # Update conformer select slider
        globals()['confId_slider_'+uid].max=mol.GetNumConformers()-1
        # Get data
        sconf = mol.GetConformer(confId)
        xyz = sconf.GetPositions()
        OwningMol = sconf.GetOwningMol()
        
    #### To Do: If the Mol object does not have any conformers????
    
    
    #### removing descriptor widget description
    globals()['moldictkeys_'+uid].description='Mol'
    globals()['prop_precalc_view_'+uid].description=''
    globals()['prop_calc_view_'+uid].description=''
    
    #### Descriptor handling chunkGetting descriptors
    
    # Get available properties and update dropdown widget 
    all_prop_from_mol=list(OwningMol.GetPropNames())
    
    if len(all_prop_from_mol)>0:
        # update widget
        globals()['prop_precalc_wg_'+uid].options=all_prop_from_mol
        
        # Get selected property name and associated value
        selected_prop=eval('prop_precalc_wg_'+uid+'.value')
        selected_prop_value=OwningMol.GetProp(selected_prop)
        
        # Update viewer
        prop_value_str = 'precalc: ' + selected_prop + ': ' + str(selected_prop_value)
        globals()['prop_precalc_view_'+uid].value = prop_value_str
        
        #
        #
    #
    #
    # Update viewer for rdkit offered descriptor with selected property name and value
    # This is real time calculation
    #
    # Descriptor calculation schema eval("Descriptors.TPSA(OwningMol)")
    # In above line descriptor_list is TPSA
    #
    prop_cmd_str="Descriptors."+calculate_prop+"(OwningMol)"
    globals()['prop_calc_view_'+uid].value='prop: ' + calculate_prop + ': ' + str(eval(prop_cmd_str))
    
    
    ##### 3Dmol.js viewer realated stuffs below
    
    #
    # Clearing previous 3Dmol objects withoiut resetting view
    #
    view.removeAllModels()
    view.removeAllSurfaces()
    view.removeAllLabels()
    
    #
    #### Adding model to viewer
    #
    if mol.GetNumAtoms()>=999 or drawAs == 'cartoon':
        # py3DMol is happier with TER and MASTER records present
        pdb = Chem.MolToPDBBlock(mol,flavor=0x20|0x10)
        view.addModel(pdb,'pdb')
    else:
        # py3Dmol does not currently support v3k mol files, so
        # we can only provide those with "smaller" molecules
        mb = Chem.MolToMolBlock(mol,confId=confId)
        view.addModel(mb,'sdf')
        
    #
    #### Making decision about model rendering style and color
    #
    
    if useDrawAs is False:
        #use from globalStyle
        view.setStyle({},{globals()['Style_'+uid]:{'colorscheme': color}})
    else:
        #update global style and use that
        globals()['Style_'+uid] = drawAs
        view.setStyle({},{globals()['Style_'+uid]:{'colorscheme': color}})
    
    #
    # This is exception for surface
    #
    if drawAs is 'surface':
        view.addSurface({}, '$3Dmol.SurfaceType.VDW');
        
    if drawAs is 'ballstick':
        view.setStyle({},{'stick':{'radius':'0.2','colorscheme': color},
                          'sphere':{'radius':'0.4', 'colorscheme': color}}
                     );
    
    #
    # Labeling conformer
    #
    if labelConf is True:
        label = molDictKeys + ':' + str(confId)
        view.addLabel(label, {'backgroundColor':'gray', 'fontColor':'white',
                              'showBackground':'true', 'alignment':'bottomCenter'})
    
    #
    # Labeling atom
    #
    if labelAtom is True:
        label_create=[OwningMol.GetAtomWithIdx(i).GetSymbol()+
                      str(OwningMol.GetAtomWithIdx(i).GetIdx()+1) 
                      for i in range(sconf.GetNumAtoms())
                     ]
        
        i = None
        for i in range(sconf.GetNumAtoms()):
            view.addLabel(label_create[i], {'inFront' : 'false', 
                                            'fontSize' : '12',
                                            'fontColor':'gray',
                                            'showBackground':'false',
                                            'position' : {'x' : xyz[i][0],
                                                          'y' : xyz[i][1],
                                                          'z' : xyz[i][2]
                                                       }
                                           })
    
    #print(drawAs)
    # zoomTo does not work well for surface and label... so, zoomTo should not be default settings
    #view.zoomTo()
    return view.update()

drawing_types_3d=['line', 'cross', 'stick', 'cartoon', 'sphere', 'surface', 'ballstick']
prop_rdkit=['MolLogP', 'MolMR', 'MolWt', 'ExactMolWt', 'HeavyAtomCount',
            'HeavyAtomMolWt', 'NHOHCount', 'NOCount', 'NumHAcceptors',
            'NumHDonors', 'NumHeteroatoms', 'NumRotatableBonds', 'NumValenceElectrons']
color_scheme_3d=['default', 'greenCarbon', 'cyanCarbon', 'magentaCarbon', 
                 'yellowCarbon', 'whiteCarbon', 'orangeCarbon', 'purpleCarbon', 
                 'blueCarbon', 'ssPyMOL', 'ssJmol', 'Jmol', 'amino', 
                 'shapely', 'nucleic', 'chain', 'chainHetatm', 'prop']

def browseMolConformers(mol_obj_list,confId=None, useDrawAs=False, drawAs=None):
    
    view_size = (400, 400)
    view_bgcolor = '0xeeeeee'
    view = py3Dmol.view(width=view_size[0],height=view_size[1])
    view.setBackgroundColor(view_bgcolor)
    view.zoomTo()
    
    containerId=str(time.time()).replace('.','')
    html=html=widgets.HTML('''<table><tr><td id="%s"></td></tr></table>'''%containerId)
    display(html)
    display(view.insert(containerId))
    
    uid=view.uniqueid
    
    # processing molecule containing object
    moldict = processMolContainingObj(mol_obj_list)
    
    # browse
    globals()['moldictkeys_'+uid] = widgets.Dropdown(options=list(moldict.keys()),
                                                     value=list(moldict.keys())[0])
    globals()['confId_slider_'+uid] = widgets.IntSlider(min=0,max=9,step=1,value=0)
    
    
    # prop handling
    globals()['prop_precalc_view_'+uid]=widgets.HTML(value='initializing...', disabled=False)
    globals()['prop_precalc_wg_'+uid]=widgets.Dropdown(options=['select'],value='select')
    globals()['prop_calc_view_'+uid]=widgets.HTML(value='initializing...', disabled=False)
    globals()['prop_calc_wg_'+uid] = widgets.Dropdown(options=globals()['prop_rdkit'],value='MolLogP')
    
    
    # Assign global Style
    globals()['Style_'+uid] = 'stick'
    useDrawAs_wg = widgets.Dropdown(options=[False, True],value=False)
    drawAs = globals()['Style_'+uid]
    
    if useDrawAs is True:
        globals()['useDrawAs_wg_'+uid] = widgets.Dropdown(options=[False, True],value=True)
        if drawAs is not None:
            globals()['Style_'+uid] = drawAs
    
    
    # style
    drawAs_wg = widgets.Dropdown(options=globals()['drawing_types_3d'],value=drawAs)
    
    # colors & labels
    colorScheme = widgets.Dropdown(options=globals()['color_scheme_3d'],value='default')
    confLabelCheckBok = widgets.Checkbox(description='confLabelCheckBok', value=False)
    atomLabelCheckBox = widgets.Checkbox(description='atomLabelCheckBox', value=False)
    
    
    # Now start interacting
    result=interact(addMolToViewForScrolling, view=fixed(view), moldict=fixed(moldict), uid=fixed(uid),
                    show_precalc_prop=eval('prop_precalc_view_'+uid),
                    show_calc_prop=eval('prop_calc_view_'+uid),
                    molDictKeys= eval('moldictkeys_'+uid), 
                    confId=eval('confId_slider_'+uid),
                    precalculated_prop=eval('prop_precalc_wg_'+uid),
                    calculate_prop=eval('prop_calc_wg_'+uid),
                    useDrawAs=useDrawAs_wg, 
                    drawAs=drawAs_wg,
                    color=colorScheme, labelConf=confLabelCheckBok, labelAtom=atomLabelCheckBox
                    ####
                   );
    return result

