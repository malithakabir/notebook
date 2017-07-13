
"""RDKit Conformer Browser
Derived from Greg Landrum's code by Malitha Humayun Kabir
As a part of GSoC 2017
Project : RDKit - 3Dmol.js integration
Mentors: Paul Czodrowski and Greg Landrum
Date: 12th July 2017
"""

import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
#from ipywidgets import interact, interactive, fixed
from ipywidgets import Layout, Label, Button, Box, HBox, VBox
import ipywidgets as widgets
from IPython.display import display
import IPython.display
import time


drawing_types_3d=['line', 'cross', 'stick', 'cartoon', 'sphere', 'surface', 'ballstick']
prop_rdkit=['MolLogP', 'MolMR', 'MolWt', 'ExactMolWt', 'HeavyAtomCount',
            'HeavyAtomMolWt', 'NHOHCount', 'NOCount', 'NumHAcceptors',
            'NumHDonors', 'NumHeteroatoms', 'NumRotatableBonds', 'NumValenceElectrons']
color_scheme_3d=['default', 'greenCarbon', 'cyanCarbon', 'magentaCarbon', 
                 'yellowCarbon', 'whiteCarbon', 'orangeCarbon', 'purpleCarbon', 
                 'blueCarbon', 'ssPyMOL', 'ssJmol', 'Jmol', 'amino', 
                 'shapely', 'nucleic', 'chain', 'chainHetatm', 'prop']


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
    
    
    

def update3D(model_id):
    
    uid=globals()['rdkit_wg_dict'][model_id]
    
    view=globals()['view_'+uid]
    molDictKey=globals()['molTupleId_'+uid].value
    confId=globals()['confId_'+uid].value
    
    calculate_prop=globals()['prop_calc_wg_'+uid].value
    
    useDrawAs=globals()['useDrawAs_wg_'+uid].value
    drawAs=globals()['drawAs_wg_'+uid].value
    
    color=globals()['colorScheme_'+uid].value
    labelConf=globals()['confLabel_'+uid].value
    labelAtom=globals()['atomLabel_'+uid].value
    
    #### Molecule and conformer selection chunk
    # 'moldict_' is dict ... NOT widget
    mol = globals()['moldict_'+uid][molDictKey]
    if mol.GetNumConformers()>0:
        # Update conformer select slider
        globals()['confId_'+uid].max=mol.GetNumConformers()-1
        # Get data
        sconf = mol.GetConformer(confId)
        xyz = sconf.GetPositions()
        OwningMol = sconf.GetOwningMol()
        
    #### To Do: If the Mol object does not have any conformers????
    
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
        prop_value_str = selected_prop + ' : ' + str(selected_prop_value)
        globals()['prop_precalc_view_'+uid].value = prop_value_str
    else:
        globals()['prop_precalc_view_'+uid].value = 'No precalculated property found!'
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
    prop_cmd_str="Descriptors."+calculate_prop + "(OwningMol)"
    globals()['prop_calc_view_'+uid].value = calculate_prop + ' : ' + str(eval(prop_cmd_str))
    
    
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
        #
        # This is exception for surface
        #
        if drawAs is 'surface':
            view.addSurface({}, '$3Dmol.SurfaceType.VDW');
        elif drawAs is 'ballstick':
            view.setStyle({},{'stick':{'radius':'0.2','colorscheme': color},
                              'sphere':{'radius':'0.4', 'colorscheme': color}}
                         );
        else:
            view.setStyle({},{globals()['Style_'+uid]:{'colorscheme': color}})
    
    
    #
    # Labeling conformer
    #
    if labelConf is True:
        label = molDictKey + ':' + str(confId)
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
    view.setBackgroundColor(globals()['background_'+uid].value)
    
    #print(drawAs)
    # zoomTo does not work well for surface and label... so, zoomTo should not be default settings
    #view.zoomTo()
    display(view.update())
    
    
def handle_change(change):
    update3D(change.owner._model_id)
    
def handle_button(b):
    b.icon='check'
    b.description="Done!"
    update3D(b._model_id)
    
def browseMolConformers2(mol_obj_list,confId=None, useDrawAs=False, drawAs=None):
    
    # Generating common suffix for all the widgets
    uid=str(time.time()).replace('.','')
    
    #########################################################################
    ### Names for interactive widgets except buttons and property viewers ###
    #########################################################################
    
    # Listing interactive widgets
    wg_int=['molTupleId_', 'confId_', 
            'prop_calc_wg_', 'prop_precalc_wg_', 
            'useDrawAs_wg_', 'drawAs_wg_',
            'colorScheme_', 'confLabel_', 'atomLabel_', 'background_']
    
    ###################################
    ### Generating required widgets ###
    ###################################
    
    # Handling molecules and conformers
    
    globals()['moldict_'+uid] = processMolContainingObj(mol_obj_list)
    
    globals()['molTupleId_'+uid] = widgets.Dropdown(description='Mol',
                                                    options=list(globals()['moldict_'+uid].keys()),
                                                    value=list(globals()['moldict_'+uid].keys())[0])
    globals()['confId_'+uid] = widgets.IntSlider(description='confId', 
                                                 min=0,max=9,step=1,value=0)
    
    # prop handling
    globals()['prop_precalc_view_'+uid] = widgets.HTML(description='', value='initializing...')
    
    globals()['prop_precalc_wg_'+uid] = widgets.Dropdown(description='precalc', 
                                                         options=['select'], value='select')
    globals()['prop_calc_view_'+uid] = widgets.HTML(description='', value='initializing...')
    
    globals()['prop_calc_wg_'+uid] = widgets.Dropdown(description='calc',
                                                      options=globals()['prop_rdkit'], value='MolLogP')
    
    # Assign global Style
    globals()['Style_'+uid] = 'stick'
    
    globals()['useDrawAs_wg_'+uid] = widgets.Dropdown(description='useDrawAs',
                                                      options=[False, True],
                                                      value=False)
    
    # resolving drawAs and useDrawAs
    if useDrawAs is True:
        globals()['useDrawAs_wg_'+uid].value = True
        if drawAs is not None:
            globals()['Style_'+uid] = drawAs
        else:
            drawAs = globals()['Style_'+uid]
    else:
        drawAs = globals()['Style_'+uid]
    
    # style
    globals()['drawAs_wg_'+uid] = widgets.Dropdown(description='drawAs',
                                                   options=globals()['drawing_types_3d'],
                                                   value=drawAs)
    
    # colors & labels
    globals()['colorScheme_'+uid] = widgets.Dropdown(description='color',
                                                     options=globals()['color_scheme_3d'],
                                                     value='default')
    globals()['confLabel_'+uid] = widgets.Checkbox(description='confLabel', value=False)
    
    globals()['atomLabel_'+uid] = widgets.Checkbox(description='atomLabel', value=False)
    
    bgcolors_3d = ['0xeeeeee', '0x000000', '0xffffff']
    
    globals()['background_'+uid] = widgets.Dropdown(description='background', 
                                                    options=bgcolors_3d, value=bgcolors_3d[1])
    
    globals()['button_'+uid] = widgets.Button(description="Start!", button_style='success')
    
    
    # Left hand panel (container holding table for viewer)
    size = (435, 430)
    vwContainerLayout=Layout(width=str(size[0]+5)+'px',
                       height=str(size[1]+5)+'px', 
                       border='solid')
    wg_html = widgets.HTML(
        '''<table><tr><td id="%s"></td></tr></table>'''%uid, layout=vwContainerLayout)
    
    globals()['wgLeftWidget_'+uid] = wg_html
    
    ###################################################################
    ### updating widgets to global dict and assigninf interactivity ###
    ###################################################################
    
    # Creating global dict for interactive widgets
    if 'rdkit_wg_dict' not in globals():
        globals()['rdkit_wg_dict'] = dict()
    
    # Adding interactive widgets in dict
    for i in wg_int:
        globals()['rdkit_wg_dict'].update({globals()[i+uid]._model_id : uid}) 
    globals()['rdkit_wg_dict'].update({globals()['button_'+uid]._model_id:uid})
    
    
    # Enabling observe functions for interactive widgets except button
    for i in wg_int:
        globals()[i+uid].observe(handle_change, names='value')
    # Enabling on_click functionality for button
    globals()['button_'+uid].on_click(handle_button)
    
    
    ##########################
    ### Creating interface ###
    ##########################
    
    # Right hand panel (widgets)
    itemLayout=Layout(display='flex',
                          flex_flow='row',
                          justify_content='space-between')
    wgListBox=list()
    wgListBox.append(Box([Label(value='calc'),globals()['prop_calc_view_'+uid]], layout=itemLayout))
    wgListBox.append(Box([Label(value='precalc'),globals()['prop_precalc_view_'+uid]], layout=itemLayout))
    for i in wg_int:
        wg=globals()[i+uid]
        wgListBox.append(Box([Label(value=wg.description),wg], layout=itemLayout))
        wg.description=''
    wgListBox.append(Box([Label(value=''), globals()['button_'+uid]], layout=itemLayout))
    
    # Right hand panel (global variable)
    globals()['wgRightBox_'+uid]=widgets.VBox(wgListBox, layout=Layout(border='solid', width='45%'))
    
    
    
    # Combining Left and Right hand panel (global variable)
    globals()['all_wg_'+uid]=widgets.HBox([globals()[i+uid] for i in ['wgLeftWidget_', 'wgRightBox_']])
    
    # Displaying combined varibale
    display(globals()['all_wg_'+uid])
    
    ###############################################################
    ### Inserting 3DMol.js viewer in existing container (table) ###
    ###############################################################
    
    view = py3Dmol.view(width=size[0],height=size[1])
    view.setBackgroundColor(globals()['background_'+uid].value)
    view.zoomTo()
    globals()['view_'+uid]=view
    display(view.insert(uid))
    
    
    