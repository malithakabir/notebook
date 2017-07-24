
"""RDKit Conformer Browser
Derived from Greg Landrum's code by Malitha Humayun Kabir
As a part of GSoC 2017
Project : RDKit - 3Dmol.js integration
Mentors: Paul Czodrowski and Greg Landrum
Date: 23th July 2017
Email# malitha12345@gmail.com
"""

import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from ipywidgets import Layout, Label, Button, Box, HBox, VBox
from ipywidgets import Dropdown, SelectMultiple, IntSlider, HTML, Checkbox, Button, Text
from IPython.display import display
import time
from collections import namedtuple


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
    
def ExtractMolFragment(mol, ResName, confId=-1,sanitize=False, includeAttachedHs=True):
    " see ExtractMolAtomsMatchingQuery() for kwargs "
    ids=list()
    for Idx in range(mol.GetNumAtoms()-1):
        atom=mol.GetAtomWithIdx(Idx)
        if atom.GetPDBResidueInfo().GetResidueName() == ResName:
            ids.append(atom.GetIdx())
    return ExtractMolAtomsMatchingQuery(mol, lambda x,y = ids : x in y, confId,sanitize, includeAttachedHs)
    
    
    
    
    
    ##### Visualization
    
bgcolors_3d = ['0xeeeeee', '0x000000', '0xffffff']

prop_rdkit=sorted([prop for prop,obj in Descriptors._descList])

drawing_types_3d=['line', 'cross', 'stick', 'cartoon', 'sphere', 'surface', 'ballstick']

color_scheme_3d=['default', 'greenCarbon', 'cyanCarbon', 'magentaCarbon', 
                 'yellowCarbon', 'whiteCarbon', 'orangeCarbon', 'purpleCarbon', 
                 'blueCarbon', 'ssPyMOL', 'ssJmol', 'Jmol', 'amino', 
                 'shapely', 'nucleic', 'chain', 'chainHetatm', 'prop']


def ProcessMolContainingObj(Mol):
    """This function checks whether the object type fits the requirements for rendering.
    If the oject doesn't have necessary attributes, it takes action to include that.
    """
    try:
        moldict = dict(Mol)
    except TypeError:
        if type(Mol) is tuple:
            moldict=list()
            moldict.append(Mol)
            moldict = dict(moldict)
        elif hasattr(Mol, '__iter__') is False:
            moldict=list()
            moldict.append(('0', Mol))
            moldict = dict(moldict)
        elif type(Mol) is list:
            Mol_keys=[str(x) for x in range(len(Mol))]
            Mol_2=[(Mol_keys[i],Mol[i]) for i in range(len(Mol))]
            moldict = dict(Mol_2)
    return moldict
    
    
    
def update3D(model_id):
    """ This function invoked whenever user interacts with widgets.
    It runs first time through handle_button() when the start button clicked.
    """
    
    uid=globals()['rdkit_wg_dict'][model_id]
    
    globals()['view_'+uid].removeAllModels()
    globals()['view_'+uid].removeAllSurfaces()
    globals()['view_'+uid].removeAllLabels()
    
    
    molDictKey=globals()['molId_'+uid].value
    
    if globals()['selectAllMols_'+uid].value is True:
        keys=list(globals()['moldict_'+uid].keys())
        for i in keys:
            globals()['rdkit_mol_selected_'+uid].append(i)
    else:
        if globals()['selectMultiMols_'+uid].value is True:
            globals()['rdkit_mol_selected_'+uid].append(molDictKey)
        else:
            globals()['rdkit_mol_selected_'+uid] = list()
            globals()['rdkit_mol_selected_'+uid].append(molDictKey)
    
    
    molNameList = list(set(globals()['rdkit_mol_selected_'+uid]))
    globals()['selected_mols_view_'+uid].value = ', '.join(molNameList)
    
    
    for i in molNameList:
        mol = globals()['moldict_'+uid][i]
        
        if mol.GetNumConformers()>1:
            
            allConfIds = range(mol.GetNumConformers()-1)
            globals()['confId_'+uid].options = allConfIds
            
            
            if globals()['selectAllConfs_'+uid].value is True:
                for i in allConfIds:
                    globals()['rdkit_conf_selected_'+uid].append(i)
            else:
                confId=globals()['confId_'+uid].value
                if globals()['selectMultiConfs_'+uid].value is True:
                    globals()['rdkit_conf_selected_'+uid].append(confId)
                else:
                    globals()['rdkit_conf_selected_'+uid] = list()
                    globals()['rdkit_conf_selected_'+uid].append(confId)
                    
            
        elif mol.GetNumConformers()==1:
            globals()['confId_'+uid].options=[0]
            confId=globals()['confId_'+uid].value
            globals()['rdkit_conf_selected_'+uid]=list()
            globals()['rdkit_conf_selected_'+uid].append(confId)
            globals()['selected_confs_view_'+uid].value = str(confId)
        
        
    confIdsList = list(set(globals()['rdkit_conf_selected_'+uid]))
    globals()['selected_confs_view_'+uid].value = ', '.join([str(x) for x in confIdsList])
    
    
    # Add molecule to viewer
    LigNumModel = -1
    for molName in molNameList:
        for confId in confIdsList:
            LigNumModel = LigNumModel + 1
            mol = globals()['moldict_'+uid][molName]
            mb = Chem.MolToMolBlock(mol,confId=confId)
            globals()['view_'+uid].addModel(mb,'sdf')
            
    
    
    if len(molNameList)==1:
        mol = globals()['moldict_'+uid][molDictKey]
        # For precalculated property
        try:
            all_prop_from_mol=list(mol.GetPropNames())
            if len(all_prop_from_mol)>0:
                # update widget
                globals()['prop_precalc_wg_'+uid].options=all_prop_from_mol
                # Get selected property
                prop_name=eval('prop_precalc_wg_'+uid+'.value')
                prop_value=mol.GetProp(prop_name)
                prop_str = prop_name + ' : ' + str(prop_value)
                # Update viewer
                globals()['prop_precalc_view_'+uid].value = prop_str
            else:
                globals()['prop_precalc_view_'+uid].value = 'No precalculated property found!'
        except:
            pass
        
        # For calculating rdkit supported property
        try:
            # Descriptor calculation schema eval("Descriptors.TPSA(mol)")
            prop_name=globals()['prop_calc_wg_'+uid].value
            prop_calc_cmd="Descriptors."+ prop_name + "(mol)"
            prop_str = prop_name + ' : ' + str(eval(prop_calc_cmd))
            # Update viewer
            globals()['prop_calc_view_'+uid].value = prop_str
        except:
            pass
    elif len(molNameList)>1:
        try:
            globals()['prop_precalc_view_'+uid].value = 'single molecule selection required!'
            globals()['prop_calc_view_'+uid].value = 'single molecule selection required!'
        except:
            pass

    
    try:
        drawAs=globals()['drawAs_wg_'+uid].value
    except:
        drawAs = globals()['drawAs_no_wg_'+uid]
        
    try:
        color=globals()['colorScheme_'+uid].value
    except:
        color = 'default'
    
    
    if drawAs is 'surface':
        globals()['view_'+uid].addSurface('SES', {});
    elif drawAs is 'ballstick':
        globals()['view_'+uid].setStyle({},{'stick':{'radius':'0.2','colorscheme': color},
                                            'sphere':{'radius':'0.4', 'colorscheme': color}});
    else:
        globals()['view_'+uid].setStyle({},{drawAs:{'colorscheme': color}})
    
    
    if len(molNameList)==1 and len(confIdsList)==1:
        mol = globals()['moldict_'+uid][molDictKey]
        confId=globals()['confId_'+uid].value
        sconf = mol.GetConformer(confId)
        xyz = sconf.GetPositions()
        OwningMol = sconf.GetOwningMol()
        
        try:
            labelConf=globals()['confLabel_'+uid].value
            labelAtom=globals()['atomLabel_'+uid].value
            if labelConf is True:
                label = molDictKey + ':' + str(confId)
                globals()['view_'+uid].addLabel(label, {'backgroundColor':'gray', 'fontColor':'white',
                                                        'showBackground':'true', 'alignment':'bottomCenter'})
            
            if labelAtom is True:
                label_create=[OwningMol.GetAtomWithIdx(i).GetSymbol()+
                              str(OwningMol.GetAtomWithIdx(i).GetIdx()+1) 
                              for i in range(sconf.GetNumAtoms())
                             ]
                i = None
                for i in range(sconf.GetNumAtoms()):
                    globals()['view_'+uid].addLabel(label_create[i], {'inFront' : 'false', 
                                                    'fontSize' : '12',
                                                    'fontColor':'gray',
                                                    'showBackground':'false',
                                                    'position' : {'x' : xyz[i][0],
                                                                  'y' : xyz[i][1],
                                                                  'z' : xyz[i][2]
                                                               }
                                                   })
        except:
            pass
    
    try:
        bg_color_selected = globals()['background_'+uid].value
        globals()['view_'+uid].setBackgroundColor(bg_color_selected)
    except:
        pass
    
    # To do : Add protein in viewer
    try:
        pVisibility = globals()['proteinVisible_'+uid].value
        
        if pVisibility is True:
            try:
                pStyle=globals()['pStyle_wg_'+uid].value
            except:
                pStyle = globals()['pStyle_no_wg_'+uid]

            if 'protein_'+uid in globals():
                pdb = Chem.MolToPDBBlock(globals()['protein_'+uid])
                globals()['view_'+uid].addModel(pdb,'pdb')
                
                if pStyle is 'surface':
                    globals()['view_'+uid].addSurface('SES', 
                                                      {'model':LigNumModel+1});
                elif pStyle is 'ballstick':
                    globals()['view_'+uid].setStyle({'model':LigNumModel+1}, 
                                                    {'stick':{'radius':'0.2','colorscheme': color}, 
                                                     'sphere':{'radius':'0.4', 'colorscheme': color}
                                                    }
                                                   );
                else:
                    try:
                        helicesAsTubes = globals()['pStyle_tube_wg_'+uid].value
                        if helicesAsTubes is True:
                            globals()['view_'+uid].setStyle({'model':LigNumModel+1},{pStyle:{'color': 'spectrum', 
                                                                                     'style': 'rectangle',
                                                                                     'arrows': 'true',
                                                                                     'tubes' : 'true'}})
                        else:
                            globals()['view_'+uid].setStyle({'model':LigNumModel+1},{pStyle:{'color': 'spectrum',
                                                                                     'arrows': 'true'}})
                    except:
                        globals()['view_'+uid].setStyle({'model':LigNumModel+1},{pStyle:{'color': 'spectrum',
                                                                                     'arrows': 'true'}})
                    
                    
    except:
        pass
    
    # zoomTo does not work well for surface and label... so, zoomTo should not be default settings
    #view.zoomTo()
    display(globals()['view_'+uid].update())
    
    
def handle_change(change):
    """This function handles all the interactive widgets except buttons and 3Dmol.js viewer"""
    update3D(change.owner._model_id)
    
def handle_start_button(b):
    """This function handles start button."""
    b.icon='check'
    b.description="Done!"
    update3D(b._model_id)
    
def handle_zoomTo_button(b):
    """This function handles zoomTo button"""
    uid=globals()['rdkit_wg_dict'][b._model_id]
    globals()['view_'+uid].zoomTo()
    display(globals()['view_'+uid].update())
    
    
def ChangeActiveLigand(uid, molId, confId, keepExistingView = False):
    """This function handles ligand select through python code"""
    if type(uid) != 'str':
        uid = str(uid)
    if type(molId) != 'str':
        molId = str(molId)
        
    globals()['molId_'+uid].value = molId
    globals()['confId_'+uid].value = confId
    
    if keepExistingView is False:
        globals()['selectMultiMols_'+uid].value=False
        globals()['selectAllMols_'+uid].value=False
        globals()['selectMultiConfs_'+uid].value=False
        globals()['selectAllConfs_'+uid].value=False
    
    update3D(globals()['molId_'+uid]._model_id)
    
    
    
def ShowConformers3D(Mol = None, protein = None,
                     useDrawAs = False, 
                     drawAs=None, pStyle=None,
                     propPanel = False, 
                     colorPanel = False, 
                     labelPanel = False):
    
    """This function initiates required widgets and 3Dmol.js viewer"""
    
    # Common suffix for widgets
    uid=str(time.time()).replace('.','')
    # print uid
    # Required global objects
    globals()['rdkit_mol_selected_'+uid] = list()
    globals()['rdkit_conf_selected_'+uid] = list()
    
    if 'rdkit_wg_dict' not in globals():
        globals()['rdkit_wg_dict'] = dict()
    
    
    # Right hand panel (widgets)
    
    itemLayout=Layout(display='flex', flex_flow='row', justify_content='space-between')
    
    wgListBox=list()
    
    wgListBox.append(Box([Label(value='uid'), HTML(description='', value=str(uid))], layout=itemLayout))
    
    # To do : Add protein
    if protein is not None:
        #print 'adding protein'
        globals()['protein_'+uid] = protein
    
    
    # molecules and conformers
    globals()['moldict_'+uid] = ProcessMolContainingObj(Mol)
    keys=list(globals()['moldict_'+uid].keys())
    
    globals()['selected_mols_view_'+uid] = HTML(description='', value=keys[0])
    globals()['selected_confs_view_'+uid] = HTML(description='', value='0')
    globals()['molId_'+uid] = Dropdown(description='', options=keys,value=keys[0])
    globals()['selectMultiMols_'+uid] = Checkbox(description='selectMultiMols', value=False)
    globals()['selectAllMols_'+uid] = Checkbox(description='selectAllMols', value=False)
    globals()['confId_'+uid] = Dropdown(description='', options=range(9), value=0)
    globals()['selectMultiConfs_'+uid] = Checkbox(description='selectMultiConfs', value=False)
    globals()['selectAllConfs_'+uid] = Checkbox(description='selectAllConfs', value=False)
    
    
    wgListBox.append(Box([Label(value='Mols'),globals()['selected_mols_view_'+uid]], layout=itemLayout))
    wgListBox.append(Box([Label(value='Confs'),globals()['selected_confs_view_'+uid]], layout=itemLayout))
    
    wgListBox.append(Box([Label(value='molId'),globals()['molId_'+uid]], layout=itemLayout))
    
    cbMolSelect=Box([Label(value=''),
                     HBox([globals()['selectMultiMols_'+uid], 
                           globals()['selectAllMols_'+uid]])
                    ], layout=itemLayout)
    
    wgListBox.append(cbMolSelect)
    
    wgListBox.append(Box([Label(value='confId'),globals()['confId_'+uid]], layout=itemLayout))
    
    cbConfSelect=Box([Label(value=''),
                      HBox([globals()['selectMultiConfs_'+uid], 
                            globals()['selectAllConfs_'+uid]])
                     ], layout=itemLayout)
    
    wgListBox.append(cbConfSelect)
    
    if protein is not None:
        globals()['proteinVisible_'+uid] = Checkbox(description='proteinVisible', value=True)
        wgListBox.append(Box([Label(value=''),globals()['proteinVisible_'+uid]], layout=itemLayout))
        globals()['rdkit_wg_dict'].update({globals()['proteinVisible_'+uid]._model_id:uid})
        globals()['proteinVisible_'+uid].observe(handle_change, names='value')
    
    
    globals()['rdkit_wg_dict'].update({globals()['molId_'+uid]._model_id:uid})
    globals()['rdkit_wg_dict'].update({globals()['selectMultiMols_'+uid]._model_id:uid})
    globals()['rdkit_wg_dict'].update({globals()['selectAllMols_'+uid]._model_id:uid})
    globals()['rdkit_wg_dict'].update({globals()['confId_'+uid]._model_id:uid})
    globals()['rdkit_wg_dict'].update({globals()['selectMultiConfs_'+uid]._model_id:uid})
    globals()['rdkit_wg_dict'].update({globals()['selectAllConfs_'+uid]._model_id:uid})
    
    
    globals()['molId_'+uid].observe(handle_change, names='value')
    globals()['selectMultiMols_'+uid].observe(handle_change, names='value')
    globals()['selectAllMols_'+uid].observe(handle_change, names='value')
    globals()['confId_'+uid].observe(handle_change, names='value')
    globals()['selectMultiConfs_'+uid].observe(handle_change, names='value')
    globals()['selectAllConfs_'+uid].observe(handle_change, names='value')
    
    
    # prop
    if propPanel is True:
        globals()['prop_precalc_view_'+uid] = HTML(description='', value='initializing...')
        globals()['prop_precalc_wg_'+uid] = Dropdown(description='',options=['select'], value='select')
        globals()['prop_calc_view_'+uid] = HTML(description='', value='initializing...')
        globals()['prop_calc_wg_'+uid] = Dropdown(description='',options=globals()['prop_rdkit'],value='MolLogP')
        wgListBox.append(Box([Label(value='calc'),globals()['prop_calc_view_'+uid]], layout=itemLayout))
        wgListBox.append(Box([Label(value='calc'),globals()['prop_calc_wg_'+uid]], layout=itemLayout))
        wgListBox.append(Box([Label(value='precalc'),globals()['prop_precalc_view_'+uid]], layout=itemLayout))
        wgListBox.append(Box([Label(value='precalc'),globals()['prop_precalc_wg_'+uid]], layout=itemLayout))
        globals()['rdkit_wg_dict'].update({globals()['prop_calc_wg_'+uid]._model_id:uid})
        globals()['rdkit_wg_dict'].update({globals()['prop_precalc_wg_'+uid]._model_id:uid})
        globals()['prop_calc_wg_'+uid].observe(handle_change, names='value')
        globals()['prop_precalc_wg_'+uid].observe(handle_change, names='value')
        
        
    # drawAs and useDrawAs
    if drawAs is None:
        drawAs = 'stick'
    
    globals()['drawAs_no_wg_'+uid] = drawAs
    
    if pStyle is None:
        pStyle = 'cartoon'
    
    globals()['pStyle_no_wg_'+uid] = pStyle
    
    if useDrawAs is True:
        globals()['drawAs_wg_'+uid] = Dropdown(description='', options=drawing_types_3d, value=drawAs)
        wgListBox.append(Box([Label(value='drawAs'),globals()['drawAs_wg_'+uid]], layout=itemLayout))
        globals()['rdkit_wg_dict'].update({globals()['drawAs_wg_'+uid]._model_id:uid})
        globals()['drawAs_wg_'+uid].observe(handle_change, names='value')
        
        globals()['pStyle_wg_'+uid] = Dropdown(description='', options=drawing_types_3d, value=pStyle)
        wgListBox.append(Box([Label(value='pStyle'),globals()['pStyle_wg_'+uid]], layout=itemLayout))
        globals()['rdkit_wg_dict'].update({globals()['pStyle_wg_'+uid]._model_id:uid})
        globals()['pStyle_wg_'+uid].observe(handle_change, names='value')
    
        globals()['pStyle_tube_wg_'+uid] = Checkbox(description='helicesAsTubes', value=False)
        wgListBox.append(Box([Label(value=''),globals()['pStyle_tube_wg_'+uid]], layout=itemLayout))
        globals()['rdkit_wg_dict'].update({globals()['pStyle_tube_wg_'+uid]._model_id:uid})
        globals()['pStyle_tube_wg_'+uid].observe(handle_change, names='value')
    
    
    # colors 
    if colorPanel is True:
        globals()['colorScheme_'+uid] = Dropdown(description='', options=color_scheme_3d, value='default')
        globals()['background_'+uid] = Dropdown(description='', options= bgcolors_3d, value=bgcolors_3d[1])
        wgListBox.append(Box([Label(value='color'),globals()['colorScheme_'+uid]], layout=itemLayout))
        wgListBox.append(Box([Label(value='background'),globals()['background_'+uid]], layout=itemLayout))
        globals()['rdkit_wg_dict'].update({globals()['colorScheme_'+uid]._model_id:uid})
        globals()['rdkit_wg_dict'].update({globals()['background_'+uid]._model_id:uid})
        globals()['colorScheme_'+uid].observe(handle_change, names='value')
        globals()['background_'+uid].observe(handle_change, names='value')
    
    
    # labels
    if labelPanel is True:
        globals()['confLabel_'+uid] = Checkbox(description='confLabel', value=False)
        globals()['atomLabel_'+uid] = Checkbox(description='atomLabel', value=False)
        
        cbLabel=Box([Label(value=''),
                     HBox([globals()['confLabel_'+uid], 
                           globals()['atomLabel_'+uid]])
                    ], layout=itemLayout)
        
        wgListBox.append(cbLabel)
        
        globals()['rdkit_wg_dict'].update({globals()['confLabel_'+uid]._model_id:uid})
        globals()['rdkit_wg_dict'].update({globals()['atomLabel_'+uid]._model_id:uid})
        globals()['confLabel_'+uid].observe(handle_change, names='value')
        globals()['atomLabel_'+uid].observe(handle_change, names='value')
    
    
    # buttons
    globals()['start_'+uid] = Button(description="Start!", button_style='success')
    globals()['zoomTo_'+uid] = Button(description="zoomTo", button_style='success')
    buttons=Box([Label(value=''),HBox([globals()['zoomTo_'+uid], globals()['start_'+uid]])],layout=itemLayout)
    globals()['rdkit_wg_dict'].update({globals()['start_'+uid]._model_id:uid})
    globals()['rdkit_wg_dict'].update({globals()['zoomTo_'+uid]._model_id:uid})
    globals()['start_'+uid].on_click(handle_start_button)
    globals()['zoomTo_'+uid].on_click(handle_zoomTo_button)
    
    wgListBox.append(buttons)
    
    
    
    # left panel (container holding table for viewer)
    size = (435, 485)
    viewerLayout=Layout(width=str(size[0]+5)+'px', height=str(size[1]+5)+'px', border='solid')
    wgLeftWidget = HTML('''<table><tr><td id="%s"></td></tr></table>'''%uid, layout=viewerLayout)
    
    # right panel
    wgRightBox=VBox(wgListBox, layout=Layout(border='solid', width='45%'))
    
    # Combining left and right panel
    all_wg = HBox([wgLeftWidget, wgRightBox])
    
    # displaying everything
    display(all_wg)
    
    
    # inserting 3DMol.js viewer in existing container (table)
    globals()['view_'+uid] = py3Dmol.view(width=size[0],height=size[1])
    globals()['view_'+uid].setBackgroundColor('0x000000')
    globals()['view_'+uid].zoomTo()
    
    display(globals()['view_'+uid].insert(uid))
    