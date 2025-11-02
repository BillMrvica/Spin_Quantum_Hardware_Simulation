# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 11:51:56 2022

@author: Francesco
"""

from core_tools.data.ds.reader import load_by_uuid, load_by_id, set_data_location
from core_tools.data.ds.ds_hdf5 import save_hdf5_uuid, save_hdf5_id
from core_tools.data.SQL.queries.dataset_gui_queries import query_for_measurement_results

from core_tools.data.SQL.connect import set_up_remote_storage
#%%



#%%
path_datasets = r'C:\Code\spin-projects\stations\vaughan\array_343\figures\hopping_spins_paper\datasets'
set_data_location(None) #loading from the database
#%%
set_data_location(path_datasets) #loading files from this folder

#%% Fig 1
uuids = [1693861160775283691,
         1692171914624283691 ]
for uuid in uuids:
    ds = load_by_uuid(uuid)
    save_hdf5_uuid(ds, path_datasets)
    my_id = ds.exp_id
    #ds = load_by_id(my_id)
    #save_hdf5_id(ds, path_datasets)


#%% Suppl Fig. 1
uuids = [1693557950033283691,
         1692020384925283691,
         1692282513469283691,
         1693415065370283691,
         1692012726066283691,
         1692171914624283691,
         1692729232124283691,
         1692013520094283691, 
         1693834872772283691, 
         1692628795405283691]

for uuid in uuids:
    ds = load_by_uuid(uuid)
    save_hdf5_uuid(ds, path_datasets)
    my_id = ds.exp_id
    
#%% Suppl. Fig. 2
uuids =   [1693558701888283691,
    1693855737612283691,
    1689846576231283691,
    1693398973179283691,
    1693847413281283691,
    1693861160775283691,
    1689780917716283691,
    1693841659702283691,
    1693836731202283691,
    1693868872715283691]

for uuid in uuids:
    ds = load_by_uuid(uuid)
    save_hdf5_uuid(ds, path_datasets)
    my_id = ds.exp_id
    
#%% Suppl. Fig. 4
uuids = [1690356693883283691, 
        1690356720952283691,
        1690356775835283691,
        1690356856203283691,
        1690356883143283691,
        1690356936770283691,
        1690356964739283691,
        1690356991668283691,
        1690357018477283691,
        1690357099743283691,
        1690357126532283691]

for uuid in uuids:
    ds = load_by_uuid(uuid)
    save_hdf5_uuid(ds, path_datasets)
    my_id = ds.exp_id
#%% Suppl. Fig. 5
uuids = [1690357153432283691]
         
for uuid in uuids:
    ds = load_by_uuid(uuid)
    save_hdf5_uuid(ds, path_datasets)
    my_id = ds.exp_id
#%% Suppl. Fig. 6

uuids = [1691762567742283691,
        1691771577064283691,
        1690556783281283691,
        1689772467030283691]

for uuid in uuids:
    ds = load_by_uuid(uuid)
    save_hdf5_uuid(ds, path_datasets)
    my_id = ds.exp_id
    
#%% Suppl. Fig. 8
uuids =  [1690976317716283691,
    1690977791498283691,
    1692778757618283691,
    1692781023443283691]

for uuid in uuids:
    ds = load_by_uuid(uuid)
    save_hdf5_uuid(ds, path_datasets)
    my_id = ds.exp_id
    
#%% utils, magnetic field calibration
uuids =  [1716540442970283691,
     1716541446593283691,
    1716542770826283691]

for uuid in uuids:
    ds = load_by_uuid(uuid)
    save_hdf5_uuid(ds, path_datasets)
    my_id = ds.exp_id