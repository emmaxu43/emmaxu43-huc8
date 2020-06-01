import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import csv
import fiona
import numpy as np
import rasterio
import copy
import time
import pickle
import sklearn.tree
import sklearn.ensemble
import sklearn.model_selection
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score


def read_nasis(path):

  data= {'name':[],'lat':[],'lon':[]}   
  fp= csv.reader(open(path,'r'), delimiter=',') 
  count= -1
  point={}
  
  for array in fp:
      count+= 1
      if count== 0:continue
      name= array[2]
      lat= array[6]
      lon= array[7]
      pedon_idx= array[0]
      if len(lat)== 0:continue
      if len(lon)== 0:continue
      #Extract pedon info
      data['name'].append(array[2])
      data['lat'].append(float(array[6]))
      data['lon'].append(float(array[7]))
      
  for var in data:
      data[var]= np.array(data[var])    
  
  data['ilat']=[]
  data['ilon']=[]
  for lat in data['lat']:
      if lat>= minlat and lat<= maxlat: data['ilat'].append(int(np.argmin(np.abs(latss- lat))))
      else: data['ilat'].append(-1)
      
  for lon in data['lon']:
      if lon>= minlon and lon<= maxlon: data['ilon'].append(int(np.argmin(np.abs(lonss- lon))))
      else: data['ilon'].append(-1)    
      
  print(len(data['lat']),' nasis points globally')
      
  return data

def read_tif(path):

  dataset= rasterio.open(path)
  data= dataset.read(1)
  dataset.close()
  
  return data
   
def csv_to_dict(input_data):
  
  new_dict={}
  count=0
  for i in input_data[0]:
      tmp= map(eval, input_data[1][count][1:-1].split("delimiter"))
      new_dict[i]=[]
      for j in tmp:
          for k in j:
              new_dict[i].append(k)
      count+=1
    
  return new_dict
  
def basin_idx_to_huc8ID(input_data):

  basin_dict={}
  tmp={}
  list1= list(final_dict.keys())
  for i in basin_idx:
      basin_dict[i]=[]
      tmp[i]=[]
      for j in neighbor_basin[i]:
        if j in list1: 
          list2= final_dict[j]
          for k in list2:
            tmp[i].append(input_data[10][k])
        else: tmp[i].append(input_data[10][j])
          
      basin_dict[i]= np.unique(np.asarray(tmp[i]))
      # print(i, len(basin_dict[i]),'(sub-basins in this cluster)')
      
  return basin_dict
  
def soil_id_globe(input_data):

  pedon_name=[]
  for num in input_data['name']:
      pedon_name.append(num)
  taxonomy_name= np.asarray(pedon_name) 
  
  tmp_list=[]
  for i in taxonomy_name:
      tmp= bytes(i, encoding="ascii")
      tmp_list.append(tmp)
  
  new_id1=[]
  idx_id1=[]
  j=0
  for i in tmp_list:
      if i in mapping[1].keys(): 
          new_id1.append(mapping[1][i])
          idx_id1.append(j)
      j+=1
  # print('points number in look up table1 = ',len(idx_id1))
  
  new_id2=[]
  idx_id2=[]
  for i in range(len(idx_id1)):
      if new_id1[i] in mapping[2].keys(): 
          new_id2.append(new_id1[i])  
          idx_id2.append(idx_id1[i])
  # print('points number in look up table2 = ',len(idx_id2))
  
  return idx_id2
  
def soil_id_conus(input_data, global_idx):

  #global coors in mapping dicts
  soil_lats= input_data['lat'][global_idx]
  soil_lons= input_data['lon'][global_idx]
  #us coors
  new_lats=[]
  new_lons=[]
  spatial_idx=[]
  for i in range(soil_lats.size):
      if soil_lats[i]> minlat and soil_lats[i]< maxlat and soil_lons[i]> minlon and soil_lons[i]< maxlon: 
          new_lons.append(soil_lons[i])
          new_lats.append(soil_lats[i])
          spatial_idx.append(global_idx[i])
  print(len(spatial_idx), 'NASIS points physically existed in CONUS')
  
  return spatial_idx
  
def taxo_name_conus(input_data, spatial_idx):

  pedon_name=[]
  for num in input_data['name']:
      pedon_name.append(num)
  taxonomy_name= np.asarray(pedon_name) 
  
  tmp_name= taxonomy_name[spatial_idx]
  taxo_name=[]
  for i in tmp_name:
      name= i.encode("utf-8")
      taxo_name.append(name)
      
  return taxo_name
  
def nasis_WI(data):

  n=0
  WI_idx=[]
  for i,j in zip(data['lat'], data['lon']):
      if i>= WI_minlat and i<= WI_maxlat and j>= WI_minlon and j<= WI_maxlon: WI_idx.append(n)
      n+=1
  #print(len(WI_idx),'NASIS points in WI (including non-physically exsited ones)\n', 'max idx: ',max(WI_idx))  #pedon csvtable idx
  
  return WI_idx
  
def rasterized_huc8(path):

  with open(path, 'rb') as handle:
      huc8_dict= pickle.load(handle)
  # print(len(huc8_dict)) #keys: lat lon; values: huc8ID
  
  return huc8_dict
  
def WI_soil(WI_idx, spatial_idx, unmask_data,data):

  ## Get grids in WI
  tmp_list= list(set(WI_idx).intersection(set(spatial_idx)))
  WI_ilat= np.asarray(data['ilat'])[tmp_list]
  WI_ilon= np.asarray(data['ilon'])[tmp_list]
  ## Get huc8ID
  wi_huc8_list= unmask_data.flatten()
  ## Get physically existed nasis poitns in WI
  final_idx={}   
  for i in range(len(WI_ilat)):
      ilat= WI_ilat[i]
      ilon= WI_ilon[i]
      final_idx[tmp_list[i]]= unmask_data.shape[1]*ilat+ilon
  print(len(final_idx), 'NASIS physically existed in WI')
  ## Get soil class ID in WI
  soil_dict={}
  for i in hrchy:
      soil_dict[i]={}  
      #keys: pedon csvtable idx; values: soil order/suborder/grtgrp/subgrp (nasis)
  n=0
  for i in spatial_idx:
      soil_dict['ord'][i]= newID['ord'][n]
      soil_dict['subord'][i]= newID['subord'][n]
      soil_dict['grtgrp'][i]= newID['grtgrp'][n]
      soil_dict['subgrp'][i]= newID['subgrp'][n]
      n+=1
  
  # Prepare soil info in WI
  WI_huc8={}
  n=0
  for i in final_WIlist:
      WI_huc8[i]=[]
      lat= np.asarray(data['lat'])[i]
      lon= np.asarray(data['lon'])[i]
      ilat= WI_ilat[n]
      ilon= WI_ilon[n]
      huc8_i= wi_huc8_list[final_idx[i]]
      
      WI_huc8[i]= np.array([huc8_i, lon, lat, soil_dict['ord'][i],soil_dict['subord'][i],\
                            soil_dict['grtgrp'][i], soil_dict['subgrp'][i], ilon, ilat])  
  
      n+=1
  #print(len(WI_huc8))
  
  return WI_huc8
  
def WI_soil_ID(WI_huc8):
  
  soil_list={}
  tmp_arr={}
  for i in hrchy:
      soil_list[i]=[]
      tmp_arr[i]=[]
      
  for i in final_WIlist:
      soil_list['ord'].append(WI_huc8[i][3])
      soil_list['subord'].append(WI_huc8[i][4])
      soil_list['grtgrp'].append(WI_huc8[i][5])
      soil_list['subgrp'].append(WI_huc8[i][6])
      
  for i in hrchy:    
    tmp_arr[i]= np.asarray(soil_list[i])

  return tmp_arr
  
def read_name_1(file_path,soil_class):

  with open(file_path, 'r') as f:
      x = f.read().split('\n')
      
  name[soil_class]=[]
  for i in x[:-1]:
      j= i.split(' ')
      name[soil_class].append(j[1])
      
  return name[soil_class]
  
def read_name_2(file_path,soil_class):

  with open(file_path, 'r') as f:
      x = f.read().split('\n')
      
  name[soil_class]=[]
  for i in x[:-1]:
      j= i.split(' ')
      if len(j)==3: name[soil_class].append(j[1]+' '+j[2])
      elif len(j)==4: name[soil_class].append(j[1]+' '+j[2]+' '+j[3])
      else: name[soil_class].append(j[1]+' '+j[2]+' '+j[3]+' '+j[4])
      
  return name[soil_class]
  
def nasis_cluster(WI_huc8, basin_idx, basin_dict):

  tmp_list=[]
  for i in list(WI_huc8.values()):
      tmp_list.append(i[0])

  assign={}
  tmp={}
  for i in basin_idx:
      tmp[i]=[]
      for j in basin_dict[i]: #basin_dict: huc8ID for each basin cluster
          idx= np.where(j== tmp_list)[0]
          for k in idx:
              tmp[i].append(k)  #idx for 26941 nasis points in wi (n= 26941)
      assign[i]= np.unique(np.asarray(tmp[i]))
      # print(i, assign[i].shape, 'nasis points in this cluster')
      
  return assign
  
def anchor_env_covars(basin_idx, unmask_data, sub_assign):

  env_anchor={}
  for i in basin_idx:
      env_anchor[i]=[]
      tmp= np.argwhere(unmask_data==sub_assign[i])
      for j in np.squeeze(np.asarray(tmp)):
          env_anchor[i].append(j)
      # print(i, len(env_anchor[i]), 'env covars in this anchor')
      
  return env_anchor

def nasis_cluster_grids(WI_idx, spatial_idx, data):

  tmp_list= list(set(WI_idx).intersection(set(spatial_idx)))
  WI_ilat= np.asarray(data['ilat'])[tmp_list]
  WI_ilon= np.asarray(data['ilon'])[tmp_list]
  
  nasis_basin_ilat={}
  nasis_basin_ilon={}
  for i in basin_idx:
      nasis_basin_ilat[i]=[]
      nasis_basin_ilon[i]=[]
      for j in assign[i]:
          nasis_basin_ilat[i].append(WI_ilat[j])
          nasis_basin_ilon[i].append(WI_ilon[j])
      
      # print(i, ": ",len(nasis_basin_ilat[i]),'nasis points in this cluster')  
      
  return nasis_basin_ilat, nasis_basin_ilon
  
def prepare_cluster_train_data(basin_idx, file_path, env_vars, nasis_basin_ilat, nasis_basin_ilon):

  input_data={}
  output={}
  datasets={}
  for k in basin_idx:
      output[k]={}
      for var in env_vars:
          output[k][var]=[]
          datasets[var]=[]
          datasets[var]= rasterio.open(file_path % var)
          input_data[var]= datasets[var].read(1)     
          for i in zip(nasis_basin_ilat[k],nasis_basin_ilon[k]): #row, col
              output[k][var].append(input_data[var][i])
          datasets[var].close()

  x={}
  y={}
  xv={}
  yv={}
  lat={}
  lon={}
  for i in basin_idx:
      x[i]= np.linspace(minlon,maxlon,input_data[env_vars[0]].shape[1])  #lon
      y[i]= np.linspace(maxlat,minlat,input_data[env_vars[0]].shape[0])  #lat        
      xv[i],yv[i]= np.meshgrid(x[i],y[i])
      lat[i]= np.asarray(yv[i])
      lon[i]= np.asarray(xv[i])
      
  for k in basin_idx:
      output[k]['lat']=[]
      output[k]['lon']=[]
      for i in zip(nasis_basin_ilat[k],nasis_basin_ilon[k]):
          output[k]['lat'].append(lat[k][i])
          output[k]['lon'].append(lon[k][i])
          
  return output
  
def reshape_train_X(basin_idx, new_vars,output):

  tmp_X={}
  X={}
  # print('Basin cluster (nasis points)')
  for k in basin_idx:
      tmp_X[k]=[]
      X[k]=[]
      for var in new_vars:
          tmp_X[k].append(output[k][var])
      X[k]= np.asarray(tmp_X[k]).reshape((len(new_vars),len(output[k][var])))
      # print(k,": ",X[k].shape)
      
  return X
  
def reshapa_train_y(tmp_arr, assign, basin_idx):

  basin_ID={}
  for i in hrchy:
      basin_ID[i]={}
      
  for i in basin_idx:
      basin_ID['ord'][i]=[]
      basin_ID['subord'][i]=[]
      basin_ID['grtgrp'][i]=[]
      basin_ID['subgrp'][i]=[]
      
      for j in assign[i]:
          basin_ID['ord'][i].append(tmp_arr['ord'][j])
          basin_ID['subord'][i].append(tmp_arr['subord'][j])
          basin_ID['grtgrp'][i].append(tmp_arr['grtgrp'][j])
          basin_ID['subgrp'][i].append(tmp_arr['subgrp'][j])
      # print(i, len(basin_ID['subgrp'][i]), 'NASIS in this cluster')
      
  return basin_ID

def RF(bain_ID, X):

  Rscore={}
  clf={}
  Xt={}
  Xv={}
  yt={}
  yv={}
  y={}
  scores={}
  # print('basin cluster')
  for i in hrchy:
      Rscore[i]={}
      clf[i]={}
      Xt[i]={}
      Xv[i]={}
      yt[i]={}
      yv[i]={}
      scores[i]=[]
      for j in basin_idx:
          y[j]= np.asarray(basin_ID[i][j])
          Xt[i][j], Xv[i][j], yt[i][j], yv[i][j] = train_test_split(X[j].T, y[j], test_size=0.3, random_state=42)
          clf[i][j] = sklearn.ensemble.RandomForestClassifier(n_estimators=100).fit(Xt[i][j],yt[i][j])
          Rscore[i][j]= clf[i][j].score(Xv[i][j], yv[i][j])
          scores[i].append(Rscore[i][j])
          # if i== hrchy[0]: print(X[j].shape, y[j].shape)
  
  return Rscore, clf, scores, yt, Xv
  
def anchor_train_data(basin_idx, file_path, env_vars, nasis_basin_ilat, nasis_basin_ilon, env_anchor, latss, lonss):

  input_data={}
  output={}
  datasets={}
  for k in basin_idx:
      output[k]={}
      for var in env_vars:
          output[k][var]=[]
          datasets[var]=[]
          datasets[var]= rasterio.open(file_path % var)
          input_data[var]= datasets[var].read(1)     
          for i in zip(nasis_basin_ilat[k],nasis_basin_ilon[k]): #row, col
              output[k][var].append(input_data[var][i])
          datasets[var].close()

  env_anchor_output={}
  for i in basin_idx:
    env_anchor_output[i]={}
    for var in env_vars:
      env_anchor_output[i][var]=[]
      for j in env_anchor[i]:  # ilat & ilon in each anchor basin
          env_anchor_output[i][var].append(input_data[var][j[0],j[1]])
      
  env_anchor_lat={}
  env_anchor_lon={}
  for i in basin_idx:
      env_anchor_lat[i]=[]
      env_anchor_lon[i]=[]
      for j in env_anchor[i]:
              env_anchor_lat[i].append(latss[j[0]])
              env_anchor_lon[i].append(lonss[j[1]])
  
  Xs={}
  env_anchor_Xs={}
  # print('Anchor basin (env covars)')
  for i in basin_idx:
      Xs[i]=[]
      env_anchor_Xs[i]=[]
      for var in env_vars:
          Xs[i].append(env_anchor_output[i][var])
      Xs[i].append(env_anchor_lat[i])
      Xs[i].append(env_anchor_lon[i])
      env_anchor_Xs[i]= np.asarray(Xs[i]).reshape(len(env_vars)+2,len(env_anchor_lat[i]))
      # print(i, env_anchor_Xs[i].shape,'env covars in this anchor basin')  
  
  return env_anchor_Xs, env_anchor_lat, env_anchor_lon
  
def sort(env_anchor_prob, unique):

  sort_prob={}
  sort_class={}
  tmp={}
  tmp_list={}
  for i in hrchy:
      sort_prob[i]={}
      sort_class[i]={}
      tmp[i]={}
      tmp_list[i]={}
      
      for j in basin_idx:
          sort_prob[i][j]= -np.sort(-env_anchor_prob[i][j]) #descending
          tmp[i][j]= np.argsort(env_anchor_prob[i][j], axis=1)
  
          tmp_array= (unique[i][j])[tmp[i][j]]
          tmp_list[i][j]=[]
          for k in tmp_array:
              tmp_list[i][j].append(k[::-1])
          sort_class[i][j]= np.asarray(tmp_list[i][j]).reshape(tmp[i][j].shape)
          
  return sort_prob, sort_class
  
def find_boundary(input_data):
  
  path= '/home/cx43/ramcharan/Jan2020/WBD/huc8_conus/HUC8_CONUS/HUC8_US.shp'
  fp= fiona.open(path, 'r')
  dataset= list(fp)  
  test={}
  vals={}
  MIN_LON={}
  MAX_LON={}
  MIN_LAT={}
  MAX_LAT={}
  for i in input_data:  
      test[i]= dataset[i]['geometry']
      vals[i]= test[i]['coordinates']
  
      tmp0=[]
      tmp1=[]
      for j in vals[i]:
          for k in j:
              if len(k)==2:
                  tmp0.append(k[0])  #lon
                  tmp1.append(k[1])  #lat
              else:
                  for k1 in k:
                      tmp0.append(k1[0])
                      tmp1.append(k1[1])
      m=0
      for item in tmp0:
          if type(item)!= float: tmp0[m]= item[0] # some coors format is tuple instead of float
          m+=1
  
      m=0
      for item in tmp1:
          if type(item)!= float: tmp1[m]= item[1]
          m+=1
  
      MAX_LON[i]= max(tmp0)
      MIN_LON[i]= min(tmp0)       
      MAX_LAT[i]= max(tmp1)
      MIN_LAT[i]= min(tmp1)

  return MAX_LAT,MIN_LAT,MAX_LON,MIN_LON

def env_grid(anchor_minlon, anchor_maxlon, anchor_minlat, anchor_maxlat, cols, rows, basin_idx):

  delta_lat= (maxlat-minlat)/cols
  delta_lon= (maxlon-minlon)/rows
  env_anchor_lats={}
  env_anchor_lons={}
  for i in basin_idx:
      env_anchor_lats[i]= np.linspace(anchor_minlat[i], anchor_maxlat[i], \
                                  int((anchor_maxlat[i]- anchor_minlat[i])/delta_lat)+1)
      env_anchor_lons[i]= np.linspace(anchor_minlon[i], anchor_maxlon[i], \
                                  int((anchor_maxlat[i]- anchor_minlat[i])/delta_lat)+1)
  anchor_x={}  # lon
  anchor_y={}  # lat
  for i in basin_idx:
      anchor_x[i], anchor_y[i]= np.meshgrid(env_anchor_lons[i],env_anchor_lats[i])
      
  return anchor_x, anchor_y, env_anchor_lats, env_anchor_lons
  
def anchor_info(env_anchor_lon, env_anchor_lat,anchor_x, basin_idx, unique, env_anchor_lons, env_anchor_lats, sort_prob, sort_class):

  anchor_prob={}
  anchor_class={}
  for i in hrchy:
      anchor_prob[i]={}
      anchor_class[i]={}
      for j in basin_idx:
          tmp_list= [-9999.0]* anchor_x[j].size * len(unique[i][j])
          anchor_prob[i][j]= np.asarray(tmp_list).reshape((len(unique[i][j]),env_anchor_lats[j].size, env_anchor_lons[j].size))
          anchor_class[i][j]= np.asarray(tmp_list).reshape((anchor_prob[i][j].shape))
          for k in range(len(unique[i][j])):
              # tmp= sort_prob[i][j][:,k]  # only select one soil class
              n=0
              for j1,j2 in zip(env_anchor_lon[j], env_anchor_lat[j]):  #nasis point data
                  k1= np.argmin(np.abs(j1-env_anchor_lons[j]))  # lons -> seamless domain in each anchor basin 
                  k2= np.argmin(np.abs(j2-env_anchor_lats[j]))  #ilat -> row -> y; ilon-> col -> x
                  anchor_prob[i][j][k,k2,k1]= sort_prob[i][j][:,k][n]  #z,y,x -> class, ilat, ilon
                  anchor_class[i][j][k,k2,k1]= sort_class[i][j][:,k][n] 
                  n+=1
                  
  return anchor_prob, anchor_class       
  


                         
# Read rasterized huc8 (CONUS map)
huc8_path2= '/home/cx43/ramcharan/Jan2020/WBD/huc8_conus/HUC8_CONUS/HUC8_US_test1.tif'
unmask_data= read_tif(huc8_path2)

# Read rasterized ssurgo (WI map)
ssurgo_WI_path= '/home/cx43/ramcharan/Jan2020/WBD/data_april/ssurgo_WI.tif'
mukey= read_tif(ssurgo_WI_path)

# Define CONUS domain
minlon= -124.7844079
minlat= 24.7391568
maxlon= -66.9538279
maxlat= 49.3457868
rows= unmask_data.shape[0]
cols= unmask_data.shape[1]
latss= np.linspace(maxlat,minlat,rows)
lonss= np.linspace(minlon,maxlon,cols)
lon,lat= np.meshgrid(lonss,latss)

# Read nasis information (table)
nasis_path= '/home/cx43/ramcharan/Jan2020/Pedons.csv'
data= read_nasis(nasis_path)

# Read cluster basin
basin_path= '/home/cx43/ramcharan/Jan2020/WBD/data_april/final_dict.csv'
tmp_dict= pd.read_csv(basin_path, header=None)
final_dict= csv_to_dict(tmp_dict)

# Find subbasins in each cluster
csv_data= pd.read_csv('/home/cx43/ramcharan/Jan2020/WBD/huc8_conus/CONUS.csv', header=None)
tmp_neighbor_basin=  pd.read_csv('/home/cx43/ramcharan/Jan2020/WBD/data_april/neighbor_basin.csv', header=None)
neighbor_basin= csv_to_dict(tmp_neighbor_basin)
basin_idx= list(neighbor_basin.keys())
basin_dict= basin_idx_to_huc8ID(csv_data)

# Read look up tables
mapping={}
file1= '/stor/soteria/hydro/private/nc153/projects/POLARIS/POLARIS/oct2017/validation/polaris_series/mapping.pck'
mapping[1]= pickle.load(open(file1,'rb'),encoding= 'latin1') 
paths= '/stor/soteria/hydro/private/nc153/projects/POLARIS/POLARIS/reclass/workspace/'
vars= [ 'series_to_soilseriesname_mapping', 'series_to_tax_order_mapping', 'series_to_subords_mapping', 'series_to_grtgrp_nms_mapping', 'series_to_subgrps_mapping']
for i in [2,3,4,5,6]:
  path= paths+ vars[i-2] + '.pck'
  mapping[i]= pickle.load(open(path,'rb'),encoding= 'latin1')

# Use taxonomy name to extract physically existed soil series id (globally)
idx_id= soil_id_globe(data)
## Down to the US
spatial_idx= soil_id_conus(data, idx_id)
## Physically exsited soil taxo in the CONUS
taxo_name= taxo_name_conus(data, spatial_idx)

# Exact ID for different soil classes
newID={}
hrchy=['ord','subord','grtgrp','subgrp']
for i in hrchy:
    newID[i]=[]

for i in taxo_name:
    tmp= mapping[1][i]
    newID['ord'].append(mapping[3][tmp])
    newID['subord'].append(mapping[4][tmp])
    newID['grtgrp'].append(mapping[5][tmp])
    newID['subgrp'].append(mapping[6][tmp])
print('Counting number of soil classes in CONUS (nasis)')
print('counting new orderID number= ', np.unique(np.asarray(newID['ord'])).size)
print('counting new suborderID number= ', np.unique(np.asarray(newID['subord'])).size)
print('counting new great group ID number= ', np.unique(np.asarray(newID['grtgrp'])).size)
print('counting new subgroupID number= ', np.unique(np.asarray(newID['subgrp'])).size)

# Define WI domain
WI_minlon= -94.72668390999999
WI_maxlon= -79.65111050199994
WI_minlat= 41.308488225000104
WI_maxlat= 49.022614884000106

# Physically exsited soil idx in WI
WI_idx= nasis_WI(data)

# Load rasterized huc8ID over the CONUS
path= '/home/cx43/ramcharan/Jan2020/WBD/huc8_conus/huc8ID.pickle'
huc8_dict= rasterized_huc8(path)

# Prepare spatial explicitly soil info in WI
## Get nasis idx in WI 
final_WIlist= list(set(WI_idx).intersection(set(spatial_idx)))
## dict keys: pedon csvtable idx; 
## dict values:huc8ID,lon,lat,orderID,subordID,grtgrpID,subgrpID,ilon(x),ilat(y)
WI_huc8= WI_soil(WI_idx, spatial_idx, unmask_data,data)

# Soil classes ID in WI
soil_arr= WI_soil_ID(WI_huc8)

# Get soil names for different taxonomies in WI
vars= ['tax_order_mapping.txt','subords_mapping.txt', 'grtgrp_nms_mapping.txt', 'subgrps_mapping.txt']
paths= '/stor/soteria/hydro/private/nc153/projects/POLARIS/POLARIS/reclass/workspace/'

name={}
count=0
for i in hrchy:
  if i== hrchy[0] or i== hrchy[1] or i==hrchy[2]: name[i]= read_name_1(paths+ vars[count],i)
  if i== hrchy[3]: name[i]= read_name_2(paths+ vars[count],i)
  count+=1
  
# Get nasis points for each basin cluster
##assign key: basin idx; 
##assign values: dict idx from WI_huc8 (26941 physically existed nasis points in WI)
assign= nasis_cluster(WI_huc8, basin_idx, basin_dict)  # cluster basin
anchor_assign= {}
for i in basin_idx:
    anchor_assign[i]= csv_data[10][i]
    
# Find ilat & ilon (conus grids) of env covars in each anchor basin
env_anchor= anchor_env_covars(basin_idx, unmask_data, anchor_assign)

# Find ilat & ilon (conus grids) of nasis points in each basin cluster
nasis_basin_ilat, nasis_basin_ilon= nasis_cluster_grids(WI_idx, spatial_idx, data)

# Generate RFs for basin cluster
## Prepare training datasets
env_vars= ['NLCD','NAMrad_K','NAMrad_U','NAMrad_Th','NED_feb26','usmag_hp500_cp','USgrv_cba_SDD_geog','USgrv_iso_SDD_geog']
file_path= "/home/cx43/ramcharan/Jan2020/data/feb_12/%s.tif"
output= prepare_cluster_train_data(basin_idx, file_path, env_vars, nasis_basin_ilat, nasis_basin_ilon)
## Reshape training datasets (X, cluster basin)
all_vars= env_vars+ ['lat','lon']
X= reshape_train_X(basin_idx, all_vars,output)
## Reshape trainning datasets (y, cluster basin)
basin_ID= reshapa_train_y(soil_arr, assign, basin_idx)
## Generate RFs
Rscore, clf, scores, yt, Xv= RF(basin_ID, X)

# Print original RF accuracy scores (before prune)
for i in hrchy:
  print(i, ' mean score: ', np.mean(scores[i]))
  for j in basin_idx:
    print(j, '   Accuracy= ',Rscore[i][j])
    
# Predict soil classes in each anchor basin
## Prepare training datasets (Xs, anchor basin)
env_anchor_Xs,env_anchor_lat, env_anchor_lon= anchor_train_data(basin_idx, file_path, env_vars, nasis_basin_ilat, nasis_basin_ilon, env_anchor, latss, lonss)
## Predict the most probable soil classes (anchor basin)
env_anchor_pre={}
for i in hrchy:
    env_anchor_pre[i]={}
    for j in basin_idx:
        env_anchor_pre[i][j]= clf[i][j].predict(np.asarray(env_anchor_Xs[j].T))
## Predict the probabilities for each class (anchor basin)
env_anchor_prob={}
unique={}
for i in hrchy:
    env_anchor_prob[i]={}
    unique[i]={}
    for j in basin_idx:
        env_anchor_prob[i][j]= clf[i][j].predict_proba(env_anchor_Xs[j].T)
        unique[i][j]= np.unique(yt[i][j])

# Sort soil classes & relevant probabilities
## probs are sorted from highest to lowest probabilities
## class ordered by probs
sort_prob, sort_class= sort(env_anchor_prob, unique) # anchor basin

# Define geographical boundary of each anchor basin
anchor_bound= find_boundary(basin_idx)
anchor_minlon= {}
anchor_maxlon= {}
anchor_minlat= {}
anchor_maxlat= {}
for i in basin_idx:
    anchor_minlon[i]= anchor_bound[3][i]
    anchor_maxlon[i]= anchor_bound[2][i]
    anchor_minlat[i]= anchor_bound[1][i]
    anchor_maxlat[i]= anchor_bound[0][i]
##Generate coors grids for each anchor basin
anchor_x, anchor_y, env_anchor_lats, env_anchor_lons= env_grid(anchor_minlon, anchor_maxlon, anchor_minlat, anchor_maxlat, cols, rows, basin_idx)

# Generate coors explicitly soil classes for each anchor basin
## z,y,x -> soil probs/ classes, lat, lon
anchor_prob, anchor_class= anchor_info(env_anchor_lon, env_anchor_lat,anchor_x, basin_idx, unique, env_anchor_lons, env_anchor_lats, sort_prob, sort_class)

# Integrating ssurgo
## keys: basin_idx; values: soil order, suborder, grtgrp, subgrp ID
with open('ssurgo_hierarchy.pickle', 'rb') as handle: 
    ssurgo_class= pickle.load(handle)
# Find soil classes from ssurgo in each cluster basin
ssurgo={}
n=0
# print('Counting number of soil classes from SSURGO')
for i in hrchy:
  ssurgo[i]= {}
  # print(i)
  for j in basin_idx:
    ssurgo[i][j]= np.unique(np.asarray(ssurgo_class[j][:,n]))
    # print(j, np.unique(ssurgo[i][j].size))
  n+=1
  
# Prune soil classes
pruned_y={}
for i in hrchy:
  pruned_y[i]={}
  for j in basin_idx:
    pruned_y[i][j]= clf[i][j].predict(Xv[i][j])
## Pruned probs and classes
pruned_prob={}
pruned_class={}
for i in hrchy:
  pruned_prob[i]={}
  pruned_class[i]={}
  for j in basin_idx:
    pruned_prob[i][j]= clf[i][j].predict_proba(Xv[i][j])
    pruned_class[i][j]= np.unique(np.asarray(yt[i][j]))

# Sort pruned arrays
sort_pruned_prob, sort_pruned_class= sort(pruned_prob, pruned_class)  # validated nasis points (cluster basin)

# Select ssurgo classes
intersect_class={}
intersect_pruned_prob={}
intersect_pruned_class={}
for i in hrchy:
  intersect_class[i]={}
  intersect_prob[i]={} # before normalizing
  intersect_pruned_class[i]={}
  for j in basin_idx:
     # soil class in both ssurgo and prediction from nasis
     intersect_class[i][j]= np.intersect1d(list(pruned_class[i][j]),list(ssurgo[i][j]))
     intersect_pruned_prob[i][j]= sort_pruned_prob[i][j][:, np.arange(intersect_class[i][j].size)]
     intersect_pruned_class[i][j]= sort_pruned_class[i][j][:, np.arange(intersect_class[i][j].size)]
     
# Accuracy scores after pruning with ssurgo
predict_class={}
pruned_scores={}
tmp={}
for i in hrchy:
  print(i)
  predict_class[i]={}
  pruned_scores[i]={}
  tmp[i]=[]
  for j in basin_idx:
    predict_class[i][j]= intersect_pruned_class[i][j][:,[0]]
    pruned_scores[i][j]= accuracy_score(yv[i][j], predict_class[i][j])
    tmp[i].append (pruned_scores[i][j])
    print(j, pruned_scores[i][j])
  print('mean score: ',np.mean(tmp[i]),'\n')
