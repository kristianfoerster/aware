#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 07:59:22 2017

@author: kristianf
"""

import rasterio
import numpy as np
import copy
import os
#import sys
#sys.path.append('/Users/kristianf/projekte/MUSICALS/AWARE/repository')
import aware
from aware.hydrographictree import HTree
#from hydrographictree import HTree

def write_gdal_file(filename, data, prj_settings):
    with rasterio.open(filename, 'w', **prj_settings) as ds:
        ds.write_band(1, data)

def get_dir(dx,dy):
    ''' Returns an ArcGIS compatible id of flow direction based on given 
    vectors in x and y direction'''
    if dx == 1 and dy == 0:
        return 1
    elif dx == 1 and dy == 1:
        return 2
    elif dx == 0 and dy == 1:
        return 4
    elif dx == -1 and dy == 1:
        return 8
    elif dx == -1 and dy == 0:
        return 16
    elif dx == -1 and dy == -1:
        return 32
    elif dx == 0 and dy == -1:
        return 64
    elif dx == 1 and dy == -1:
        return 128
    return 0

def get_vec(fd):
    '''Transforms ArcGIS flow direction ids to vectors. '''
    if fd == 1:
        return [1,0]
    elif fd == 2:
        return [1,1]
    elif fd == 4:
        return [0,1]
    elif fd == 8:
        return [-1,1]
    elif fd == 16:
        return [-1,0]
    elif fd == 32:
        return [-1,-1]
    elif fd == 64:
        return [0,-1]
    elif fd == 128:
        return [1,-1]
    return [0,0]

def compute_flow_dir(dtm):
    ''' Computes a flow direction map based on digital elevation map'''
    flowdir = np.zeros(dtm.shape)
    rows = dtm.shape[0]
    cols = dtm.shape[1]
    for i in range(1,rows-1):
        for j in range(1,cols-1):
            slope_max = 0
            direction = 0
            for dx in range(-1,2):
                for dy in range(-1,2):
                    if dx == 0 and dy == 0:
                        continue
                    dz = dtm[i,j] - dtm[i+dy,j+dx]
                    dl = 1
                    if abs(dx) == abs(dy):
                        dl = np.sqrt(2)
                    slope = dz / dl
                    if slope > slope_max:
                        slope_max = slope
                        direction = get_dir(dx,dy)
                flowdir[i,j] = direction
    return flowdir

def calc_flow_accumulation_single_iteration(flowdir, accu_grid, i, j):
    '''Calculates a single step of flow accumulation based on flow direction'''
    accu = 1
    rows = flowdir.shape[0]
    cols = flowdir.shape[1]
    for dx in range(-1,2):
        for dy in range(-1,2):
            if dx == 0 and dy == 0:
                continue
            if not (i+dy < 0 or j+dx < 0 or i+dy >= rows or j+dx >= cols):
                fd_next = flowdir[i+dy,j+dx]
                vec = get_vec(fd_next)
                if dx == -vec[0] and dy == -vec[1]:
                    accu+=calc_flow_accumulation_single_iteration(flowdir, accu_grid, i+dy,j+dx)
    accu_grid[i,j] = accu
    return accu

def calc_flow_accumulation(dtm, flowdir):
    '''Calculates the flow accumulation based on a flow direction map'''
    # create a copy of the dtm. this copy includes blank areas already processed
    dtm2 = copy.deepcopy(dtm)

    flow_accu2 = np.zeros(dtm.shape)
    flow_accu2[:] = np.nan
          
    loop = True
    li = 1
    last_min_i = -1
    last_min_j = -1
    max_iter = dtm2.size
    while loop:
        if li > max_iter:
            loop=False
        #print('Iteration %i' % li)
        li+=1
        pos_min = np.nanargmin(dtm2)
        #print(pos_min)
        pos_min2 = np.unravel_index(pos_min, dtm2.shape)  
        pos_min_i = pos_min2[0]
        pos_min_j = pos_min2[1]
        #print(pos_min_i,pos_min_j)
        if last_min_i == pos_min_i and last_min_j == pos_min_j:
            dtm2[pos_min_i,pos_min_j] = np.nan
            continue
    
        last_min_i = pos_min_i
        last_min_j = pos_min_j
    
        calc_flow_accumulation_single_iteration(flowdir,flow_accu2, pos_min_i, pos_min_j)
        n_nan = len(flow_accu2[np.isnan(flow_accu2)])
        n_calc = n_nan / flow_accu2.size
        dtm2[~np.isnan(flow_accu2)] = np.nan
        #print(n_calc)
        if n_calc < 0.001 or n_calc == 1.:
            loop = False

    print('# of iterations for flow accumulation = %i' % li)
    return flow_accu2

def get_upstream_area(flowdir,accu_grid,cid_grid,i,j,id, pplist, downstream_pp):

    copy_of_pplist = copy.deepcopy(pplist)
    
    rows = flowdir.shape[0]
    cols = flowdir.shape[1]

    #if i==90 and j==100:
    #    print('Test')

    accu = 1
    del_p = -1
    pp_added=False
    current_pp = downstream_pp
    updated_pp = downstream_pp
    for k,pi in enumerate(copy_of_pplist):
        if i == pi.x and j == pi.y:
            id += 1
            print('Create catchment for %s at i=%i, j=%i, set to id=%i' % (pi.name,i,j,id))
            
            if downstream_pp is not None:
                id = downstream_pp.add_new_tributary_id(id)
                pp_added=True
                print('node %i: new node added %i' % (downstream_pp.id, id))
            else:
                
                print('root node added %i' % id)
                
            # current_pp = HTree(id=id, name=pi.name)
            # now the new_instance method is capable of selecting different types of nodes
            current_pp = HTree.new_instance(new_id=id, new_name=pi.name, node_type_in=pi.type)
            updated_pp = current_pp
            
            del_p = k
            break
    if del_p >= 0:
        print(del_p)
        del copy_of_pplist[del_p]
    #outgrid[i,j] = id
    
    for dx in range(-1,2):
        for dy in range(-1,2):
            if dx == 0 and dy == 0:
                continue
            if not (i+dy < 0 or j+dx < 0 or i+dy >= rows or j+dx >= cols):
                fd_next = flowdir[i+dy,j+dx]
                vec = get_vec(fd_next)
                #print(vec)
                if dx == -vec[0] and dy == -vec[1]:
                    # go upstream!
                    accu_up_i, id, copy_of_pplist, updated_pp = get_upstream_area(flowdir,accu_grid,cid_grid,i+dy,j+dx,id, copy_of_pplist, current_pp)
                    accu += accu_up_i

    
    if pp_added:
        updated_pp.set_area(accu)
        grid_id = updated_pp.id
        downstream_pp.add_tributary(updated_pp)
        updated_pp = downstream_pp

    if updated_pp is not None:
        if not pp_added:
            grid_id = updated_pp.id
        cid_grid[i,j] = grid_id
    
    accu_grid[i,j] = accu    
    
    return accu, id, copy_of_pplist, updated_pp

def delineate_watersheds(dtm,pplist):
    id = 0
    treelist = []
    n_nodes = 0

    accu[:] = np.zeros(dtm.shape)
    accu[:] = np.nan

    cid = np.zeros(dtm.shape)
    cid[:] = np.nan
    
    ii = 0
    while len(pplist) > 0:
        p = pplist[ii]
        last_length = len(pplist)
        area,id, pplist, tree = get_upstream_area(flowdir,accu,cid,p.x,p.y,id, pplist, None)
        tree.area = area
        n_nodes_new, max_depth = tree.postprocess_tree(n_nodes)
        n_nodes += n_nodes_new
        tree.to_str()
        tree.list_tree()
        treelist.append(tree)
        id += 1
        # pplist.pop(ii)
        if last_length == len(pplist):
            break
    
    return accu,cid,treelist


class Point:
   'Class representing point data within a flow structure'

   def __init__(self, name, x, y, ptype=HTree.STR_HTREE_STD):
      self.name = name 
      self.x    = x
      self.y    = y
      self.accu = 0
      self.type = ptype

   def assign_accu(self,accu):
       self.accu = accu[self.x,self.y]
      
def sort_points_in_flow_direction(list_of_points, accu_grid):
    for pi in list_of_points:
        pi.assign_accu(accu_grid)
    return sorted(list_of_points, key=lambda p: -p.accu)

# Testing
import matplotlib.pyplot as plt
import pickle

# Load DTM
dtm_file = '../awaretestdata/data/filled_dtm_tirol_1k.tif'
with rasterio.open(dtm_file) as ds:
    dtm = ds.read(1)
    crs = ds.crs
    prj_settings = ds.meta
    center_lon, center_lat = ds.lnglat()

# modify Achensee pixels
dtm[38,158] = 1011.0
dtm[38,158] = 1011.1
dtm[39,158] = 1011.2
dtm[40,158] = 1011.3
dtm[41,158] = 1011.4
dtm[42,158] = 1011.5
dtm[43,158] = 1011.6
dtm[44,158] = 1011.7
dtm[45,159] = 1011.8
dtm[45,160] = 1011.9
dtm[45,161] = 1097.0

# compute flowdir
flowdir = compute_flow_dir(dtm)

# compute flow accumulation
accu = calc_flow_accumulation(dtm,flowdir)

# include data
#p1 = [173,125] # Etsch
#p2 = [34,186]  # Kirchbichl
#p3 = [90,100]  # Mini-Testgebiet
#p4 = [43,12]   # outside
#p5 = [92,54]   # trisanna / rosanna
#p6 = [137,90]  # obere etsch
#p7 = [62,141]  # innsbruck
#p8 = [71,95]   # oetztal

p2 = Point('Kirchbichl', 34, 186)
p3 = Point('Mini-Testgebiet', 90, 100)
p4 = Point('Bregenzer Ache', 43, 12)
p5 = Point('Trisanna', 92, 54)
p6 = Point('Adige', 155, 132)
p7 = Point('Innsbruck', 62, 141)
p8 = Point('Oetztal', 71, 95)
p9 = Point('Speicher Gepatsch', 100, 86, ptype=HTree.STR_RESERVOIR_CONST)
p10 = Point('Gepatschalm', 106, 86)
p11 = Point('Taschachbach-Fassung', 101, 95)
p12 = Point('Pitzbach-Fassung', 100, 97)
p13 = Point('Isar', 20, 148)
p14 = Point('Achensee', 34, 157)
p15 = Point('Baechental', 32, 145)
p16 = Point('Bieltalbach', 110, 38)
p17 = Point('Jambachtal', 106, 43)
p18 = Point('Larainbachtal', 99, 47)
p19 = Point('Fimbatal', 100, 51)
p20 = Point('Obere Rosana', 90, 42)
p21 = Point('Fasulbach', 90, 45)
p22 = Point('Tiroler Vermuntbach', 106, 37)
p23 = Point('Idbach', 98, 53)
p24 = Point('Sanna', 80, 71)
p25 = Point('Martina', 112, 62)

list_pp = [p16,p6,p8,p3,p17,p5,p4,p7,p18,p2,p9,p10,p11,p15,p12,p13,p19,p14,p20,p21,p22,p23,p24,p25]
# warum geht p6 nicht?

# sort list according to flow accumulation
# list_pp2 = sort_pplist(accu,list_pp)

newlist = sort_points_in_flow_direction(list_pp,accu)

#for pi in list_pp:
#    print('%s: %i' % (pi.name, pi.accu))
#print('')
#for pi in newlist:
#    print('%s: %i' % (pi.name, pi.accu))

# delineate watersheds
flow_accu,basin,treelist = delineate_watersheds(dtm,newlist)

basin_structure_file = '/tmp/treelist.dump'
pickle.dump(treelist, open(basin_structure_file,'wb'), protocol=2 ) # python2 compatibility


plt.close()
plt.figure
plt.imshow(flowdir)
plt.colorbar()
plt.show()


plt.close()
plt.figure
plt.imshow(accu)
plt.colorbar()
plt.show()

plt.close()
plt.figure
plt.imshow(basin)
plt.colorbar()
plt.savefig('/tmp/sub-catchments.pdf')
plt.show()

os.system('python /Users/kristianf/projekte/MUSICALS/AWARE/treevis.py')
