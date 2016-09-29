# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 08:15:26 2016

@author: kristianf
"""

import numpy as np

class HTree:
    '''Hydrographic Tree class for AWARE'''
    
    def __init__(self, id=0, tributaries=None):
        self.id = id
        self.order=0
        self.level=0
        self.area = 0
        self.subarea = 0
        self.tributaries = []
        self.downstream_node = None
        if tributaries is not None:
            for tributary in tributaries:
                tributary.downstream_node = self.id
                self.add_tributary(tributary)
        
    def add_tributary(self,tributary):
        assert isinstance(tributary, HTree)
        tributary.downstream_node = self.id
        self.tributaries.append(tributary)

    def add_new_tributary_id(self, suggested_id=1):
        # get unique ids that are in use
        if len(self.tributaries) == 0:
            if suggested_id == self.id:
                new_id=suggested_id+1
            else:
                new_id = suggested_id
        else:
            list_ids = []
            list_ids.append(self.id)
            new_id = suggested_id
            for tributary in self.tributaries:
                new_id = tributary.id
                list_ids.append(new_id)
            # generate new unique id that is different from the existing ones
            # modifiy number if already existing
            if suggested_id > new_id:
                new_id = suggested_id
            else:
                new_id += 1
            loop = True
            while loop:
                loop = False
                for existing_id in list_ids:
                    if existing_id == new_id:
                        new_id+=1
                        loop = True
                        break
        return new_id

    def get_tree_size(self,order=0):
        size = 0
        n = len(self.tributaries)
        if n > 0:
            maxii = np.zeros(n, dtype=int)
            for i,tributary in enumerate(self.tributaries):
                si,maxii[i]=tributary.get_tree_size(order+1)
                size+=si
            maxi=max(maxii)
        else:
            maxi=order+1
        
        return size+1,maxi
          
    def calculate_computation_level(self,level=0):
        self.level=level
        for i,tributary in enumerate(self.tributaries):
            tributary.calculate_computation_level(level+1)

    def calculate_computation_priority(self,lowest):
        self.order = lowest
        for ti in self.tributaries:
            ti.assign_priority(self)
        self.order = lowest

    def assign_priority(self, start_node):
        start_node.order -= 1
        self.order=start_node.order
        for i,tributary in enumerate(self.tributaries):
            tributary.assign_priority(start_node)

    def calculate_sub_areas(self):
        if len(self.tributaries) == 0:
            self.subarea = self.area
        defined_areas = self.area
        for i,tributary in enumerate(self.tributaries):
            defined_areas -= tributary.calculate_sub_areas()
        self.subarea = defined_areas
        return defined_areas

    def postprocess_tree(self, offset=0):
        # get number of nodes and number of levels (depth)
        n_nodes,max_depth = self.get_tree_size()      
        # compute computation level 
        self.calculate_computation_level(max_depth)
        # compute computation priority
        self.calculate_computation_priority(n_nodes+offset)
        # define sub-cacthment areas without upstream areas
        self.calculate_sub_areas()
       
        return n_nodes, max_depth

    def to_str(self):
         for i,tributary in enumerate(self.tributaries):
            print('%i<-%i' % (self.id,tributary.id))
            tributary.to_str()
    
    def list_tree(self):
         for i,tributary in enumerate(self.tributaries):
             print('%i' % (self.id))
             tributary.list_tree()

    def set_area(self,area):
        self.area = area        

    def get_tree_by_id(self,id):
        target_node = None
        if self.id == id:
            return self
        for i,tributary in enumerate(self.tributaries):
            target_node_i = tributary.get_tree_by_id(id)
            if target_node_i is not None:
                target_node = target_node_i
        return target_node

    
    def get_upstream_areas(self):
        list_ids = list()
        list_areas = list()
        for i,tributary in enumerate(self.tributaries):
            list_ids.append(tributary.id)
            list_areas.append(tributary.subarea)
            ids_up, areas_up = tributary.get_upstream_areas()
            list_ids += ids_up
            list_areas +=areas_up
        return list_ids, list_areas