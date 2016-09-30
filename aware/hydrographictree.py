# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 08:15:26 2016

@author: kristianf
"""

import numpy as np
import copy

class HTree:
    '''Hydrographic Tree class for AWARE
    
    HTree represnts the structure of a basin and its sub-catchments. Nodes
    represent sub-catchments and segements (i.e. the way how these nodes are
    connected) can be viewed as rivers. The class and its methods provides
    the most important features to generate a dendritic structure of sub-
    catchments and the correct interlinkage of tributaries. This helps in
    assigning model parameters and evaluating results in hydrologic modelling.
    A hydrographic tree can be generated through analyses of digital elevation
    models.
    '''    
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
        '''Adds a node to the list of tributaries.
        
        Parameters
        ----
        tributary: a HTree object that will be treated as upstream node (tributary)
    
        Returns
        ----    

        '''       
        assert isinstance(tributary, HTree)
        tributary.downstream_node = self.id
        self.tributaries.append(tributary)

    def add_new_tributary_id(self, suggested_id=1):
        '''Checks the structure of the tree and returns a new id inteded to be used
        for a new node. A suggestion can be passed to the function as argument. In
        case, the function checks if this id is available yet.
        
        Parameters
        ----
        suggested_id: an id that should be used in case it is not yet defined in the tree
    
        Returns
        ----    
        new id integer
        '''       
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
        '''This function returns the number of nodes in tree. Moreover, the depth
        of the tree is returned, i.e. the maximum number of levels prevailing in
        each branch.
        
        Parameters
        ----
        order: optional argument that is viewed as offset
    
        Returns
        ----    
        number of nodes and depth of tree (integer)
        '''
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
        '''Returns the level of computation (each branching increases this number
        by one).
        
        Parameters
        ----
        level: optional argument that is viewed as offset
    
        '''
        self.level=level
        for i,tributary in enumerate(self.tributaries):
            tributary.calculate_computation_level(level+1)

    def calculate_computation_priority(self,lowest):
        '''This function computes the computation order / priority of all nodes
        in the tree. It warrants that upstreams nodes are generally computed 
        prior to their corresponding following nodes - even in complex trees.
        
        Parameters
        ----
        lowest: the minimum number to be assigned
    
        '''
        self.order = lowest
        for ti in self.tributaries:
            ti.assign_priority(self)
        self.order = lowest

    def assign_priority(self, start_node):
        '''This function assigns the priority of nodes in a tree. It is uesed by
        its calling function calculate_computation_priority() only. The priority
        decreases in order to force root nodes beeing computed at first, which is
        why a positive number greater than the tree size has to be passed as
        argument.
        
        Parameters
        ----
        start_node: the first priority to be assigned.
    
        '''

        start_node.order -= 1
        self.order=start_node.order
        for i,tributary in enumerate(self.tributaries):
            tributary.assign_priority(start_node)

    def calculate_sub_areas(self):
        '''The area of the sub-catchment upstream of each node is calculated
        taking it into account that delineated upstream areas are defined
        seprately (i.e., only the are betwenn two nodes is returned).

        Returns:
        upstream area that is represented by other nodes
        '''
        if len(self.tributaries) == 0:
            self.subarea = self.area
        defined_areas = self.area
        for i,tributary in enumerate(self.tributaries):
            defined_areas -= tributary.calculate_sub_areas()
        self.subarea = defined_areas
        return defined_areas

    def postprocess_tree(self, offset=0):
        '''Calling this function is a convenient way to poceed with the most
        important steps to link all nodes of tree in order to generate valid
        tree properties.
        
        
        Parameters 
        ----
        offest: optional argument used as offest for calculate_computation_priority
    
        '''
        # get number of nodes and number of levels (depth)
        n_nodes,max_depth = self.get_tree_size()      
        # compute computation level 
        self.calculate_computation_level(max_depth)
        # compute computation priority
        self.calculate_computation_priority(n_nodes+offset)
        # define sub-catchment areas without upstream areas
        self.calculate_sub_areas()
       
        return n_nodes, max_depth

    def to_str(self):
        '''to_str() simply prints the tree structure to the terminal. This 
        might be helpful to check tree consistency.
  
        '''
        for i,tributary in enumerate(self.tributaries):
            print('%i<-%i' % (self.id,tributary.id))
            tributary.to_str()
    
    def list_tree(self):
        '''list_tree() lists all nodes including their ids to the screen.
        '''
        for i,tributary in enumerate(self.tributaries):
             print('%i' % (self.id))
             tributary.list_tree()

    def set_area(self,area):
        '''Adds area information the object.
        
        
        Parameters 
        ----
        area: the area which will be assigned to the HTree object
    
        '''
        self.area = area        

    def get_tree_by_id(self,id):
        '''Searches for the node associated to id and returns the node and its
        tree.
        
        Parameters
        ----
        id: id to be looked up
    
        Returns
        ----    
        HTree object. The root node has the id passed as argument to the function
        '''
        target_node = None
        if self.id == id:
            return self
        for i,tributary in enumerate(self.tributaries):
            target_node_i = tributary.get_tree_by_id(id)
            if target_node_i is not None:
                target_node = target_node_i
        return target_node

    
    def get_upstream_areas(self):
        '''Lists ids and areas of all upstream nodes of a tree.

        Returns
        ----    
        two lists including the ids and sub-areas, respectively
        '''
        list_ids = list()
        list_areas = list()
        for i,tributary in enumerate(self.tributaries):
            list_ids.append(tributary.id)
            list_areas.append(tributary.subarea)
            ids_up, areas_up = tributary.get_upstream_areas()
            list_ids += ids_up
            list_areas +=areas_up
        return list_ids, list_areas
        
    def get_downstream_list(self, id, idlist=None):
        '''Generates a list including all nodes downstream of the corresponding
        HTree object. It is called by get_downstream_path only.
        
        Parameters
        ----
        id: id in the tree to be found
        idlist: a list including the current results of the search algorithm

        Returns
        ----    
        updated list including the downstream nodes
        '''
        if idlist is None:
            new_list = list()
        else:
            new_list = copy.deepcopy(idlist)
        new_list.append(self.id)
        if id in new_list:
            return new_list
        for i,tributary in enumerate(self.tributaries):
            li = tributary.get_downstream_list(id, new_list)
            if id in li:
                return li
            else:
                del li
        
        return new_list
        
    def get_downstream_path(self, id):
        '''Computes the downstream path in a tree given that the path start at
        the id specified as argument.
        
        Parameters
        ----
        id: id which represents the sart node

        Returns
        ----    
        list of nodes representing the path downstream of the start node
        '''
        reverse_list = self.get_downstream_list(id)
        return reverse_list[::-1]
        if id in reverse_list:
            return reverse_list[::-1]
        else:
            return [id]
        