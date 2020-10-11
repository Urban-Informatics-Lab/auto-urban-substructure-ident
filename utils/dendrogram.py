import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from .polygon_based.cluster import Cluster

def hierarchy2dendrogram(cluster, reverse=False, block_leaves=True, return_root=False):
    """
    @param reverse: edges point toward parent (default False: edges point toward child)
    @param block_leaves: leaves are blocks (default True)
    @param return_root: the root_node_id should always be 0 but return just in case...
    """
        
    def helper(cluster, G):
        root_node_id = G.number_of_nodes()
        assert root_node_id not in G
        G.add_node(root_node_id, size=len(cluster.parcels), cluster=cluster)
        
        positions = {}
#         first_x = 0
        next_x = 0
        new_next_x = next_x
        temp = sorted(cluster.sub_clusters, key=lambda x: x.birth, reverse=True)
        scs = temp[len(temp)%2::2] + temp[::-2]
        for sc in scs:
            persistence = cluster.birth - np.maximum(sc.birth, 0)  # np.maximum is for backwards 
                                                                   # compatability for when leafs have -inf birth
            if cluster.birth > 0 or (not block_leaves and sc.sub_clusters is not None):
                node_id, desc_poss = helper(sc, G)
                G.add_node(node_id, size=len(sc.parcels), cluster=sc)
                if reverse:
                    G.add_edge(node_id, root_node_id, persistence=persistence)
                else:
                    G.add_edge(root_node_id, node_id, persistence=persistence)
                    
                new_next_x = 0
                for k, v in desc_poss.items():
                    positions[k] = (v[0] + next_x, v[1])
                    new_next_x = max(new_next_x, v[0] + next_x)
                next_x = new_next_x + 1                
                
            elif not block_leaves:
                node_id = G.number_of_nodes()
                G.add_node(node_id, size=len(sc.parcels), cluster=sc)
                if reverse:
                    G.add_edge(node_id, root_node_id, persistence=persistence)
                else:
                    G.add_edge(root_node_id, node_id, persistence=persistence)
                positions[node_id] = (next_x, np.maximum(sc.birth, 0))
                next_x += 1
            else:
                desc_poss = None
                
        if next_x == 0:  # if node has no subtree
            root_x = 0
        else:
            root_x = (next_x - 1) / 2
        
        positions[root_node_id] = (root_x, np.maximum(cluster.birth, 0))
        return root_node_id, positions
        
    G = nx.DiGraph()
    root_node_id, positions = helper(cluster, G)
    nx.set_node_attributes(G, positions, name='position')
    if return_root:
        return G, root_node_id
    return G

def plot_dendrogram(G, ax=None, block_leaves=True, include_labels=True, **kwargs):
    if type(G) is Cluster:
        G = hierarchy2dendrogram(G, block_leaves=block_leaves)        
    positions = dict(nx.get_node_attributes(G, 'position'))
    sizes = list(nx.get_node_attributes(G, 'size').values())
    if ax is None:
        f, ax = plt.subplots(figsize=(50,10), dpi=200)
    
    if 'node_size' not in kwargs:
        kwargs['node_size'] = 2
    if 'width' not in kwargs:
        kwargs['width'] = .2
    if 'node_color' not in kwargs:
        kwargs['node_color'] = np.log(sizes)
    if 'cmap' not in kwargs:
        kwargs['cmap'] = plt.cm.viridis_r

    im = nx.draw(G=G, pos=positions, ax=ax, **kwargs)
    if include_labels:
        ax.axis('on')
        ax.set_ylabel('Birth $\epsilon$')
        ax.yaxis.set_tick_params(labelleft=True)