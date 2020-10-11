import graph_tool.topology
import graph_tool as gt
import numpy as np
from sklearn import neighbors
import time
from ..dendrogram import hierarchy2dendrogram


class Cluster:
    def __init__(self, death, stat, parcels=None, sub_clusters=None, threshold=5):
        assert (parcels is None) is not (sub_clusters is None)
        self.death = death
        self._parcels = parcels
        self.stat = stat
        self._sub_clusters = sub_clusters
        
        self._major_sub_clusters, self._minor_sub_clusters = [], []
        self.hidden = False
        self.threshold = threshold
        self.re_thresh(threshold)
            
    def re_thresh(self, threshold):
        if self._sub_clusters is not None:
            for sc in self._sub_clusters:
                sc.re_thresh(threshold)
        
        self._major_sub_clusters, self._minor_sub_clusters = [], []
        self.hidden = False
        # https://stackoverflow.com/a/12135169/7587232 for a fun time
        if self._sub_clusters is not None:
            for sc in self._sub_clusters:
                if sc.points.shape[0] >= threshold:
                    self._major_sub_clusters.append(sc)
                else:
                    self._minor_sub_clusters.append(sc)
                    
            self.hidden = len(self._major_sub_clusters) == 1  # TODO (should this be <= or ==)
            
        self.hidden = self.hidden or self.points.shape[0] < threshold
        
        if not self._minor_sub_clusters:
            self._minor_sub_clusters = None
        if not self._major_sub_clusters:
            self._major_sub_clusters = None
            
    def add_sub_clusters(scs):
        if type(scs) is not list:
            scs = [scs]
        self._sub_clusters.extend(scs)
        # self.re_thresh(self.threshold)

    @property
    def parcels(self):
        if self._parcels is None:
            def get_subcluster_parcels(cluster):
                if cluster._parcels is not None:
                    return cluster._parcels
                return [parcel for sc in cluster._sub_clusters 
                        for parcel in get_subcluster_parcels(sc)]
            return get_subcluster_parcels(self)
        return self._parcels
    
    @property
    def _points(self):
        return np.concatenate([p.np_road_points for p in self.parcels])
    
    @property
    def points(self):
        if self._parcels is None:
            def get_subcluster_parcels(cluster):
                if cluster._parcels is not None:
                    return cluster._parcels[0].np_road_points
                return np.concatenate([get_subcluster_parcels(sc)
                                       for sc in cluster._sub_clusters])
            return get_subcluster_parcels(self)
        return self._parcels[0].np_road_points

    @property
    def sub_clusters(self):
        if self._sub_clusters is None:
            return None
        
        def collapse_hidden_chain(cluster):
            """returns first non-hidden cluster"""
            if not cluster.hidden:
                return cluster
            
            if cluster._major_sub_clusters is None:
                return None
            
            assert len(cluster._major_sub_clusters) == 1
            return collapse_hidden_chain(cluster._major_sub_clusters[0])
                    
        def recurse(cluster):
            """returns all non-hidden sub-clusters"""
            if cluster._major_sub_clusters is None:
                return None
            elif not cluster.hidden:
                scs = []
                for sc in cluster._major_sub_clusters:
                    new_sc = collapse_hidden_chain(sc) 
                    if new_sc is not None:
                        scs.append(new_sc)
                if not scs:
                    return None
                return scs
            else:
                assert cluster.hidden
                assert len(cluster._major_sub_clusters) == 1
                return [collapse_hidden_chain(cluster)]
        
        return recurse(self)

    @property
    def _birth(self):
        if self._parcels is None:
            assert len(set([s.death for s in self._sub_clusters])) == 1
            return self._sub_clusters[0].death
        return None
    
    @property
    def birth(self):
        if self._parcels is None and self.sub_clusters is not None:
            return max(sc.death for sc in self.sub_clusters if sc is not None)
        elif self.sub_clusters is None and self._sub_clusters is not None:
            return self._sub_clusters[0].death
        return -np.inf

    @property
    def num_branches(self):
        if self._sub_clusters is None:
            return 1
        return sum([sc.num_branches for sc in self._sub_clusters])

    @property
    def num_np_branches(self):
        """gets the number of non singleton branches"""
        if self._parcels is not None:
            return int(len(self._parcels) != 1)

        sub = sum([sc.num_np_branches for sc in self._sub_clusters])
        if sub == 0:
            return 1
        return sub
    
    def __str__(self):
        sc_string = ""
        if self.sub_clusters is not None:
            sc_string += '{'
            for sc in self.sub_clusters:
                sc_string += '\n\t' + '\n\t'.join(str(sc).split('\n')) + ('\n\t' + '-' * 20)
            sc_string += '\n}'
        else:
            sc_string = "None"
        return ("death: %s\nbirth: %s\nlen(points): %s\nstat: %s\nsub_clusters: %s"
                % (self.death, self.birth, len(self.points), self.stat, sc_string))


def cluster_finder(parcels_series, final_eps, threshold=5, verbose=False, point_groups=None):
    points = np.concatenate([parcel.np_road_points for parcel in parcels_series if parcel.road_points])

    tree = neighbors.BallTree(points)
    neigh = tree.query_radius(points, final_eps, return_distance=True)
    G = gt.Graph()
    G.add_vertex(points.shape[0])

    if verbose: print('created ball tree, neighbors, init graph')
    
    ## create initial cluster map and prune grouped points from neighbors to create idx_dists
    cluster_map = np.array([None] * points.shape[0])
    idx_dists = [np.vstack([x, y.astype(int)]) for x, y in zip(neigh[0], neigh[1])]
    base_idx = 0
    for parcel in parcels_series:
        cluster = Cluster(np.inf, None, parcels=[parcel], sub_clusters=None, threshold=5)
        group = np.arange(base_idx, base_idx + parcel.np_road_points.shape[0])
        for idx in group:
            if idx != base_idx:
                G.add_edge(base_idx, idx)
                G.add_edge(idx, base_idx)
            cluster_map[idx] = cluster
            mask = ~np.isin(idx_dists[idx][0], group)
            idx_dists[idx] = idx_dists[idx][:, mask]
        
        base_idx += group.shape[0]
            
    if verbose: print('created cluster_map array thing, initialized graph and pruned neighbors distances from point groups')
    
    ## create dendogram
    if verbose: print(len(set(cluster_map)), len(cluster_map))
    comp, hist = gt.topology.label_components(G)
    for eps in range(1, final_eps + 1):
        try:
            start = time.time()
            new_edge_set = set()
            for i, x in enumerate(idx_dists):
                mask = (x[1] <= eps) & (x[1] > (eps - 1))
                idxs = list(x[0, mask])

                for idx in idxs:
                    if comp.a[i] != comp.a[idx]:
                        G.add_edge(i, idx)
                        new_edge_set.add(idx)
            if verbose: print(time.time() - start)

            if len(new_edge_set) > 0:
                comp, hist = gt.topology.label_components(G)
                for component in set(comp.a[list(new_edge_set)]):
                    scs = []
                    for sc in set(cluster_map[comp.a==component]):
                        sc.death = eps  # TODO: should this be eps or eps - 1
                        scs.append(sc)

                    cluster_map[comp.a==component] = Cluster(np.inf, None, parcels=None,
                                                             sub_clusters=scs, threshold=threshold)

            if verbose: print(time.time() - start)
            curr_num_clusters = len(set(cluster_map))
            if curr_num_clusters == 1:
                break
            if verbose:    
                print('%s clusters at eps=%s in %ss' % (curr_num_clusters, eps, time.time() - start))
                print('---')
        except KeyboardInterrupt:
            print('interupted during graph part @ eps=%s. you should probably ignore clusters born at this eps' % eps)
            break

            
    clusters = set(cluster_map)
    if len(clusters) == 1:
        return clusters.pop()
    return clusters


def convert2Z(cluster):
    """converts dendrogram to scipy Z"""
    
    G = hierarchy2dendrogram(cluster)
    
    # below from https://stackoverflow.com/questions/35490371/how-to-visualize-dendrogram-a-dictionary-of-hierarchical-items
    nodes       = G.nodes()
    leaves      = set( n for n in nodes if G.out_degree(n) == 0 )
    inner_nodes = [ n for n in nodes if G.out_degree(n) > 0 ]

    # Compute the size of each subtree
    subtree = dict( (n, [n]) for n in leaves )
    for u in inner_nodes:
        children = set()
        node_list = list(G.neighbors(u))
        while len(node_list) > 0:
            v = node_list.pop(0)
            children.add( v )
            node_list += list(G.neighbors(v))

        subtree[u] = sorted(children & leaves)

    inner_nodes.sort(key=lambda n: len(subtree[n])) # <-- order inner nodes ascending by subtree size, root is last

    # Construct the linkage matrix
    leaves = sorted(leaves)
    index  = dict( (tuple([n]), i) for i, n in enumerate(leaves) )
    Z = []
    k = len(leaves)
    for i, n in enumerate(inner_nodes):
        children = list(G.neighbors(n))
        x = children[0]
        for y in children[1:]:
            z = tuple(subtree[x] + subtree[y])
            i, j = index[tuple(subtree[x])], index[tuple(subtree[y])]
            Z.append([i, j, float(len(subtree[n])), len(z)]) # <-- float is required by the dendrogram function
            index[z] = k
            subtree[z] = list(z)
            x = z
            k += 1
    return Z