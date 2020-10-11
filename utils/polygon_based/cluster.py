import numpy as np
import fastcluster
from shapely.ops import cascaded_union
from shapely.geometry import MultiPolygon
from .misc import cache


class Cluster:
    def __init__(self, death, parcels=None, parent=None, sub_clusters=None):
        assert (parcels is None) is not (sub_clusters is None)
        self.death = death
        self._parcels = parcels
        self.sub_clusters = sub_clusters 
        self.parent = parent
    
    @property
    def parcels(self):
        if self._parcels is None:
            def get_subcluster_parcels(cluster):
                if cluster._parcels is not None:
                    return cluster._parcels
                return [parcel for sc in cluster.sub_clusters 
                        for parcel in get_subcluster_parcels(sc)]
            return get_subcluster_parcels(self)
        return self._parcels
    
    def get_polygon_list(self):
        """
        in case I ever decide to implement parcels differently/remove them
        might contain MultiPolygons if that's what the original parcel was
        """
        return [p.polygon for p in self.parcels]
    
    def as_MultiPolygon(self):
        """returns a multipolygon of all contained parcels"""
        from .parcel import proper_MultiPolygonizer
        return proper_MultiPolygonizer(self.parcels)
                
    @property
    def birth(self):
        if self._parcels is None:
            assert len(set([s.death for s in self.sub_clusters])) == 1
            return self.sub_clusters[0].death
        return 0
    
    def __str__(self):
        sc_string = ""
        if self.sub_clusters is not None:
            sc_string += '{'
            for sc in self.sub_clusters:
                sc_string += '\n\t' + '\n\t'.join(str(sc).split('\n')) + ('\n\t' + '-' * 20)
            sc_string += '\n}'
        else:
            sc_string = "None"
        return ("death: %s\nbirth: %s\nlen(parcels): %s\nsub_clusters: %s"
                % (self.death, self.birth, len(self.parcels), sc_string))
    
    @cache
    def merged_polygon(self):
        eps = self.birth
        polygon = cascaded_union([parcel.polygon.buffer(eps) for parcel in self.parcels])
        polygon = polygon.buffer(-eps/2)
        if type(polygon) is MultiPolygon:
            debug = len(list(polygon))
            for poly in list(polygon):
                if poly.area > .01:
                    polygon = poly
        return polygon


def cluster_finder(polygons):
    """returns a matrix Z as described in scipy.cluster.hierarchy.linkage"""
    def distfunc(u, v):
        return polygons[int(u)].distance(polygons[int(v)])
    
    X = np.arange(len(polygons))[:, np.newaxis]
    return fastcluster.linkage_vector(X, method='single', metric=distfunc)


def cluster_wrapper(Z, parcel_obj_list, discretize=False, verbose=True):
    """
    :param Z: a matrix Z as described in scipy.cluster.hierarchy.linkage
    :param parcel_obj_list an iterable (e.g. Series/List/1D array) of parcels
    :param discretize: if True, ceilings all birth/deaths. The tree will be different
                       as a result (default False)
    :returns a Cluster object
    """
    n = parcel_obj_list.shape[0]
    
    cluster_map = np.array([Cluster(np.inf, parcels = [p], sub_clusters=None) 
                            for p in parcel_obj_list])
    cluster_dict = {}

    for i in range(Z.shape[0]):
        if i % 1000 == 0 and verbose: print(i)
        sc1_index, sc2_index, death, _ = Z[i]
        sc1_index, sc2_index = int(sc1_index), int(sc2_index)

        leaf_nodes = []
        if sc1_index < n:
            leaf_nodes.append(sc1_index)
        else:
            leaf_nodes.extend(cluster_dict[sc1_index])
            del cluster_dict[sc1_index]

        if sc2_index < n:
            leaf_nodes.append(sc2_index)
        else:
            leaf_nodes.extend(cluster_dict[sc2_index])
            del cluster_dict[sc2_index]

        cluster_dict[n + i] = leaf_nodes

        scs = []
        new_cluster = Cluster(np.inf, parcels=None, parent=None, sub_clusters=scs)
        for sc in set(cluster_map[leaf_nodes]):
            sc.death = death
            sc.parent = new_cluster
            scs.append(sc)

        cluster_map[leaf_nodes] = new_cluster

    final_cluster = list(set(cluster_map))[0]
    
    if discretize:
        discretize_eps(final_cluster)
        
    return final_cluster


def discretize_eps(cluster):
    """
    Ceilings all birth/deaths for sub clusters in cluster.
    The tree will be modified in place.
    """
    if cluster.death < np.inf:
        cluster.death = int(np.ceil(cluster.death))
        
    if cluster.sub_clusters is None:
        return [cluster], cluster.death
    
    both_scs = []
    for sc in cluster.sub_clusters:
        scs, sc_death = discretize_eps(sc)
        both_scs.extend(scs)
    
    if cluster.death == sc_death:
        return both_scs, sc_death
    else:
        for sc in both_scs:
            sc.parent = cluster
        cluster.sub_clusters = both_scs
        if cluster in cluster.sub_clusters:
            raise Exception
        return [cluster], cluster.death
