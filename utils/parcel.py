import functools
from shapely.geometry import *
from multiprocessing import Pool
import numpy as np
import geopandas as gpd

class Parcel(GeometryCollection):  
    def __init__(self, polygon, road_points, crs, add_attr={}):  
        """
        @param road_points: list of shapely Point objects denoting parcel/road intersection points
        """
        self.polygon = polygon
        self.road_points = MultiPoint(road_points)
        if road_points: 
            self.np_road_points = np.stack([np.array((point.x, point.y)) for point in road_points])
        else:
            self.np_road_points = np.array([])
        self.crs = crs
        
        self.__add_attr_keys = add_attr.keys()
        for name, value in add_attr.items():
            setattr(self, name, value)
        
        super().__init__([self.road_points] + [polygon])
        
    # functions needed to make pickling work---see how shapely objects are pickled 
    def __new__(cls, polygon, road_points, crs, add_attr={}):
        instance = super(Parcel, cls).__new__(cls)
        instance.__init__(polygon, road_points, crs, add_attr)
        return instance
        
    def __reduce__(self):
        add_attr = {key: getattr(self, key) for key in self.__add_attr_keys}
        init_params = (self.polygon, self.road_points, self.crs, add_attr)
        return (self.__class__, init_params, self.wkb)
            

def proper_MultiPolygonizer(parcel_list):
    polygon_list = []
    for parcel in parcel_list:
        geom = parcel.polygon
        if type(geom) is Polygon:
            polygon_list.append(geom)
        else:
            polygon_list.extend(list(geom))
    return MultiPolygon(polygon_list)


def get_parcel_obj_list(parcel_polys, all_block_polys, crs, min_dist=2, **kwargs):
    """
    @param parcel_polys: a series of polygons representing parcels
    @param all_block_polys: a shapely MultiPolygon representing unique blocks
    @param crs: the cordinate reference system
    @param min_dist: the minimum distance between points sample on block perimeter
    kwargs must be series objects with indices that correspond to parcel_polys
    """
    
    func = functools.partial(sp_road_points, all_block_polys=all_block_polys,
                             min_dist=min_dist)
    with Pool() as p: 
        points_lists = p.map(func, parcel_polys)
    
    parcel_list = []
    rp_geometry_list = []
    for i, points_list in enumerate(points_lists):
        pointobjs = [Point(point) for point in points_list]
        rp_geometry_list.append
        parcel_attr = {field: series.iloc[i] for field, series in kwargs.items()}
        parcel_list.append(Parcel(parcel_polys.iloc[i], pointobjs, crs, parcel_attr))
    
    gdf = gpd.GeoDataFrame({'parcels': parcel_list}, crs=crs, geometry='parcels')
    gdf.set_index(parcel_polys.index, inplace=True)
    return gdf 


def sp_road_points(parcel_poly, all_block_polys, min_dist):
    frontage = parcel_poly.intersection(all_block_polys.boundary)
    points = []
    if type(frontage) == MultiLineString or type(frontage) == GeometryCollection:
        for line in list(frontage):
            assert type(line) == LineString or type(line) == Point
            if type(line) == Point:
                curr_points = [line.coords[0]]
            else:
                num_splits = int(line.length // min_dist)
                if num_splits > 1:
                    curr_points = [line.interpolate((i / num_splits), normalized=True).coords[0]
                                   for i in range(1, num_splits)]
                else:
                    curr_points = [line.centroid.coords[0]]
                    
            if len(curr_points) == 1:  
                if len(points) > 0:
                    dist = (np.array(points) - np.array(curr_points[0]))**2
                    dist = np.sum(dist, axis=1)
                    dist = np.sqrt(dist)
                    if dist.min() < min_dist:
                        continue
            points.extend(curr_points)
        
    elif type(frontage) == LineString:
        num_splits = int(frontage.length // min_dist)
        if num_splits > 1:
            points.extend([frontage.interpolate((i / num_splits), normalized=True).coords[0]
                           for i in range(1, num_splits)])
        else:
            try:
                points.append(frontage.centroid.coords[0])
            except IndexError:
                print('error!', frontage.centroid, type(parcel_poly), type(all_block_polys.boundary))
    
    return points