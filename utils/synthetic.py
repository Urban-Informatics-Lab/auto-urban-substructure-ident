import random
import numpy as np
import pandas as pd
import geopandas as gpd
import networkx as nx
import matplotlib.pyplot as plt

from shapely.geometry import Polygon, LineString, MultiPolygon, MultiPoint, Point, MultiLineString, GeometryCollection, JOIN_STYLE
from shapely.affinity import translate, rotate, scale
from shapely.ops import linemerge, unary_union, polygonize, cascaded_union, split
from shapely.geos import TopologicalError

class block_map:
    def __init__(self, block_list):
        self.blocks = block_list
        self.update()
    
    def update(self):
        """updates parcel information from downstream + map border"""
        self.parcels = []
        self.__blocks = []
        points = []
        for block in self.blocks:
            self.parcels += list(block.parcels)
            self.__blocks.append(block.polygon)
            try:
                points += list(block.polygon.exterior.coords)
            except AttributeError:
                points += [point for poly in list(block.polygon) 
                           for point in poly.exterior.coords]
        self.parcels = MultiPolygon(self.parcels)
        self.__blocks = MultiPolygon(self.__blocks)
            
        # See https://gist.github.com/dwyerk/10561690 for posible improvement (concave hull)
        self.polygon = MultiPoint(points).convex_hull
    
    
    def random_parcel_merge(self, num_merges, **kwargs):
        for _ in range(num_merges):
            successful = False
            while not successful:
                try:
                    random.choice(self.blocks).parcel_merge(**kwargs)
                    successful = True
                except IndexError:
                    pass
                except Warning:
                    pass
        self.update()
    
    
    def random_merge_all(self, num_merges):
        for _ in range(num_merges):
            blocks_left = [block for block in self.blocks if len(list(block.parcels)) > 1]
            random.choice(blocks_left).merge_all()
        self.update()
    
    
    def random_parcel_subdivide(self, num_splits, **kwargs):
        for _ in range(num_splits):
            successful = False
            while not successful:
                try:
                    random.choice(self.blocks).parcel_subdivide(**kwargs)
                    successful = True
                except ValueError:
                    pass
        self.update()

    
    def merge_with_map(self, other, road_width, x_offset, y_offset, min_area=5,
                       cap_style=2, join_style=2, return_copy=False):
        # use erase, but draw a buffer around other.polygon representing dividing road width
        other.translate(x_offset, y_offset)
        self.erase(other.polygon.buffer(road_width, cap_style=cap_style,
                                        join_style=join_style),
                   min_area=min_area)
        self.blocks += other.blocks
        self.update()
        if return_copy:
            return block_map(self.blocks.copy())
        
    def erase(self, other_polygon, min_area=5):
        for i, block in reversed(list(enumerate(self.blocks))):
            try:
                block.erase(other_polygon, min_area=min_area)
            except Warning as w:
                if 'Empty Block' not in str(w):
                    raise w
                del self.blocks[i]
        self.update()
        
    
    def rotate(self, theta, origin='centroid'):
        assert origin=='centroid' or (type(origin)==tuple and len(origin)==2)
        if origin == 'centroid':
            origin_point = self.__blocks.centroid
        else:
            origin_point = Point(origin)
        
        for block in self.blocks:
            block.rotate(theta, origin_point)
        self.update()
        
    def translate(self, x_offset, y_offset):
        for block in self.blocks:
            block.translate(x_offset, y_offset)
        self.update()
    
    
    def characterize_centroid(self):
        index_list = [] #[[x] for x in range(len(self.parcels))]
        coords = np.array([parcel.centroid.coords[0] for parcel in self.parcels])
        return index_list, coords
    
    
    def characterize_roadpoint(self, oneroad=False):
        if oneroad:
            raise NotImplementedError('block.frontage would have to be implemented ' +
                                      'differently (return a multilinestring for each parcel)')
        points_list = []
        index_list = []
        replace_map = {}
        last_i = 0
        for block in self.blocks:
            for frontage in block.get_frontage():
                index_inner = []
                if type(frontage) == MultiLineString or type(frontage) == GeometryCollection:
                    for line in list(frontage):
                        points_list.append(line.centroid.coords[0])
                        index_inner.append(len(points_list)-1)
                        replace_map[last_i] = len(index_list)
                        last_i += 1
                    index_list.append(index_inner)
                elif type(frontage) == LineString:
                    points_list.append(frontage.centroid.coords[0])
                    index_inner.append(len(points_list)-1)
                    replace_map[last_i] = len(index_list)
                    last_i += 1
                    index_list.append(index_inner)
                else:
                    print(type(frontage))
                    raise Exception('something weird is happening')
        return index_list, replace_map, np.array(points_list)
            
            
class grid_block_map(block_map):
    """
    Produces a grid of rectangular blocks as a block_map
    @param grid_h: number of blocks of height
    @param grid_w: number of blocks of width
    @param block_h: height of each block
    @param block_w: width of each block
    @param vert_road_w: widths of roads running vertically
    @param horz_road_w: widths of roads running horizontally
    @param frontage_w: controls the width of parcel frontage
    """
    def __init__(self, grid_h=10, grid_w=10, block_h=80, block_w=274,
                 vert_road_w=50, horz_road_w=20, frontage_w=20, theta=0):
        block_list = []
        for x in range(grid_w):
            for y in range(grid_h):
                block = rectangle_block(x * (block_w + vert_road_w),
                                        y * (block_h + horz_road_w), 
                                        height=block_h,
                                        width=block_w,
                                        frontage_w=frontage_w,
                                        theta=theta)
                
                block_list.append(block)
        super().__init__(block_list)


class block:
    def __init__(self, corners, parcel_borders=[], x=0, y=0, theta=0):
        self.__polygon = translate(Polygon(corners), x, y)
        merged = linemerge([self.__polygon.boundary] + 
                           [LineString(_) for _ in parcel_borders])
        borders = unary_union(merged)
        self.__parcels = MultiPolygon(list(polygonize(borders)))
        self.polygon = self.__polygon
        self.parcels = self.__parcels
        self.__transform_list = []
        self.rotate(theta)  
        # TODO: raise a warning if one of the parcels has no road access    
    
    def parcel_merge(self, parcel1=None, parcel2=None, max_sides=6):
        """
        raise a warning if the two parcels are not neighbors
        randomly choose a neighboring parcel of none given
        """
        
        # randomly select parcels if none gives
        parcel_list = list(self.__parcels)
        if parcel1 is None and parcel2 is None:
            random.shuffle(parcel_list)
            parcel1 = parcel_list.pop()
            merged_parcel = None
            for i, parcel in enumerate(parcel_list):
                # find adjacent parcel
                if type(parcel1.intersection(parcel)) in (LineString, MultiLineString):
                    parcel2 = parcel
                    merged_parcel = parcel1.union(parcel2).simplify(0)
                    if len(merged_parcel.boundary.coords) <= max_sides + 1: 
                        del parcel_list[i]
                        break
                    else:
                        parcel2 = None
                        merged_parcel = None
            if merged_parcel is None:
                merged_parcel = parcel1
                raise Warning('a merge did not occur')
            parcel_list.append(merged_parcel)
            self.__parcels = MultiPolygon(parcel_list)
        else:
            raise NotImplementedError
        self.__polygon = cascaded_union(self.__parcels)
        self.__update()
        
        
    def merge_all(self):
        self.__parcels = MultiPolygon([self.__polygon])
        self.__update()
        return self.parcels
    
    
    def parcel_subdivide(self, parcel=None, num_subdivisions=1, min_frontage=3, eps=0.001):
        """DO THIS FIRST (or parcel merge) before any 
        transformations/map merges just in case"""
        if num_subdivisions != 1:
            raise NotImplementedError
        
        parcel_list = list(self.__parcels)
        if parcel==None:
            parcel = random.randint(0, len(parcel_list)-1)
        intersection = parcel_list[parcel].intersection(self.__polygon.exterior)
        if type(intersection) == MultiLineString:
            line = intersection[0]
            for elem in intersection[1:]:
                if elem.length > line.length:
                    line = elem
        elif type(intersection) == LineString:
            line = intersection
        else:
            print(intersection)
            raise Exception('The parcels dont seem to line up with their block')
        
        if line.length < min_frontage * 2:
            raise ValueError('frontage too small')
        
        perp_line = rotate(line, 90, origin='centroid')

        scaling_factor = 2 * parcel_list[parcel].length / line.length
        perp_line = scale(perp_line, scaling_factor, scaling_factor, origin='centroid')

        merged = linemerge([parcel_list[parcel].boundary, perp_line])
        borders = unary_union(merged)
        new_parcels = list(polygonize(borders))
        
        if len(new_parcels) == 2:
            del parcel_list[parcel]
            parcel_list += new_parcels
        self.__parcels = MultiPolygon(parcel_list)
        self.__polygon = cascaded_union(self.__parcels)
        self.__update(eps=eps)
        
        if type(self.__polygon) in (MultiPolygon, GeometryCollection):
            print('...')
            gpd.GeoDataFrame({'parcels': [self.polygon.boundary, self.polygon.boundary]},
                             geometry='parcels').plot()
            raise Exception('weird rotation float errors probably')
        return self.parcels
    
    
    def erase(self, other_polygon, min_area=5):
        """delete the part of the block that intersects with the polygon"""
        # anti rotate/translate:
        orig_other_polygon = other_polygon
        for trans_type, params in reversed(self.__transform_list):
            if trans_type == 'rotation':
                theta, origin = params
                other_polygon = rotate(other_polygon, -theta, origin)
            elif trans_type == 'translation':
                x_offset, y_offset = params
                other_polygon = translate(other_polygon, -x_offset, -y_offset)

        remaining_parcels = []
        for parcel in self.__parcels:
            remaining_parcel = parcel.difference(other_polygon)
            if remaining_parcel.area > min_area:
                remaining_parcels.append(remaining_parcel)
        self.__parcels = MultiPolygon(remaining_parcels)
        self.__polygon = cascaded_union(self.__parcels)
        self.__update()
            
        
    def __update(self, eps=0.001):
        if len(list(self.__parcels)) == 0:
            raise Warning('Empty Block')

        self.parcels = self.__parcels
        self.polygon = self.__polygon
        for trans_type, params in self.__transform_list:
            if trans_type == 'rotation':
                theta, origin = params
                self.parcels = rotate(self.parcels, theta, origin=origin)
                self.polygon = rotate(self.polygon, theta, origin=origin)
            elif trans_type == 'translation':
                x_offset, y_offset = params
                self.parcels = translate(self.parcels, x_offset, y_offset)
                self.polygon = translate(self.polygon, x_offset, y_offset)
            elif trans_type == 'scaling':
                xfact, yfact, origin = params
                self.parcels = scale(self.parcels, xfact, yfact, origin=origin)
                self.polygon = scale(self.polygon, xfact, yfact, origin=origin)
    
    
    def affine_transform(self, transformation, inplace, args):
        if inplace:
            self.__transform_list.append((transformation, args))
            self.__update()
        else:
            import copy
            new = copy.deepcopy(self)
            new.affine_transform(transformation, True, args)
            return new            
    
    
    def rotate(self, theta, origin='centroid', inplace=True):
        return self.affine_transform('rotation', inplace, (theta, origin))
        
        
    def translate(self, x_offset, y_offset, inplace=True):
        return self.affine_transform('translation', inplace, (x_offset, y_offset))
        
        
    def scale(self, xfact=1, yfact=1, origin='center', inplace=True):
        if type(origin) is int:
            sides = split(self.polygon.boundary, MultiPoint(self.polygon.exterior.coords))
            origin = sides[origin].centroid
        return self.affine_transform('scaling', inplace, (xfact, yfact, origin)) 
        
        
    def get_frontage(self):
        frontage_list = []
        for parcel in self.__parcels:
            frontage = parcel.intersection(self.__polygon.boundary)
            for trans_type, params in self.__transform_list:
                if trans_type == 'rotation':
                    theta, origin = params
                    frontage = rotate(frontage, theta, origin=origin)
                elif trans_type == 'translation':
                    x_offset, y_offset = params
                    frontage = translate(frontage, x_offset, y_offset)
            frontage_list.append(frontage)
        return frontage_list


class rectangle_block(block):
    def __init__(self, x, y, parcel_borders=None, height=80, width=274, frontage_w=None, theta=0):
        if type(parcel_borders) is not list:  # create default partition
            parcel_borders = []
            if height > width:
                parcel_borders.append([(x + (width / 2), y), (x + (width / 2), y + height)])
                num_divisions = height // frontage_w
                parcel_dimension = height / num_divisions
                for i in range(1, num_divisions):
                    coord = parcel_dimension * i
                    parcel_borders.append([(x, y + coord), (x + width, y + coord)])
            else:
                parcel_borders.append([(x, y + (height / 2)), (x + width, y + (height / 2))])
                num_divisions = width // frontage_w
                parcel_dimension = width / num_divisions
                for i in range(1, num_divisions):
                    coord = parcel_dimension * i
                    parcel_borders.append([(x + coord, y), (x + coord, y + height)])
        
        corners = [(x, y), (x, y + height), 
                   (x + width, y + height), (x + width, y)]
        super().__init__(corners, parcel_borders, theta=theta)
        
        
def plot_polygon(poly, colors=None, ax=None):
    if type(poly) is not list:
        poly = [poly] * 2
    if ax is None:
        f, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect('equal')
    if not colors:
        gpd.GeoDataFrame({'parcels': poly}, geometry='parcels').plot(ax=ax, color='white',
                                                                         edgecolor='black')
    else:
        gpd.GeoDataFrame({'parcels': poly, 'colors': colors}, geometry='parcels').plot(column='colors')