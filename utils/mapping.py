import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
from .queryeps import get_exists_at, get_born_at
from shapely.ops import cascaded_union
from shapely.geometry import LineString, Point


def plot_clusters_at_eps(cluster, eps, ax=None):
    """convinience function to plots clusters at certain eps"""
    cl = get_exists_at(cluster, eps)
    plot_clusters(cl, ax=ax)

    
def plot_clusters(cluster_list, color=None, cmap='tab20', ax=None, **kwargs):
    """
    plots a list of clusters onto the same map
    """
    if not ax:
        f, ax = plt.subplots(figsize=(10,10))
    
    ax.axis('equal')
    
    if color is None:
        iterable = enumerate(cluster_list)
    else:
        iterable = zip(color, cluster_list)
    plot_data = {'geometry': [], 'cluster_number': []}
    for i, sc in iterable:
        plot_data['geometry'].extend([parcel.polygon for parcel in sc.parcels])
        plot_data['cluster_number'].extend([i % 20] * len(sc.parcels))
    
    gpd.GeoDataFrame(plot_data).plot(
        ax=ax, 
        column='cluster_number', 
        edgecolor='white', 
        cmap=cmap, 
        linewidth=.1,
        **kwargs
    )


def plot_with_interal_edge(cluster, eps, ax=None, border_color='k', seperator_color='black', border_width=1, seperator_width=2,
                           debug=False, return_df=False):
    """plots all clusters born at eps and demarcates the barriers of size eps"""
    if ax is None and not return_df:
        f, ax = plt.subplots(figsize=(10,10))
    if ax is not None:
        ax.axis('equal')
    cmap = plt.cm.get_cmap('tab20')

    background = gpd.GeoDataFrame({  # plots grey background for all parcels
        'geometry': [parcel.polygon for parcel in cluster.parcels]
    })
    
    plot_data = {'geometry': [], 'cluster_number': []}
    sub_data = {'geometry': []}
    for i, sc in enumerate(get_born_at(cluster, eps)):
        plot_data['geometry'].extend([parcel.polygon for parcel in sc.parcels])
        plot_data['cluster_number'].extend([i % 20] * len(sc.parcels))
        for j, ssc in enumerate(sc.sub_clusters):
            ssc_polygon = cascaded_union([parcel.polygon.buffer(eps) 
                                          for parcel in ssc.parcels])
            try:
                ssc_polygon = ssc_polygon.buffer(-eps/2)
                ssc_boundary = LineString(ssc_polygon.exterior)
            except AttributeError as e:
                if debug: print('debug', len(list(ssc_polygon)))
                for poly in list(ssc_polygon):
                    if poly.area > .01:
                        ssc_boundary = LineString(poly.exterior)
                        sub_data['geometry'].append(ssc_boundary)
                # plot_clusters([ssc])
                # print(i, j)
                # return ssc
                # raise e
            sub_data['geometry'].append(ssc_boundary)
    
    
    parcels = gpd.GeoDataFrame(plot_data)
    outlines = gpd.GeoSeries(sub_data['geometry'])
    
    if return_df:
        return background, parcels, outlines
    
    background.plot(ax=ax, color='#f0f3f7')
    
    parcels.plot(
        ax=ax, 
        column='cluster_number',
        edgecolor='white', 
        cmap='tab20', 
        linewidth=.05
    )
    
    outlines.plot(ax=ax, edgecolor=seperator_color, linewidth=seperator_width)
    
    for side in ['bottom', 'top', 'right', 'left']:
        ax.spines[side].set_color(border_color)
        ax.spines[side].set_linewidth(border_width)
    ax.tick_params(
        axis='both',
        which='both',
        bottom=False,
        left=False,
        labelbottom=False,
        labelleft=False
    )
    return ax


def utm_scale_bar_guide(ax, pos=(.1, .9), size=.2, **kwargs):
    xl0, xl1 = ax.get_xlim()
    yl0, yl1 = ax.get_ylim() 
    xl_size = xl1 - xl0
    yl_size = yl1 - yl0
    
    length = (xl_size) * size
    
    rounded_length = 10 ** np.floor(np.log10(length))
    print(rounded_length, 'm')
    
    p1 = (xl_size * pos[0] + xl0, yl_size * pos[1] + yl0)
    p2 = (p1[0] + rounded_length, p1[1])
    
    ax.plot(*zip(p1, p2), c='k', **kwargs)
    return ax