import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from .plotting_helpers import axline


def get_persistent_global_scales(transition_mat, tau=np.inf, ignore0=True,
                                 verbose=False, firstok=1):
    colsums = transition_mat.sum(axis=0, keepdims=True)
    f = np.hstack([colsums[:, :firstok].sum(axis=1, keepdims=True),
                   colsums[:, firstok:]]).sum(axis=0)
    return get_persistent_f(f, tau, ignore0, verbose, firstok-1)
    

def get_persistent_f(f, tau=np.inf, ignore0=False, verbose=False, offset=0):
    pers_clusters = {}
    C = np.zeros(len(f), dtype=int) - 1  # -1 means no assignment
    PD_pairs = {}
    
    # "We then process the vertices in decreasing value of f"
    xs = f.argsort()[::-1]
        
    for x in xs:  # x is an eps        
        # "we first determine if it is a local maximum in the mesh by comparing
        # f(x) with f(y) for all y in a one-ring neighbor hood of x"
        # "If x is a local maximum, a new component is born and the vertex is
        # assigned to itself in the segmentation, C(x)=x."
        C[x] = x  # base case
        max_f = f[x]
        
        not_at_right_edge = (x + 1 < len(f))
        if not_at_right_edge:
            right_bigger = (f[x + 1] > max_f)
            right_equal_but_part_cluster = (f[x + 1] == max_f and C[x + 1] != x and C[x + 1] >= 0)
            if right_bigger or right_equal_but_part_cluster:
                C[x] = C[x + 1]  # assign x to the root of x + 1
                max_f = f[x + 1]

        not_at_left_edge = (x - 1 >= 0)
        if not_at_left_edge:
            left_bigger = (f[x - 1] > max_f)
            left_same_but_part_cluster = (f[x - 1] == max_f and C[x - 1] != x and C[x - 1] >= 0)
            if left_bigger or left_same_but_part_cluster:
                C[x] = C[x - 1]  # assign x to the root of x - 1
                max_f = f[x - 1]
        # -------
           
        # "If the vertex is adjacent to two or more existing components"
        if 0 < x < (len(f) - 1) and C[x - 1] >= 0 and C[x + 1] >= 0:
            small_neigh, big_neigh = sorted([C[x-1], C[x+1]], key=lambda neigh: f[neigh])
            if verbose: 
                print(C[x-1], x, C[x+1])
                print(C)
                print(small_neigh, big_neigh)
            # "we check the persistence of the components and merge them only if
            # they are not τ-persistent"
            if np.abs(f[C[x - 1]] - f[C[x + 1]]) < tau:
                C[C==small_neigh] = big_neigh
                PD_pairs[(f[x] + offset, f[small_neigh] + offset)] = small_neigh + offset
    
    if not ignore0:
        PD_pairs[(0, f[C[0]])] = C[0]
    if verbose:
        print(C)
    return PD_pairs


def plot_PD(PD_points, tau=None, ax=None, labels=False, **kwargs):
    if ax is None:
        f, ax = plt.subplots(figsize=(7, 7))
        
    if 'alpha' not in kwargs:
        kwargs['alpha'] = .2
    if 'size' not in kwargs:
        kwargs['s'] = 10
    if 'c' not in kwargs:
        kwargs['c'] = 'k'
    plt.scatter(*zip(*PD_points), **kwargs)
    
    axline(ax, 1, -.01)

    if tau is not None:
        tau_points = {(b, d): eps for (b, d), eps in PD_points.items() if (d - b >= tau)}
        kwargs['c'] = 'r'
        plt.scatter(*zip(*tau_points), **kwargs)
        axline(ax, 1, tau-.01, linestyle='--', )
        
        if labels:
            for (b, d), eps in tau_points.items():
                ax.annotate(eps, (b, d))

    
    plt.axis('square')
    ax.set_xlim(left=-.25)
    ax.set_ylim(bottom=-.25)
    
    if tau is not None:
        return tau_points


def composite_persistence_plot(transition_mat, tau, firstok=1):
    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(13,4), gridspec_kw={'width_ratios': [2, 1]})
    PD = get_persistent_global_scales(transition_mat, firstok=firstok)
    persistent_scales = plot_PD(PD, tau=tau, alpha=1, ax=ax2, labels=True)
    ax2.set_title('Scale Persistence Diagram')
    ax2.set_xlabel('Scale Birth Probsum')
    ax2.set_ylabel('Scale Death Probsum')
    
    sns.kdeplot([d - b for b, d in PD], bw=.05, ax=ax1)
    ax1.axvline(tau, linestyle='--')
    ax1.set_title('density of scales with given persistence')
    ax1.set_ylabel('scale density')
    ax1.set_xlabel('persistence')

    plt.show()
    return persistent_scales
