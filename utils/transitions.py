import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar
from matplotlib.patches import Rectangle


def desc_probs(cluster, normalize=0):
    """
    normalize = -1: unnormalized
                 0: rows add to 1 (right stochastic)
                 1: columns add to 1 (left stochastic)
                 2: matrix adds to 1 (doubly stochastic)
    """
    transition_mat = np.zeros((cluster.birth + 1, cluster.birth + 1))
    
    def recurse(cluster):
        if cluster.sub_clusters is None:
            return
        for sc in cluster.sub_clusters:
            if sc.birth < 0:  # TODO fix cluster birth property
                sc_birth = 0
            else:
                sc_birth = sc.birth
            transition_mat[sc.death, sc_birth] += 1
            recurse(sc)
    
    recurse(cluster)
    
    if normalize == 0:
        transition_mat = transition_mat / transition_mat.sum(axis=1)[:, np.newaxis]
    elif normalize == 1:
        transition_mat = transition_mat / transition_mat.sum(axis=0)[np.newaxis, :]
    transition_mat[np.isnan(transition_mat)] = 0
    
    return transition_mat


def plot_transitions(transition_mat, ax=None, firstok=11, maxdim=50, highlights=[], labels=True,
                     highlight_colors=['b', 'g', 'r', 'c', 'm', 'y']):
    """only works for desc_probs normalize=0 for now"""

    colsums = transition_mat.sum(axis=0, keepdims=True)
    colsums = np.hstack([colsums[:, :firstok].sum(axis=1, keepdims=True), 
                         colsums[:, firstok:maxdim]])
    print(colsums.shape)
    transition_mat = transition_mat[firstok:maxdim]
    transition_mat = np.hstack([transition_mat[:, :firstok].sum(axis=1, keepdims=True),
                                transition_mat[:, firstok:maxdim]])

    if not ax:
        f, ax = plt.subplots(figsize=(20, 15), dpi=200)
    im = ax.imshow(transition_mat * 100, cmap='Blues')
    
    ax.set_yticks(range(transition_mat.shape[0]))
    ax.set_xticks(range(transition_mat.shape[1]))
    ax.xaxis.tick_top()
    ax.axvline(0.5, color='black')

    if labels:
        ax.set_yticklabels(range(firstok, transition_mat.shape[1]))
        ax.set_xticklabels(['"0"'] + list(range(firstok, transition_mat.shape[1])))
        ax.set_xlabel('Bounding Epsilon (m)')
        ax.set_ylabel('Subdividing Epsilon (m)')
    else:
        ax.set_yticklabels([])
        ax.set_xticklabels([])

        
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

        
    ax_divider = make_axes_locatable(ax)
    cbax = ax_divider.append_axes("right", size="10%", pad="3%")
    cb1 = colorbar(im, cax=cbax, format='%.0f%%')
    if labels:
        cb1.set_label_text('Subdivision Probability')
    
    tax = ax_divider.append_axes("top", size="4%", pad="4%", sharex=ax)
    tax.set_yticks([])
    tax.xaxis.set_tick_params(labelbottom=False)
    colsums[:, 0] = np.nan # just to make the other values more distinct
    tax.imshow(colsums, cmap='Blues')#, aspect='auto')
    tax.axvline(0.5, color='black')
#     tax.set_xlim(1, tax.get_xlim()[1])
    if labels:
        tax.set_title('Subdivison Matrix and Disprortionality Vector')
    tax.spines['top'].set_visible(False)
    tax.spines['right'].set_visible(False)
    tax.spines['bottom'].set_visible(False)
    tax.spines['left'].set_visible(False)

    for i, highlight in enumerate(highlights):
        size = transition_mat.shape[1]
        xlocation = highlight - 0.5 - firstok + 1
        tax.add_patch(Rectangle((xlocation, -0.5), 1, 1, fill=False, 
                                edgecolor=highlight_colors[i], lw=2.5))
        ax.add_patch(Rectangle((xlocation, -0.5), 1, size-1, fill=False, 
                               edgecolor=highlight_colors[i], lw=2))