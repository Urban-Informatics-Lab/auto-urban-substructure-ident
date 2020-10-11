def axline(ax, slope, intercept, **kwargs):
    '''
    originally described (July '18) in ../pluto/hole_areas.ipynb
    see https://github.com/dstansby/matplotlib/blob/49da75e46ccc4714009b124a58784eaeaea53970/lib/matplotlib/axes/_axes.py
    and
    https://github.com/matplotlib/matplotlib/pull/9321
    '''
    import matplotlib.transforms as mtransforms
    import matplotlib.lines as mlines

    if "transform" in kwargs:
        raise ValueError("'transform' is not allowed as a kwarg; "
                         "axline generates its own transform.")

    xtrans = mtransforms.BboxTransformTo(ax.viewLim)
    viewLimT = mtransforms.TransformedBbox(
        ax.viewLim,
        mtransforms.Affine2D().rotate_deg(90).scale(-1, 1))
    ytrans = (mtransforms.BboxTransformTo(viewLimT) +
              mtransforms.Affine2D().scale(slope).translate(0, intercept))
    trans = mtransforms.blended_transform_factory(xtrans, ytrans)

    line = mlines.Line2D([0, 1], [0, 1],
                         transform=trans + ax.transData,
                         **kwargs)
    ax.add_line(line)
    return line