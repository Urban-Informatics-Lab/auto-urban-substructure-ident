def get_dead_after(cluster, eps):
    """
    return clusters that were dead after a certain epsilon 
    (e.g. all clusters are the first split below eps apart)
    # TODO: should this be changed to "all clusters are at least eps apart"?
    """
    if cluster.sub_clusters is None:
        return []
    if cluster.birth <= eps:
        return cluster.sub_clusters
    candidates = []
    for sc in cluster.sub_clusters:
        candidates.extend(get_dead_after(sc, eps))
    return candidates


def get_born_after(cluster, eps):
    """
    return clusters that were born before a certain epsilon 
    """
    if cluster.sub_clusters is None:
        return []
    if cluster.birth <= eps:
        return [cluster]
    candidates = []
    for sc in cluster.sub_clusters:
        candidates.extend(get_born_after(sc, eps))
    return candidates


def get_born_before(cluster, eps):
    """
    return clusters that were born before a certain epsilon 
    """
    candidates = []
    if cluster.birth >= eps:
        candidates.append(cluster)
    if cluster.sub_clusters is None:
        return candidates
    for sc in cluster.sub_clusters:
        candidates.extend(get_born_before(sc, eps))
    return candidates


def get_born_at(cluster, eps):
    """
    return clusters that were born at a specific epsilon
    """
    if cluster.sub_clusters is None:
        return []
    if cluster.birth == eps:
        return [cluster]
    candidates = []
    for sc in cluster.sub_clusters:
        candidates.extend(get_born_at(sc, eps))
    return candidates


def get_exists_at(cluster, eps):
    """return clusters that exist at certain epsilon"""
    candidates = []
    if cluster.sub_clusters is None:  
        return candidates
    for sc in cluster.sub_clusters:  
        if sc.birth <= eps < cluster.birth: 
            candidates.append(sc)
        else:
            candidates.extend(get_exists_at(sc, eps))
    return candidates

def get_exists_at_d(dendrogram, eps, root=0):
    """return clusters that exist at certain epsilon"""
    candidates = []
    birth = dendrogram.nodes[root]['position'][1]
    if not dendrogram.successors(root):
        return candidates
    for child in dendrogram.successors(root):
        child_birth = dendrogram.nodes[child]['position'][1]
        if child_birth <= eps < birth: 
            candidates.append(child)
        else:
            candidates.extend(get_exists_at_d(dendrogram, eps, child))
    return candidates