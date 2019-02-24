import numpy as np


def dpm(vector):
    n = len(vector)
    dpm = np.std(vector, ddof=1) * np.sqrt(n)
    return dpm

def tukey_biweight(vector, c=5):
    m = np.median(vector)
    s = np.median(np.abs(vector-m))
    u = (vector-m) / ((c*s)+1e-4)
    mask = (np.abs(u) > 1)
    w = (1-u**2) ** 2
    w[mask] = 0
    tbi = np.sum(w*vector) / np.sum(w)
    return tbi

def entropy(vector):
    n = len(vector)
    if not np.any(vector):
        return np.log2(n)
    else:
        p = vector / np.sum(vector)
        p = p[np.nonzero(p)[0]]
        h = -1 * np.dot(p, np.log2(p))
        return h

def roku(vector):
    tbi = tukey_biweight(vector)
    vector_p = np.abs(vector-tbi)
    roku = entropy(vector_p)
    return roku

def js_distance(p, q):
    p = p / np.sum(p)
    q = q / np.sum(q)
    left = entropy((p+q)/2)
    right = (entropy(p)+entropy(q)) / 2
    js = left - right
    jsd = np.sqrt(js)
    return jsd
