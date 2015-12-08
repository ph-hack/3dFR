from numpy import array, zeros, argmin, inf, gradient, ones
from numpy.linalg import norm


def dtw(x, y, dist=lambda x, y: norm(x - y, ord=1)):
    """ Computes the DTW of two sequences.
    :param array x: N1*M array
    :param array y: N2*M array
    :param func dist: distance used as cost measure (default L1 norm)
    Returns the minimum distance, the accumulated cost matrix and the wrap path.
    """
    x = array(x)
    if len(x.shape) == 1:
        x = x.reshape(-1, 1)
    y = array(y)
    if len(y.shape) == 1:
        y = y.reshape(-1, 1)

    r, c = len(x), len(y)

    D = zeros((r + 1, c + 1))
    D[0, 1:] = inf
    D[1:, 0] = inf

    for i in range(r):
        for j in range(c):
            D[i+1, j+1] = dist(x[i], y[j])

    for i in range(r):
        for j in range(c):
            D[i+1, j+1] += min(D[i, j], D[i, j+1], D[i+1, j])

    D = D[1:, 1:]

    dist = D[-1, -1] / sum(D.shape)

    return dist, D, _trackeback(D)


def _trackeback(D):
    i, j = array(D.shape) - 1
    p, q = [i], [j]
    while (i > 0 and j > 0):
        tb = argmin((D[i-1, j-1], D[i-1, j], D[i, j-1]))

        if (tb == 0):
            i = i - 1
            j = j - 1
        elif (tb == 1):
            i = i - 1
        elif (tb == 2):
            j = j - 1

        p.insert(0, i)
        q.insert(0, j)

    p.insert(0, 0)
    q.insert(0, 0)
    return (array(p), array(q))


def dtw_gradient(x, y, n=10):
    """ Computes the DTW of two sequences.
    :param array x: N1*M array
    :param array y: N2*M array
    :param func dist: distance used as cost measure (default L1 norm)
    Returns the minimum distance, the accumulated cost matrix and the wrap path.
    """
    x = array(x)
    y = array(y)

    gx = gradient(x)
    gy = gradient(y)

    r, c = len(x), len(y)

    D = zeros((r + 1, c + 1))
    D[0, 1:] = inf
    D[1:, 0] = inf

    for i in range(r):
        for j in range(c):
            D[i+1, j+1] = gradient_dist(gx, gy, i, j, n)

    for i in range(r):
        for j in range(c):
            D[i+1, j+1] += min(D[i, j], D[i, j+1], D[i+1, j])

    D = D[1:, 1:]

    dist = D[-1, -1] / sum(D.shape)

    return dist, D, _trackeback(D)


def gradient_dist(g1, g2, i1, i2, n=5):

    l1 = i1 if i1 - n < 0 else n
    l2 = i2 if i2 - n < 0 else n
    r1 = len(g1) - i1 if i1 + n > len(g1) - 1 else n
    r2 = len(g2) - i2 if i2 + n > len(g2) - 1 else n

    l = min(l1, l2)
    r = min(r1, r2)

    w = ones((len(range(i1-l,i1+r)),)) * 1/(2*n)

    w[l] = 1

    d = abs(g1[slice(i1-l,i1+r)] - g2[slice(i2-l,i2+r)]) * w

    return sum(d)/sum(w)
