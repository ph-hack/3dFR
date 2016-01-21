import numpy as np
import re
import scipy.spatial.distance as dist
import scipy.stats as stats
import cv2
import os
import matplotlib.pyplot as plt
import copy as cp
from distances import dtw, dtw_gradient, p2p_dist
from unittest.case import TestCase
from unittest.loader import TestLoader
from unittest.runner import TextTestRunner

def load_data_split(train_path, test_path):

    train_file = open(train_path, 'r')
    test_file = open(test_path, 'r')

    train = train_file.read()
    test = test_file.read()

    train = train.split('\n')
    test = test.split('\n')

    train = train[:-1]
    test = test[:-1]

    train = [get_sample_id(t) for t in train]
    test = [get_sample_id(t) for t in test]

    return train, test

def read_curves(file_path, **kwargs):

    curve_idx = 0 if 'curve_idx' not in kwargs else kwargs['curve_idx']

    try:
        with open(file_path, 'r') as f:

            for i in range(0, curve_idx+1):
                line = f.readline()

            sts = line.split(' ')
            curve = np.float32(sts)

            return curve

    except Exception as e:

        print 'Error reading the curves file: ', file_path, ':\n', e.message

def read_saliency(folder, sample_id, curve_idx=0):

    saliencies = np.loadtxt(''.join([folder, 'saliency_', str(sample_id), '.txt']))

    return saliencies[curve_idx]

def get_person_id(file_path):

    name = get_sample_id(file_path)

    id = name.split('d')[0]

    return id

def get_sample_id(file_path):

    m = re.search('[0-9]+d[0-9]+', file_path)

    if m is None:

        print 'Error in file ', file_path
        name = 'error'

    else:
        name = m.group(0)

    return name

def representant_curve(curves, distance='cosine'):

    distances = {

        'cosine': dist.cosine,
        'euclidean': dist.euclidean,
        'dtw': dtw,
        'dtwg': dtw_gradient
    }

    m = np.mean(curves, 0)

    ds = np.array([distances[distance](c, m) for c in curves[:]])

    rep_idx = ds.argmin()

    return curves[rep_idx,:]

def smooth(x, n):

    x = cp.deepcopy(x)

    for i in range(0, len(x)):

        l = 0 if i - n <= 0 else i - n
        r = len(x) if i + n >= len(x) else i + n

        w = slice(i,r)

        x[i] = np.mean(x[w])

    return x

def histogram(x, n_bins=20):

    if type(x) is not np.ndarray:

        x = np.array(x)

    if len(x.shape) > 1:

        raise Exception('x should be 1-D')

    x_max = max(x)
    x_min = min(x)

    delta = (x_max - x_min)/n_bins

    H = np.zeros(n_bins)
    X = np.zeros(n_bins)
    Y = np.zeros(n_bins)

    for i in window(x_min, delta, n_bins):

        H[i[2]] = len(np.where((x >= i[0]).__and__(x < i[1]))[0])
        X[i[2]] = i[0]
        Y[i[2]] = i[1]

    return H, X, Y

def window(m, d, n):

    ds = np.ones(n+1) * d
    ds[0] = m
    ds = np.cumsum(ds)

    for i in range(0,n):

        yield (ds[i], ds[i+1], i)

def saliency(from_folder, distance, to_folder, curves_idx=range(0,11)):

    file_list = os.listdir(from_folder)
    files = [f for f in file_list if f.endswith('.lines')]

    curves = []

    for c in curves_idx:

        curve_of_files = []
        curves.append(curve_of_files)

        for f in files:

            curve_of_files.append(read_curves(''.join([from_folder, f]), curve_idx=c))

    f_classes = np.array([get_person_id(f) for f in files])
    classes = np.array(set(f_classes))

    for f in range(0,len(files)):

        S = np.zeros((len(curves_idx), 1))

        for c in curves_idx:

            curve = curves[c][f]

            this_class = f_classes[f]

            other_classes_idx = np.where(f_classes != this_class)[0]

            dists= []

            for f2 in other_classes_idx:

                other = curves[c][f2]
                d = distance(curve, other)
                dists.append(d)

            h, x = np.histogram(dists, 20)

            weib = stats.exponweib.fit(dists, 1, 1, floc=0)

            weib_curve = stats.exponweib.pdf(x[:20], *weib)

            x_max = weib_curve.argmax()
            S[c] = (x[x_max] + x[x_max+1])/2

            # plt.plot(x[:20], h)
            # plt.plot(x[:20], weib_curve)
            # plt.show()

        np.savetxt(''.join([to_folder, 'saliency_', get_sample_id(files[f]), '.txt']), S)

        print 'file ', f


def remove_non_face(curve):

    c_sum = np.cumsum(curve)

    begin = np.where(c_sum == 0)
    end = np.where(c_sum == c_sum[-1])

    # print 'begin=', begin, ' end=', end

    if len(begin[0]) > 0:

        begin = begin[0][-1]
    else:

        begin = -1

    if len(end[0]) > 0:

        end = end[0][0]
    else:

        end = len(curve)

    return curve[begin+1:end-1]


# UNIT TESTS ##########################################################################################################

class PreProcessorsTests(TestCase):

    def est_01_read_curves(self):

        curve_file = '/home/hick/Documents/Mestrado/Research/Code/Experiments5/training/02463d452.lines'

        curve = read_curves(curve_file, curve_idx=0)

        self.assertEqual(len(curve), 130)

    def est_02_person_id(self):

        id = get_person_id('alsdf/asdfad/0200d345.lines')

        self.assertEqual(id, '0200')

        id = get_sample_id('alsdf/asdfad/0200d345.lines')

        self.assertEqual(id, '0200d345')

    def est_03_dtw(self):

        f1 = '/home/hick/Documents/Mestrado/Research/Code/Experiments5/training/02463d452.lines'
        f2 = '/home/hick/Documents/Mestrado/Research/Code/Experiments5/training/02463d654.lines'
        f3 = '/home/hick/Documents/Mestrado/Research/Code/Experiments5/training/04202d350.lines'

        c1 = read_curves(f1)
        c2 = read_curves(f2)
        c3 = read_curves(f3)

        dist = dtw_gradient(c1,c2)
        print 'dtwg: c1 x c2 = ', dist

        dist = dtw(c1,c2)
        print 'dtw: c1 x c2 = ', dist

        dist = dtw_gradient(c1,c3)
        print 'dtwg: c1 x c3 = ', dist

        dist = dtw(c1,c3)
        print 'dtw: c1 x c3 = ', dist

        # plt.imshow(D)
        # plt.show()

    def est_04_smooth(self):

        f1 = '/home/hick/Documents/Mestrado/Research/Code/Experiments5/training/02463d452.lines'
        f2 = '/home/hick/Documents/Mestrado/Research/Code/Experiments5/training/02463d654.lines'
        f3 = '/home/hick/Documents/Mestrado/Research/Code/Experiments5/training/04202d350.lines'

        c1 = read_curves(f1)
        c2 = read_curves(f2)
        c3 = read_curves(f3)

        plt.plot(c1, 'r-')

        c1 = smooth(c1, 5)
        c2 = smooth(c2, 5)
        c3 = smooth(c3, 5)

        plt.plot(c1, 'b-')
        plt.show()

        dist = dtw_gradient(c1,c2)
        print '\ndtwg: c1 x c2 = ', dist

        dist = dtw(c1,c2)
        print 'dtw: c1 x c2 = ', dist

        dist = dtw_gradient(c1,c3)
        print 'dtwg: c1 x c3 = ', dist

        dist = dtw(c1,c3)
        print 'dtw: c1 x c3 = ', dist

    def est_05_representant(self):

        f1 = '/home/hick/Documents/Mestrado/Research/Code/Experiments5/training/02463d452.lines'
        f2 = '/home/hick/Documents/Mestrado/Research/Code/Experiments5/training/02463d654.lines'
        f3 = '/home/hick/Documents/Mestrado/Research/Code/Experiments5/training/02463d456.lines'
        f4 = '/home/hick/Documents/Mestrado/Research/Code/Experiments5/training/02463d550.lines'

        curves = []
        curves.append(read_curves(f1))
        curves.append(read_curves(f2))
        curves.append(read_curves(f3))
        curves.append(read_curves(f4))

        # ax = plt.subplot(1,4,0)
        # ax.plot(curves[0], 'k-')
        # ax = plt.subplot(1,4,1)
        # ax.plot(curves[1], 'k-')
        # ax = plt.subplot(1,4,2)
        # ax.plot(curves[2], 'k-')
        # ax = plt.subplot(1,4,3)
        # ax.plot(curves[3], 'k-')
        # plt.show()

        curves = np.array(curves)

        rep = representant_curve(curves, 'cosine')
        # ax = plt.subplot(1,4,0)
        # ax.plot(rep, 'k-')
        # plt.plot(rep, 'k-')
        #
        # rep = representant_curve(curves, 'euclidean')
        # ax = plt.subplot(1,4,1)
        # ax.plot(rep, 'g-')
        #
        # rep = representant_curve(curves, 'dtw')
        # ax = plt.subplot(1,4,2)
        # ax.plot(rep, 'r-')
        #
        # rep = representant_curve(curves, 'dtwg')
        # ax = plt.subplot(1,4,3)
        # ax.plot(rep, 'b-')
        #
        # plt.show()

    def est_00_saliency(self):

        def d(x1, x2):

            # x1 = remove_non_face(x1)
            # x2 = remove_non_face(x2)

            x1 = smooth(x1, 5)
            x2 = smooth(x2, 5)

            return dist.cosine(x1, x2)

        saliency('/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/', d,
                 '/home/hick/Documents/Mestrado/Research/Code/Experiments5/saliency/')

        plt.plot(range(0,5))
        plt.show()

    def est_07_histogram(self):

        x = [5,7,2,2,4,6,7,8,54,4,3,2,3,4,6,7,7,4,3,3,2,2,3,4,5,6,7,8,9.34,5,3,2,2,45,6,7,65,4,3,22,3,4,5,6]

        h, x, y = histogram(x)

        plt.plot(x,h)
        plt.show()

    def est_08_weibull(self):

        data = np.array([22,12,3,4,43,23,56,43,23,58,29,39,47,9,89,78,98,78,83,37,29,34,33,25,37,24,27])
        h, x = np.histogram(data, 20)

        weib = stats.exponweib.fit(data, 1, 1)
        weib_curve = stats.exponweib.pdf(x[:20],  *weib)

        plt.plot(h)
        plt.plot(weib_curve)
        plt.show()

    def est_09_p2p(self):

        f1 = '/home/hick/Documents/Mestrado/Research/Code/Experiments5/training/02463d452.lines'
        f2 = '/home/hick/Documents/Mestrado/Research/Code/Experiments5/training/02463d654.lines'
        f3 = '/home/hick/Documents/Mestrado/Research/Code/Experiments5/training/04202d350.lines'

        c1 = read_curves(f1,curve_idx=9)
        c2 = read_curves(f2, curve_idx=9)
        c3 = read_curves(f3, curve_idx=9)

        dist = p2p_dist(c1,c2)
        print 'p2p: c1 x c2 = ', dist

        dist = dtw_gradient(c1,c2)
        print 'dtwg: c1 x c2 = ', dist

        dist = p2p_dist(c1,c3)
        print 'p2p: c1 x c3 = ', dist

        dist = dtw_gradient(c1,c3)
        print 'dtwg: c1 x c3 = ', dist

        print '\nwith non face removal\n'

        c1 = remove_non_face(c1)
        c2 = remove_non_face(c2)
        c3 = remove_non_face(c3)

        dist = p2p_dist(c1,c2)
        print 'p2p: c1 x c2 = ', dist

        dist = dtw_gradient(c1,c2)
        print 'dtwg: c1 x c2 = ', dist

        dist = p2p_dist(c1,c3)
        print 'p2p: c1 x c3 = ', dist

        dist = dtw_gradient(c1,c3)
        print 'dtwg: c1 x c3 = ', dist

    def est_10_remove_non_face(self):

        f1 = '/home/hick/Documents/Mestrado/Research/Code/Experiments5/training/02463d452.lines'

        c9 = read_curves(f1, curve_idx=9)
        c = remove_non_face(c9)

        plt.plot(c9, 'r-')
        plt.plot(c, 'g-')
        plt.show()

        c9 = read_curves(f1)
        c = remove_non_face(c9)

        plt.plot(c9, 'r-')
        plt.plot(c, 'g-')
        plt.show()

    def test_00_get_data_split(self):

        folder = '/home/hick/Documents/Mestrado/Research/Code/Experiments5/'
        train, test = load_data_split(''.join([folder, 'trainingFiles.txt']), ''.join([folder, 'testFiles.txt']))


if __name__ == '__main__':

    suite = TestLoader().loadTestsFromTestCase(PreProcessorsTests)
    TextTestRunner(verbosity=2).run(suite)