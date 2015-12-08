import numpy as np
import re
import scipy.spatial.distance as dist
import cv2
import matplotlib.pyplot as plt
from distances import dtw, dtw_gradient
from unittest.case import TestCase
from unittest.loader import TestLoader
from unittest.runner import TextTestRunner


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

def get_person_id(file_path):

    m = re.search('[0-9]+d[0-9]+', file_path)
    name = m.group(0)

    id = name.split('d')[0]

    return id

def representant_curve(curves, distance='cosine'):

    distances = {

        'cosine': dist.cosine,
        'euclidean': dist.euclidean,
        'dtw': dtw,
        'dtwg': dtw_gradient
    }

def smooth(x, n):

    for i in range(0, len(x)):

        l = 0 if i - n <= 0 else i - n
        r = len(x) if i + n >= len(x) else i + n

        w = slice(i,r)

        x[i] = np.mean(x[w])

    return x


# UNIT TESTS ##########################################################################################################

class PreProcessorsTests(TestCase):

    def test_1_read_curves(self):

        curve_file = '/home/hick/Documents/Mestrado/Research/Code/Experiments5/training/02463d452.lines'

        curve = read_curves(curve_file, curve_idx=0)

        self.assertEqual(len(curve), 130)

    def test_2_person_id(self):

        id = get_person_id('alsdf/asdfad/0200d345.lines')

        self.assertEqual(id, '0200')

    def test_0_dtw(self):

        f1 = '/home/hick/Documents/Mestrado/Research/Code/Experiments5/training/02463d452.lines'
        f2 = '/home/hick/Documents/Mestrado/Research/Code/Experiments5/training/02463d654.lines'
        f3 = '/home/hick/Documents/Mestrado/Research/Code/Experiments5/training/04202d350.lines'

        c1 = read_curves(f1)
        c2 = read_curves(f2)
        c3 = read_curves(f3)

        dist, D, _ = dtw_gradient(c1,c2)
        print 'dtwg: c1 x c2 = ', dist

        dist, D, _ = dtw(c1,c2)
        print 'dtw: c1 x c2 = ', dist

        dist, D, _ = dtw_gradient(c1,c3)
        print 'dtwg: c1 x c3 = ', dist

        dist, D, _ = dtw(c1,c3)
        print 'dtw: c1 x c3 = ', dist

        # plt.imshow(D)
        # plt.show()

    def test_4_smooth(self):

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

        dist, D, _ = dtw_gradient(c1,c2)
        print '\ndtwg: c1 x c2 = ', dist

        dist, D, _ = dtw(c1,c2)
        print 'dtw: c1 x c2 = ', dist

        dist, D, _ = dtw_gradient(c1,c3)
        print 'dtwg: c1 x c3 = ', dist

        dist, D, _ = dtw(c1,c3)
        print 'dtw: c1 x c3 = ', dist

if __name__ == '__main__':

    suite = TestLoader().loadTestsFromTestCase(PreProcessorsTests)
    TextTestRunner(verbosity=2).run(suite)