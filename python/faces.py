from unittest import TestLoader, TextTestRunner
from unittest.case import TestCase
import math
from pre_processors import get_sample_id, get_person_id, read_curves, read_saliency, smooth
from distances import dtw
import numpy as np
import matplotlib.pyplot as plt


class Face:

    def __init__(self, fullpath, n=11, metric=dtw, saliency_folder=None, top=11):

        self.curves = []
        self.saliencies = []
        self.n_curves = 11 if n > 11 or n < 1 else n
        self.id = get_sample_id(fullpath)
        self.person = get_person_id(fullpath)
        self.metric = metric
        self.top = top

        for i in range(0, self.n_curves):

            self.curves.append(read_curves(fullpath, curve_idx=i))

        if saliency_folder is not None:
            for i in range(0, self.n_curves):
                self.saliencies.append(read_saliency(saliency_folder, sample_id=self.id, curve_idx=i))

    def __getitem__(self, item):

        return self.curves[item]

    def __sub__(self, other):

        if other.n_curves != self.n_curves:

            raise Exception('The number of curves should be the same!')

        errors = []
        for c in range(0, self.n_curves):

            errors.append(self.metric(self.curves[c], other.curves[c]))

        errors = np.array(errors)
        saliencies = np.array(self.saliencies)

        if saliencies.shape[0] > 0:
            errors = errors/(errors + saliencies)

        errors.sort()

        return np.mean(errors[:self.top])

    def __repr__(self):

        return ''.join(['face ', self.id, ' with ', str(self.n_curves), ' curves'])

    def __str__(self):

        return self.id

    def show(self, n=-1):

        if n == -1:
            n = len(self.curves)

        r = int(math.ceil(n/2.))

        for c in range(n):
            ax = plt.subplot(2,r,c+1)
            ax.plot(self.curves[c], 'k-')

        plt.show()

    def compare_show(self, other, n=-1, measure=False):

        if n == -1:
            n = len(self.curves)

        r = int(math.ceil(n/2.))

        for c in range(n):
            ax = plt.subplot(2,r,c+1)
            ax.plot(smooth(self.curves[c], 5), 'b-')
            ax.plot(smooth(other.curves[c], 5), 'r-')

            if measure:

                error = self.metric(self.curves[c], other.curves[c])
                ax.set_title(str(error))

        plt.show()


def stats_template():

    return {

        'max': 0,
        'mean': 0,
        'var': 0
    }

def elect(candidates):

    keys = candidates.keys()
    values = np.array(candidates.values())

    v = values.argmin()

    return (keys[v], values[v])


class FaceTests(TestCase):

    def est_01_faces(self):

        f1 = Face('/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/02463d452.lines',
                  saliency_folder='/home/hick/Documents/Mestrado/Research/Code/Experiments5/saliency/', top=5)

        f2 = Face('/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/02463d456.lines',
                  saliency_folder='/home/hick/Documents/Mestrado/Research/Code/Experiments5/saliency/', top=5)

        f3 = Face('/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04202d350.lines',
                  saliency_folder='/home/hick/Documents/Mestrado/Research/Code/Experiments5/saliency/', top=5)

        self.assertEqual(f1.id, '02463d452')
        self.assertEqual(f3.person, '04202')
        self.assertEqual(len(f3[0]), 130)

        error = f1 - f2
        print 'f1 - f2 = ', error

        error = f1 - f3
        print 'f1 - f3 = ', error

        f1 = Face('/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/02463d452.lines', top=4)
        f2 = Face('/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/02463d456.lines', top=4)
        f3 = Face('/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04202d350.lines', top=4)

        error = f1 - f2
        print 'f1 - f2 = ', error
        error = f1 - f3
        print 'f1 - f3 = ', error

    def test_00_show(self):

        f1 = Face('/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/02463d452.lines', top=5)
        f2 = Face('/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/02463d456.lines', top=5)
        f3 = Face('/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04202d350.lines', top=5)

        f1.show()
        f1.compare_show(f2)
        f1.compare_show(f3)
        f1.compare_show(f3, 1)


if __name__ == '__main__':

    suite = TestLoader().loadTestsFromTestCase(FaceTests)
    TextTestRunner(verbosity=2).run(suite)