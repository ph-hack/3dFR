from unittest import TestLoader, TextTestRunner
from unittest.case import TestCase
import math
from pre_processors import get_sample_id, get_person_id, read_curves, read_saliency, get_augmentation_id
from distances import dtw
import numpy as np
import matplotlib.pyplot as plt
import copy as cp


class Face:

    def __init__(self, fullpath, n=11, metric=dtw, saliency_folder=None, top=11):

        self.curves = []
        self.saliencies = []
        self.n_curves = 11 if n > 11 or n < 1 else n
        self.id = get_sample_id(fullpath)
        self.person = get_person_id(fullpath)
        self.metric = metric
        self.top = top
        self.augmentation = get_augmentation_id(fullpath)

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
        errors[np.where(errors < 0.)] = 0.

        saliencies = np.array(self.saliencies)

        if saliencies.shape[0] > 0:
            errors = errors/(errors + saliencies)

        errors.sort()

        return np.mean(errors[:self.top])

    def __repr__(self):

        return ''.join(['face ', self.id, '(', self.augmentation, ') with ', str(self.n_curves), ' curves'])

    def __str__(self):

        return '{}({})'.format(self.id, self.augmentation)

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
            # ax.plot(smooth(self.curves[c], 5), 'b-')
            ax.plot(self.curves[c], 'b-')
            # ax.plot(smooth(other.curves[c], 5), 'r-')
            ax.plot(other.curves[c], 'r-')

            if measure:

                error = self.metric(self.curves[c], other.curves[c])
                ax.set_title(str(error))

        plt.show()

    def save_to_file(self, file_path):

        try:
            with open(file_path, 'w') as f:

                curves = self._curves_to_strings()

                f.writelines(curves)

        except Exception as e:

            print 'Error writing the curves into the file: ', file_path, ':\n', e.message

    def apply(self, F, **kwargs):

        for c in range(len(self.curves)):

            self.curves[c] = F(self.curves[c], **kwargs)

    def _curves_to_strings(self):

        strings = []

        curves = []

        for i in range(self.n_curves):

            curves.append([str(c) for c in self.curves[i]])

        for i in range(self.n_curves-1):

            s = str(curves[i]).replace(', ', ' ')
            s = s.replace('[', '')
            s = s.replace(']', '')
            s = s.replace("'", '')

            s += '\n'

            strings.append(s)

        s = str(curves[-1]).replace(', ', ' ')
        s = s.replace('[', '')
        s = s.replace(']', '')
        s = s.replace("'", '')

        strings.append(s)

        return strings

    @staticmethod
    def mean(faces):

        M = cp.deepcopy(faces[0])

        M.id = '{}_mean'.format(M.person)

        for c in range(len(faces[0].curves)):

            curves = np.array([f.curves[c] for f in faces])

            M.curves[c] = curves.mean(axis=0)

        distances = np.array([M - f for f in faces])

        chosen = distances.argmin()

        return cp.deepcopy(faces[chosen])


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