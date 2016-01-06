from unittest import TestLoader, TextTestRunner
from unittest.case import TestCase
from pre_processors import get_sample_id, get_person_id, read_curves, read_saliency, smooth
from distances import dtw, p2p_dist, dtw_gradient
from faces import Face
import numpy as np
import scipy.spatial.distance as dist
import copy as cp
import matplotlib.pyplot as plt


def augment_faces(faces, transformations):

    aug_faces = []

    for f in faces:

        aug_faces.append(f)

        for t in transformations:

            face = cp.deepcopy(f)

            for c in range(0, face.n_curves):

                face.curves[c] = t.apply(face.curves[c])

            face.id = ''.join([face.id, ': ', str(t)])
            aug_faces.append(face)

    return aug_faces


class Transformation:

    def __init__(self):

        self.functions = []
        self.params = []
        self.names = []

    def __str__(self):

        text = ''

        for f in range(0, len(self)):

            text = ''.join([text, self.names[f], ' ', str(self.params[f]), ' '])

        return text

    def __len__(self):

        return len(self.functions)

    def add_function(self, f, **kwargs):

        transforms = {

            'noise': Transformation.noise,
            'rotation': Transformation.rotation
        }

        self.functions.append(Transformation._make_function(transforms[f], kwargs))
        self.params.append(kwargs)
        self.names.append(f)

    def apply(self, x):

        result = cp.deepcopy(x)

        for f in self.functions:

            result = f(result)

        return result

    @staticmethod
    def _make_function(f, param):

        def t(x):

            return f(x, **param)

        return t

    @staticmethod
    def noise(x, sigma=0.3, seed=1):

        np.random.seed(seed)
        noise = np.random.normal(0, sigma, len(x))
        return x + noise

    @staticmethod
    def rotation(x, theta=5, unit='degree', center=-1):

        if unit == 'degree':

            theta = np.pi * theta/180.

        rot_matrix = np.array([[np.cos(theta), -np.sin(theta)],
                               [np.sin(theta),  np.cos(theta)]])

        if len(x.shape) != 2:

            x = np.vstack((x, np.array(range(0,len(x))))).T

        if center == -1:

            center = x.shape[0]/2

        x[:,1] = x[:,1] - center

        rotated_x = np.dot(x, rot_matrix)

        return rotated_x[:,0]


class AugmentationTests(TestCase):

    def test_01_noise(self):

        c1 = read_curves('/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04202d350.lines')
        c2 = Transformation.noise(c1)
        c3 = Transformation.noise(c1, seed=2)
        c4 = Transformation.noise(c1, seed=3)

        plt.plot(c1, 'r-')
        plt.plot(c2, 'b-')
        plt.plot(c3, 'g-')
        plt.plot(c4, 'k-')
        plt.show()

    def test_02_rotation(self):

        c1 = read_curves('/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04202d350.lines')
        c2 = Transformation.rotation(c1, 3, 'degree', 30)

        plt.plot(c1, 'r-')
        plt.plot(c2, 'b-')
        plt.show()

    def test_03_transformations(self):

        c1 = read_curves('/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04202d350.lines')

        t1 = Transformation()
        t1.add_function('noise')
        t1.add_function('rotation')

        t2 = Transformation()
        t2.add_function('noise')
        t2.add_function('rotation', theta=4, center=30)

        t3 = Transformation()
        t3.add_function('noise')
        t3.add_function('rotation', theta=-4, center=70)

        c2 = t1.apply(c1)
        c3 = t2.apply(c1)
        c4 = t3.apply(c1)

        plt.plot(c1, 'r-')
        plt.plot(c2, 'b-')
        plt.plot(c3, 'g-')
        plt.plot(c4, 'k-')
        plt.show()

    def test_04_faces_aug(self):

        f1 = Face('/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04202d350.lines')

        faces = [f1]

        t1 = Transformation()
        t1.add_function('noise', sigma=1)
        t1.add_function('rotation')

        t2 = Transformation()
        t2.add_function('noise')
        t2.add_function('rotation', theta=4, center=30)

        t3 = Transformation()
        t3.add_function('noise', sigma=1.3)
        t3.add_function('rotation', theta=-4, center=70)

        transformations = [t1, t2, t3]

        aug_faces = augment_faces(faces, transformations)

        print 'faces = ', len(faces), '\naug_faces = ', len(aug_faces)

        plt.plot(smooth(aug_faces[0].curves[3], 5), 'r-')
        plt.plot(smooth(aug_faces[1].curves[3], 5), 'b-')
        plt.plot(smooth(aug_faces[2].curves[3], 5), 'g-')
        plt.plot(smooth(aug_faces[3].curves[3], 5), 'k-')
        plt.show()

    def test_00_str(self):

        f1 = Face('/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04202d350.lines')

        faces = [f1]

        t1 = Transformation()
        t1.add_function('noise', sigma=1)
        t1.add_function('rotation')

        t2 = Transformation()
        t2.add_function('noise')
        t2.add_function('rotation', theta=4, center=30)

        print 't1 = ', str(t1)
        print 't2 = ', str(t2)


if __name__ == '__main__':

    suite = TestLoader().loadTestsFromTestCase(AugmentationTests)
    TextTestRunner(verbosity=2).run(suite)