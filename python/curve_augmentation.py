from unittest import TestLoader, TextTestRunner
from unittest.case import TestCase
from pre_processors import read_curves, smooth
from faces import Face
import numpy as np
import copy as cp
import matplotlib.pyplot as plt
import random
import math


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

RANGES = {

    'noise': {

        'sigma': [0.3, 0.4, 0.5, 0.6, 0.75], #, 0.82, 1., 1.3],
        'seed': range(10)
    },
    'rotation': {

        'theta': [3, 4, 5, 6, 7],# 8, 9, 10],
        'center': [0.25, 0.5, 0.75]
    }
}

def augment_training_set(train_x, train_y, min_count=-1, transformations_range=RANGES):

    classes = set(train_y)
    new_train_x = []

    if min_count == -1:

        min_count = np.min([train_y.count(c) for c in classes])

    train_y = np.array(train_y)
    train_x = np.array(train_x)

    for c in classes:

        x = np.where(train_y == c)[0]

        if len(x) >= min_count:

            random.seed(5)
            filtered = [t for t in train_x[random.sample(x, min_count)]]
            new_train_x.extend(filtered)

        else:

            diff_count = min_count - len(x)

            transformations = get_transformations(transformations_range, diff_count)

            for i in range(diff_count):

                s = x[random.randint(0, len(x)-1)]
                face = Face(train_x[s], 11)

                new_face = augment_faces([face], [transformations[i]])

                print transformations[i]

                new_face_file = '{}{}_{}.lines'.format(train_x[s].replace('{}.lines'.format(face.id), ''), face.id, i)
                new_face[-1].save_to_file(new_face_file)

                new_train_x.append(new_face_file)

            for i in range(len(x)):

                new_train_x.append(train_x[x[i]])

    return new_train_x, min_count

def get_transformations(d, n):

    sigmas = d['noise']['sigma']
    seeds = d['noise']['seed']
    thetas = d['rotation']['theta']
    centers = d['rotation']['center']

    max_transf = len(sigmas) * len(seeds) * len(thetas) * len(centers)

    # print 'max transformations = ', max_transf

    if n > max_transf:

        raise Exception('There is only {} possible combinations and {} were asked!'.format(max_transf, n))

    transformations = []
    random.seed(None)
    chosen = random.sample(range(max_transf), n)
    # print 'chosen = ', chosen

    for c in chosen:

        sigma = int(math.ceil(c/(len(centers) * len(thetas) * len(seeds)))) % len(sigmas)
        seed = int(math.ceil(c/(len(centers) * len(thetas)))) % len(seeds)
        theta = int(math.ceil(c/len(centers))) % len(thetas)
        center = int(c % len(centers))

        # print 'c =', c
        # print 'sigma =', sigma, ' seed =', seed, ' theta =', theta, ' center =', center

        t = Transformation()
        t.add_function('noise', sigma=sigmas[sigma], seed=seeds[seed])
        t.add_function('rotation', theta=thetas[theta], center=centers[center])

        transformations.append(t)

    return transformations


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

    def __repr__(self):

        return 'Transformation <{}>'.format(str(self))

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
    def rotation(x, theta=5, unit='degree', center=0.5):

        if unit == 'degree':

            theta = np.pi * theta/180.

        rot_matrix = np.array([[np.cos(theta), -np.sin(theta)],
                               [np.sin(theta),  np.cos(theta)]])

        if len(x.shape) != 2:

            x = np.vstack((x, np.array(range(0,len(x))))).T

        # if center == -1:
        #
        #     center = x.shape[0]/2

        center = x.shape[0] * center

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
        c2 = Transformation.rotation(c1, 3, 'degree', 0.3)

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
        t2.add_function('rotation', theta=4, center=0.3)

        t3 = Transformation()
        t3.add_function('noise')
        t3.add_function('rotation', theta=-4, center=0.7)

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
        t2.add_function('rotation', theta=4, center=0.3)

        t3 = Transformation()
        t3.add_function('noise', sigma=1.3)
        t3.add_function('rotation', theta=-4, center=0.7)

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
        t2.add_function('rotation', theta=4, center=0.3)

        print 't1 = ', str(t1)
        print 't2 = ', str(t2)


if __name__ == '__main__':

    suite = TestLoader().loadTestsFromTestCase(AugmentationTests)
    TextTestRunner(verbosity=2).run(suite)