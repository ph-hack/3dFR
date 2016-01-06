from unittest import TestLoader, TextTestRunner
from unittest.case import TestCase
from pre_processors import get_sample_id, get_person_id, read_curves, read_saliency, smooth, remove_non_face
from distances import dtw, p2p_dist, dtw_gradient
from faces import Face
import numpy as np
import scipy.spatial.distance as dist
import curve_augmentation as aug


class HierarchicalClassifier:

    def __init__(self, level=1, n_curves=11, distance=dtw, saliency_folder=None, top=11, augmentations=None):

        self.level = level
        self.root = HierarchicalNode()
        self.n_curves = n_curves
        self.saliency_folder = saliency_folder
        self.top = top
        self.distance = distance
        self.augmentations = augmentations

    def fit(self, X, y=None):

        if y is None:

            y = [get_person_id(s) for s in X]

        # makes all samples nodes
        leafs = self._make_nodes(X)
        last_level = leafs

        for l in range(1, self.level):

            pass

        for n in last_level:

            self.root.add_child(n)

    def predict(self, X):
        pass

    def decision_function(self, X):

        decision = []

        for x in X:

            face = Face(x, self.n_curves)

            node_queue = [self.root]

            candidates = {}

            while(len(node_queue) > 0):

                current_node = node_queue.pop(0)

                if current_node.isleaf:

                    error = current_node.sample - face

                    if current_node.sample.person not in candidates:

                        candidates[current_node.sample.person] = error

                    else:

                        candidates[current_node.sample.person] = error \
                            if error < candidates[current_node.sample.person] \
                            else candidates[current_node.sample.person]

                elif current_node.test(face):

                    node_queue.extend(current_node.children)

            print '\ncandidates: ', candidates,

            chosen = elect(candidates)

            decision.append(chosen[0])

        return decision

    def _make_nodes(self, X):

        nodes = []
        for x in X:

            face = Face(x, self.n_curves, saliency_folder=self.saliency_folder, top=self.top, metric=self.distance)

            #applies data augmentation and creates the new nodes
            if self.augmentations is not None:

                faces = aug.augment_faces([face], self.augmentations)

                for f in faces:

                    nodes.append(HierarchicalNode(f))

            else:
                nodes.append(HierarchicalNode(face))

        return nodes


class HierarchicalNode:

    def __init__(self, sample=None):

        self.sample = sample
        self.children = []
        self.parent = None
        self.stats = stats_template()

    def __str__(self):

        return str(self.sample)

    def __repr__(self):

        root = str(self.sample)
        parent = str(self.parent)
        children = str([s.sample for s in self.children])
        brothers = str(self.get_brothers())

        return ''.join(['root = ', root, '\nparent = ', parent, '\nbrothers = ', brothers, '\nchildren = ', children,
                        '\nstats:\n\tmax = ', str(self.stats['max']), '\n\tmean = ', str(self.stats['mean']),
                        '\n\tvar = ', str(self.stats['var'])])

    def __len__(self):

        return len(self.children)

    def get_brothers(self):

        if self.parent is None:

            return []

        candidates_brothers = self.parent.children
        brothers = []

        for b in candidates_brothers:

            if b.sample != self.sample:

                brothers.append(b.sample)

        return brothers

    def add_child(self, child):

        self.children.append(child)
        child.parent = self

    def test(self, face):

        if self.stats['max'] == 0 and self.stats['mean'] == 0 and self.stats['var'] == 0:

            return True

        error = self.sample - face

        return error <= self.stats['max']

    @property
    def isleaf(self):

        return len(self.children) == 0

    @property
    def isroot(self):

        return self.parent is None


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


class HierarchicalTests(TestCase):

    def est_01_node(self):

        n1 = HierarchicalNode([1,2,3,4])
        n2 = HierarchicalNode([5,6,7])
        n3 = HierarchicalNode([8,9,10])

        n1.add_child(n2)
        n1.add_child(n3)

        reprn1 = 'root = [1, 2, 3, 4]\nparent = None\nbrothers = []\nchildren = [[5, 6, 7], [8, 9, 10]]\nstats:\n\tmax = 0\n\tmean = 0\n\tvar = 0'
        reprn2 = 'root = [5, 6, 7]\nparent = [1, 2, 3, 4]\nbrothers = [[8, 9, 10]]\nchildren = []\nstats:\n\tmax = 0\n\tmean = 0\n\tvar = 0'

        # print 'n1:\n', repr(n1), '\n'
        # print 'n2:\n', repr(n2), '\n'

        self.assertEqual(repr(n1), reprn1)
        self.assertEqual(repr(n2), reprn2)

        self.assertTrue(n2.isleaf)
        self.assertFalse(n1.isleaf)
        self.assertTrue(n1.isroot)

    def est_02_classifier(self):

        x_train = [
            '/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04202d350.lines',
            '/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/02463d452.lines',
            '/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04201d302.lines',
            '/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04203d346.lines',
            '/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04211d341.lines',
            '/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04213d344.lines'
        ]
        y_train = ['04202', '02463', '04201', '04203', '04211', '04213']

        def smoothed_cosine(x,y):

            x = smooth(x,5)
            y = smooth(y,5)

            return dist.cosine(x,y)

        classifier = HierarchicalClassifier(top=4, distance=smoothed_cosine,
                                            saliency_folder='/home/hick/Documents/Mestrado/Research/Code/Experiments5/saliency/')

        classifier.fit(x_train, y_train)

        print repr(classifier.root)

        x_test = [

            '/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04202d348.lines',
            '/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/02463d456.lines',
            '/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04201d304.lines',
            '/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04203d348.lines',
            '/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04211d343.lines',
            '/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04213d348.lines'
        ]

        decision = classifier.decision_function(x_test)

        print 'decision = ', decision

    def test_00_classifier_augmentation(self):

        x_train = [
            '/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04202d350.lines',
            '/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/02463d452.lines',
            '/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04201d302.lines',
            '/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04203d346.lines',
            '/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04211d341.lines',
            '/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04213d344.lines'
        ]
        y_train = ['04202', '02463', '04201', '04203', '04211', '04213']

        def smoothed_cosine(x,y):

            x = remove_non_face(x)
            y = remove_non_face(y)

            x = smooth(x,5)
            y = smooth(y,5)

            #return dist.cosine(x,y)
            return dtw_gradient(x,y)

        augmentations = []
        augmentations.append(aug.Transformation())
        augmentations[0].add_function('noise', seed=2)
        augmentations.append(aug.Transformation())
        augmentations[1].add_function('noise', sigma=1)
        augmentations.append(aug.Transformation())
        augmentations[2].add_function('noise', sigma=1.5, seed=3)
        augmentations.append(aug.Transformation())
        augmentations[3].add_function('noise', seed=4)

        classifier = HierarchicalClassifier(top=6, distance=smoothed_cosine
                                            #, saliency_folder='/home/hick/Documents/Mestrado/Research/Code/Experiments5/saliency/'
                                            )#,augmentations=augmentations)

        classifier.fit(x_train, y_train)

        print 'root lenght: ', len(classifier.root)
        print repr(classifier.root)

        x_test = [

            '/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04202d348.lines',
            '/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/02463d456.lines',
            '/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04201d304.lines',
            '/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04203d348.lines',
            '/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04211d343.lines',
            '/home/hick/Documents/Mestrado/Research/Code/Experiments5/from_saliency/04213d348.lines'
        ]

        decision = classifier.decision_function(x_test)
        print 'hahaha'
        print 'decision = ', decision


if __name__ == '__main__':

    suite = TestLoader().loadTestsFromTestCase(HierarchicalTests)
    TextTestRunner(verbosity=2).run(suite)