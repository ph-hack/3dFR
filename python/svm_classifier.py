from unittest import TestCase, TestLoader, TextTestRunner

__author__ = 'hick'

import sklearn.svm as svm
from pre_processors import get_person_id, smooth
from faces import Face
import numpy as np
import curve_augmentation as aug
import copy as cp

class SVM_Classifier:

    def __init__(self, augmentations=None):

        self.classifiers = []
        self.classes = []
        self.augmentations = augmentations

    def fit(self, X, y):

        if y is None:

            y = [get_person_id(s) for s in X]

        print "getting the data"
        faces = [Face(x, 11) for x in X]

        if self.augmentations is not None:

            faces = aug.augment_faces(faces, self.augmentations)
            y = [f.person for f in faces]

        data = []

        for f in faces:

            s = []
            for c in f.curves:

                # s.extend(smooth(c, 5).tolist())
                s.extend(c)


            data.append(s)

        self.classes = list(set(y))

        for c in self.classes:

            yc = [1 if i == c else 0 for i in y]

            print {0:0.5, 1:len(yc)/np.sum(yc)}

            classifier = svm.SVC(class_weight={0:0.5, 1:len(yc)/np.sum(yc)}, kernel='linear')
            classifier.fit(data, yc)

            self.classifiers.append(classifier)

            print "classifier ", c, " trained"

    def decision_function(self, X):

        faces = [Face(x, 11) for x in X]

        data = []

        for f in faces:

            s = []
            for c in f.curves:

                # s.extend(smooth(c, 5).tolist())
                s.extend(c)

            data.append(s)

        pres = []

        for c in self.classifiers:

            pres.append(c.decision_function(data))

        pres = np.array(pres)
        des = np.argmax(pres, 1)

        return self.get_classes(des)

    def get_classes(self, x):

        classes = []

        for i in x:

            classes.append(cp.copy(self.classes[i]))

        return classes


class SvmTests(TestCase):

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

        classifier = SVM_Classifier()
        classifier.fit(x_train, y_train)

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

        augmentations = []
        augmentations.append(aug.Transformation())
        augmentations[0].add_function('noise', seed=2)
        augmentations.append(aug.Transformation())
        augmentations[1].add_function('noise', sigma=1)
        augmentations[1].add_function('rotation', theta=4)
        augmentations.append(aug.Transformation())
        augmentations[2].add_function('noise', sigma=1.5, seed=3)
        augmentations.append(aug.Transformation())
        augmentations[3].add_function('noise', seed=4)

        classifier = SVM_Classifier(augmentations=augmentations)
        classifier.fit(x_train, y_train)

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

    suite = TestLoader().loadTestsFromTestCase(SvmTests)
    TextTestRunner(verbosity=2).run(suite)