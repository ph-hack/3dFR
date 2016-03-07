import pre_processors as pre
import hierarchical_classifier as cla
import sklearn.metrics as metrics
from scipy.spatial.distance import cosine
from curve_augmentation import augment_training_set
from faces import Face

f1 = 'slfda/adsfaf/adfaf/0001d00_0.lines'
f2 = 'slfda/adsfaf/adfaf/0001d00.lines'

print 'f1 = ', pre.get_augmentation_id(f1), ' should be 0'
print 'f2 = ', pre.get_augmentation_id(f2), ' should be orig'

d = {
    'a': 10.,
    'b': 1.5,
    'c': 4.5,
    'd': 0.5,
    'e': 12.
}

print cla.top_candidates(d, 3)
print cla.top_candidates(d, 3, 'desc')