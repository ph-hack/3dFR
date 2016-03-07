import pre_processors as pre
import hierarchical_classifier as cla
import sklearn.metrics as metrics
from scipy.spatial.distance import cosine
from curve_augmentation import augment_training_set

pwd = '../../Experiments6/'

train_x, test_x = pre.load_data_split('{}trainingFiles.txt'.format(pwd), '{}testFiles.txt'.format(pwd))

train_y = [pre.get_person_id(t) for t in train_x]
test_y = [pre.get_person_id(t) for t in test_x]

print 'train size before: ', len(train_x)

train_x = ['{}lines/{}.lines'.format(pwd, t) for t in train_x]
# train_x, N = pre.filter_training_set(train_x, train_y)
train_x, N = augment_training_set(train_x, train_y, 2)

train_y = [pre.get_person_id(t) for t in train_x]

# test_x = pre.filter_testing_set(test_x, set(train_y))
test_x = ['{}lines/{}.lines'.format(pwd, t) for t in test_x]
# test_y = [pre.get_person_id(t) for t in test_x]

print 'samples per class =', N
print 'total samples =', len(train_x), ' = ', len(train_y) 

def dtw_gradient(x,y):

  # x = pre.remove_non_face(x)
  # y = pre.remove_non_face(y)

  x = pre.smooth(x,5)
  y = pre.smooth(y,5)

  return cosine(x,y) #dist.dtw_gradient(x,y)

classifier = cla.HierarchicalClassifier(top=6, distance=dtw_gradient, log_file='classifier2.log')
             #saliency_folder='{}saliency/'.format(pwd))

classifier.fit(train_x, train_y)

decision = classifier.decision_function(test_x)

print metrics.classification_report(test_y, decision)
