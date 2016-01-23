import pre_processors as pre
import distances as dist
import hierarchical_classifier as cla

pwd = '../../Experiments5/'

train_x, test_x = pre.load_data_split('{}/trainingFiles.txt'.format(pwd), '{}/testFiles.txt'.format(pwd))

train_x = ['{}/training/{}.lines'.format(pwd, t) for t in train_x]
test_x = ['{}/test/{}.lines'.format(pwd, t) for t in test_x]

train_y = [pre.get_person_id(t) for t in train_x]
test_y = [pre.get_person_id(t) for t in test_x]

def dtw_gradient(x,y):

  x = pre.remove_non_face(x)
  y = pre.remove_non_face(y)

  x = pre.smooth(x,5)
  y = pre.smooth(y,5)

  return dist.dtw_gradient(x,y)

classifier = cla.HierarchicalClassifier(top=6, distance=dtw_gradient)

classifier.fit(train_x[:200], train_y[:200])

decision = classifier.decision_function(test_x[33:35])

print 'decision = \n', decision
