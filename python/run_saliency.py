import pre_processors as pre
import distances as dist
from scipy.spatial.distance import cosine

pwd = '../../Experiments6/'
lines_path = '../../Experiments6/lines/'
saliency_path = '../../Experiments6/saliency/'

train_x, test_x = pre.load_data_split('{}trainingFiles.txt'.format(pwd), '{}testFiles.txt'.format(pwd))

train_y = [pre.get_person_id(t) for t in train_x]
test_y = [pre.get_person_id(t) for t in test_x]

train_x, N = pre.filter_training_set(train_x, train_y)
print 'train =', train_x[:2]
train_x = ['{}lines/{}.lines'.format(pwd, t) for t in train_x]
print 'train =', train_x[:2]
train_y = [pre.get_person_id(t) for t in train_x]

# test_x = pre.filter_testing_set(test_x, set(train_y))
test_x = ['{}lines/{}.lines'.format(pwd, t) for t in test_x]

print 'samples per class =', N
print 'total samples =', len(train_x), ' = ', len(train_y)

def d(x,y):

  # x = pre.remove_non_face(x)
  # y = pre.remove_non_face(x)
  
  x = pre.smooth(x, 5)
  y = pre.smooth(y, 5)

  # return dist.dtw_gradient(x,y) 
  return cosine(x,y)

pre.saliency(train_x, d, saliency_path)
