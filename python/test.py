import pre_processors as pre
import hierarchical_classifier as cla
import sklearn.metrics as metrics
from scipy.spatial.distance import cosine
from curve_augmentation import augment_training_set
from faces import Face

pwd = './ut/'

# train_x, test_x = pre.load_data_split('{}trainingFiles.txt'.format(pwd), '{}testFiles.txt'.format(pwd))

train_x = [
    '/home/phack/Documents/Mestrado/ut/02463d546.lines',
    '/home/phack/Documents/Mestrado/ut/02463d548.lines',
    '/home/phack/Documents/Mestrado/ut/02463d550.lines',
    '/home/phack/Documents/Mestrado/ut/02463d552.lines',
    '/home/phack/Documents/Mestrado/ut/02463d666.lines',
    '/home/phack/Documents/Mestrado/ut/04213d241.lines',
    '/home/phack/Documents/Mestrado/ut/04217d331.lines',
    '/home/phack/Documents/Mestrado/ut/04217d333.lines',
    '/home/phack/Documents/Mestrado/ut/04217d335.lines'
]

train_y = [pre.get_person_id(t) for t in train_x]

train_x, N = augment_training_set(train_x, train_y, 3)

print 'train =', train_x

# train_x = ['{}lines/{}.lines'.format(pwd, t) for t in train_x]

# print 'train =', train_x

train_y = [pre.get_person_id(t) for t in train_x]

print 'classes =', train_y

print 'samples per class =', N
print 'total samples =', len(train_x), ' = ', len(train_y)

def dtw_gradient(x,y):

  # x = pre.remove_non_face(x)
  # y = pre.remove_non_face(y)

  x = pre.smooth(x,5)
  y = pre.smooth(y,5)

  return cosine(x,y) #dist.dtw_gradient(x,y)

f1 = Face(train_x[4], 11, metric=dtw_gradient)
f2 = Face(train_x[5], 11, metric=dtw_gradient)
f3 = Face(train_x[0], 11, metric=dtw_gradient)

f1.compare_show(f2, measure=True)
f1.compare_show(f3, measure=True)