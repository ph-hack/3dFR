
from skimage import data, io, filters
from skimage.feature import CENSURE
import matplotlib.pyplot as plt

image = io.imread("/home/hick/Documents/Mestrado/Research/Code/04201d444.unholed.jpg")
detector = CENSURE(1, 10, 'DoB', 0.005, 100)

fig, ax = plt.subplots(nrows=1, ncols=1)
plt.gray()

detector.detect(image)

ax.imshow(image)
ax.axis('off')
ax.scatter(detector.keypoints[:,1], detector.keypoints[:,0],
              2 ** detector.scales, facecolors='none', edgecolors='r')

plt.show()
