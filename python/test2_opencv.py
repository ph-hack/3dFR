# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 18:24:14 2015

@author: hick
"""

import cv2, os
import matplotlib.pyplot as plt
import numpy as np

path = '/home/hick/Documents/Mestrado/Research/Code/'
toPath = '/home/hick/Documents/Mestrado/Research/Code/Experiments6/data/'

surf = cv2.xfeatures2d.SURF_create(50)

for root, dirs, files in os.walk(path):
    
    for file in files:
        
        if file.endswith('.jpg'):
            
            img = cv2.imread(path + file)
            
            kp, desc = surf.detectAndCompute(img, None)
            
            s = len(kp)
            coords = np.zeros((s,2), np.int)
            
            np.savetxt(topath + file + '_kp.txt', coords, '%d')
            np.savetxt(topath + file + '_desc.txt', desc, '%f')
            

#img = np.loadtxt(path + '04201d434.unholed.jpg.dat')
#img2 = cv2.imread(path + '04201d434.unholed.jpg')

#i = 0
#
#for p in kp:
#    
#    x = np.int(p.pt[0])
#    y = np.int(p.pt[1])
#    coords[i,:] = (x,y)
#    
#    cv2.circle(img2, (x,y), 4, (255,0,0), 1)
#    
#    i = i + 1
#    
#np.savetxt(path + 'keyPoints.txt', coords, '%d')
#
#plt.imshow(img2)
#plt.show()