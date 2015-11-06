# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 18:24:14 2015

@author: hick
"""

import cv2, os
import matplotlib.pyplot as plt
import numpy as np
import sys



#path = '/home/hick/Documents/Mestrado/Research/Code/Experiments6/meanFaces/images/'
#toPath = '/home/hick/Documents/Mestrado/Research/Code/Experiments6/meanFaces/'
path = sys.argv[1]
toPath = sys.argv[2]

surf = cv2.xfeatures2d.SURF_create(50)

for root, dirs, files in os.walk(path):
    
    print len(files), ' files'

    for file in files:
        
        if file.endswith('.jpg'):
            
            img = cv2.imread(path + file)
            
            kp, desc = surf.detectAndCompute(img, None)
            
            s = len(kp)
            coords = np.zeros((s,2), np.int)

            i = 0

            for p in kp:

                x = np.int(p.pt[0])
                y = np.int(p.pt[1])
                coords[i,:] = (x,y)

                i += 1
            
            np.savetxt(toPath + 'surf_kps/' + file + '_kp.txt', coords, '%d')
            np.savetxt(toPath + 'surf_features/' + file + '_desc.txt', desc, '%f')

            print 'saving ', file, ' data'
