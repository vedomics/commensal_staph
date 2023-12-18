#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Created on Thu Feb 16 11:56:35 2023
# This script was made to measure the Euclidean distance between two sets of x-y coordinates which represent the boundaries of the uppermost layer of skin in a representative image of a histology section
# For samples in which serocellular crust was present, the width of the serocellular crust was included between the "top" and "bottom" lines which delineated the uppermost skin layer
# When crust was not present, the uppermost layer of skin included only the epidermis. Samples included in analysis had either epidermis, serocellular crust or a combination of the two which comprised the uppermost skin layer

#@author: Laura Markey

import os
import pandas as pd
import numpy as np
from math import dist

#%%read in files and get list of images
# each sample is represented by two txt files:
# "sample_t.txt": containing the x-y coordinates for a line delineating the top boundary of the skin layer
# "sample_b.txt": containing the x-y coordinates for a line indicating the bottom boundary 
# code was written for x-y coordinate files output from ImageJ after manual drawing of lines as follows:
# 1) open TIF image of histology section
# 2) draw boundary line of top or bottom of skin layer using "freehand line" tool
# 3) get x-y coordinates by running the ImageJ macro getselectioncoordinates.ijm
# 4) save x-y coordinates to file

#navigate to directory that holds your coordinate files for top and bottom lines
os.chdir('/Users/liebermanlab/Dropbox (MIT)/Lieberman Lab/Projects/Commensal_memory_mouse/data_and_code/Epidermal_thickness/coordinates')#iterate through contents of directory and make a list of samples assuming all files that end with txt contain coordinates
dir=os.getcwd()
files=[]
for file in os.listdir(dir):
    if file.endswith(".txt"):
        files.append(file)
samples=np.unique([item[:-5] for item in files])
#measure shortest distance between two lines
thickness=[]
sample=[]
for m in samples:
    sample.append(m)
    toppoints=pd.read_csv(str(m)+"t.txt", sep=" ", header=None, index_col=0)
    bpoints=pd.read_csv(str(m)+"b.txt", sep=" ", header=None, index_col=0)
    #print(toppoints)
    #print(bpoints)
    top=toppoints.apply(lambda x: [x[1],x[2]],axis=1) #convert x and y columns into a list of coordinates
    bottom=bpoints.apply(lambda x: [x[1],x[2]],axis=1)
    points=[]#gets you a list of all distances between each point in top and each point in bottom in order
    mins=[] #shortest distance from a point in top to one of the points in bottom
    for i in top:
        upper=[]
        for t in bottom:
            d1=dist(i,t)
            upper.append(d1)
            points.append(d1)
        m=min(upper)
        mins.append(m)
    ave=np.mean(mins)
    #print(ave)
    thickness.append(ave)
    os.chdir(dir)
summary=pd.DataFrame({"sample":sample, "ave_thickness":thickness})
summary.to_csv('/Users/liebermanlab/Dropbox (MIT)/Lieberman Lab/Projects/Commensal_memory_mouse/data_and_code/Epidermal_thickness/output/testset_thickness.csv')