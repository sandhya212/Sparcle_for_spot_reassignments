## Code author SP

###########
### Python libraries
############

# uncomment below two rows when running code for the first time. This is to create these folders in the cwd.
#os.mkdir(path_python+'/count_data') 
#os.mkdir(path_python+'/figures_python')

import sys

import numpy as np
np.set_printoptions(threshold=sys.maxsize)
import cv2

import time

from matplotlib.patches import Polygon

from PIL import Image, ImageDraw
from skimage import measure #package is scikit-image


import pandas as pd

import seaborn as sns
import fastcluster

import umap

import os

from numpy.lib.stride_tricks import as_strided

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
plt.close('all')

import random
from scipy.stats import zscore

import scipy.stats as stats

from joblib import Parallel, delayed
import multiprocessing
import colorsys

from scipy.spatial import distance

from scipy.stats import multivariate_normal
import phenograph as pg

from math import copysign, hypot

from sklearn.mixture import BayesianGaussianMixture
from sklearn.manifold import TSNE

##############function definitions################
##################################################

def connectpoints(x,y,p1,p2):
    x1, x2 = x[p1], x[p2]
    y1, y2 = y[p1], y[p2]
    plt.plot([x1,x2],[y1,y2],'k-',linewidth=0.5)
    


def setdiff2d(lst1,lst2):
    aset = set([tuple(x) for x in lst1])
    bset = set([tuple(x) for x in lst2])
    l = np.array([x for x in aset - bset])
    return l

def is_nan(x):
    return (x is np.nan or x != x)

  
### D.Bartel: https://github.com/danbar/qr_decomposition

"""Module for Givens rotation."""


def givens_rotation(A):
    """Perform QR decomposition of matrix A using Givens rotation."""
    (num_rows, num_cols) = np.shape(A)

    # Initialize orthogonal matrix Q and upper triangular matrix R.
    Q = np.identity(num_rows)
    R = np.copy(A)

    # Iterate over lower triangular matrix.
    (rows, cols) = np.tril_indices(num_rows, -1, num_cols)
    for (row, col) in zip(rows, cols):

        # Compute Givens rotation matrix and
        # zero-out lower triangular matrix entries.
        if R[row, col] != 0:
            (c, s) = _givens_rotation_matrix_entries(R[col, col], R[row, col])

            G = np.identity(num_rows)
            G[[col, row], [col, row]] = c
            G[row, col] = s
            G[col, row] = -s

            R = np.dot(G, R)
            Q = np.dot(Q, G.T)

    return (Q, R)


def _givens_rotation_matrix_entries(a, b):
    """Compute matrix entries for Givens rotation."""
    r = hypot(a, b)
    c = a/r
    s = -b/r

    return (c, s)
  

# perspective projection of stage coord to image coord


pts_src = np.array([[-3357.4,-4078.2], [-3357.4,-3855.3], [-3134.5,-4078.2], [-3134.5,-3855.3]])
# corresponding points from image 2 (i.e. (154, 174) matches (212, 80))
pts_dst = np.array([[0.0,0.0], [0.0,2048.0], [2048.0,0.0],[2048.0,2048.0]])


# calculate matrix H
h, status = cv2.findHomography(pts_src, pts_dst)




###############Initialise variables #####################
#########################################################

colours=[ "darkviolet", "red", "orange", "limegreen", "blue", "purple", "seagreen","gold","lightpink","thistle","mistyrose","saddlebrown","slategrey",
            "palevioletred","mediumvioletred","yellowgreen","darkolivegreen","lemonchiffon","chocolate","lightsalmon",
            "lightcyan","lightblue", "black","cyan","yellow","burlywood","khaki","goldenrod","darkkhaki","slateblue","lightcoral","tomato","sienna","rosybrown",
            "tan","deepskyblue","mediumaquamarine","mediumturquoise","mediumorchid","violet","olivedrab","teal","lawngreen","orangered"]

bayes_MM = 0 # to run the DPMM else set to 0 for phenograph

num_comp = 10 #number of expected clusters if using Bayesian finite mixture model

neigh_size = 15 #  number of neighbours for the mock cell 

dist_ub = 150 #number of neighbours lie within this upper bound

dpi_set = 100
num_pg_k = 50 #number of PG clusters, optional
num_fov = 2 #change to total number of FoVs 

#random.seed(9001)

mRNA_asgn_counter = np.zeros([num_fov,1])

s_size = 20

num_colour = 25
win_size = 150 #vary between 70 - 150 to construct 'mock' cell

cm = colours

f_size=10

fig2p = 0

z = 2
num_cell_fov = np.zeros(num_fov, dtype=int)

asgn_stats_list=[] 
asgn_stats_listoflists = [] 
mRNA_assign_listoflists = [] #1st July 2020

iterations = 4 # was 2
num_cores = multiprocessing.cpu_count() - 1
inputs = range(num_fov)

Merfish_genes = 140

givens = True
wgt_mock_cell = True

exec(open(os.path.join(path_python+"/Sparcle.py")).read())
