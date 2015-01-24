## richness_est.py
"""
This program is a python port of Chris Miller's richness estimator make_optical_rev1.pro.
It uses a multi-variable approach to estimate members of a galaxy cluster and determine its richness.
"""

## Import Modules ##
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as mp
import kcorrect as kc

## Flags ##


## Constants ##
c = 2.998e8		# Speed of Light, m/s
H0 = 100.0		# Hubble parameter, km/s / Mpc
h = H0 * 1.0		# small h


clus_num = 100		# Number of Clusters to Work Over


## Program ##

## Load Data






