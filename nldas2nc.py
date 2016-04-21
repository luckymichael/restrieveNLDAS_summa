# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 21:24:12 2016

@author: mgou
"""

import pygrib
import numpy as np

grbs = pygrib.open('/home/mgou/summaProj/summaData/summa_columbia/annual_forcing/grib/NLDAS_FORA0125_H.A20130501.0000.002.grb')

grbs[1]['values']
grbs[1]['getNumberOfValues']