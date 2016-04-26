# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 21:24:12 2016

@author: mgou
"""

import pygrib
import numpy as np
import sys

grbs = pygrib.open('/home/mgou/uwhydro/summaProj/summaData/summa_columbia/annual_forcing/grib/NLDAS_FORA0125_H.A20130501.0000.002.grb')

a = [[1,2,3,4],[2,3,4,5]]
a[1][2]



dir(file)
dir(grbs[1])
grbs.name
grbs[1].latlons()
grbs[1].data()[0][10,100]
len(grbs[1]['values'])
grbs[1]['shortName']

cc = 'abcdefg'
cc[4]

cc=['abcde','fghij']
cc[1]