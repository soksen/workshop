'''
import matplotlib.pyplot as plt
from scipy.fftpack import fft

import re
import os
import math
import random
from sklearn import svm

from scipy.fftpack import fftn
import numpy as np

from TEMPy.StructureParser import PDBParser
from TEMPy.StructureBlurrer import StructureBlurrer
from TEMPy.MapParser import MapParser

import csv
'''

print 'Importing...',
import prepare
import learn
print 'done'
working_directory = 'D:/university/biology'
prepare.prepare(working_directory)
learn.learn(working_directory)