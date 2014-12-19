from itertools import izip
import os

   
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3, venn3_circles

import seaborn as sns

import screed

from peasoup import configtools
from peasoup.plot import FigManager
from peasoup.plot import plot_dendro

from collections import defaultdict as ddict

import pandas as pd
from pandas import Categorical
pd.set_option('display.max_rows', 200)

import numpy as np
from numpy.random import rand

from scipy import stats 
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
import scipy.cluster.hierarchy as sch

from sklearn.decomposition import PCA, KernelPCA
from sklearn.cluster import KMeans

from IPython.display import SVG
from IPython.display import FileLink
from IPython.display import HTML

import mygene
from bioservices import UniProt

from pprint import pprint

def wdir(fn):
    return os.path.join('../../work/', fn)

outfmt6 = ['sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', \
           'qend', 'sstart', 'send', 'evalue', 'bitscore']

metadata = configtools.get_cfg('../../metadata.ini', '../../metadata.spec.ini')

import warnings
warnings.filterwarnings("ignore",category=DeprecationWarning)
warnings.filterwarnings("ignore",category=UserWarning)
