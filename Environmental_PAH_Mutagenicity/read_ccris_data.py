# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 09:27:56 2019

@author: twsle
"""
#acknowledgements
#https://stackoverflow.com/a/54932071/10226776 for the CIR convert function
import os
import pandas as pd

from urllib.request import urlopen

def CIRconvert(ids):
    try:
        url = 'http://cactus.nci.nih.gov/chemical/structure/' + ids + '/smiles'
        ans = urlopen(url).read().decode('utf8')
        return ans
    except:
        return 'Did not work'
