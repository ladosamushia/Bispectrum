using PyCall

py""" 
import imp
readsnap = imp.load_source('readsnap', './readsnap.py')
readgadget = imp.load_source('readgadget', './readgadget.py')
import numpy as np

def read_quijote_snapshot(infile):
    return readgadget.read_block(infile, "POS ", [1])
"""
