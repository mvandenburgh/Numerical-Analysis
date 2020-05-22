import numpy as np


r = lambda theta: np.sin(2*theta)

def toss_needle(L):
    center = (np.random.uniform(-1,1), np.random.uniform(-1, 1))
    np.arctan2