
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import scipy.integrate as integrate
from scipy.optimize import brentq as root
import math
import numpy as np
import scipy.special as scp
from scipy.special import iv


# In[2]:


def dvonmises(x, mu, kappa, log = False):
    if (type(x) == int):
        x = [x]
    if (type(x) == float):
        x = [x]
    vm = np.zeros(len(x))
    if (log):
        if (kappa == 0):
            vm = np.log(np.repreat(1/(2*pi), len(x)))
        elif (kappa < 100000):
            vm = -(np.log(2*math.pi)+np.log(scp.ive(0, kappa)) + kappa) + kappa*(np.cos(np.subtract(x - mu)))
        else:
            if (((x-mu)%(2*math.pi))==0):
                vm = math.inf
            else:
                vm = -math.inf
    else:
        if (kappa == 0):
            vm = np.repeat(1/(2*np.pi), len(x))
        elif (kappa < 100000):
            vm = 1/(2 * np.pi * scp.ive(0, kappa)) * (np.exp(np.subtract(np.cos(np.subtract(x, mu)), 1)))**kappa
        else:
            if (np.mod(np.subtract(x, mu),(2*np.pi))==0):
                vm = math.inf
            else:
                vm = 0
    return(vm)

