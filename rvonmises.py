
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


def rvonmises(n, mu, kappa):
    vm = np.zeros(n)
    a = 1 + (1 + 4 * (kappa**2))**0.5
    b = (a - (2 * a)**0.5)/(2 * kappa)
    r = (1 + b**2)/(2 * b)
    obs = 0
    while (obs < n):
        U1 = np.random.uniform(0, 1, 1)
        z = np.cos(np.pi * U1)
        f = (1 + r * z)/(r + z)
        c = kappa * (r - f)
        U2 = np.random.uniform(0, 1, 1)
        if (c * (2 - c) - U2 > 0):
            U3 = np.random.uniform(0, 1, 1)
            vm[obs] = np.sign(U3 - 0.5) * math.acos(f) + mu
            vm[obs] = vm[obs] % (2 * np.pi)
            obs = obs + 1
        else:
            if (math.log(c/U2) + 1 - c >= 0):
                U3 = np.random.uniform(0, 1, 1)
                vm[obs] = np.sign(U3 - 0.5) * math.acos(f) + mu
                vm[obs] = vm[obs] % (2 * math.pi)
                obs = obs + 1
    return(vm)

