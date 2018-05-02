
# coding: utf-8

# In[9]:


import numpy as np
import pandas as pd
import scipy.integrate as integrate
from scipy.optimize import brentq as root
import math
import numpy as np
import scipy.special as scp
from scipy.special import iv
from vonMi


# In[11]:


def qvonmises(p, mu = 0 ,  kappa = None, from_ = None, tol = np.finfo(float).eps**0.6):
    epsilon = 10 * np.finfo(float).eps     ##epsilon is Python equivalent of .Machine$double.eps
    if (type(p) == int):
        p = np.array([p])
    elif (type(p) == float):
        p = np.array([p])
    else:
        p = np.array(p)
    if (np.any(p > 1)): 
        raise ValueError("p must be in [0,1]")
    elif (np.any(p < 0)):
        raise ValueError("p must be in [0,1]")

    if (pd.isnull(from_)): ##from is a keyword
        from_ = mu - np.pi
   
    n = p.size
    mu = (mu - from_)%(2 * np.pi)      ## from is a keyword    
    if (pd.isnull(kappa)): 
        raise ValueError("kappa must be provided")   
        
    def zeroPvonmisesRad(x, p, mu, kappa):
        if (np.isnan(x)):         
            y = np.nan              
        else: 
            integration = integrate.quad(lambda x: dvonmises(x, mu, kappa), 0, x)
            y = integration[0] - p         ##integration[0] will give the value
        return(y);
    
    value = np.repeat(np.nan, p.size)
    for i in range(p.size):
        try:
            value[i] = root(lambda x: zeroPvonmisesRad(x, p[i], mu, kappa), 0, 2 * np.pi - epsilon)
        except:
            pass
            if(p[i] < (10 * epsilon)):
                value[i] = 0
            elif (p[i] > (1 - 10 * epsilon)):
                value[i] = 2 * np.pi - epsilon         
    value += from_
    return(value)

