
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


def pvonmises(q, mu, kappa, tol = 1e-020):
    from_ = mu - np.pi
    mu = (mu - from_) % (2 * np.pi)
    if (type(q) == int):
        q = [q]
    if(type(q) == float):
        q =[q]
    q = np.mod(np.subtract(q, from_), (2 * np.pi))
    q = np.mod(q,(2 * np.pi))
    n = len(q)
    mu = mu % (2 * np.pi)
    def pvm_mu0(q, kappa, tol):
        flag = True
        p = 1
        sum_ = 0
        while (flag):
            term = (iv(p, kappa) * np.sin(np.multiply(q, p)))/p
            sum_ = sum_ + term
            p = p + 1
            if (abs(term) < tol):
                flag = False
        return(np.divide(q,(2 * np.pi)) + sum_/(np.pi * iv(0, kappa)))

    result = np.repeat(np.nan, n)
    if (mu == 0):
        for i in range(0,n):
            result[i] = pvm_mu0(q[i], kappa, tol)
    else:
        for i in range(0,n):
            if (q[i] <= mu):
                upper = (q[i] - mu) % (2 * np.pi)
                if (upper == 0):
                    upper = 2 * np.pi
                lower = (-mu) % (2 * np.pi)
                result[i] = pvm_mu0(upper, kappa, tol) - pvm_mu0(lower, kappa, tol)
            else:
                upper = q[i] - mu
                lower = mu % (2 * np.pi)
                result[i] = pvm_mu0(upper, kappa, tol) + pvm_mu0(lower, kappa, tol)
    return(result)

