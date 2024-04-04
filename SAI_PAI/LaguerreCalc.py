#!/usr/bin/env python
# coding: utf-8

# In[1]:


import math
import numpy as np


# In[2]:


# # define Laguerre polynomial of order j with time intervals n
def phi_j(n, j, alpha):
    phi=[]
    sum_= 0
    for i in range(0, j+1):
        sum_ += (-1)**i*math.comb(n,i)*math.comb(j, i)*alpha**(j-i)*(1-alpha)**i
        phi = alpha**((n-j)/2)*(1-alpha)**(1/2)*sum_
    return np.array(phi) #.T


# In[3]:


# Calculating Laguerre filter output using laguerre polynomials and RR interval at time t in the given window (w_n)
def l_j(J,t, RR_clean, W_n, alpha):
    phi_mat = np.array([[phi_j(n,j, alpha)  for j in range (J)] for n in range(1, W_n)])
    l_j_t = [np.sum([phi_mat[n-1, j]*RR_clean['int_RR'].values[t-n] for n in range(1, W_n) if (t-n)>=0 ]) for j in range (J)]
    return np.array(l_j_t).T


# In[4]:


# function to calculate lagerre filter output at each time point
def lagl_j(J, RR_clean, W_n, alpha):
    lagl_j = [l_j(J, t, RR_clean, W_n, alpha) for t in range(1, len(RR_clean)+1)]
    return lagl_j


# In[ ]:




