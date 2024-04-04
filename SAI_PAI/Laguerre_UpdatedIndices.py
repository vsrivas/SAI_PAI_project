#!/usr/bin/env python
# coding: utf-8

# In[7]:


import numpy as np
import math


# In[8]:


# # define Laguerre function of order j
def phi_j(n, j, alpha):
    phi=[]
    sum_= 0
    for i in range(0, j+1):
        sum_ += (-1)**i*math.comb(n,i)*math.comb(j, i)*alpha**(j-i)*(1-alpha)**i
        phi = alpha**((n-j)/2)*(1-alpha)**(1/2)*sum_
    return np.array(phi) #.T


# In[9]:


# Define laguerre filter output at time t
def l_j(J, t, RR_clean, W_n, alpha):  
    phi_mat = np.array([[phi_j(n,j, alpha)  for j in range (J+1)] for n in range(W_n)])
    l_j_t = [np.sum([phi_mat[n, j]*RR_clean['RR'].values[t-n-2] for n in range(W_n)]) for j in range (J+1)]
    return np.array(l_j_t).T


# In[10]:


# Define laguerre filter output for the duration of ECG excluding the window W_n
def lagl_j(J, RR_clean, W_n, alpha):
    lagl_j = [l_j(J, t, RR_clean, W_n, alpha) for t in range(W_n, len(RR_clean)+1)]
    return lagl_j


# In[ ]:




