import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Gaussian1DKernel
import csv

ch = {}

wl_min = 0.4
wl_max = 20.0
delta_wl = 0.005

eta_dichroic = 0.9
eta_grating = 0.5
qe = 0.55

file_suffix = '_v0.csv'

ch[0] = {'wl_min':0.55, 'wl_max': 1}
ch[1] = {'wl_min':1.95, 'wl_max': 3.9}
ch[2] = {'wl_min':3.9, 'wl_max': 7.8}

wl = np.arange(wl_min, wl_max, delta_wl)

d0_r = np.zeros(wl.size)
d0_t = np.zeros(wl.size)
d1_r = np.zeros(wl.size)
d1_t = np.zeros(wl.size)

idx = np.logical_and(wl >= ch[0]['wl_min'], wl <= ch[0]['wl_max'])
d0_r[idx] = eta_dichroic

idx = np.where(wl >= ch[0]['wl_max'])
d0_t[idx] = eta_dichroic

idx = np.logical_and(wl >= ch[1]['wl_min'], wl <= ch[1]['wl_max'])
d1_r[idx] = eta_dichroic

idx = np.logical_and(wl >= ch[2]['wl_min'], wl <= ch[2]['wl_max'])
d1_t[idx] = eta_dichroic

gaussian_kernel = Gaussian1DKernel(10)
d0_t = convolve(d0_t, gaussian_kernel)
d1_t = convolve(d1_t, gaussian_kernel)
d0_r = convolve(d0_r, gaussian_kernel)
d1_r = convolve(d1_r, gaussian_kernel)



plt.ion()
plt.subplot(211)
plt.plot(wl, d0_r)
plt.plot(wl, d1_r)
plt.plot(wl, d1_t)
plt.subplot(212)
plt.plot(wl, d0_t)
plt.plot(wl, d1_t)

with open('d0_re'+file_suffix, 'wb') as csvfile:
  writer = csv.writer(csvfile, delimiter=',')
  writer.writerows([ ['{:.3f}'.format(wl[i]), '{:.3f}'.format(d0_r[i]), '{:.3f}'.format(0.03)] for i in xrange(wl.size)])
with open('d0_te'+file_suffix, 'wb') as csvfile:
  writer = csv.writer(csvfile, delimiter=',')
  writer.writerows([ ['{:.3f}'.format(wl[i]), '{:.3f}'.format(d0_t[i]), '{:.3f}'.format(0.03)] for i in xrange(wl.size)])
with open('d1_re'+file_suffix, 'wb') as csvfile:
  writer = csv.writer(csvfile, delimiter=',')
  writer.writerows([ ['{:.3f}'.format(wl[i]), '{:.3f}'.format(d1_r[i]), '{:.3f}'.format(0.03)] for i in xrange(wl.size)])
with open('d1_te'+file_suffix, 'wb') as csvfile:
  writer = csv.writer(csvfile, delimiter=',')
  writer.writerows([ ['{:.3f}'.format(wl[i]), '{:.3f}'.format(d1_t[i]), '{:.3f}'.format(0.03)] for i in xrange(wl.size)])

with open('qe'+file_suffix, 'wb') as csvfile:
  writer = csv.writer(csvfile, delimiter=',')
  writer.writerows([ ['{:.3f}'.format(wl[i]),  '{:.3f}'.format(qe)] for i in xrange(wl.size)])
