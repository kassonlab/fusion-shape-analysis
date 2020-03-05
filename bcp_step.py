#!/usr/bin/python

import math
import numpy

def bayesian_step(data):
  """Bayesian changepoint detection.  Based on code by Frank Kelly."""
  num_samples = len(data)
  dbar = numpy.mean(data)
  dsbar = numpy.mean(numpy.square(data))
  fac = dsbar - numpy.square(dbar)

  summ = 0
  summup = []
  y = []

  for z in range(num_samples):
      summ += data[z]
      summup.append(summ)
  for m in range(num_samples - 1):
    pos = m + 1
    mscale = 4 * (pos) * (num_samples - pos)
    Q = summup[m] - (summ - summup[m])
    U = -numpy.square(dbar*(num_samples - 2*pos) + Q)/float(mscale) + fac
    y.append(-(num_samples / float(2) -1)
             * math.log(num_samples * U / 2)
             - 0.5 * math.log((pos * (num_samples - pos))))
  z, zz = numpy.max(y), numpy.argmax(y)
  mean1 = sum(data[:zz+ 1 ]) / float(len(data[:zz + 1]))
  mean2=sum(data[(zz + 1):num_samples]) / float(num_samples - 1 - zz)
  return y, zz, mean1, mean2


