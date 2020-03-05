#!/usr/bin/python
# Code 2017 by Peter Kasson

"""Analyze alpha-shapes connected components to get stalk formation times."""

import gflags
import json
import numpy
import sys
import bcp_step

def rolling_window(a, window):
    """Code from Erik Rigtorp for rolling window"""
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return numpy.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

def stalk_times(component_dict):
  """Calculate stalk formation times or stalk-free interval.
  Args:
    component_dict:  dictionary where values are [time, # components]
  Rets:
    stalktimes:  times after which stalks are formed
    nostalktimes:  intervals of stalk-free simulation
  """ 
  # require stalk_interval consecutive single components
  stalk_interval = '1 1 1 1 1'
  stalktimes = []
  nostalktimes = []
  for listdat in component_dict.values():
    dat = numpy.array(listdat)  # JSON array has to be as list, so convert
    firststalk = numpy.where(numpy.mean(rolling_window(dat[:, -1], 5), -1) == 1)[0]
    if not firststalk.any():
      nostalktimes.append(dat[-1, 0])
    else:
      # find stalk
      bcp_likelihood = bcp_step.bayesian_step(dat[:, -1])
      stalktimes.append(dat[numpy.argmax(bcp_likelihood[0]), 0])
  return stalktimes, nostalktimes


if __name__ == '__main__':
  FLAGS = gflags.FLAGS
  gflags.DEFINE_string('infile', '', 'input json')
  gflags.DEFINE_string('outfile', 'out.json', 'output file')
  argv = FLAGS(sys.argv)
  infile = open(FLAGS.infile, 'r')
  component_dat = json.load(infile)
  infile.close()
  (stalktimes, nostalktimes) = stalk_times(component_dat)
  outfile = open(FLAGS.outfile, 'w')
  json.dump({'stalktimes' : stalktimes, 'nostalktimes' : nostalktimes}, outfile)
  outfile.close()
