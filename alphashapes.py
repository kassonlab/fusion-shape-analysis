#!/usr/bin/python
# Code 2017 by Peter Kasson

"""Compute alpha-shapes connected components."""

import gflags
import glob
import json
import numpy
import os
import sys

def compute_alpha(trajname, outdir, preprocess=True, ndx=None,
                  tpr='CLONE0/frame0.tpr', center_res=1040,
                  skip_existing=False,
                  ndxgroup='Tails',
                  gmxdir='/usr/local/gromacs5/bin',
                  connbin='/home/kasson/conncomp/connected_comp'):
  """Compute alpha-shapes components.
  Args:
    trajname:  xtc to process
    outdir:  path to write files
    preprocess:  whether to preprocess traj
    ndx:  input index if preprocessing
    tpr:  input tpr if preprocessing
    center_res:  residue to center on if preprocessing
    skip_existing:  skip existing conncomp output files
    gmxdir:  path to gromacs binaries if preprocessing
    connbin: binary for connected components
  Rets:
  """
  # construct based on gromacs versions
  if os.path.exists('%s/gmx' % gmxdir):
    selbin = '%s/gmx select' % gmxdir
    convbin = '%s/gmx trjconv' % gmxdir
  else:
    selbin = '%s/g_select' % gmxdir
    convbin = '%s/trjconv_axial' % gmxdir

  file_ext = trajname[-4:]
  outname = '%s/%s' % (outdir, trajname.replace(file_ext, '.conncomp.txt'))
  if skip_existing and os.path.exists(outname):
    return numpy.loadtxt(outname, skiprows=1)
  # preprocess data
  if preprocess:
    newtraj = '%s/%s' % (outdir, trajname.replace(file_ext, '_center_tails.xtc'))
    os.system('%s -s %s -n %s -select "group %s; '
              'residue %d and group %s" -on  %s/centertails.ndx'
              % (selbin, tpr, ndx, ndxgroup, center_res, ndxgroup, outdir))
    os.system('echo 1 0 | %s -s %s -pbc mol -xycenter -f %s '
              '-n %s/centertails.ndx -o %s'
              % (convbin, tpr, trajname, outdir, newtraj))
    os.unlink('%s/centertails.ndx' % outdir)
    traj = newtraj
  os.system('%s -xtcfile %s > %s' % (connbin, traj, outname))
  # keep preprocessed traj for now
  return numpy.loadtxt(outname, skiprows=1)


if __name__ == '__main__':
  FLAGS = gflags.FLAGS
  gflags.DEFINE_string('inspec', '',
                       'Pattern matching trajectory files')
  gflags.DEFINE_string('outfile', 'out.txt',
                       'Output filename')
  gflags.DEFINE_string('ndxfile', 'new_tails.ndx',
                       'Indexfile for tails')
  gflags.DEFINE_string('ndxgroup', 'Tails',
                       'Index group to use for tails')
  gflags.DEFINE_string('tprfile', 'CLONE0/frame0.tpr',
                       'TPR file to use')
  gflags.DEFINE_string('connbin', '/home/kasson/conncomp/connected_comp',
                       'Binary for connected components code')
  gflags.DEFINE_boolean('preprocess', True,
                       'Center trajectories')
  gflags.DEFINE_integer('center_res', 1040,
                       'Center around this residue')
  gflags.DEFINE_boolean('skip_existing', False,
                       'Skip trajectories for which output exists')
  argv = FLAGS(sys.argv)
  trajdict = {}
  for trajfile in glob.glob(FLAGS.inspec):
    trajdict[trajfile] = compute_alpha(trajfile,
                         'traj', preprocess=FLAGS.preprocess,
                          ndx=FLAGS.ndxfile, connbin=FLAGS.connbin,
                          tpr=FLAGS.tprfile, skip_existing=FLAGS.skip_existing,
                          ndxgroup=FLAGS.ndxgroup,
                          center_res=FLAGS.center_res).tolist()
  outfile = open(FLAGS.outfile, 'w')
  json.dump(trajdict, outfile)
  outfile.close()
