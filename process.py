#!/usr/bin/python
#
# Code 2019 by Peter Kasson

"""Preprocess trajectories, then run fusion pore analysis."""

import glob
import os
import sys
import gflags

def preprocess(fstring, tprname, outprefix,
               gmxbin='/home/kasson/gromacs/build-2018/bin/gmx'):
  """Preprocess everything matching globstring."""
  os.system('%s select -s %s -on arg_sel -select \'System;'
            ' resid 216 and resname ARG\'' % (gmxbin, tprname))
  for fil in glob.glob(fstring):
    os.system('echo 1 1 0 | %s trjproc -pbc cluster -s %s -n arg_sel.ndx'
              ' -f %s -o temp_proc.xtc -center ' % (gmxbin, tprname, fil))
    os.system('echo 0 | %s trjconv -s %s -f temp_proc.xtc -o %s'
              ' -pbc mol -ur compact' % (gmxbin, tprname, outprefix+fil))
    os.system('rm tmp_proc.xtc')
  os.system('rm arg_sel.ndx')

def make_fpore_ndx(tprfile, atomistic=False,
                   gmxbin='/home/kasson/gromacs/build-2018/bin/gmx'):
  """Make indices for fusion pore analysis."""
  if atomistic:
    selstring = ('"Heads" (resname POPE or resname POPC) and (name N or name '
                 'C12 or name C13 or name C14 or name C15 or name C11 or name P'
                 ' or name O11 or name O12 or name O13 or name O14 or name C1 '
                 'or name C3)')
  else:
    selstring = ('"Heads" (resname POPE or resname POPC) and (name NC3 or '
                 'name NH3 or name PO4 or name GL1 or name GL2)')
  cmdstr = ('%s select -s %s -on %s_heads.ndx -select \'%s\' '
            % (gmxbin, tprfile, tprfile.split('.')[0], selstring))
  os.system(cmdstr)
  cmdstr = ('%s select -s %s -on centerheads.ndx -n %s_heads '
            '-select "group Heads; residue 2711 and group Heads"'
            % (gmxbin, tprfile, tprfile.split('.')[0]))
  os.system(cmdstr)

def fpore(xtcfile, tprfile=None, make_ndx=True, atomistic=False,
          alpha=0.3, manual_ndx=None,
          gmxbin='/home/kasson/gromacs/build-2018/bin/gmx',
          connbin='/home/kasson/analysis-code/conncomp/connected_comp'):
  """Fusion pore analysis."""
  if make_ndx:
    make_fpore_ndx(tprfile, atomistic, gmxbin=gmxbin)
  if not manual_ndx:
    ndxfile = 'centerheads.ndx'
  else:
    ndxfile = manual_ndx
  cmdstr = ('echo 0 | %s trjconv -s %s -pbc mol -f %s -n %s '
            '-o %s_heads.xtc'
            % (gmxbin, tprfile, xtcfile, ndxfile, xtcfile.split('.')[0]))
  os.system(cmdstr)
  cmdstr = ('%s --xtcfile=%s_heads.xtc --alpha=%f > %s.conncomp%.2f_heads.txt'
            % (connbin, xtcfile.split('.')[0], alpha, xtcfile.split('.')[0],
               alpha))
  os.system(cmdstr)

if __name__ == '__main__':
  FLAGS = gflags.FLAGS
  gflags.DEFINE_string('xtcs', '',
                       'Pattern matching trajectory files')
  gflags.DEFINE_string('outprefix', 'traj/',
                       'Prefix output files with this')
  gflags.DEFINE_string('tprfile', '',
                       'TPR filename')
  gflags.DEFINE_string('headndx', '',
                       'Optional ndx file')
  gflags.DEFINE_float('alpha', 0.3,
                      'alpha value in nm')
  gflags.DEFINE_boolean('atomistic', False, 'Atomistic trajectories')
  gflags.DEFINE_boolean('skip_preprocess', False, 'Skip preprocessing')
  argv = FLAGS(sys.argv)
  if not FLAGS.skip_preprocess:
    preprocess(FLAGS.xtcs, FLAGS.tprfile, FLAGS.outprefix)
  make_fpore_ndx(FLAGS.tprfile, FLAGS.atomistic)
  for xtc in glob.glob(FLAGS.xtcs):
    fpore(FLAGS.outprefix+xtc, tprfile=FLAGS.tprfile, make_ndx=False,
          manual_ndx=FLAGS.headndx, alpha=FLAGS.alpha)
