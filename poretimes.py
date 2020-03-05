#!/usr/bin/python
# Code to calculate pore-formation times
# Peter Kasson, 2017

import commands
import gflags
import glob
import numpy
import os
import sys

def make_conncomp(trajname, overwrite=False, gmxbin='/usr/local/gromacs5/bin',
                  connbin='/home/kasson/analysis-code/conncomp'):
  """Calculate connected components.
  """
  trjconv_bin = '%s/gmx trjconv' % gmxbin
  sel_bin = '%s/gmx select' % gmxbin

  sel_cmd = ('%s -s %s -on %s_head2.ndx -select \'"Heads" (resname POPE or '
             'resname POPC) and (name N or name C12 or name C13 or name C14 '
             'or name C15 or name C11 or name P or name O11 or name O12 or '
             'name O13 or name O14 or name C1 or name C2)\''
             % (sel_bin, trajname[:-4], trajname[:-4]))
  sel2_cmd = ('%s -s %s -on tmp_centertail2.ndx -n %s_head2.ndx -select '
              '"group Heads; residue 2711 and group Heads" '
              % (sel_bin, trajname[:-4], trajname[:-4]))
  trjconv_cmd = ('echo 0 | %s -s %s -pbc mol -f %s -n tmp_centertail2.ndx -o '
                 'traj/%s_center3_tails.xtc -trans -5 5 0'
                 % (trjconv_bin, trajname[:-4], trajname[:-4], trajname[:-4]))
  conn_cmd = ('%s/connected_comp --xtcfile=traj/%s_center3_tails.xtc --alpha=0.3'
              % (connbin, trajname[:-4]))

  if not overwrite and not os.path.exists('traj/%s_center3_tails.xtc'
                                          % trajname[:-4]):
    os.system(sel_cmd)
    os.system(sel2_cmd)
    os.system(trjconv_cmd)
  conn_data = commands.getoutput(conn_cmd).splitlines()
  return numpy.genfromtxt(conn_data[1:])

def find_poretime(conn_data, debug=False, porestr='22222'):
  """Find approximate pore formation time.
  Args:  conn_data:  numpy array with time and number of components.
  Rets:  first time 2 components stably formed; -1 * time simulated if none.
  """
  firstidx = ''.join(str(int(x)) for x in conn_data[:, 1]).find(porestr)
  if debug:
    print conn_data
  return conn_data[firstidx, 0] if firstidx >=0 else -1*conn_data[firstidx, 0]

if __name__ == '__main__':
  FLAGS = gflags.FLAGS
  gflags.DEFINE_string('trajfiles', '', 'glob string for trajectories')
  gflags.DEFINE_string('outfile', 'out.txt', 'output file')
  gflags.DEFINE_string('porestr', '22222', 'Connected components string for pore')
  gflags.DEFINE_boolean('read_conncomp', False, 'read conncomp instead of traj')
  argv = FLAGS(sys.argv)
  outfile = open(FLAGS.outfile, 'w')
  for trajfile in glob.glob(FLAGS.trajfiles):
    if not FLAGS.read_conncomp:
      compdat = make_conncomp(trajfile)
    else:
      compdat = numpy.loadtxt(trajfile, skiprows=1)
    poretime = find_poretime(compdat, porestr=FLAGS.porestr)
    outfile.write('%s %d\n' % (trajfile, poretime))
  outfile.close()
