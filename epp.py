#
# Copyright (c) 2013, Daniel Guterding <guterding@itp.uni-frankfurt.de>
#
# This file is part of vsmd.
#
# vsmd is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# vsmd is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with vsmd. If not, see <http://www.gnu.org/licenses/>.
#

import os
import shutil
import subprocess
import numpy as np

os.system("taskset -p 0xff %d" % os.getpid()) #importing numpy sets cpu affinity of ALL SUBPROCESSES to one specific core

def call_vsmd(processid, T_des, P_des):
  mode = 5
  potential = 1
  gridsize = 6
  nsteps = 50000
  boxlength = 20
  rcutoff = 2.1
  deltat = 3E-3
  t_thermo = 10
  t_baro = 150
  
  command = './vsmd %i %i %i %i %i %f %f %f %f %f %f %f' % (processid, mode, potential, gridsize, nsteps, boxlength, rcutoff, deltat,
                                                           t_thermo, t_baro, T_des, P_des)
  p = subprocess.Popen(command.split())       
  return p
  
def call_procs(nruns, T_des, P_des):
  processes = []
  for i in range(nruns):
    p = call_vsmd(i, T_des, P_des)
    processes.append(p)
  for p in processes:
    p.wait()
  
def get_results(data_path, nruns, P):
  results = []
  for i in range(nruns):
    input_path = data_path + "P=%i.%i.dat" % (P, i)
    input_file = open(input_path, 'r')
    line = input_file.readline()
    input_file.close()
    results.append(str(line).split(';'))
    
  return results
    
def get_average(results):
  res_col = 2
  average = 0
  for res in results:
    average += float(res[res_col])
  average /= float(len(results))
  
  return average
  
def plot_results(Pressure, data_path, datfile_path):
  pyxfile_path = 'data/plot.pyx'
  pyxfile = open(pyxfile_path, 'w')
  
  script = """ set multiplot\n
               set nodisplay\n
               \n
               set title "Energy per Particle in Lennard-Jones Potential"\n
               set xlabel "T [MD]"\n
               set ylabel "E [MD]"\n
               set key outside\n
               set terminal pdf\n
               \n
               set output "%sP=%i.pdf"\n
               plot "%s" using 1:3 title '$P=%i$'\n
               set display ; refresh  """ % (data_path, Pressure, datfile_path, Pressure)

  pyxfile.write(script)
  pyxfile.close()
  
  command = 'pyxplot %s' % pyxfile_path
  p = subprocess.Popen(command.split())
  p.wait()
  
def epp():
  nruns = 5 #number of runs per p-T-point
  Pressure = 5.0
  Temperature_list = list(np.arange(1.055, 1.20, 0.005))
  
  data_path = 'data/';
  logfile_path = 'data/log.txt'
  datfile_path = 'data/plot.dat'
  
  if (not os.path.exists(data_path)):
    os.mkdir(data_path)
    
  for path in [logfile_path, datfile_path]:
    if (os.path.isfile(path)):
      shutil.copyfile(path, path + '.bak')
  
  outfile = open(datfile_path, 'w')
  outfile.write('# Temperature, Pressure, Energy per Particle\n')
  outfile.close()
  
  outfile = open(logfile_path, 'a')
  outfile.write('# Temperature, Pressure, Energy per Particle\n')
  outfile.close()
  
  for Temperature in Temperature_list:
    T_des = Temperature
    P_des = Pressure
    
    call_procs(nruns, T_des, P_des)
      
    results = get_results(data_path, nruns, Pressure)
    average = get_average(results)
    
    logfile = open(logfile_path, 'a')
    datfile = open(datfile_path, 'a')
    
    for res in results:
      logfile.write('%f %f %f\n' % (float(res[0]), float(res[1]), float(res[2])))
    logfile.write('Average: %f\n\n' % average)
    logfile.close()
    
    datfile.write('%f %f %f\n' % (T_des, P_des, average))
    datfile.close()
    
  plot_results(P_des, data_path, datfile_path)
    
epp()
