/*
* Copyright (c) 2013, Daniel Guterding <guterding@itp.uni-frankfurt.de>
*
* This file is part of vsmd.
*
* vsmd is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* vsmd is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with vsmd. If not, see <http://www.gnu.org/licenses/>.
*/

//settings.hpp
#ifndef __SETTINGS_H_INCLUDED__ 
#define __SETTINGS_H_INCLUDED__ 

#include <eigen3/Eigen/Dense>

#include "typedefs.hpp"

struct GlobalSettings {
  uint processid; //process number for bias in random seed
  uint mode; // "1" for verletNVE, "2" for verletNVTscaling, "3" for leapfrogNVTberendsen, "4" for leapfrogNPTberendsen 
  uint potential; //choose between lennard jones ("1") and buckingham potential ("2")
  uint gridsize; //number of conventional fcc cells in each spatial direction, i.e. N=4*gridsize**3
  uint n; //number of particles
  uint nsteps;
  m33f h; //simulation box h = (a,b,c)
  fptype rcutoff; //cutoff radius in units of lennard jones sigma
  fptype deltat; //time step width in lennard jones reduced units
  fptype t_thermo; //thermostat rise time in lennard jones reduced units
  fptype t_baro; //barostat rise time in lennard jones reduced units
  fptype T_des; //desired temperature in lennard jones reduced units
  m33f P_des; //desired pressure tensor in lennard jones reduced units
};

#endif