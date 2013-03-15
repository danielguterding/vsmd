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

//init.cpp
#include "init.hpp"

m3Xf initialize_positions_fcc(GlobalSettings& settings){
  
  int i=0, j=0, k=0, l=0;
  fptype d = settings.h(0,0)/sqrt(2.0)/settings.gridsize; //cubic initial h is assumed 
  
  m3Xf positions(3,settings.n);
  positions = positions.Zero(positions.innerSize(),positions.outerSize());
  v3f temp, bias1, bias2, bias3, bias4;
  temp << 0, 0, 0;
  bias1 << 0.5, 0.5, 0.5; //base vectors of the fcc lattice
  bias2 << 0.5, 1.0, 1.0;
  bias3 << 1.0, 0.5, 1.0;
  bias4 << 1.0, 1.0, 0.5;
  
  for(i=0;i<(int)settings.gridsize;i++){
    for(j=0;j<(int)settings.gridsize;j++){
      for(k=0;k<(int)settings.gridsize;k++){
	temp << i, j, k;
	positions.col(l++) = ((bias1 + temp)*d);
	positions.col(l++) = ((bias2 + temp)*d);
	positions.col(l++) = ((bias3 + temp)*d);
	positions.col(l++) = ((bias4 + temp)*d);
      }
    }
  }
  
  return positions;
}

m3Xf maxwell_boltzmann_velocities(GlobalSettings& settings){
  
  int i=0, j=0;
  default_random_engine generator((uint)time(0) + settings.processid*1234); //normal distributed random number generator
  normal_distribution<fptype> normal(0,1); //normal distribution with zero mean and unity standard deviation
  
  m3Xf velocities(3,settings.n);
  velocities = velocities.Zero(velocities.innerSize(),velocities.outerSize());
  
  for(i=0;i<(int)settings.n;i++){
    for(j=0;j<3;j++){
      velocities(j,i) = normal(generator);
    }
  }
  
  velocities /= sqrt(settings.T_des);
  
  return velocities;
}

m3Xf correct_momentum_to_zero(GlobalSettings& settings, m3Xf velocities){
  
  int i=0;
  v3f momentum(3);
  momentum << 0, 0, 0;
  
  for(i=0;i<velocities.innerSize();i++){
    momentum += velocities.col(i);
  }
  
  momentum /= settings.n;
  
  for(i=0;i<velocities.innerSize();i++){
    velocities.col(i) -= momentum;
  }
  
  return velocities;
}

m3Xf initialize_velocities_maxwell_boltzmann(GlobalSettings& settings){
  
  m3Xf velocities(3,settings.n);
  
  velocities = maxwell_boltzmann_velocities(settings);
  velocities = correct_momentum_to_zero(settings, velocities);
  
  return velocities;
}

void initialize(GlobalSettings& settings, Results& res){ //initialize positions, velocities and accelerations with total momentum constrained to zero
  
  m3Xf positions(3,settings.n), velocities(3,settings.n);
  positions = initialize_positions_fcc(settings);
  res.set_positions(positions);
  velocities = initialize_velocities_maxwell_boltzmann(settings);
  velocities = correct_momentum_to_zero(settings, velocities);
  res.set_velocities(velocities);
  forces(settings, res);
  res.set_accelerations(res.get_accelerations_future());
  res.set_accelerations_future(velocities.Zero(velocities.innerSize(), velocities.outerSize()));
  res.update_kinetic_tensor();
  res.update_pressure_tensor(settings);
}
 
