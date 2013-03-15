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

//lfNVTberendsen.cpp

#include "lfNVTberendsen.hpp"

void leapfrogNVTberendsen(GlobalSettings& settings, Results& res){
  
  fptype lambda;
  m3Xf positions_future(3,settings.n), velocities_future(3,settings.n), accelerations_future(3,settings.n);
  positions_future = positions_future.Zero(positions_future.innerSize(), positions_future.outerSize());
  velocities_future = velocities_future.Zero(velocities_future.innerSize(), velocities_future.outerSize());
  accelerations_future = accelerations_future.Zero(accelerations_future.innerSize(), accelerations_future.outerSize());
  
  lambda = sqrt(1 + settings.deltat/settings.t_thermo*(settings.T_des/res.get_temperature() - 1));
  if(lambda < 0.9){
    lambda = 0.9;
  }
  if(lambda > 1.1){
    lambda = 1.1;
  }
  
  positions_future = res.get_positions() + res.get_velocities() * settings.deltat + 0.5 * res.get_accelerations() * pow(settings.deltat,2);
  velocities_future = res.get_velocities() * lambda + 0.5 * res.get_accelerations() * settings.deltat;
  
  positions_future = correct_positions(settings, positions_future);
  
  forces(settings, res);
  
  velocities_future.noalias() += 0.5 * res.get_accelerations_future() * settings.deltat;
  
  res.set_positions_past(res.get_positions());
  res.set_velocities_past(res.get_velocities());
  res.set_accelerations_past(res.get_accelerations());
  
  res.set_positions(positions_future);
  res.set_velocities(velocities_future);
  res.set_accelerations(res.get_accelerations_future());
  
  res.set_positions_future(positions_future.Zero(positions_future.innerSize(), positions_future.outerSize()));
  res.set_velocities_future(velocities_future.Zero(velocities_future.innerSize(), velocities_future.outerSize()));
  res.set_accelerations_future(accelerations_future.Zero(accelerations_future.innerSize(), accelerations_future.outerSize()));
  
  res.update_kinetic_tensor();
  res.update_pressure_tensor(settings);
}
