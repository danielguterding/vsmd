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

//lfNPTberendsen.cpp

#include "lfNPTberendsen.hpp"

void leapfrogNPTberendsen(GlobalSettings& settings, Results& res){
  
  fptype lambda=0;
  m33f mu;
  mu = mu.Identity();
  m3Xf positions_future(3,settings.n), velocities_future(3,settings.n), accelerations_future(3,settings.n);
  positions_future = positions_future.Zero(positions_future.innerSize(), positions_future.outerSize());
  velocities_future = velocities_future.Zero(velocities_future.innerSize(), velocities_future.outerSize());
  accelerations_future = accelerations_future.Zero(accelerations_future.innerSize(), accelerations_future.outerSize());
  
  forces(settings, res);
  res.update_kinetic_tensor();
  res.update_pressure_tensor(settings);
  
  mu << res.get_pressure_tensor()(0,0), 0, 0,
        0, res.get_pressure_tensor()(1,1), 0, 
        0, 0, res.get_pressure_tensor()(2,2);
  mu = mu.Identity() + settings.deltat/settings.t_baro/3.0 * (mu - settings.P_des);
  
  lambda = sqrt(1 + settings.deltat/settings.t_thermo*(settings.T_des/res.get_temperature() - 1));
  if(lambda < 0.9){
    lambda = 0.9;
  }
  if(lambda > 1.1){
    lambda = 1.1;
  }
  
  velocities_future = res.get_velocities() + res.get_accelerations_future() * settings.deltat;
  velocities_future *= lambda;
  
  positions_future = res.get_positions() + velocities_future * settings.deltat;
  positions_future = correct_positions(settings, positions_future);
  
  positions_future = mu * positions_future;
  settings.h = mu * settings.h;
  
  
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
