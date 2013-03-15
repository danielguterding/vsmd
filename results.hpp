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

//results.hpp
#ifndef __RESULTS_H_INCLUDED__ 
#define __RESULTS_H_INCLUDED__ 

#include <eigen3/Eigen/Dense>

#include "typedefs.hpp"
#include "settings.hpp"
using namespace Eigen;
 
class Results{
  
    m3Xf positions;
    m3Xf positions_past;
    m3Xf positions_future;
    m3Xf velocities;
    m3Xf velocities_past;
    m3Xf velocities_future;
    m3Xf accelerations;
    m3Xf accelerations_past;
    m3Xf accelerations_future;
  
    fptype temperature;
    fptype pressure; //hydrostatic pressure
    m33f pressure_tensor; //tensorial pressure
    fptype virial; //scalar virial
    m33f virial_tensor; //tensorial virial
    fptype kinetic_energy; //scalar kinetic energy
    m33f kinetic_tensor; //tensorial kinetic energy
    fptype potential_energy;
    fptype total_energy;
  public:
    Results(GlobalSettings& settings);
    void set_positions(const m3Xf&);
    void set_positions_past(const m3Xf&);
    void set_positions_future(const m3Xf&);
    void set_velocities(const m3Xf&);
    void set_velocities_past(const m3Xf&);
    void set_velocities_future(const m3Xf&);
    void set_accelerations(const m3Xf&);
    void set_accelerations_past(const m3Xf&);
    void set_accelerations_future(const m3Xf&);
    
    void set_virial_tensor(const m33f&);
    void set_potential_energy(const fptype&);
    void update_kinetic_tensor();
    void update_pressure_tensor(const GlobalSettings&);
    
    m3Xf get_positions();
    m3Xf get_positions_past();
    m3Xf get_positions_future();
    m3Xf get_velocities();
    m3Xf get_velocities_past();
    m3Xf get_velocities_future();
    m3Xf get_accelerations();
    m3Xf get_accelerations_past();
    m3Xf get_accelerations_future();
    
    m33f get_pressure_tensor();
    m33f get_virial_tensor();
    m33f get_kinetic_tensor();
    fptype get_temperature();
    fptype get_pressure();
    fptype get_virial();
    fptype get_kinetic_energy();
    fptype get_potential_energy();
    fptype get_total_energy();
  private:
    void update_total_energy();
    void update_temperature();
    
}; 

#endif
