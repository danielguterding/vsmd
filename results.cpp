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

//results.cpp
#include "results.hpp"

Results::Results(GlobalSettings& settings){
  
  positions.resize(3,settings.n);
  positions_past.resize(3,settings.n);
  positions_future.resize(3,settings.n);
  velocities.resize(3,settings.n);
  velocities_past.resize(3,settings.n);
  velocities_future.resize(3,settings.n);
  accelerations.resize(3,settings.n);
  accelerations_past.resize(3,settings.n);
  accelerations_future.resize(3,settings.n);
}

void Results::set_positions(const m3Xf& input){
  
  positions = input;
}

void Results::set_positions_past(const m3Xf& input){
  
  positions_past = input;
}
void Results::set_positions_future(const m3Xf& input){
  
  positions_future = input;
}
void Results::set_velocities(const m3Xf& input){
  
  velocities = input;
}
void Results::set_velocities_past(const m3Xf& input){
  
  velocities_past = input;
}
void Results::set_velocities_future(const m3Xf& input){
  
  velocities_future = input;
}
void Results::set_accelerations(const m3Xf& input){
  
  accelerations = input;
}
void Results::set_accelerations_past(const m3Xf& input){
  
  accelerations_past = input;
}
void Results::set_accelerations_future(const m3Xf& input){
  
  accelerations_future = input;
}

void Results::set_virial_tensor(const m33f& input){
  
  virial_tensor = input;
  virial = virial_tensor.trace();
}

void Results::set_potential_energy(const fptype& input){
  
  potential_energy = input;
  Results::update_total_energy();
}

void Results::update_kinetic_tensor(){
  
  int i=0;
  
  kinetic_tensor.noalias() = kinetic_tensor.Zero();
  
  for(i=0;i<velocities.outerSize();i++){
    kinetic_tensor.noalias() += 0.5 * velocities.col(i) * velocities.col(i).transpose();
  }
  
  kinetic_energy = kinetic_tensor.trace();
  Results::update_total_energy();
  Results::update_temperature();
}

void Results::update_pressure_tensor(const GlobalSettings& settings){
  
  fptype volume=0;
  pressure_tensor = pressure_tensor.Zero();
  
  volume = settings.h.determinant();
  pressure_tensor = 1.0/volume * (2 * kinetic_tensor + virial_tensor);
  pressure = pressure_tensor.trace();
}

void Results::update_total_energy(){
  
  total_energy = kinetic_energy + potential_energy;
}

void Results::update_temperature(){
  
  temperature = 2.0*kinetic_energy/3.0/(fptype)velocities.outerSize();
}

m3Xf Results::get_positions(){
  
  return positions;
}

m3Xf Results::get_positions_past(){
  
  return positions_past;
}

m3Xf Results::get_positions_future(){
  
  return positions_future;
}

m3Xf Results::get_velocities(){
  
  return velocities;
}

m3Xf Results::get_velocities_past(){
  
  return velocities_past;
}

m3Xf Results::get_velocities_future(){
  
  return velocities_future;
}

m3Xf Results::get_accelerations(){
  
  return accelerations;
}

m3Xf Results::get_accelerations_past(){
  
  return accelerations_past;
}

m3Xf Results::get_accelerations_future(){
  
  return accelerations_future;
}
    
m33f Results::get_pressure_tensor(){
  
  return pressure_tensor;
}

m33f Results::get_virial_tensor(){
  
  return virial_tensor;
}

m33f Results::get_kinetic_tensor(){
  
  return kinetic_tensor;
}

fptype Results::get_temperature(){
  
  return temperature;
}

fptype Results::get_pressure(){
  
  return pressure;
}

fptype Results::get_virial(){
  
  return virial;
}

fptype Results::get_kinetic_energy(){
  
  return kinetic_energy;
}

fptype Results::get_potential_energy(){
  
  return potential_energy;
}

fptype Results::get_total_energy(){
  
  return total_energy;
} 