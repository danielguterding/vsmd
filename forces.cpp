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

//forces.cpp
#include "forces.hpp"

//-----Functions for ClosestMirrorParticle Class-----

ClosestMirrorParticle::ClosestMirrorParticle(const GlobalSettings& settings){

  h = settings.h;
  h_inv = h.inverse();
}

void ClosestMirrorParticle::set_posi(const v3f& vec){
  
  posi = vec;
}

void ClosestMirrorParticle::set_posj(const v3f& vec){
  
  posj = vec;
}

void ClosestMirrorParticle::update_distance(){
  
  distvec = posj - posi;
  dist2 = distvec.array().pow(2);
  dist = sqrt(dist2.sum());
}

v3f ClosestMirrorParticle::get_distance_vector(){
  
  return distvec;
}

fptype ClosestMirrorParticle::get_distance_scalar(){
  
  return dist;
}

void ClosestMirrorParticle::update_closest_mirror(){
  
  int i=0;
  
  solution = h_inv * distvec;
  
  for(i=0;i<3;i++){
    if(fabs(solution(i)) > 0.5){
      posj.noalias() -= h.col(i) * copysignf(1.0,solution(i));
    }
  }
  
  update_distance();
}

//-----Main Functions of Force Module-----

fptype scalar_force(const GlobalSettings& settings, const fptype& dist){
  
  fptype force=0;
  
  switch(settings.potential){
    case 1: 
      force = 24*(2*pow(dist, -13) - pow(dist, -7));
      break;
    case 2: 
      force = 10.692*exp(12*(1-0.891*dist)) - 24*pow(dist, -7);
      break;
  }
  
  return force;
}


v3f vector_force(const GlobalSettings& settings, ClosestMirrorParticle& mirror){
  
  fptype quot=0;
  quot = scalar_force(settings, mirror.get_distance_scalar())/mirror.get_distance_scalar();
  
  v3f forcevec;
  forcevec = -quot*mirror.get_distance_vector();
  
  return forcevec;
}

fptype cutoff_energy(const GlobalSettings& settings){
  
  fptype energy=0;
  
  switch(settings.potential){
    case 1: 
      energy = 4*(pow(settings.rcutoff, -12) - pow(settings.rcutoff, -6));
      break;
    case 2: 
      energy = exp(12*(1-0.891*settings.rcutoff)) - 4*pow(settings.rcutoff, -6);  
      break;
  }
  
  return energy;
}

fptype potential_energy_pairwise(const GlobalSettings& settings, const fptype dist, const fptype cutoff_energy){
  
  fptype energy=0;
  
  switch(settings.potential){
    case 1: 
      energy = 4*(pow(dist, -12) - pow(dist, -6)) - cutoff_energy;
      break;
    case 2: 
      energy = exp(12*(1-0.891*dist)) - 4*pow(dist, -6) - cutoff_energy;
      break;
  }
  
  return energy;
}

void forces(GlobalSettings& settings, Results& res){
  
  v3f force(3);
  force = force.Zero();
  m33f virial_tensor;
  m3Xf accelerations(3,settings.n), positions(3,settings.n);
  virial_tensor = virial_tensor.Zero();
  positions = res.get_positions();
  accelerations = accelerations.Zero(accelerations.innerSize(), accelerations.outerSize());
  fptype cenerg=0, penerg=0;
  int i=0,j=0;
  ClosestMirrorParticle mirror(settings);
  
  
  cenerg = cutoff_energy(settings);
  
  for(i=0;i<accelerations.outerSize()-1;i++){
    mirror.set_posi(positions.col(i));
    for(j=i+1;j<accelerations.outerSize();j++){
      mirror.set_posj(positions.col(j));
      mirror.update_distance();
      
      if(mirror.get_distance_scalar() > settings.rcutoff){
	mirror.update_closest_mirror();
      }
      
      if(mirror.get_distance_scalar() < settings.rcutoff){
	force = vector_force(settings, mirror);
	accelerations.col(i).noalias() += force;
	accelerations.col(j).noalias() -= force;
	
	virial_tensor.noalias() -= force * mirror.get_distance_vector().transpose();
	
	penerg += potential_energy_pairwise(settings, mirror.get_distance_scalar(), cenerg);
      }
    }
  }
  res.set_accelerations_future(accelerations);
  res.set_virial_tensor(virial_tensor);
  res.set_potential_energy(penerg);
}