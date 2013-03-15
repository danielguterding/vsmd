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

//lf.cpp

#include "lf.hpp"

void leapfrog_start(GlobalSettings& settings, Results& res){
  
  res.set_velocities(res.get_velocities() - 0.5 * res.get_accelerations() * settings.deltat);
}