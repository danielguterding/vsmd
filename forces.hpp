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

//forces.hpp

#include <math.h>

#include "results.hpp"
using namespace Eigen;

//-----Class Definitions-----

#ifndef __MIRROR_H_INCLUDED__ 
#define __MIRROR_H_INCLUDED__ 

class ClosestMirrorParticle{
    v3f posi, posj, distvec, dist2, solution;
    m33f h, h_inv;
    fptype dist;
  public:
    ClosestMirrorParticle(const GlobalSettings& set);
    void set_posi(const v3f& vec);
    void set_posj(const v3f& vec);
    void update_distance();
    v3f get_distance_vector();
    fptype get_distance_scalar();
    void update_closest_mirror();
};
#endif

//-----Function Declarations-----

void forces(GlobalSettings& settings, Results& results);


