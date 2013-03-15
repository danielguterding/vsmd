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

//poscorrect.cpp

#include "poscorrect.hpp"

m3Xf correct_positions(GlobalSettings& settings, m3Xf& positions_in){
  
  int i=0, j=0, factor=0;
  m3Xf positions_out(3,settings.n);
  positions_out = positions_in;
  m33f h_inv;
  h_inv = settings.h.inverse();
  v3f solution;
  solution << 0,0,0;
 
  for(i=0;i<positions_in.outerSize();i++){
    solution = h_inv * positions_in.col(i);
    for(j=0;j<3;j++){
      if(solution(j)<0){
	factor = abs((int)floor(solution(j)));
	positions_out.col(i).noalias() += factor * settings.h.col(j); 
      }
      if(solution(j)>1){
	factor = (int)floor(solution(j));
	positions_out.col(i).noalias() -= factor * settings.h.col(j);
      }
    }
  }
 return positions_out;
}
