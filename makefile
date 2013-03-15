#
# Copyright (c) 2013, Daniel Guterding <guterding@itp.uni-frankfurt.de>
#
# This file is part of vsmd.
#
# vsmd is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# vsmd is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with vsmd. If not, see <http://www.gnu.org/licenses/>.
#

CXX      = g++
CXXFLAGS = -Wall -march=native -O3 -flto -fuse-linker-plugin -std=c++0x
LDFLAGS  = -lm -lboost_system -lboost_filesystem -lgomp

OBJECTS = main.o init.o forces.o results.o poscorrect.o verletNVE.o verletNVTscaling.o lf.o lfNVTberendsen.o lfNPTberendsen.o 
DEFINES =

vsmd : $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(DEFINES) $(OBJECTS) $(LDFLAGS) -o vsmd

main.o : main.cpp typedefs.hpp results.hpp settings.hpp init.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c main.cpp -o main.o

init.o : init.hpp init.cpp typedefs.hpp results.hpp settings.hpp forces.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c init.cpp -o init.o
	
forces.o : forces.hpp forces.cpp typedefs.hpp results.hpp settings.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c forces.cpp -o forces.o
	
results.o : results.hpp results.cpp settings.hpp typedefs.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c results.cpp -o results.o
	
poscorrect.o : poscorrect.hpp poscorrect.cpp settings.hpp typedefs.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c poscorrect.cpp -o poscorrect.o
	
verletNVE.o : verletNVE.hpp verletNVE.cpp results.hpp settings.hpp typedefs.hpp forces.hpp poscorrect.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c verletNVE.cpp -o verletNVE.o
	
lf.o : lf.hpp lf.cpp settings.hpp results.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c lf.cpp -o lf.o
	
verletNVTscaling.o : verletNVTscaling.hpp verletNVTscaling.cpp results.hpp settings.hpp typedefs.hpp forces.hpp poscorrect.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c verletNVTscaling.cpp -o verletNVTscaling.o
	
lfNVTberendsen.o : lfNVTberendsen.hpp lfNVTberendsen.cpp results.hpp settings.hpp typedefs.hpp forces.hpp poscorrect.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c lfNVTberendsen.cpp -o lfNVTberendsen.o
	
lfNPTberendsen.o : lfNPTberendsen.hpp lfNPTberendsen.cpp results.hpp settings.hpp typedefs.hpp forces.hpp poscorrect.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c lfNPTberendsen.cpp -o lfNPTberendsen.o

clean:
	rm vsmd $(OBJECTS)
	rm -R data
