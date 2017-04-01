//alel.h
//
// Copyright 2015 Jessen T. Havill, Chinmoy Bhatiya, Steven M. Johnson, 
// Mitchell J. Keller, J. D. Sheets, and Jeffrey S. Thompson.
// 
// This file is part of DRD.
// 
// DRD is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// DRD is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with DRD.  If not, see <http://www.gnu.org/licenses/>.


#include <string>
using namespace std;

#ifndef ALEL_H
#define ALEL_H

// A structure representing a motif alignment

class AlEl
{
  public:
  
	int score;                // alignment score
	unsigned int position;    // index in the first string where alignment begins
	string str;               // string representing the alignment
	int index;                // index of the match in adjacency matrix of multipartite graph
};

#endif
