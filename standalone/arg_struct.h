// arg_struct.h
// struct for passing args to finder() thread
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


#ifndef ARG_STRUCT_H
#define ARG_STRUCT_H

#include <string>
#include <iostream>
#include "ribodef.h"

using namespace std;

struct arg_struct 
{
	string forwSeq;         // plus sequence to search
	string backSeq;         // minus sequence to search
	int seed;               // global index of first nt in forwSeq
	RiboDef *riboswitch;    // riboswitch to search for
	string foldername;      // directory in user directory to place results
	int* orfS;              // array of start indices of ORFs
	int* orfF;              // array of finish indices of ORFs
	int numOrfs;            // number of forward ORFs represented in above arrays
	int numROrfs;           // number of reverse ORFs represented in above arrays
	string fullSeq;         // the original whole sequence
	ofstream* testfile;      // file to which to write test output (if opened)
};

#endif
