//ribodef.h
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


#ifndef RIBODEF_H
#define RIBODEF_H

#include <string>
#include <iostream>
using namespace std;

class RiboDef
{
  public:

	unsigned int snippet_size;
	unsigned int overlap_size;
	string * motifs;
	unsigned int * minscores;
	unsigned int nummotifs;
	unsigned int * distances;
	unsigned int * mindistances;
	unsigned int lowestTotal;
	unsigned int maxLength;
	string vienna;
	int minViennaScore;
	string motifs_vienna;
	
	RiboDef();
	RiboDef(ifstream& infile);
	RiboDef(const RiboDef &r);
	~RiboDef();
	RiboDef& operator=(const RiboDef &r);

  private:
	void copy(const RiboDef &r);
};

ostream& operator<<(ostream& os, const RiboDef& r);

#endif
