//ribodef.cpp
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


#include<string>
#include<iostream>
#include<fstream>
#include "ribodef.h"

using namespace std;

RiboDef::RiboDef()
{
	snippet_size = 0;
	overlap_size = 0;
	motifs = NULL;
	minscores = NULL;
	nummotifs = 0;
	distances = NULL;
	mindistances = NULL;
	lowestTotal = 0;
	maxLength = 0;
	vienna = "";
	minViennaScore = 0;
	motifs_vienna = "";
}

RiboDef::RiboDef(ifstream& infile)
{
	infile >> snippet_size;
	if (snippet_size <= 0)
		snippet_size = 100;
		
	infile >> overlap_size;
	if (overlap_size < 0)
		overlap_size = 0;
	if (overlap_size >= snippet_size)
	{
		if (snippet_size > 10)
			overlap_size = snippet_size - 10;
		else
			overlap_size = 0;
	}
	infile >> nummotifs;

	motifs = new string[nummotifs];
	minscores = new unsigned int[nummotifs];
	distances = new unsigned int[nummotifs + 1];
	mindistances = new unsigned int[nummotifs - 1];

	for(unsigned int i = 0; i < nummotifs ; i++)
	{
		infile >> motifs[i];
		infile >> minscores[i];
	}
	
	for(unsigned int i = 0; i < nummotifs + 1 ; i++)
		infile >> distances[i];
		
	for(unsigned int i = 0; i < nummotifs - 1 ; i++)
		infile >> mindistances[i];

	infile >> lowestTotal;       // min motif identities
	infile >> minViennaScore;    // min Vienna identities
	infile >> maxLength;         // maximum length of found riboswitch
	infile >> vienna;            // consensus Vienna string
	infile >> motifs_vienna;     // consensus Vienna string with motifs inserted
	infile.close();
}

void RiboDef::copy(const RiboDef &r)
{
	snippet_size = r.snippet_size;
	overlap_size = r.overlap_size;
	nummotifs = r.nummotifs;

	motifs = new string[nummotifs];
	minscores = new unsigned int[nummotifs];
	distances = new unsigned int[nummotifs + 1];
	mindistances = new unsigned int[nummotifs - 1];

	for(unsigned int i = 0; i < nummotifs ; i++)
	{
		motifs[i] = r.motifs[i];
		minscores[i] = r.minscores[i];
	}
	
	for(unsigned int i = 0; i < nummotifs + 1 ; i++)
		distances[i] = r.distances[i];
		
	for(unsigned int i = 0; i < nummotifs - 1 ; i++)
		mindistances[i] = r.mindistances[i];

	lowestTotal= r.lowestTotal;
	vienna = r.vienna;
	minViennaScore = r.minViennaScore;
	maxLength = r.maxLength;
	motifs_vienna = r.motifs_vienna;
}

RiboDef::RiboDef(const RiboDef &r)
{
	copy(r);
}

RiboDef& RiboDef::operator=(const RiboDef &r)
{
        if (&r != this)
        {
                delete this;
                copy(r);
        }
        return *this;
}

RiboDef::~RiboDef()
{
	if (motifs != NULL)
		delete [] motifs;
	if (minscores != NULL)
		delete [] minscores;
	if (distances != NULL)
		delete [] distances;
	if (mindistances != NULL)
		delete [] mindistances;
}

ostream& operator<<(ostream& os, const RiboDef& r)
{
//	os << "Snippet/overlap size: " << r.snippet_size << "/" << r.overlap_size << endl;
	
	os << "5'- &lt;&lt; " << r.distances[0] << " &gt;&gt; ";
	for(unsigned int i = 0; i < r.nummotifs; i++)
	{	
		if (i > 0)
			os << " &lt;&lt; [" << r.mindistances[i-1] << "," << r.distances[i] << "] &gt;&gt; ";
		os << r.motifs[i] << "(&ge;" << r.minscores[i] << ")";
	}
	
	os << " &lt;&lt; " << r.distances[r.nummotifs] << " &gt;&gt; -3'" << endl;
	
	return os;
}

