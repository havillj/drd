// util.cpp
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


#include "drd.h"

using namespace std;

// Scramble a string by performing num passes of random swaps of each char.
string scrambleString(string str, int num)
{
     int n = str.length();
     srand(time(NULL));
     for (int i = 0; i < num; i++)
		 for(int j = n; j > 0; j--)
		      swap(str[rand() % n], str[j - 1]);
     return str;
}

// Uniformly randomly choose a character from nts.
char choose(const string& nts)
{
	return nts[rand() % nts.length()];
}

// Replace ambiguous bases in dna with ACGT.
string noAmbiguous(const string& dna)
{
	string newDNA = dna;
	for (unsigned int i = 0; i < dna.length(); i++)
	{	
		switch (dna[i])
		{
			case 'R': newDNA[i] = choose(string("AG")); break;
			case 'Y': newDNA[i] = choose(string("CT")); break;
			case 'M': newDNA[i] = choose(string("AC")); break;
			case 'K': newDNA[i] = choose(string("GT")); break;
			case 'S': newDNA[i] = choose(string("CG")); break;
			case 'W': newDNA[i] = choose(string("AT")); break;
			case 'B': newDNA[i] = choose(string("CGT")); break;
			case 'D': newDNA[i] = choose(string("AGT")); break;
			case 'H': newDNA[i] = choose(string("ACT")); break;
			case 'V': newDNA[i] = choose(string("ACG")); break;
			case 'N': newDNA[i] = choose(string("ACGT"));
		}
	}
	return newDNA;
}

// Return reverse of string.
string stringrev(string str)
{
	int n = str.length();
	for(int i = n; i >= (n / 2) + 1; i--)
		swap(str[n - i], str[i - 1]);
	return str;
}

//Return standardized version of dna (No U's, all caps).
string standardize(const string& dna)
{
	unsigned int n = dna.length();
	char a[n + 1];
	unsigned int j = 0;
	for (unsigned int i = 0; i < n; i++)
	{
		if (isalpha(dna[i]))
		{
			a[j] = toupper(dna[i]);
			if (a[j] == 'U')
				a[j] = 'T';
			j++;
		}
	}
	a[j] = 0;
	return string(a);
}

// Return complement of dna.
string complement(string dna)
{
	for(unsigned int i = 0; i < dna.length(); i++)
	{
		if (dna[i] == 'A')
			dna[i] = 'T';
		else if (dna[i] == 'T')
			dna[i] = 'A';	
		else if (dna[i] == 'G')
			dna[i] = 'C';	
		else if (dna[i] == 'C')
			dna[i] = 'G';	
	}
    return dna;	
}
