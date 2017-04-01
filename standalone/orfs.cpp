// orfs.cpp
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

// Returns all significant open reading frames (length > threshold) across all forward reading frames.
// orfS[i] and orfF[i] contain the start and finish index of ORF i
int findOrfs(const string& seq, int *orfS, int *orfF, int threshold)
{
	int s = 0;
	int seqLength = seq.length();
	
	for (int i = 0; i < 3; i++)  // start of reading frame
	{
		int start = -1, stop;
		for (int j = i; j < seqLength - 2; j += 3)
		{
			if (start == -1 and seq.substr(j, 3) == "ATG")
			{
				start = j;
			}
			else if (start >= 0 and (seq.substr(j, 3) == "TAA" or seq.substr(j, 3) == "TAG" or seq.substr(j, 3) == "TGA"))
			{
				stop = j;
				if (stop - start + 3 >= threshold)
				{
					orfS[s] = start;
					orfF[s] = stop + 2;
					s++;
				}
				start = -1;
			}
		}
	}
	return s;
}

// Returns all significant open reading frames (length > threshold) across all reverse reading frames.
// revSeq is reversed prior to being passed in.
// orfS[i] and orfF[i] contain the start and finish index of ORF i (i starts after forward ORFs)
int findRevOrfs(const string& revSeq, int *orfS, int *orfF, int threshold, int nextIndex)
{
	int s = nextIndex;
	int seqLength = revSeq.length();
	for (int i = 0; i < 3; i++)
	{
		int startR = -1, stopR;
		for (int j = i; j < seqLength - 2; j += 3)
		{
			if (startR == -1 and revSeq.substr(j, 3) == "ATG")
			{
				startR = j;
			}
			else if (startR >=0 and (revSeq.substr(j, 3) == "TAA" or revSeq.substr(j, 3) == "TAG" or revSeq.substr(j, 3) == "TGA"))
			{
				stopR = j;
				if (stopR - startR + 3 >= threshold)
				{
					orfS[s] = seqLength - stopR - 3;
					orfF[s] = seqLength - startR - 1;
					s++;
				}
				startR = -1;
			}
		}
	}
	return s - nextIndex;
}
