// maxalign.cpp
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

int gaps(string s)
{
	int count = 0;
	for (size_t i = 0; i < s.length(); i++)
		if (s[i] == '-')
			count++;
	return count;
}

// Finds max total alignment for a sequence (snippet) with the given motif sequences.
string max_alignment(const string& in, int *st, int *en, int *tScore, const RiboDef& ribo, list<AlEl> &motifs, bool rev, int seed, ofstream *testfile)
{
	int start = 0;
	int end = 0;
 	string input = in;
	stringstream stream1;
	list<AlEl> seqlists[ribo.nummotifs];     // list of alignment elements matches for each motif
	list<AlEl> completeList;                 // list of all motif matches
	int matrix_size = 0;                     // total number of motifs found

	for (unsigned int i = 0; i < ribo.nummotifs; i++)    // for each motif...
	{
		seqlists[i] = aligner(input, ribo.motifs[i], ribo.minscores[i]);  // find list of matches for motif in input
		if (seqlists[i].size() == 0)
		{
			//cerr << "Missing motif sequence " << i + 1 << endl;
			
			if ((testfile != NULL) && (testfile->is_open()))
			{
				for (unsigned int i = 0; i < ribo.nummotifs + 2; i++)  // motifs + scores + length
					(*testfile) << " \t";
			}
		
			return string("");
		}
		matrix_size = matrix_size + seqlists[i].size();
	}

	// Aggregate all matches.
	
	unsigned int *all_positions = new unsigned int [matrix_size];
	for (int x = 0; x < matrix_size; x++)
	{
		all_positions[x] = 0;
	}
	list<AlEl>::iterator p;
	int n = 0;
	int *scoresarray = new int[matrix_size];
	for (unsigned int i=0; i < ribo.nummotifs; i++)
	{
		p = seqlists[i].begin();
		while (p != seqlists[i].end())
		{
			p->index = n;
			scoresarray[n] = p->score;
			completeList.push_back(*p);
			all_positions[n] = p->position;
			p++;
			n++;
		}
	}
	
	// Set up the adjacency matrix for heaviest path algorithm.
	
	int **adjacency_matrix = new int*[matrix_size];
	for (int i = 0; i < matrix_size; i++)
	{
		adjacency_matrix[i] = new int[matrix_size];
	}
	for (int i = 0; i < matrix_size; i++)
		for(int j = 0; j < matrix_size; j++)
			adjacency_matrix[i][j] = 0;
	  
	list<AlEl> fillLists[ribo.nummotifs - 1];
	fillLists[0] = seqlists[0];
	list<AlEl>::iterator q;
	for (unsigned int i=0; i < ribo.nummotifs - 1; i++)
	{
		p = fillLists[i].begin();
		while (p != fillLists[i].end())
		{
			q = seqlists[i+1].begin();
			while (q != seqlists[i+1].end())
			{
				if ((p->position < q->position) && (q->position - (p->position + (p->str).length() / 2 - 1) <= ribo.distances[i + 1]) && (q->position - (p->position + (p->str).length() / 2 - 1) >= ribo.mindistances[i]))
				{
					if (i + 1 < ribo.nummotifs - 1)
					{
						fillLists[i + 1].push_back(*q);
					}
					adjacency_matrix[p->index][q->index] = 1;
					adjacency_matrix[q->index][p->index] = 1;
				}
				q++;
			}
			p++;
		}
	}

	unsigned int maxscore = 0;
	for (unsigned int i=0; i < ribo.nummotifs; i++)
	{
		maxscore += ribo.motifs[i].size();
	}
	
	// Find the best sequence of motifs.
	
 	list<int> rev_motif_locations =  find_heaviest(adjacency_matrix, matrix_size, scoresarray, seqlists, all_positions, ribo); 
	
	// Report results
	if (rev_motif_locations.size() < ribo.nummotifs)
	{
		//cerr << "No heaviest path was found" << endl;
		delete [] all_positions;
		delete [] scoresarray;
		for(int i = 0; i < matrix_size;i++)
			delete [] adjacency_matrix[i];
		delete [] adjacency_matrix;
		
		if ((testfile != NULL) && (testfile->is_open()))
		{
			for (unsigned int i = 0; i < ribo.nummotifs + 2; i++)
				(*testfile) << " \t";
		}
		
		return string("");
	}

	// Create a string containing the alignment between the sequence and the motifs.
	
	int M[16][16]= {{1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0},   // A
			    {0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0},   // C    (Candidate strand)
			    {0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0},   // G
			    {0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0},   // T
			    {1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0},   // R
			    {0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0},   // Y    (Candidate strand)
			    {0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0},   // K
			    {1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0},   // M
			    {0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0},   // S
			    {1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0},   // W   (Candidate strand)
			    {0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0},   // B
			    {1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0},   // D
			    {1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0},   // H
			    {1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0},   // V    (Candidate strand)
			    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0},   // N
			    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} }; // -
			    
	map<char,int> nucToInt;
	nucToInt['A'] = 0;
	nucToInt['C'] = 1;
	nucToInt['G'] = 2;
	nucToInt['T'] = 3;
	nucToInt['R'] = 4;
	nucToInt['Y'] = 5;
	nucToInt['K'] = 6;
	nucToInt['M'] = 7;
	nucToInt['S'] = 8;
	nucToInt['W'] = 9;
	nucToInt['B'] = 10;
	nucToInt['D'] = 11;
	nucToInt['H'] = 12;
	nucToInt['V'] = 13;
	nucToInt['N'] = 14;
	nucToInt['-'] = 15;
	
	list<AlEl>::iterator vit;
	list<int>::iterator kit = rev_motif_locations.end();
	unsigned int totalScore = 0;
	int numMotif = 1;	
	int pos;
	string motifstr, matches, rseq;
	int lastpos = 0, firstpos;
	int lastlen = 0;
	int indels = 0;
						
	while (kit != rev_motif_locations.begin())
	{
		kit--;
		vit = completeList.begin();
		while (vit != completeList.end())
		{
			if (vit->index == *kit)
			{
				pos = vit->position;
				if (numMotif > 1)
				{
					if (pos - (lastpos + lastlen - indels) > 0)
					{
						rseq += in.substr(lastpos + lastlen - indels, pos - (lastpos + lastlen - indels));
						matches.append(pos - (lastpos + lastlen - indels), ' ');
						motifstr.append(pos - (lastpos + lastlen - indels), ' ');
					}
					
					if ((testfile != NULL) && (testfile->is_open()))
						(*testfile) << pos - (lastpos + lastlen - indels) << "\t";
				}
				else
					firstpos = vit->position;
				motifstr += vit->str.substr(((vit->str).length()/2));
				motifstr.erase(motifstr.length() - 1);
				rseq += vit->str.substr(0,((vit->str).length()/2)-1);
				indels = gaps(vit->str.substr(0,((vit->str).length()/2)-1));
				for (unsigned int i = 0; i < (vit->str).length()/2 - 1; i++)
				{
					if (M[nucToInt[vit->str.substr(((vit->str).length()/2))[i]]][nucToInt[vit->str.substr(0,((vit->str).length()/2)-1)[i]]] == 1)
						matches += '|';
					else
						matches += ' ';
				}
				lastpos = pos;
				lastlen = (vit->str).length()/2 - 1;  // (vit->str.substr(((vit->str).length()/2))).length() - 1;
				totalScore = totalScore+vit->score;
				motifs.push_back(*vit);     // Populates a list of motifs in AlEl format
				numMotif++;
			}
			vit++;	
		}
	}
	
	*tScore = totalScore;
	
	if (totalScore < ribo.lowestTotal)
	{	
		//cerr << "Total score too low" << endl;
		delete [] all_positions;
		delete [] scoresarray;
		for(int i = 0; i < matrix_size;i++)
			delete [] adjacency_matrix[i];
		delete [] adjacency_matrix;
		return string("");
	}
				
	// Set putative start and end positions of riboswitch
	
 	int indeX = rev_motif_locations.front();
 	int lastmotifcap = lastlen - indels + ribo.distances[ribo.nummotifs];
	if (all_positions[indeX] + lastmotifcap - 1 > input.length() - 1)
	{
		end = input.length() - 1;
	}
	else
	{
		end = all_positions[indeX] + lastmotifcap - 1;
	}
	
	indeX = rev_motif_locations.back();
	if (all_positions[indeX] < ribo.distances[0])
	{
		start = 0;
	}
	else
	{
		start = all_positions[indeX] - ribo.distances[0];
	}
	
	*st = start;
	*en = end;
	
	motifstr.insert(0, firstpos - start, ' ');
	rseq = in.substr(start, firstpos - start) + rseq + in.substr(lastpos + lastlen - indels, end - (lastpos + lastlen - indels) + 1);
	matches.insert(0, firstpos - start, ' ');
		
	delete [] all_positions;
	delete [] scoresarray;
	for(int i = 0; i < matrix_size;i++)
		delete [] adjacency_matrix[i];
	delete [] adjacency_matrix;
	
	return rseq + "\n" + matches + "\n    " + motifstr;
}
