//drd.h
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


#ifndef DRD_H
#define DRD_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <list>
#include <iomanip>
#include <stdio.h>
#include <map>
#include <cmath>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/dir.h>
#include <dirent.h>
#include <pthread.h>
#include "arg_struct.h"
#include "ribodef.h"
#include "alel.h"

using namespace std;

// const string default_inputfile = "sequence.fasta";
// const string prefix = "user/";
// const string query_output_filename = "query.html";
// const string query_update_filename = "query_update.html";

const char default_inputfile[] = "sequence.fasta";
const char prefix[] = "";
const char query_output_filename[] = "query.html";
const char query_update_filename[] = "query_update.html";

const int MAX_RIBOS = 500;

const int NUM_THREADS = 8;

// util.cpp

// Scramble a string by performing num passes of random swaps of each char.
string scrambleString(string str, int num);

// Uniformly randomly choose a character from nts.
char choose(const string& nts);

// Return reverse of string.
string stringrev(string str);

//Return standardized version of dna (No U's, all caps).
string standardize(const string& dna);

// Return complement of dna.
string complement(string dna);

// orfs.cpp

// Returns all significant open reading frames (length > threshold) across all forward reading frames.
// orfS[i] and orfF[i] contain the start and finish index of ORF i
int findOrfs(const string& seq, int *orfS, int *orfF, int threshold);

// Returns all significant open reading frames (length > threshold) across all reverse reading frames.
// revSeq is reversed prior to being passed in.
// orfS[i] and orfF[i] contain the start and finish index of ORF i (i starts after forward ORFs)
int findRevOrfs(const string& revSeq, int *orfS, int *orfF, int threshold, int nextIndex);

// align.cpp

const double v_mb = 10;        // match base
const double v_mp = 1;         // match ().
const double v_mmb = 5;	       // mismatch base
const double v_mmbp = -3;      // mismatch base with ().
const double v_mmp = 0;        // mismatch (). with ().
const double v_extendGap = 0;  // Extend Gap Penalty	    
const double v_openGap = -1;   // Open Gap Penalty
const double v_deletion = -10000;    // "infinite" penalty for matching deletions	

double similarity_score(char a , char b, double match, double mu);

double find_array_max(double array[], int length, int& index);

// Return all local alignments of two sequences with scores at least lowerbound.
list<AlEl> aligner(const string& sequence1, const string& sequence2, int lowerbound);

// Global alignment algorithm between Vienna strings
bool vienna_alignment(const string& seq1, const string& seq2, int lowerbound, string& alignment, double& score, int& identities, int& length);

//Global Alignment between shape strings
bool shape_alignment(const string& seq1, const string& seq2, int lowerbound, string& shapealignment, double& score, int& identities, int& length);

// Algorithm to convert Vienna string into shape sequence
// From Lorenz, W., Ponty, Y. and Clote, P. (2008). "Asymtotics of RNA shapes".
string Viennatoshape (const string& vienna, const string& shapetype);

// maxalign.cpp

string max_alignment(const string& in, int * st, int * en, int * tScore, const RiboDef& ribo, list<AlEl> &motifs, bool rev, int seed, ofstream * testfile);

// heaviestpath.cpp

list<int> find_heaviest(int ** adj_matrix, int matrix_size, int * score_array, list<AlEl> seqlists[], unsigned int all_positions[], const RiboDef& ribo);

// finder.cpp

void* finder(void* arguments);

#endif
