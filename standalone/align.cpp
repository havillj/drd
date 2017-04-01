// align.cpp
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

double similarity_score(char a , char b, double match, double mu)
{
	if (a == b)
		return match;
	return mu;
}

double find_array_max(double array[], int length, int& index)
{
	double max = array[0];            // start with max = first element
	index = 0;
	for(int i = 1; i < length; i++)   // finds the max value (and index) in an array
	{
		if(array[i] >= max)
		{
			max = array[i];
			index = i; 
		}
	}
	return max;                    // return highest value in array
}

int computeIdentities(const string& seq1, const string& seq2, string& matches)
{
	//               A  C  G  T  R  Y  K  M  S  W  B  D  H  V  N  (  )  .  -
	int M[19][19]= {{1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0},   // A
			        {0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0},   // C    (Candidate strand)
			        {0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0},   // G
			        {0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0},   // T
			        {1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0},   // R
			        {0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0},   // Y    (Candidate strand)
			        {0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0},   // K
			        {1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0},   // M
			        {0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0},   // S
			        {1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0},   // W   (Candidate strand)
			        {0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0},   // B
			        {1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0},   // D
			        {1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0},   // H
			        {1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0},   // V    (Candidate strand)
			        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0},   // N
			        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},   // (
			        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},   // )
			        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0},   // .
			        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};  // -
			        
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
	nucToInt['('] = 15;
	nucToInt[')'] = 16;
	nucToInt['.'] = 17;
	nucToInt['-'] = 18;
	
	int n = seq1.length();
	int score = 0;
	matches = "";
	for (int i = 0; i < n; i++)
	{
		score += M[nucToInt[seq1[i]]][nucToInt[seq2[i]]];
		matches += (M[nucToInt[seq1[i]]][nucToInt[seq2[i]]] == 1 ? "|" : "&nbsp;");
	}
	return score;
}


// Return all local alignments of two sequences with scores at least lowerbound.
list<AlEl> aligner(const string& sequence1, const string& sequence2, int lowerbound)
{
	double delta = 10;            // indel penalty
	
	string seq_a = sequence1;
	string seq_b = sequence2;
	
	int N_a = seq_a.length();     // get the actual lengths of the sequences
	int N_b = seq_b.length();
	
	double H[N_a + 1][N_b + 1];   // Initialize dynamic programming matrix H
	int I_i[N_a + 1][N_b + 1], 
	    I_j[N_a + 1][N_b + 1];    // Index matrices to remember the 'path' for backtracking
	    
	H[0][0] = 0;
	I_i[0][0] = I_j[0][0] = 0;
	
	for (int i=1; i <= N_a; i++)
	{	
		H[i][0] = 0;
		I_i[i][0] = 0;
	}
	for (int j=1; j <= N_b; j++)
	{	
		H[0][j] = 0;
		I_i[0][j] = 0;
	}

	//Scoring Matrix
	//                    A   C   G   T   R   Y   K   M   S   W   B   D   H   V   N 
	int nuc44[15][15]= {{10,  1,  1,  1, 10,  1,  1, 10,  1, 10,  1, 10, 10, 10, 10},   // A
			    { 1, 10,  1,  1,  1, 10,  1, 10, 10,  1, 10,  1, 10, 10, 10},   // C    (Candidate strand)
			    { 1,  1, 10,  1, 10,  1, 10,  1, 10,  1, 10, 10,  1, 10, 10},   // G
			    { 1,  1,  1, 10,  1, 10, 10,  1,  1, 10, 10, 10, 10,  1, 10},   // T		    
			    {10,  1, 10,  1, 10,  1,  1,  1,  1,  1,  1, 10,  1, 10, 10},   // R
			    { 1, 10,  1,  1,  1, 10,  1,  1,  1,  1, 10,  1, 10,  1, 10},   // Y    (Candidate strand)
			    { 1,  1, 10,  1,  1,  1, 10,  1,  1,  1, 10, 10,  1,  1, 10},   // K
			    { 1,  1,  1, 10,  1,  1,  1, 10,  1,  1,  1,  1, 10, 10, 10},   // M
			    {10,  1,  1,  1,  1,  1,  1,  1, 10,  1, 10,  1,  1, 10, 10},   // S
			    { 1, 10,  1,  1,  1,  1,  1,  1,  1, 10,  1, 10, 10,  1, 10},   // W    (Candidate strand)	
			    { 1, 10, 10, 10,  1, 10, 10,  1, 10,  1, 10,  1,  1,  1, 10},   // B
			    {10,  1, 10, 10, 10,  1, 10,  1,  1, 10,  1, 10,  1,  1, 10},   // D
			    {10, 10,  1, 10,  1, 10, 10,  1,  1, 10,  1,  1, 10,  1, 10},   // H
			    {10, 10, 10,  1, 10,  1,  1, 10, 10,  1,  1,  1,  1, 10, 10},   // V    (Candidate strand)
			    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0} }; // N

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

	double temp[4];                                      // holds the four possible values for each matrix position
	
	// Smith-Waterman algorithm
	
	for(int i = 1; i <= N_a; i++)
	{
		for(int j = 1; j <= N_b; j++)
		{
			
			temp[0] = H[i-1][j-1] + nuc44[ nucToInt[ seq_a[i-1] ] ][ nucToInt[ seq_b[j-1] ] ];
			temp[1] = H[i-1][j] - delta;                  
			temp[2] = H[i][j-1] - delta;                 
			temp[3] = 0;
			
			int index;
			H[i][j] = find_array_max(temp, 4, index);
			
			switch(index)
			{
				case 0:                                  // score in (i,j) stems from a match/mismatch
					I_i[i][j] = i-1;
					I_j[i][j] = j-1;
					break;
				case 1:                                  // score in (i,j) stems from a deletion in sequence A
					I_i[i][j] = i-1;
					I_j[i][j] = j;
					break;
				case 2:                                  // score in (i,j) stems from a deletion in sequence B
					I_i[i][j] = i;
					I_j[i][j] = j-1;
					break;
				case 3:                                  // (i,j) is the beginning of a subsequence
					I_i[i][j] = i;
					I_j[i][j] = j;	
					break;
			}
		}
	}

	// Get all nonzero alignments

	stringstream stream;             // for recording the alignment as a string
	list<AlEl> listing;              // list of alignments
	AlEl element;                    // individual alignment element
	list<AlEl>::iterator it, it2;
	
	int minscore = 10 * lowerbound;

	for(int i = 1; i <= N_a; i++)
	{
		for(int j = 1; j <= N_b; j++)
		{
			if (H[i][j] >= minscore)
			{
				stream.str("");	
		
				int current_i = i, 
				    current_j = j;
				int next_i = I_i[current_i][current_j];
				int next_j = I_j[current_i][current_j];
				int tick = 0;
		
				char consensus_a[N_a + N_b + 2],     // two strings hold the alignment
				     consensus_b[N_a + N_b + 2];
		
				while(((current_i != next_i) || (current_j != next_j)) && (current_i != 0) && (current_j != 0))
				{
					if(next_i == current_i)  
						consensus_a[tick] = '-';                    // deletion in A
					else                   
						consensus_a[tick] = seq_a[current_i - 1];   // match/mismatch in A
			
					if(next_j == current_j)  
						consensus_b[tick] = '-';                    // deletion in B
					else                   
						consensus_b[tick] = seq_b[current_j - 1];   // match/mismatch in B
			
			
					H[current_i][current_j] = 0;  // reset position in H so not found again
	
					current_i = next_i;
					current_j = next_j;
					next_i = I_i[current_i][current_j];
					next_j = I_j[current_i][current_j];
					tick++;
				}
	
				element.position = current_i;
		
				string con_a;
				string con_b;
		
				for(int i = tick - 1;i >= 0; i--)
				{
					stream << consensus_a[i];
					con_a = con_a + consensus_a[i];
				}
				stream << endl;
		
				for(int j = tick - 1; j >= 0; j--)
				{
					stream << consensus_b[j];
					con_b = con_b + consensus_b[j];
				}
				stream<<endl;

				string matches;
				element.score = computeIdentities(con_a, con_b, matches);
				element.str = stream.str();
				
// 				ofstream f;
// 				string fn = prefix + "blah";
// 				f.open(fn.c_str(), ios::app);
// 				f << "-----BEGIN-----" << endl;
// 				f << element.str << endl;
// 				f << element.score << "/" << lowerbound << " " << element.position << endl;
// 				f << "-----END-----" << endl;
// 				f.close();

				if (element.score >= lowerbound) // and con_b.length() == (unsigned int)N_b)   // if high enough score, add to listing
				{	
					listing.push_back(element);
				}
			}
		}
	}
	return listing;
}

const string tab = "      ";

// Global alignment algorithm between Vienna strings
bool vienna_alignment(const string& seq1, const string& seq2, int lowerbound, string& alignment, double& score, int& identities, int& length)
{
	int N_1 = seq1.length();                       // Get the lengths of the strings to be aligned
	int N_2 = seq2.length();
	
	double H[N_1 + 1][N_2 + 1];                    // Dynamic programming matrix
	int I_i[N_1 + 1][N_2 + 1], 
	    I_j[N_1 + 1][N_2 + 1];     // Index matrices to remember the 'path' for backtracking
	    
	H[0][0] = 0;
	I_i[0][0] = I_j[0][0] = 0;

	//Mapping of characters to integers (used for scoring matrix)
 	map<char,int> charToInt;
	charToInt['A'] = 0;
	charToInt['C'] = 1;
	charToInt['G'] = 2;
	charToInt['T'] = 3;
	charToInt['('] = 4;
	charToInt[')'] = 5;
	charToInt['.'] = 6;
	charToInt['-'] = 7;
	charToInt['R'] = 8;
	charToInt['Y'] = 9;
	charToInt['K'] = 10;
	charToInt['M'] = 11;
	charToInt['S'] = 12;
	charToInt['W'] = 13;
	charToInt['B'] = 14;
	charToInt['D'] = 15;
	charToInt['H'] = 16;
	charToInt['V'] = 17;
	charToInt['N'] = 18;
	
	// values defined in drd.h
	double p = v_mb;                // match base
	double q = v_mp;                // match ().
	double o = v_mmb;	            // mismatch base
	double n = v_mmbp;              // mismatch base with ().
	double m = v_mmp;               // mismatch (). with ().
	double extendGap = v_extendGap; // Extend Gap Penalty	    
	double openGap = v_openGap;     // Open Gap Penalty
	double z = v_deletion;          // "infinite" penalty for matching deletions	

	//Scoring Matrix       	   A  C  G  T  (  )  .  -  R  Y  K  M  S  W  B  D  H  V  N
	double Scoring[19][19]= {{ p, o, o, o, n, n, n, n, p, o, o, p, o, p, o, p, p, p, p},   // A
			    			 { o, p, o, o, n, n, n, n, o, p, o, p, p, o, p, o, p, p, p},   // C    
			    	 		 { o, o, p, o, n, n, n, n, p, o, p, o, p, o, p, p, o, p, p},   // G
			    			 { o, o, o, p, n, n, n, n, o, p, p, o, o, p, p, p, p, o, p},   // T   (Candidate strand)
			    	 		 { n, n, n, n, q, m, m, n, n, n, n, n, n, n, n, n, n, n, n},   // (
			    	 		 { n, n, n, n, m, q, m, n, n, n, n, n, n, n, n, n, n, n, n},   // )
			    			 { n, n, n, n, m, m, q, n, n, n, n, n, n, n, n, n, n, n, n},   // .
				 			 { n, n, n, n, n, n, n, z, n, n, n, n, n, n, n, n, n, n, n},   // -
				 			 { p, o, p, o, n, n, n, n, p, o, o, o, o, o, o, p, o, p, p},   // R
			    			 { o, p, o, p, n, n, n, n, o, p, o, o, o, o, p, o, p, o, p},   // Y    
			    			 { o, o, p, p, n, n, n, n, o, o, p, o, o, o, p, p, o, o, p},   // K
			    			 { p, p, o, o, n, n, n, n, o, o, o, p, o, o, o, o, p, p, p},   // M   (Candidate strand)
			    			 { o, p, p, o, n, n, n, n, o, o, o, o, p, o, p, o, o, p, p},   // S
			    			 { p, o, o, p, n, n, n, n, o, o, o, o, o, p, o, p, p, o, p},   // W
			    			 { o, p, p, p, n, n, n, n, o, p, p, o, p, o, p, o, o, o, n},   // B
				 			 { p, o, p, p, n, n, n, n, p, o, p, o, o, p, o, p, o, o, n},   // D
				 	 		 { p, p, o, p, n, n, n, n, o, p, p, o, o, p, o, o, p, o, p},   // H
			    	 		 { p, p, p, o, n, n, n, n, p, o, o, p, p, o, o, o, o, p, p},   // V    
			    	 		 { p, p, p, p, n, n, n, n, p, p, p, p, p, p, p, p, p, p, p} }; // N

	for (int i=1; i <= N_1; i++)
	{	
		H[i][0] = H[i-1][0] + n;                // Initialize matrix values for global alignment
		I_i[i][0] = i - 1;
		I_j[i][0] = 0;
	}
	for (int j=1; j <= N_2; j++)
	{	
		H[0][j] = H[0][j-1] + n;                // Initialize matrix values for global alignment
		I_i[0][j] = 0;
		I_j[0][j] = j - 1;
	}
	
 	// Global Alignment algorithm
 	
	double temp[3];                      // holds the three possible values for each matrix position

	for(int i = 1; i <= N_1; i++)                        //
	{                                                    // For each matrix cell...
		for(int j = 1; j <= N_2; j++)                    //
		{
			char tempchar1 = seq1[i-1];
			char tempchar2 = seq2[j-1];

			temp[0] = H[i-1][j-1] + Scoring[ charToInt[tempchar1] ][ charToInt[tempchar2] ]; // Match/mismatch
			
			if (i >= 1 and I_i[i-1][j] == i-2 and I_j[i-1][j] == j)    //Del in Seq 1 with extend gap
				temp[1] = H[i-1][j] + extendGap;
			else temp[1] = H[i-1][j] + openGap;			 //Del in Seq 1 with open gap

			if (j >=1 and I_i[i][j-1] == i and I_j[i][j-1] == j-2)      //Del in Seq 1 with extend gap
				temp[2] = H[i][j-1] + extendGap;
			else temp[2] = H[i][j-1] + openGap;		         //Del in Seq 1 with open gap


			int index;
			H[i][j] = find_array_max(temp, 3, index);    // Find the maximum of the above cases
		
			if (index == 0)            // score in (i,j) stems from a match/mismatch
			{
				I_i[i][j] = i-1;
				I_j[i][j] = j-1;
			}
			else if (index == 1)       // score in (i,j) stems from a deletion in sequence 1
			{
				I_i[i][j] = i-1;
				I_j[i][j] = j;
			}
			else                       // score in (i,j) stems from a deletion in sequence 2
			{
				I_i[i][j] = i;
				I_j[i][j] = j-1;
			}
		}
	}

	// Backtrack from end of both strings to get the maximal alignment
	int current_i = N_1,
		current_j = N_2;
	int next_i = I_i[current_i][current_j];
	int next_j = I_j[current_i][current_j];
	int tick = 0;

	char alignment_1[N_1 + N_2 + 2],     // two strings hold the alignment
		 alignment_2[N_1 + N_2 + 2];
	
	while((current_i != 0) || (current_j != 0))
	{	

		if (next_i == current_i)  
			alignment_1[tick] = '-';                    // deletion in 1
		else
			alignment_1[tick] = seq1[current_i - 1];    // match/mismatch in 1
		
		if (next_j == current_j)  
			alignment_2[tick] = '-';                    // deletion in 2
		else
			alignment_2[tick] = seq2[current_j - 1];    // match/mismatch in 2
		
		current_i = next_i;
		current_j = next_j;
		next_i = I_i[current_i][current_j];
		next_j = I_j[current_i][current_j];
		tick++;
	}

	alignment_1[tick] = '\0';                           // "ends" the array of chars
	alignment_2[tick] = '\0';
	
	//Reverse the alignment strings
	int size1 = (string(alignment_1)).size();           
	int size2 = (string(alignment_2)).size();           
	                                                    
	int halfsize1 = int(size1/2);                       
	int halfsize2 = int(size2/2);                       
	
	for (int i=0; i < halfsize1 ; i++)                 
	{                                                      
		char temp = alignment_1[i];                         
		alignment_1[i] = alignment_1[size1 - 1 - i];    
		alignment_1[size1 - 1 - i] = temp;              
	}                                                   
	                                                    
	for (int j=0; j < halfsize2 ; j++)                 
	{                                                   
		char temp2 = alignment_2[j];                    
		alignment_2[j] = alignment_2[size2 - 1 - j];    
		alignment_2[size2 - 1 - j] = temp2;             
	}	                                                
	
	//Assign alignment, score
	
	string matches;
	identities = computeIdentities(string(alignment_1), string(alignment_2), matches);
	
	// alignment string for web display
	alignment = "<table><tr><td><tt>" + string(alignment_1) + "</tt>&nbsp;&nbsp; (candidate)</td></tr><tr><td><tt>" + matches + "</tt></td></tr><tr><td><tt>" + string(alignment_2) + "</tt>&nbsp;&nbsp; (consensus)</td></tr></table>";
	
	length = strlen(alignment_1);
	
	score = H[N_1][N_2];
				
	return (score >= lowerbound);
}


//Global Alignment between shape strings
bool shape_alignment(const string& seq1, const string& seq2, int lowerbound, string& shapealignment, double& score, int& identities, int& length)
{
	double match = +3;             // Match     
	double mu    = -2;             // Mismatch  
	double in    =  0;             // Insertion  
	double del   =  0;             // Deletion   

	int N_1 = seq1.length();       // Get the sequence lengths
	int N_2 = seq2.length();
	
	double H[N_1 + 1][N_2 + 1];
	
	int I_i[N_1 + 1][N_2 + 1], 
	    I_j[N_1 + 1][N_2 + 1];     // Index matrices to remember the 'path' for backtracking
	    
	H[0][0] = 0;
	I_i[0][0] = I_j[0][0] = 0;
	
	//Fill in first row and column
	for (int i=1; i <= N_1; i++)
	{	
		H[i][0] = -3*i;
		I_i[i][0] = i-1;
		I_j[i][0] = 0;
	}
	for (int j=1; j <= N_2; j++)
	{	
		H[0][j] = -3*j;
		I_i[0][j] = 0;
		I_j[0][j] = j-1;
	}
	
	double temp[3];          //Holds three possible values for recursive definition

	//Create 2D table of results
	for(int i = 1; i <= N_1; i++)
	{
		for(int j = 1; j <= N_2; j++)
		{
			temp[0] = H[i-1][j-1] + similarity_score(seq1[i-1], seq2[j-1], match, mu); 
			temp[1] = H[i-1][j] + in;                  
			temp[2] = H[i][j-1] + del;                 
			
			int index;
			H[i][j] = find_array_max(temp, 3, index);

			if (index == 0)                   // score in (i,j) stems from a match/mismatch
			{
				I_i[i][j] = i-1;
				I_j[i][j] = j-1;
			}
			else if (index == 1)              // score in (i,j) stems from a deletion in sequence A
			{
				I_i[i][j] = i-1;
				I_j[i][j] = j;
			}
			else if (index == 2)              // score in (i,j) stems from a deletion in sequence B
			{
				I_i[i][j] = i;
				I_j[i][j] = j-1;
			}
		}
	}
	
	double H_max = H[N_1][N_2];
	stringstream stream;
	stream.str("");	

	// Backtrack from H_max to get the maximal alignment
	
	int current_i =  N_1, 
		current_j =  N_2;
	int next_i = I_i[current_i][current_j];
	int next_j = I_j[current_i][current_j];
	int tick = 0;
	
	char consensus_a[N_1 + N_2 + 3],     // two strings hold the alignment
		 consensus_b[N_1 + N_2 + 3];
	
	while (current_i != 0 || current_j !=0)
	{
		if(next_i == current_i)  
			consensus_a[tick] = '-';                    // deletion in A
		else                   
			consensus_a[tick] = seq1[current_i - 1];    // match/mismatch in A
		
		if(next_j == current_j)  
			consensus_b[tick] = '-';                    // deletion in B
		else                   
			consensus_b[tick] = seq2[current_j - 1];    // match/mismatch in B
		
		current_i = next_i;
		current_j = next_j;

		next_i = I_i[current_i][current_j];
		next_j = I_j[current_i][current_j];

		tick++;
	}
	consensus_a[tick] = '\0';
	consensus_b[tick] = '\0';
	
	string con_a;                             // Initialize strings to hold alignments
	string con_b;
	
	stream << tab;

	stream << "Candidate: ";
	
	for(int i = tick - 1;i >= 0; i--)
	{
		stream << consensus_a[i];
		con_a = con_a + consensus_a[i];       // Populate strings with characters
	}
	stream << endl << tab << "Consensus: ";
	
	for(int j = tick - 1; j >= 0; j--)
	{
		stream << consensus_b[j];
		con_b = con_b + consensus_b[j];
	}
	stream<<endl;
	
	score = H_max;                            // Score of alignment

	string matches;
	identities = computeIdentities(con_a, con_b, matches);
	
	shapealignment = "<table><tr><td><tt>" + con_a + "</tt>&nbsp;&nbsp; (candidate)</td></tr><tr><td><tt>" + matches + "</tt></td></tr><tr><td><tt>" + con_b + "</tt>&nbsp;&nbsp; (consensus)</td></tr></table>";
	
	length = con_a.length();
	
	return (score >= lowerbound);             // Only return True if score is at least lowerbound
}

// Algorithm to convert Vienna string into shape sequence
// From Lorenz, W., Ponty, Y. and Clote, P. (2008). "Asymtotics of RNA shapes".
string Viennatoshape (const string& vienna, const string& shapetype)
{
	string viennanew;
	
	if (shapetype == "short")                   // "short" means dots are omitted: [[][][]]
	{
		for (unsigned int i=0; i<vienna.size(); i++)
		{
			if (vienna[i] == '(')               // Convert parenthesis into brackets and ignore dots
				viennanew += "[";
			else if (vienna[i] == ')')
				viennanew += "]";
		}
	}
	else if (shapetype == "long")               // "long" means dots are condensed: [.[.][.].[.]].
	{
		for (unsigned int i=0; i<vienna.size(); i++)
		{
			if ((vienna[i] == '.') && (vienna[i+1] != '.'))
				viennanew += ".";
			else if (vienna[i] == '(')          // Convert parenthesis into brackets and condense dots
				viennanew += "[";
			else if (vienna[i] == ')')
				viennanew += "]";
		}
	}
	int n = viennanew.size();
	if (n<=2)
	{
		viennanew = ".";
		return viennanew;
	}
	list<unsigned int> stackpairs;
	int A[n+1];
	for (int i=1; i<=n; i++)
	{
		A[i] = 0;
		if (viennanew[i-1] == '[')             // Use a stack to condense groups of brackets
			stackpairs.push_back(i);           // e.g. "[[[[][]]]]" becomes "[[][]]"
		else if (viennanew[i-1] == ']')
		{
			A[stackpairs.back()] = i;
			stackpairs.pop_back();
		}
	}
	
	int x = 0;
	int y = 0;
	int j;

	for (int i=1; i<=n; i++)
	{
		if (A[i] > 0)
		{
			j = A[i];
			if ((i == x+1) && (A[i] == y-1))
			{
				x = i;
				y = j;
				viennanew[i-1] = '-';          // Mark redundant brackets for deletion
				viennanew[j-1] = '-';
			}
			else
			{
				x = i;
				y = j;
			}
		}
		else
		{
			x = 0;
			y = 0;
		}
	}
	string viennanewer;
	for (int i=0; i<n; i++)
	{
		if (viennanew[i] != '-')
			viennanewer += viennanew[i];       // Add only nonredundant brackets to viennanewer
	}
	return viennanewer;
}
