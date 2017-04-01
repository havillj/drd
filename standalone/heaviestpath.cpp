// heaviestpath.cpp
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

// Convert node number in adjacency matrix to (row, col) in layered graph.
void nodetoboard(int node_num, int * row, int *col, int num_columns, int *column_sizes)
{
	int sum = 0;
	for(int i = 0; i < num_columns; i++)
	{	
		sum += column_sizes[i];
		if (sum >= node_num + 1)
		{
			*row = i;
			*col = node_num - (sum - column_sizes[i]);
			return;
		}
	}
}

// Convert (row, col) in layered graph to node number in adjacency matrix.
int boardtonode(int row, int col, int num_columns, int *column_sizes)
{	
	int sum = 0;
	for(int i = 0; i < row; i++)
	{	
		sum += column_sizes[i];
	}
	return sum + col;
}

// Get the shortest path from among all optimal paths.
void getBestPath(int i, int j, list<int> **pathboard, int num_columns, int column_sizes[], unsigned int all_positions[], list<int> temp_path, list<int>& winning_path, int& min_length)
{
	temp_path.push_back(boardtonode(i, j, num_columns, column_sizes));    // find node corresponding to (i, j) and add it to the list
	if (!pathboard[i][j].empty())
	{
		list<int>::iterator it = pathboard[i][j].begin();
		while (it != pathboard[i][j].end())
		{
 			int pred = *it;                                          // get a predecessor node
 			int pred_i, pred_j;
 			nodetoboard(pred, &pred_i, &pred_j, num_columns, column_sizes);    // convert node to coordinates
 			getBestPath(pred_i, pred_j, pathboard, num_columns, column_sizes, all_positions, temp_path, winning_path, min_length);
 			it++;
 		}
 	}
 	else
 	{
 		if (abs(all_positions[temp_path.front()] - all_positions[temp_path.back()]) < min_length)
 		{
 			winning_path = temp_path;
			min_length = abs(all_positions[temp_path.front()] - all_positions[temp_path.back()]);
 		}
 	}
}

// heaviest path algorithm
list<int> find_heaviest(int ** adj_matrix, int matrix_size, int * score_array, list<AlEl> seqlists[], unsigned int all_positions[], const RiboDef& ribo)
{
	int **matrix = adj_matrix;
	int *scores = score_array;
	
	// Create and initialize scoreboard for the dynamic programming algorithm.
	// Column size is determined by the number of motifs.
	
	int num_columns = ribo.nummotifs;
	int *column_sizes = new int [num_columns];
	for (unsigned int i=0; i < ribo.nummotifs; i++)
		column_sizes[i] = seqlists[i].size();
	
	int **scoreboard = new int * [num_columns];
	list<int> **pathboard = new list<int> * [num_columns];	
	for (int i = 0; i < num_columns; i++)
	{
		scoreboard[i] = new int[column_sizes[i]];
		pathboard[i] = new list<int>[column_sizes[i]];
	}
	
	for (int i = 0; i < num_columns; i++)
	{
		for (int j = 0; j< column_sizes[i]; j++)
		{
			scoreboard[i][j] = 0;
//			pathboard[i][j].push_back(-1);
		}
	}
	
	// Here comes the algorithm... first score the matrix
	
	//Set up first column
	for(int i = 0; i < column_sizes[0]; i++)
		scoreboard[0][i] = scores[i];
	
	int dest = column_sizes[0];  // current node to which to find heaviest path
	int p_row, p_col;            // coordinates of current predecessor node
	int d_row, d_col;            // coordinates of current dest node
	for (int i = 1; i < num_columns; i++)
	{	
		for (int j = 0; j < column_sizes[i]; j++)
		{
			// Find best way to get to node dest.
			for (int pred = 0; pred < dest; pred++)
			{
				if (matrix[pred][dest] == 1)   // if there is an edge between nodes pred and dest
				{	
					nodetoboard(pred, &p_row, &p_col, num_columns, column_sizes);   // (x_row, x_col) = coordinates of node pred
					nodetoboard(dest, &d_row, &d_col, num_columns, column_sizes);   // (n_row, n_col) = coordinates of node dest
					
					if (scores[dest] + scoreboard[p_row][p_col] >= scoreboard[d_row][d_col])   // if we can get to dest via pred better than current best
					{	
						if (scores[dest] + scoreboard[p_row][p_col] > scoreboard[d_row][d_col])
						{
							scoreboard[d_row][d_col] = scores[dest] + scoreboard[p_row][p_col];                   // update score
							pathboard[d_row][d_col].clear();
						}
						pathboard[d_row][d_col].push_back(boardtonode(p_row, p_col, num_columns, column_sizes));  // add an optimal path
					}
				}
			}
			dest++;	
		}
	}
	
	//Find the maximum value in last column of scoreboard.
	int max = 0;
	for (int j = 0; j < column_sizes[num_columns - 1]; j++)
	{
		if (scoreboard[num_columns-1][j] > max)
			max = scoreboard[num_columns - 1][j];
	}
	
	// Find the shortest path with the maximum value.
	list<int> winning_path, temp_path;
	int min_length = 0x7FFFFFFF;
	for (int j = 0; j < column_sizes[num_columns - 1]; j++)
	{
 		if (scoreboard[num_columns-1][j] == max)
 		{
 			getBestPath(num_columns - 1, j, pathboard, num_columns, column_sizes, all_positions, temp_path, winning_path, min_length);
 		}
 	}

	delete [] column_sizes;
	
	for (int i = 0; i < num_columns; i++)
	{
		delete [] scoreboard[i];
		delete [] pathboard[i];
	}
	delete [] scoreboard;
	delete [] pathboard;
	
	return winning_path;
}
