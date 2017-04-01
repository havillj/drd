// finder.cpp
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


int numRibo = 0;  // number of riboswitches found so far
struct Riboswitch
{
	int start;
	int end;
	int score;
	int viennaid;
	bool reverse;         // False = forward strand, True = reverse strand
};
Riboswitch ribos[MAX_RIBOS];   // array of riboswitch boundaries for duplicate checks

pthread_mutex_t numLock = PTHREAD_MUTEX_INITIALIZER;   // guards numRibo global variable
pthread_mutex_t riboLock = PTHREAD_MUTEX_INITIALIZER;  // guards ribos global variable

// Find all matches in a sequence containing many snippets.
// Results are written to output files.
void* finder(void* arguments)
{
	arg_struct* args = (arg_struct*) arguments;
	string forwardSeq = args->forwSeq;
	string revSeq = args->backSeq;
	int curseed = args->seed;
	RiboDef ribo(*(args->riboswitch));
	string foldername = args->foldername;
	string path = foldername;	
	int* orfS = args->orfS;
	int* orfF = args->orfF;
	int numOrfs = args->numOrfs;
	int numROrfs = args->numROrfs;
	string fullSeq = args->fullSeq;
	int seqLength = fullSeq.length();
	ofstream *testfile = args->testfile;
	
	bool rev = 0;
	int totalSnippets;
	string snippet;                 // current snippet
	unsigned int starts = 0, ends;  // relative start and end of current snippet
	int seed;                       // absolute start of current snippet
	
	//ofstream queryUpdateFile;
	//char queryupdatefname[256]; // = prefix + foldername + "/" + query_update_filename;
	//char path[256];
	//strncpy(path, prefix, 255);
	//strncat(path, foldername.c_str(), 255);
	//strncat(path, "/", 255);
	//strncpy(queryupdatefname, path, 255);
	//strncat(queryupdatefname, query_update_filename, 255);
	char resultfname[256];  
	ofstream resultFile;
	char cmd[256];
	char seqfoldername[256];
	char imgfoldername[256]; // = prefix + foldername + "/Ribo/";  // folder in which to place fold images
	strncpy(imgfoldername, path.c_str(), 255);
	strncat(imgfoldername, "/images/", 255);
	char motif_clean[256];
	string motif;
	string seq = forwardSeq;
	string found;
	string riboswitch;
	char riboFilename[256];
	char seqfname[256], seqfname_short[256];
	char viennafname[256];
	string vienna_string;
	string vienna_string_motifs;
	string alignment1;
	string alignment;
	double viennascore, viennascore1;
	int viennaid, viennalen, viennaid1, viennalen1;
	string shape_string;   // make shape string (Candidate)
	string lit_shape;        // make shape string (Literature)
	string shape_align;
	double shape_score;
	int shape_id, shape_len;
	//char uniqueDirname[256];
	
	int found_something = 0;

	//Calculate the total number of snippets for the thread
	
	if (forwardSeq.length() == 0)
		return NULL;
	else if (forwardSeq.length() <= ribo.snippet_size)
		totalSnippets = 1;
	else if ((forwardSeq.length() - ribo.overlap_size) % (ribo.snippet_size - ribo.overlap_size) == 0)
		totalSnippets = (forwardSeq.length() - ribo.overlap_size) / (ribo.snippet_size - ribo.overlap_size);
	else
		totalSnippets = ((forwardSeq.length() - ribo.overlap_size) / (ribo.snippet_size - ribo.overlap_size)) + 1;
		
	unsigned int maxscore = 0;
	for (unsigned int i = 0; i < ribo.nummotifs; i++)
	{
		maxscore += ribo.motifs[i].size();
	}
	
	//Once for the forward strand, once for the backwards strand
	for (int back = 0; back < 2; back++)
	{
		if (back == 1)
		{
			if (revSeq.length() == 0)
				break;
			rev = 1;
			seq = revSeq;
			starts = 0;
		}
	
		for (int count = 1; count <= totalSnippets; count++)        // Iterate over each snippet
		{                                                           // Example: SS=700, O=200
			 ends = starts + ribo.snippet_size;                     // ends = starts + 700
			 
			if (ends > seq.length())                                // if ends > length of sequence...
				snippet = seq.substr(starts);                       //    snippet = substring [starts =>]
			else                                                    // else...
				snippet = seq.substr(starts, ribo.snippet_size);    //    snippet = substring [starts => starts + 700]
			
			if (rev)
				seed = seqLength - (curseed + starts) - 1;
			else
				seed = curseed + starts;

			int startSwitch = 0, endSwitch = 0;   // local bounds of found riboswitch returned by alignment
			int score = 0;                        // score returned by alignment
	
			// Search for a riboswitch in snippet.
			
			list<AlEl> motif_locations;
			found = max_alignment(snippet, &startSwitch, &endSwitch, &score, ribo, motif_locations, rev, seed, testfile);
			
//			if (found != "") cout << found << endl;
			
			if ((found != "") and (endSwitch - startSwitch + 1) <= (int)(ribo.maxLength))  // if we found something
			{
				//Calculate absolute start and ending points.
				
				int absStart, absEnd;
				if (rev)
				{
					absEnd = seed - startSwitch;
					absStart = seed - endSwitch;
				}
				else
				{	
					absStart = seed + startSwitch;
					absEnd = seed + endSwitch;
				}

				// Check for overlap with ORFs.
				
				bool done = false;
				if (!rev)
				{
					for (int i = 0; i < numOrfs; i++)
						if (((absStart >= orfS[i]) && (absStart <= orfF[i])) || ((absEnd >= orfS[i]) && (absStart <= orfF[i])))
							  done = true;
				}
				else
				{
					for (int i = numOrfs; i < numOrfs + numROrfs; i++)
						if (((absStart >= orfS[i]) && (absStart <= orfF[i])) || ((absEnd >= orfS[i]) && (absStart <= orfF[i])))
							  done = true;
				}

				if (!done)
				{ 
					if (rev)
					{
						riboswitch = complement(fullSeq.substr(absStart, absEnd - absStart));
						riboswitch = stringrev(riboswitch);
					}
					else
						riboswitch = fullSeq.substr(absStart, absEnd - absStart);
					
					// fold potential riboswitch
					
					strncpy(seqfoldername, "/tmp/drd_XXXXXX", 255);  // create unique folder for mfold output files
					mktemp(seqfoldername);
					strncpy(cmd, "mkdir ", 255);
					strncat(cmd, seqfoldername, 255);
					//cmd = string("mkdir ") + seqfoldername;
					system(cmd);
					strncat(seqfoldername, "/", 255);
					
					snprintf(riboFilename, 255, "%d_XXXXXX", absStart + 1);
					mktemp(riboFilename);
					strncpy(seqfname, riboFilename, 255);
					strncat(seqfname, "_sequence.fasta", 255);
					strncpy(seqfname_short, riboFilename, 255);
					strncat(seqfname_short, "_sequence", 255);
					char seqfname_long[256];
					strncpy(seqfname_long, seqfoldername, 255);
					strncat(seqfname_long, seqfname, 255);
					ofstream seqfile(seqfname_long, ios::trunc);
					seqfile << riboswitch;
					seqfile.close();
					
					// run mfold
					//cmd = "export PATH=$PATH:/usr/local/bin; cd " + seqfoldername + "; mfold SEQ=" + seqfname + " NA=DNA MODE=bases > /dev/null 2> /dev/null";
					strncpy(cmd, "cd ", 255);
					strncat(cmd, seqfoldername, 255);
					strncat(cmd, "; mfold SEQ=", 255);
					strncat(cmd, seqfname, 255);
					strncat(cmd, " NA=DNA MODE=bases > /dev/null 2> /dev/null", 255);
					system(cmd);
					
					// convert first result to vienna format
					//viennafname = seqfoldername + seqfname_short + "_1.b";
					strncpy(viennafname, seqfoldername, 255);
					strncat(viennafname, seqfname_short, 255);
					strncat(viennafname, "_1.b", 255);
					//cmd = "./Ct2B.pl " + seqfoldername + seqfname_short + "_1.ct > " + viennafname;
					strncpy(cmd, "./Ct2B.pl ", 255);
					strncat(cmd, seqfoldername, 255);
					strncat(cmd, seqfname_short, 255);
					strncat(cmd, "_1.ct > ", 255);
					strncat(cmd, viennafname, 255);
					system(cmd);
					
					// read vienna string
					ifstream viennafile(viennafname);
					vienna_string = "";
					viennafile >> vienna_string;  // read sequence off the top
					viennafile >> vienna_string;
					viennafile.close();
	
					// Align Vienna strings.
					
					vienna_alignment(vienna_string, ribo.vienna, ribo.minViennaScore, alignment1, viennascore1, viennaid1, viennalen1);
		
					// Replace sections of the Vienna string with explicit motifs.
					
					vienna_string_motifs = vienna_string;
					list<AlEl>::iterator a = motif_locations.begin();
					while (a != motif_locations.end())
					{                                                              
						unsigned int pos = a->position - startSwitch;			   
						motif = a->str.substr(0, a->str.find('\n', 0));
						size_t j = 0;
						for (size_t i = 0; i < motif.length(); i++)
							if (motif[i] != '-')
							{
								motif_clean[j] = motif[i];
								j++;
							}
						motif_clean[j] = '\0';
						vienna_string_motifs.replace(pos, j, motif_clean);
						a++;
					}
		
					// Align Vienna strings with motifs inserted.
					
					vienna_alignment(vienna_string_motifs, ribo.motifs_vienna, ribo.minViennaScore, alignment, viennascore, viennaid, viennalen);
									
					// Check if Vienna identity score surpasses threshold.				
					if (viennaid < ribo.minViennaScore)
						done = true;
					
					if ((testfile != NULL) && (testfile->is_open()))
					{
						(*testfile) << score << "\t" << viennaid << "\t" << found.find('\n') << "\t";
					}

					if (!done)
					{
						// Check if this is a duplicate.
						bool dup = false;
						bool replace = false;
						pthread_mutex_lock(&numLock);
						for (int i = 0; i < numRibo; i++)
						{
							if (((absStart >= ribos[i].start) && (absStart <= ribos[i].end)) || ((absEnd >= ribos[i].start) && (absStart <= ribos[i].end)))
							{
								if (score > ribos[i].score)
								{
									// replace
									replace = true;
									
									char replaceRiboFilename[256];
									snprintf(replaceRiboFilename, 255, "%d-%d-%d", ribos[i].score, ribos[i].viennaid, ribos[i].start);
									
									char number[2] = "0";
									for (int count = 1; count < 10; count++)
									{
										number[0] = (char)(count + '0');
										strncpy(cmd, "[ -e ", 255);
										strncat(cmd, imgfoldername, 255);
										strncat(cmd, replaceRiboFilename, 255);
										strncat(cmd, "_", 255);
										strncat(cmd, number, 255);
										strncat(cmd, ".png ]", 255);
										if (system(cmd) != 0)
											break;
										strncpy(cmd, "rm ", 255);
										strncat(cmd, imgfoldername, 255);
										strncat(cmd, replaceRiboFilename, 255);
										strncat(cmd, "_", 255);
										strncat(cmd, number, 255);
										strncat(cmd, ".png", 255);
										system(cmd);
									}
									
									strncpy(cmd, "rm ", 255);
									strncat(cmd, path.c_str(), 255);
									strncat(cmd, "/output/", 255);
									strncat(cmd, replaceRiboFilename, 255);
									strncat(cmd, ".html", 255);
									system(cmd);
									
									pthread_mutex_lock(&riboLock);
									ribos[i].start = absStart;   // add riboswitch to the list
									ribos[i].end = absEnd;
									ribos[i].score = score;
									ribos[i].viennaid = viennaid;
									ribos[i].reverse = rev;
									pthread_mutex_unlock(&riboLock);
								}
								else
								{
									// ignore
									dup = true;
								}
								break;
							}
						}
						
						if (dup)
						{
							pthread_mutex_unlock(&numLock);
						}
						else
						{
							found_something = 1;
							if (!replace)
							{
								pthread_mutex_lock(&riboLock);
								ribos[numRibo].start = absStart;   // add riboswitch to the list
								ribos[numRibo].end = absEnd;
								ribos[numRibo].score = score;
								ribos[numRibo].viennaid = viennaid;
								ribos[numRibo].reverse = rev;
								pthread_mutex_unlock(&riboLock);
								numRibo = numRibo + 1;
							}
							pthread_mutex_unlock(&numLock);
							
							// Align shape strings.
					
							shape_string = Viennatoshape(vienna_string, "short");   // make shape string (Candidate)
							lit_shape = Viennatoshape(ribo.vienna, "short");        // make shape string (Literature)
							shape_alignment(shape_string, lit_shape, 0, shape_align, shape_score, shape_id, shape_len);
						
							// Update query update file.
						
							//if (!replace)
							//{
							//	queryUpdateFile.open(queryupdatefname, ios::trunc);
							//	if (numRibo == 1)
							//		queryUpdateFile << numRibo << " potential riboswitch has been found so far. <br /><br />" << endl;
	//							else
	//								queryUpdateFile << numRibo << " potential riboswitches have been found so far. <br /><br />" << endl;
	//							queryUpdateFile.close();
	//						}
			
							snprintf(riboFilename, 255, "%d-%d-%d", score, viennaid, absStart + 1);

							// Move mfold image files for riboswitch to image folder.
							
							char number[2] = "0";
							for (int count = 1; count < 10; count++)
							{
								//cmd = "[ -e " + seqfoldername + seqfname_short + "_" + (char)(count + '0') + ".png ]";
								number[0] = (char)(count + '0');
								strncpy(cmd, "[ -e ", 255);
								strncat(cmd, seqfoldername, 255);
								strncat(cmd, seqfname_short, 255);
								strncat(cmd, "_", 255);
								strncat(cmd, number, 255);
								strncat(cmd, ".png ]", 255);
								if (system(cmd) != 0)
									break;
								//cmd = "mv " + seqfoldername + seqfname_short + "_" + (char)(count + '0') + ".png " + imgfoldername + riboFilename + "_" + (char)(count + '0') + ".png";
								strncpy(cmd, "mv ", 255);
								strncat(cmd, seqfoldername, 255);
								strncat(cmd, seqfname_short, 255);
								strncat(cmd, "_", 255);
								strncat(cmd, number, 255);
								strncat(cmd, ".png ", 255);
								strncat(cmd, imgfoldername, 255);
								strncat(cmd, riboFilename, 255);
								strncat(cmd, "_", 255);
								strncat(cmd, number, 255);
								strncat(cmd, ".png", 255);
								system(cmd);
							}
							
							// Write info about this riboswitch to a result file.
							
							//resultfname = prefix + foldername + "/output/" + riboFilename + ".html";
							strncpy(resultfname, path.c_str(), 255);
							strncat(resultfname, "/output/", 255);
							strncat(resultfname, riboFilename, 255);
							strncat(resultfname, ".html", 255);
							resultFile.open(resultfname, ios::app);
							resultFile << "position ";
							if (rev) 
								resultFile << absEnd + 1 << " &ndash; " << absStart + 1;
							else
								resultFile << absStart + 1 << " &ndash; " << absEnd + 1;
							resultFile <<" (" << score << " / " << maxscore << " motif identities)<br />" <<  endl;
							resultFile << "<br />" << endl;
							resultFile << "<div style='margin-left:50px'>" << endl;
							resultFile << "<pre>" << endl;
							int t = found.find('\n');   // find end of riboswitch sequence in alignment
							if (rev)
								resultFile << "3'- <a href=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&DATABASE=nr&LAYOUT=TwoWindows&AUTO_FORMAT=Semiauto&PAGE=Nucleotides&HITLIST_SIZE=100&QUERY=" << riboswitch << "\" target=\"_blank\">" << found.substr(0, t) << "</a> -5'" << endl;
							else
								resultFile << "5'- <a href=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&DATABASE=nr&LAYOUT=TwoWindows&AUTO_FORMAT=Semiauto&PAGE=Nucleotides&HITLIST_SIZE=100&QUERY=" << riboswitch << "\" target=\"_blank\">" << found.substr(0, t) << "</a> -3'" << endl;
							resultFile << "    " << found.substr(t + 1) << endl;
														
							resultFile << "</pre>" << endl;
							resultFile << "<br />" << endl;
							resultFile << "Vienna global alignment (" << viennaid1 << " / " << ribo.vienna.length() << " identities)" << endl;    //   Vienna align
							resultFile << alignment1 << "<br /><br />" << endl;
							resultFile << "Vienna global alignment with motifs (" << viennaid << " / " << ribo.motifs_vienna.length() << " identities)" << endl;    //   Vienna align
							resultFile << alignment << "<br /><br />" << endl;
							resultFile << "Shape global alignment (" << shape_id << " / " << lit_shape.length() << " identities)" << endl;    //   Shape align
							resultFile << shape_align << "<br /><br />" << endl;
							resultFile << "</div>" << endl;
							resultFile.close();		
						}
					}
					//cmd = string("rm -rf ") + seqfoldername;
					strncpy(cmd, "rm -rf ", 255);
					strncat(cmd, seqfoldername, 255);
					system(cmd);
				}
			}
			starts = ends - ribo.overlap_size;   
		}
	}
	return (void*)found_something;
}
