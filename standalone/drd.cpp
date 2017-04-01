// drd.cpp
// Denison Riboswitch Detector: main program
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
extern int numRibo;  // defined in finder.cpp


string getSafeName(string str)
{

/*Description: Generate a "terminal safe" directory name based on string str.
  Replaces all special characters in str with '_'.
*/

	string newstr = "";
	string::iterator str_it;

	str_it = str.begin();
	while (str_it != str.end())
	{
		if ( (*str_it == '/') or
		(*str_it == '|') or
		(*str_it == ' ') or
		(*str_it == '\'') or
		(*str_it == '>') or
		(*str_it == '<') or
		(*str_it == '[') or
		(*str_it == ']') or
		(*str_it == '*') or
		(*str_it == '$') or
		(*str_it == '?') or
		(*str_it == '@') or
		(*str_it == '#') or
		(*str_it == '{') or
		(*str_it == '}') or
		(*str_it == '~') or
		(*str_it == '!') or
		(*str_it == '^') )
		{
			newstr = newstr + "_";
		}

		else
		{
			newstr = newstr + *str_it;
		}
		
		++str_it; 
	}

	return newstr;	
}




bool inList(list<string> lst, string str)
{

/*Description: A function that returns true if the string str is in a the list of strings lst, and false
  Otherwise. 
*/

	list<string>::iterator it;
	it = lst.begin();

	while (it != lst.end())
	{
		if (*it == str)
		{
			return true;
		}
		
		++it;
	}
	
	return false;
}


string ridPrefix(string str)
{

/*Description: A function used to strip the "--" or '-' prefixes of the command line arguments. 
*/

  if (str.substr(0,2) == "--")
    return str.substr(2);

  else 
    return str.substr(1);
}

string formatTime(string time)
{
  //Rearrange string time so that it produces the current time in format: JUN-00-2015_HH:MM:SS for uniqueID 
  string newTime = "", temp;

  newTime = newTime + time.substr(4,3) + "-";

  temp = time.substr(8,2);
  if (temp.substr(0,1) == " ")
    {
      temp = "0" + temp.substr(1,1);
    }

  newTime = newTime + temp + "-";
  newTime = newTime + time.substr(20,4) + "_";
  newTime = newTime + time.substr(11,8);
  
  return newTime;
}

string get_rimg(string riboname)
{

/*Description: A function used to get the location of the Rfam rimg that corresponds to the 
  riboswitch type riboname. If the riboswitch type is user-defined, the rimg returned is an
  empty string.
*/

	string rimg;

	//map that converts riboname to an int
	map<string, int> riboToInt;                         
		    riboToInt["Cobalamin"] = 0;
		    riboToInt["FMN"] = 1;
		    riboToInt["glmS"] = 2;
		    riboToInt["Glycine"] = 3;
		    riboToInt["Lysine"] = 4;
		    riboToInt["PreQ1"] = 5;
		    riboToInt["Purine"] = 6;
		    riboToInt["SAM-II"] = 7;
		    riboToInt["SAM-I"] = 8;	
		    riboToInt["SAM-IV"] = 9;
		    riboToInt["TPP"] = 10;
		    riboToInt["ykkCyxkD"] = 11;
		    riboToInt["yybPykoY"] = 12;
		
		    //switch table to get corresponding rimg 
		    switch(riboToInt[riboname])              
		      {

		      case 0:
			rimg = "'http://rfam.xfam.org/family/RF00174/image/cons'";
			break;

		      case 1:
			rimg = "'http://rfam.xfam.org/family/RF00050/image/cons'";
			break;

		      case 2:
			rimg = "'http://rfam.xfam.org/family/RF00234/image/cons'";
			break;

		      case 3:
			rimg = "'http://rfam.xfam.org/family/RF00504/image/cons'";
			break;

		      case 4:
			rimg = "'http://rfam.xfam.org/family/RF00168/image/cons'";
			break;

		      case 5:
			rimg = "'http://rfam.xfam.org/family/RF00522/image/cons'";
			break;

		      case 6:
			rimg = "'http://rfam.xfam.org/family/RF00167/image/cons'";
			break;

		      case 7:
			rimg = "'http://rfam.xfam.org/family/RF00521/image/cons'";
			break;

		      case 8:
			rimg = "'http://rfam.xfam.org/family/RF00162/image/cons'";
			break;

		      case 9:
			rimg = "'http://rfam.xfam.org/family/RF00634/image/cons'";
			break;

		      case 10:
			rimg = "'http://rfam.xfam.org/family/RF00059/image/cons'";
			break;

		      case 11:
			rimg = "'http://rfam.xfam.org/family/RF00442/image/cons'";
			break;

		      case 12:
			rimg = "'http://rfam.xfam.org/family/RF00080/image/cons'";
			break;

		      default:
			rimg = "";

		      }
	return rimg;
}



list<string> getFiles(string path)
{
/*Description: A function used to create a list of strings representing all of the files within the directory 
  "path."
*/
	__dirstream *dir;
	dirent *file;
	string file_name;
	list<string> files;

	dir = opendir(path.c_str());
	if (dir == NULL)
	{
		cout << "Could not open directory: " << path << endl;
	}
	
	else
	{
		file = readdir(dir);
		while (file != NULL)
		{
			file_name = file->d_name;
			if (file_name.substr(0,1) != ".")
			{
				files.push_back(file_name);
			}

			file = readdir(dir);
		}	
	}

	closedir(dir);
	return files;
}


// Read input sequence from file.
string readSeqFile(ifstream& infile)
{
	char cap;           // '>' character at beginning of input file
	string gen_name;    // sequence name from input file
	string curSeq;
	string inSeq = "";
	
	cap = infile.peek();
	if (cap == '>')
		getline(infile, gen_name);
		
	while (infile.good())
	{
		getline(infile, curSeq);
		if (curSeq[0] == '>')      // only read first sequence
			break;
		inSeq += curSeq.substr(0, curSeq.length());
	}
	infile.close();
	
	return standardize(inSeq);  // make all caps, U -> T, remove non-alpha;
}

int main(int argc, char **argv)
{
	
	string mainArg;

	//If minimum number of arguments (6) is not reached upon execution of the drd 
	if (argc < 6)                                  
	{
		cout << endl;

		//Check for the --help command line opiton
		if (argc == 2)                        
		{
			mainArg = argv[1];
			mainArg = ridPrefix(mainArg);
		
			if (mainArg == "help")
			{
				//--help, output Usage and command line options, and exit.

				cout << "Usage: drd [--orfs orf_len] --ouput output_directory --ribo ribo_type1 [ribo_type2 ... ribo_typen] sequence_filename" << endl << endl;

				cout << "Command line options:" << endl;
				cout << "--orfs orf_len" << endl;
	    			cout << "    indicates that only riboswitches that are not within an open reading" << endl;
	    			cout << "    frame with minimum length orf_len (measured in codons) should be" << endl;
	    			cout << "    returned." << endl << endl; 

				cout << "--output output_directory" << endl;
	    			cout << "    indicates that output_directory is the directory to which results" << endl;
	    			cout << "    will be written.  If the directory does not already exist, it will" << endl; 
	    			cout << "    be created." << endl << endl;

				cout << "--ribo ribo_type1 [ribo_type2 ... ribo_typen]" << endl;
	    			cout << "    indicates the riboswitch types to search for.  The types are" << endl; 
	    			cout << "    specified by the names of riboswitch definition files that" << endl; 
	    			cout << "    are located in the defFiles folder.  (Do not include the .def" << endl; 
	    			cout << "    extension.)" << endl << endl;

				cout << "sequence_filename" << endl;
	    			cout << "    the name of the FASTA file containing the input sequence(s)." << endl; 
	    			cout << "    The sequence file may contain multiple sequences." << endl << endl;
	    
				cout << "Example: drd --orfs 120 --output DRD_RESULTS --ribo TPP SAM-I seq_file.fasta" << endl << endl;

				cout << "For more information on the DRD algorithm and instructions" << endl; 
				cout << "on how to create your own riboswitch definition file, go to" << endl; 
				cout << "drd.denison.edu and select the Help option." << endl << endl;

				exit(0);		
			}

		}

		//Output corresponding error statement and exit
		if (argc < 5)                                      
			cerr << "Missing arguments in";
	  
	  	else
	   		cerr << "Missing an argument in";

		for (int i=0; i < argc; i++)
		{
			cerr << ' ' << argv[i]; 
		}
	 	
		cerr << endl << "Run " << argv[0] << " --help to see a list of command line options." << endl << endl;
	  	exit(EXIT_FAILURE);
	}

	//initialize variables to read in command line arguments
	string sequence_file, riboswitch_def_file, results_folder, arg;                   
	bool doOrfs, ribos_com, output_com, ribos_found, debug;
	int minORFLength;
	list<string> rbs_def_files, ribo_types;
	list<string>::iterator it;
	
	//set debug to true if attempting to fix an error	
	debug = false;

	//bools to check if the orfs, ribo, and output commands were used upon execution.             
	doOrfs = false;		//optional    
	ribos_com = false;	//required
	output_com = false;	//required

	//get list of riboswitch definition files in defFiles
	rbs_def_files = getFiles("defFiles");
	

	//Read in command line arguments

	mainArg = argv[1];
	if ( mainArg.substr(0,1) != "-" )
	  {
	    cerr << "Command " << argv[1] << " can not be recognized." << endl;
	    cerr << "Run " << argv[0] << " --help to see Usage and a whole list of command line options." << endl;
	    exit(EXIT_FAILURE);
	  }
	
	mainArg = ridPrefix(mainArg);
	
	for(int i=2; i< (argc-1); i++)
	  {
	    	
	    if ( mainArg.substr(0,1) == "-" )
	      {
		cerr << "Command " << argv[i-1] << " can not be recognized." << endl;
		cerr << "Run " << argv[0] << " --help to see Usage and a whole list of command line options." << endl;
		exit(EXIT_FAILURE);
	      }

	    arg = argv[i];	      
	    
	    if ( arg.substr(0,1) == "-" )
	      {
		mainArg = ridPrefix(arg);
		continue;
	      }
		
	    
	    else if ( (mainArg == "ribo") or (mainArg == "r") )
	      {
		ribos_com = true;
		if ( inList(rbs_def_files, arg + ".def") )
		  {
		    ribo_types.push_back(arg);
		  }

		else
		  {
		    cerr << "Riboswitch Definition File " << arg << " can not be found." << endl;
		    exit(EXIT_FAILURE);
		  }
	      }

	    else if (mainArg == "orfs")
	      {
		doOrfs = true;
		minORFLength = atoi(argv[i]);
	      }

	    else if (mainArg == "output")
	      {
		output_com = true;
		results_folder = arg;
	      }
	    
	    else
	      {
		cerr << "Command " << argv[i-1] << " can not be recognized." << endl;
		cerr << "Run " << argv[0] << " --help to see Usage and a whole list of command line options." << endl;
		exit(EXIT_FAILURE);
	      }
	  }
	
	if ( (!ribos_com) or (!output_com) )
	{
		cerr << endl;		
		cerr << "Missing an argument in";

		for (int i=0; i < argc; i++)
		{
			cerr << ' ' << argv[i]; 
		}
	 	
		cerr << "Run " << argv[0] << " --help to see Usage and a whole list of command line options." << endl;
	  	exit(EXIT_FAILURE);
	}

	sequence_file = argv[(argc-1)];


	//intialize variables to develop directory configuration and Path
	string uniqueID, seq_name, acc_num, newFile, line, riboname, path, resultsFile, tercmd;  
	ifstream infile, riboswitch, sequences;
	ofstream outfile, results, overall, overall_trunc;
	int char_found, dir_found;
	
	
	//uniqueID = name of sequence file without the .fasta extension
	uniqueID = sequence_file;
	char_found = uniqueID.find_last_of('.');

	if (char_found != -1)
	{
		uniqueID = uniqueID.erase(char_found, -1); 
	}
	uniqueID = getSafeName(uniqueID);
	

	//build path variable, results_folder/uniqueID, and check if it already exists
	path = results_folder + '/' + uniqueID;
	tercmd = "[ -d " + path + " ]";
	dir_found = system(tercmd.c_str());


	//If results_folder/uniqueID directory exists, prompt user if they would like to overwrite the existing directory
	if (dir_found == 0)
	{
		cout << "The Directory " + path + " already exists, would you like to overwrite it? (Y/N) ";
		cin >> tercmd;
		
		if ( (tercmd == "Y") or (tercmd == "y") )
		{
			tercmd = "rm -rf " + path;
			system(tercmd.c_str());	
		}

		else if ( (tercmd == "N") or (tercmd == "n") )
		{
			cout << "The directory " + path + " was not overwritten." << endl;
			exit(0);
		}

		else 
		{
			cerr << "The command \"" + tercmd + "\" was not recognized." << endl;
			exit(EXIT_FAILURE);
		}
	}
	
	
	//Create the results_folder/uniqueID directory 
	tercmd = "mkdir -p " + path;
	system(tercmd.c_str());

	//Create the index.html file (Overall results page)
	newFile = path + "/index.html";
	overall.open(newFile.c_str());

	//Create the overall_trunc.html file (Truncated overall results page)
	newFile = path + "/overall_trunc.html";
	overall_trunc.open(newFile.c_str());

	
	// Create Overall results web page
	overall << "<html><head><title>DRD Results: Overall</title>" << endl;
	overall << "<link rel=\"stylesheet\" type=\"text/css\" href=\"http://drd.denison.edu/drd.css\">" << endl;
	overall << "</head>" << endl; 
	overall << "<body>" << endl;

	// Create Truncated Overall results web page 
	overall_trunc << "<html><head><title>DRD Results: Overall</title>" << endl;
	overall_trunc << "<link rel=\"stylesheet\" type=\"text/css\" href=\"http://drd.denison.edu/drd.css\">" << endl;
	overall_trunc << "</head>" << endl; 
	overall_trunc << "<body>" << endl;
	
	// Overall results page header 
	overall <<  "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">\n\n<!-- Header Table Begin -->\n<table width=\"100%\" border=\"2\">\n<tr>\n<td>\n\n<table width=\"100%\" border=\"0\" class=\"header\">\n<tr border=\"0\">\n<td class=\"drd\" border=\"0\" width=\"80\">\n<a href=\"http://drd.denison.edu\"><img alt=\"DRD-logo\" src= \"http://drd.denison.edu/images/DRD-logo.png\" align=\"left\"></a>\n</td>\n<td>\n<h1>&nbsp; <a href=\"http://drd.denison.edu\">Denison Riboswitch Detector</a></h1>\n<h5>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  A tool for finding riboswitches within a DNA sequence</h5>\n</td>\n<td class=\"drd\" align=\"right\" valign=\"bottom\">\n<table width=\"50\" border=\"0\" class=\"drd\">\n<tr border=\"0\">\n<td class=\"drd\" valign=\"top\" width=\"25\">\n<a href=\"http://drd.denison.edu/help.php\">\n<input type=\"submit\" value=\"Help\">\n</a>\n </td>\n<td class=\"drd\" valign=\"top\" width=\"25\">\n<form method='GET' action='mailto:havill@denison.edu'>\n<input type=\"submit\" value=\"Contact\">\n</form>\n</td>\n</tr>\n</table>\n\n</td>\n</tr>\n</table>\n\n</td>\n</tr>\n</table>\n<!-- Header Table End -->\n\n<!-- Message Table Begin -->\n<table width=\"100%\" border=2 cellpadding=10>\n<tr>\n<td>\n" << endl;

	// Truncated overall results page header
	overall_trunc <<  "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">\n\n<!-- Header Table Begin -->\n<table width=\"100%\" border=\"2\">\n<tr>\n<td>\n\n<table width=\"100%\" border=\"0\" class=\"header\">\n<tr border=\"0\">\n<td class=\"drd\" border=\"0\" width=\"80\">\n<a href=\"http://drd.denison.edu\"><img alt=\"DRD-logo\" src= \"http://drd.denison.edu/images/DRD-logo.png\" align=\"left\"></a>\n</td>\n<td>\n<h1>&nbsp; <a href=\"http://drd.denison.edu\">Denison Riboswitch Detector</a></h1>\n<h5>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  A tool for finding riboswitches within a DNA sequence</h5>\n</td>\n<td class=\"drd\" align=\"right\" valign=\"bottom\">\n<table width=\"50\" border=\"0\" class=\"drd\">\n<tr border=\"0\">\n<td class=\"drd\" valign=\"top\" width=\"25\">\n<a href=\"http://drd.denison.edu/help.php\">\n<input type=\"submit\" value=\"Help\">\n</a>\n </td>\n<td class=\"drd\" valign=\"top\" width=\"25\">\n<form method='GET' action='mailto:havill@denison.edu'>\n<input type=\"submit\" value=\"Contact\">\n</form>\n</td>\n</tr>\n</table>\n\n</td>\n</tr>\n</table>\n\n</td>\n</tr>\n</table>\n<!-- Header Table End -->\n\n<!-- Message Table Begin -->\n<table width=\"100%\" border=2 cellpadding=10>\n<tr>\n<td>\n" << endl;


	//Overall/Truncated results page titles 	
	overall << "<h2>Overall Results - "+ uniqueID + "</h2>" << endl;
	overall_trunc << "<h2>Overall Results (Truncated) - "+ uniqueID + "</h2>" << endl;

	//Link to overall/truncated results pages
	overall << "&nbsp;&nbsp;<tt><font color=\'blue\'><a href=\'overall_trunc.html\'>For truncated results, click here.</a></font></tt>";
	overall_trunc << "&nbsp;&nbsp;<tt><font color=\'blue\'><a href=\'index.html\'>For elongated results, click here.</a></font></tt>";  

	//Check that the sequence file can be found, if so, open it.
	sequences.open(sequence_file.c_str());
	if (!sequences.good())
	  {
	    sequences.close();
	    cerr << "Input sequence file could not be read." << endl;
	    exit(EXIT_FAILURE);
	  }
	
	getline(sequences, acc_num, '>');
	while(!sequences.eof())                  //============================================ Multiple Sequences Loop =============================================================
	{
		
		//Read in the sequence accession number
		getline(sequences, acc_num);

		//Add sequence to Overall/Truncated results pages		
		overall << "<h3>" + acc_num + "</h3>" << endl;
		overall << "<table>" << endl;
		overall_trunc << "<h3>" + acc_num + "</h3>" << endl;
		overall_trunc << "<table>" << endl;

		//Generate terminal-friendly directory name
		seq_name = getSafeName(acc_num.substr(0, 20));
		
		if (debug)       
			cout << "======================== " + seq_name + " ===========================" << endl;
		
		
		//Add seq_name to path and create directory for the new sequence
		path = results_folder + '/' + uniqueID + '/' + seq_name;

		tercmd = "mkdir -p " + path;
		system(tercmd.c_str());
	
		//Create a copy of the sequence file to be placed in the seq_name folder 
		newFile = path + '/' + seq_name + ".fasta";                                 
		outfile.open(newFile.c_str());
		getline(sequences, line, '>');	
		outfile << acc_num << endl;
		outfile << line;
		outfile.close();

		//Read in sequence from new sequence file
		string inSeq = "";                         // input sequence
		infile.open(newFile.c_str());  
		inSeq = readSeqFile(infile);
		infile.close();
		
		//Reset to no riboswitches found and start of riboswitch types list
		ribos_found = false;
		it = ribo_types.begin();


		while(it != ribo_types.end())    //========================================== Multiple Riboswitches Loop ============================================================
		{
			//Read in the ribo_type, and reset number of ribo_type riboswitches found to 0.
			numRibo = 0;
			riboname = *it;
			++it;
			
			//Build name of the riboswitch definition file
			riboswitch_def_file = "defFiles/" + riboname + ".def"; 
		
			if (debug)
			{
				cout << endl << "START RIBOSWITCH: " << riboname << endl; 
				cout << "   BUILD RIBOSWITCH" << endl; 
			}
			
			//Create new directory for riboswitch type being searched for in input sequence
			//Also add riboswitch to path variable
			path = results_folder + '/' + uniqueID + '/' + seq_name + '/' + riboname;
			tercmd = "mkdir -p " + path;
			system(tercmd.c_str());
	
			if (debug)
				cout << "   BUILD IMAGES" << endl; 
			
			//Create images directory to store mfold images
			tercmd = "mkdir -p " + path + "/images";                                    
			system(tercmd.c_str());

			//Create temporary output directory for mfold output
			tercmd = "mkdir -p " + path + "/output";                                    
			system(tercmd.c_str());

			if (debug)
			cout << "   BUILD RESULTS" << endl;
			

			//Create results file for above riboswitch type
			newFile = path + '/' + "results.html";                                     
			resultsFile = newFile;
			results.open(newFile.c_str());
			
			// reverse complement sequence
			string rev_comp_seq;  
	
			int totalSnippets;
			int numOrfs, numROrfs;
			int *orfS;
			int *orfF;
			
			if (debug)
				cout << "   READ RIBODEF " << endl; 

			//Check that the riboswitch definition file can be found, if so, open it.
			infile.open(riboswitch_def_file.c_str(), ifstream::binary);
			if (!infile.good())
			{
				cerr << "Riboswitch definition file was not found." << endl;
				results << endl << "Riboswitch definition file was not found." << endl;
				results.close();
				return EXIT_FAILURE;
			}
			RiboDef *ribo = new RiboDef(infile);
			
	     
			// Calculate the total number of snippets
			if ( (inSeq.length() <= ribo->overlap_size) or (inSeq.length() <= ribo->snippet_size) )
			  totalSnippets = 1;

			else if ((inSeq.length() - ribo->overlap_size) % (ribo->snippet_size - ribo->overlap_size) == 0)
			  {
			    totalSnippets = (inSeq.length() - ribo->overlap_size) / (ribo->snippet_size - ribo->overlap_size);		
			  }
			
			else
				totalSnippets = ((inSeq.length() - ribo->overlap_size) / (ribo->snippet_size - ribo->overlap_size)) + 1;
	
			rev_comp_seq = stringrev(complement(inSeq));       

	
			//Create results page for the individual riboswitch type 
			results << "<html><head><title>DRD Results: " + riboname + "</title>" << endl;
			results << "<link rel=\"stylesheet\" type=\"text/css\" href=\"http://drd.denison.edu/drd.css\">" << endl;
			results << "</head>" << endl; 
			results << "<body>" << endl;

			//Results page header
			results <<  "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">\n\n<!-- Header Table Begin -->\n<table width=\"100%\" border=\"2\">\n<tr>\n<td>\n\n<table width=\"100%\" border=\"0\" class=\"header\">\n<tr border=\"0\">\n<td class=\"drd\" border=\"0\" width=\"80\">\n<a href=\"http://drd.denison.edu\"><img alt=\"DRD-logo\" src= \"http://drd.denison.edu/images/DRD-logo.png\" align=\"left\"></a>\n</td>\n<td>\n<h1>&nbsp; <a href=\"http://drd.denison.edu\">Denison Riboswitch Detector</a></h1>\n<h5>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  A tool for finding riboswitches within a DNA sequence</h5>\n</td>\n<td class=\"drd\" align=\"right\" valign=\"bottom\">\n<table width=\"50\" border=\"0\" class=\"drd\">\n<tr border=\"0\">\n<td class=\"drd\" valign=\"top\" width=\"25\">\n<a href=\"http://drd.denison.edu/help.php\">\n<input type=\"submit\" value=\"Help\">\n</a>\n </td>\n<td class=\"drd\" valign=\"top\" width=\"25\">\n<form method='GET' action='mailto:havill@denison.edu'>\n<input type=\"submit\" value=\"Contact\">\n</form>\n</td>\n</tr>\n</table>\n\n</td>\n</tr>\n</table>\n\n</td>\n</tr>\n</table>\n<!-- Header Table End -->\n\n<!-- Message Table Begin -->\n<table width=\"100%\" border=2 cellpadding=10>\n<tr>\n<td>\n" << endl;


			//Output query information to results.html for inclusion on results page. 
			results << "<h3>Query Parameters</h3>" << endl;
			results << "<table>" << endl;
			results << "<tr><td>&nbsp;</td><td><b>Input file</b> </td><td><tt>" << sequence_file << "</tt> (" << inSeq.length() << " bases, " << totalSnippets << " segments, segment/overlap length = " << ribo->snippet_size << "/" << ribo->overlap_size << ")</tt></td></tr>" << endl;
			results  << "<tr><td>&nbsp;</td><td><b>Riboswitch target&nbsp;&nbsp;</b></td><td>" << riboname << " : <tt>" << *ribo << "</tt></td></tr>" << endl;
			results << "<tr><td>&nbsp;</td><td><b>Threshold scores</b></td><td> motif identities &geq; " << ribo->lowestTotal << ", Vienna identities &geq; " << ribo->minViennaScore << ", maximum length = " << ribo->maxLength << "</td></tr>" << endl;
			if (doOrfs)
				results << "<tr><td>&nbsp;</td><td><b>ORF search</b></td><td> on, threshold = " << minORFLength / 3 << " amino acids</td></tr>" << endl;
			else
				results << "<tr><td>&nbsp;</td><td><b>ORF search</b></td><td> off</td></tr>" << endl;
			results << "</table>" << endl;
			results.close();	
	

			// Find ORFs if option is checked.
			if (doOrfs)
			{
				int maxOrfs = inSeq.length() / minORFLength;
				orfS = new int[maxOrfs];
				orfF = new int[maxOrfs];
				numOrfs = findOrfs(inSeq,orfS,orfF,minORFLength);
				numROrfs = findRevOrfs(rev_comp_seq,orfS,orfF,minORFLength,numOrfs);
			}
			else
			{
				numOrfs = 0;
				numROrfs = 0;
				orfS = NULL;
				orfF = NULL;
			}

			//Determine number of threads based on the number of snippets
			int numThreads;
			if (totalSnippets < NUM_THREADS)
			{
				numThreads = totalSnippets;
			}
			else
			{
				numThreads = NUM_THREADS;
			}			
			int per = totalSnippets / numThreads;
			int lper = per * (ribo->snippet_size - ribo->overlap_size);

			//Create a list of thread objects and build the list of arg_struct objects 
			pthread_t thread[numThreads];
			arg_struct **args = new arg_struct*[numThreads];
			args[0] = new arg_struct;	
			args[0]->forwSeq = inSeq.substr(0, lper + ribo->overlap_size);
			args[0]->backSeq = rev_comp_seq.substr(0, lper + ribo->overlap_size);				
			args[0]->seed = 0;
			args[0]->riboswitch = ribo;
			args[0]->foldername = path;      //holds path to riboswitch directory instead of results_folder
			args[0]->orfS = orfS;
			args[0]->orfF = orfF;
			args[0]->numOrfs = numOrfs;
			args[0]->numROrfs = numROrfs;
			args[0]->fullSeq = inSeq;
			args[0]->testfile = NULL;

			if (debug)
				cout << "   START SEARCH" << endl;
			
			//Continue to build the list of arg_struct objects
			for (int i = 1; i < numThreads; i++)
			{
				args[i] = new arg_struct;
				*args[i] = *args[0];
				if (i == numThreads - 1)
				{
					args[i]->forwSeq = inSeq.substr(i * lper);
					args[i]->backSeq = rev_comp_seq.substr(i * lper);	
				}
				else 
				{		
								
				  args[i]->forwSeq = inSeq.substr(i * lper, lper + ribo->overlap_size);
				  args[i]->backSeq = rev_comp_seq.substr(i * lper, lper + ribo->overlap_size);
				}			
				args[i]->seed = i * lper;
			}
		
			//Generate the threads and search for riboswitches in finder.cpp
			for (int i = 0; i < numThreads; i++)
			{
			  if (pthread_create(&thread[i], NULL, finder, (void *) args[i]) != 0)
			    {
			      results.open(resultsFile.c_str(), ios::app);
			      results << endl << "Thread Create Error"<<i<<": "<< EXIT_FAILURE << endl;
			      results.close();
			      return EXIT_FAILURE;
			    }
			}

			//Join threads
			for (int i = 0; i < numThreads; i++)
			{
				if (pthread_join(thread[i], NULL) != 0)
				{
				  results.open(resultsFile.c_str(), ios::app);
				  results << endl << "Join Thread Error"<<i<<": "<< EXIT_FAILURE <<endl;
				  results.close();
				  return EXIT_FAILURE;
				}
			}
			
			if (debug)
			{
				cout << "   SEARCH COMPLETE" << endl;
				cout << "   WRITE RESULTS START" << endl;
			}
			
			//Reopen results page
			results.open(resultsFile.c_str(), ios::app);			
			results << "<h3>Query Results</h3>" << endl;
			
			//Add riboswitch type to Overall results page
			overall <<"<tr><td>&nbsp;</td><td><b>" + riboname + "&nbsp;&nbsp;</b></td><td>"; 
			

			//No riboswitches found
			if (numRibo == 0)
			{
				//Add "No potential riboswitches were found" to Overall and individual results pages
			 	overall << "<tt><font color=\'red\'>No potential riboswitches were found</font></tt></td></tr>";
				results << endl << "&nbsp;<font color=\'red\'>No potential riboswitches were found</font></b><br /><br /><br />";
				
				//results page footer
				results << "<!-- Footer Table Begin -->\n<table width=\"100%\" border=\"0\" class=\"footer\" cellpadding=\"4px\">\n<tr style=\"vertical-align:bottom\">\n<td style=\"text-align:left\">\nIf you find this software useful, please cite:\n<small>\n<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; J.T. Havill, C. Bhatiya, S.M. Johnson, J.D. Sheets, and J.S. Thompson.  <a href=\"http://bioinformatics.oxfordjournals.org/content/30/21/3012\">A new approach for detecting riboswitches in DNA sequences</a>.  <i>Bioinformatics</i> 30(21):3012-3019, 2014.\n</small>\n</td>\n<td style=\"text-align:right\">\n<p><a href='http://denison.edu/academics/computer-science'>Department of Mathematics and Computer Science</a>\n<br />Denison University</a>\n<br />Granville, Ohio, USA\n</td>\n</tr>\n</table>\n<!-- Footer Table End -->\n" << endl;
				results.close();
			}

			//1 riboswitch found	
			else if (numRibo == 1)
			{
				//Add "1 potential riboswitch was found" to the Overall, Truncated, and individual results pages
				overall << "<tt><font color=\'blue\'><a href=\'" + seq_name + '/' + riboname + "/results.html" + "\'>1 potential riboswitch was found</font></a></tt></td></tr>";
				overall_trunc <<"<tr><td>&nbsp;</td><td><b>" + riboname + "&nbsp;&nbsp;</b></td><td>"; 				
				overall_trunc << "&nbsp;&nbsp;<tt><font color=\'blue\'><a href=\'" + seq_name + '/' + riboname + "/results.html" + "\'>1 potential riboswitch was found</font></a></tt></td></tr>";
				results << endl << "&nbsp;<font color=\'lime\'>1 potential riboswitch was found</font></b><br /><br /><br />";
				ribos_found = true;
			}

			//Multiple riboswitches found	
			else
			{
				//Add "numRibo potential riboswitches were found" to the Overall, Truncated, and individual results pages
				overall << "<tt><font color=\'blue\'><a href=\'" + seq_name + '/' + riboname + "/results.html" + "\'>" << numRibo << " potential riboswitches were found</font></a></tt></td></tr>";
				overall_trunc <<"<tr><td>&nbsp;</td><td><b>" + riboname + "&nbsp;&nbsp;</b></td><td>"; 				
				overall_trunc << "&nbsp;&nbsp;<tt><font color=\'blue\'><a href=\'" + seq_name + '/' + riboname + "/results.html" + "\'>" << numRibo << " potential riboswitches were found</font></a></tt></td></tr>";
				results << endl << "&nbsp;<font color=\'lime\'>" << numRibo << " potential riboswitches were found</font></b><br /><br /><br />";
				ribos_found = true;	
			}

			//Initialize variables for copying over the mfold images and DRD output for inclusion in the results page
			list<string> outputFiles, imageFiles, temp;
			list<string>::iterator output_it;			
			string output_file, output_file_path, image_file, image_file_path;	
			int imageLen, count;

			if (debug)
				cout << "   TRANSFER OUTPUT START" << endl;
			
			
			//Get list of filenames in the output directory
			outputFiles = getFiles(path + "/output");
			

			//Iterate over the output files and take first line of each for results page
			output_it = outputFiles.begin();
			count = 0;
			while (output_it != outputFiles.end())
			{
				output_file = *output_it;
		
				output_file_path = path + "/output/" + output_file;
				infile.open(output_file_path.c_str());
				getline(infile, line);

				results << "<strong><a href=\'#result" << count+1 << "\'>Result " << count+1 << ":&nbsp;&nbsp;</strong>" + line + "</a>" << endl;
				infile.close();
				
				count++;
				++output_it;
			}
			
			char numTochar[256];
			string extension;

			//Iterate over the output files and get full query results for all riboswitches found for inclusion on individual results page
			output_it = outputFiles.begin();
			results << "<h3>Query Result Details</h3>" << endl;
			count = 0;
			while (output_it != outputFiles.end())
			{
				output_file = *output_it;
				++output_it;

				output_file_path = path + "/output/" + output_file;
				infile.open(output_file_path.c_str());

				results << "<a name=\'result" << count+1 << "\'> <strong>Result " << count+1 << ":&nbsp;&nbsp;</strong></a>" << endl;
		
				//Read in output file until EOF
				while(!infile.eof())
				{
						
					getline(infile, line);
					results << line << endl;
				}
				infile.close();
				output_file.erase(output_file.find_last_of('.'));

				//Reference to mfold and Rfam for individual results page images
				results << "<table style=\'margin-left:50px\' border=\'0\'>" << endl;
				results << "<tr><td style=\'text-align:center\'>Candidate sequence folded with <a href=\'http://mfold.rna.albany.edu/?q=mfold\'>mfold</a></td><td style=\'width:30px\'></td><td style=\'text-align:center\'>Consensus image courtesy of <a href=$rimg>Rfam</a></td></tr>" << endl;	
		
				if (debug)
					cout << "      TRANSFER IMAGES START" << endl;
				

				//Include mfold images in the individual results page
				imageFiles = getFiles(path + "/images/");
				imageLen = imageFiles.size();	
				for (int k=0; k < imageLen; k++)
				{	
					snprintf(numTochar, 255, "_%d.png", k+1);
					extension = numTochar;
					image_file = output_file + extension;
					
					if (inList(imageFiles, image_file))
					{
					
						results << "<tr height=\'400\' style=\'vertical-align:top\'><td style=\'text-align:center\'><img src = images/" + image_file + " height = \'400\' style=\'border:2px solid black;\'></td>" << endl;
						results << "<td style=\'width:30px\'></td>" << endl;
						if (riboname != "Userdefined")
						{
							results << "<td style=\'text-align:center\'><img src = " + get_rimg(riboname) + " height = \'400\' style=\'border:2px solid black;v-align:top\'></td></tr>";
						}
						
					}
				}
				
				if (debug)
					cout << "      TRANSFER IMAGES END" << endl;
				
				results << "</table><br />" << endl;
				count++;
			}

			//Individual results page footer
			results << "<p></p>" << endl;
			results << "<!-- Footer Table Begin -->\n<table width=\"100%\" border=\"0\" class=\"footer\" cellpadding=\"4px\">\n<tr style=\"vertical-align:bottom\">\n<td style=\"text-align:left\">\nIf you find this software useful, please cite:\n<small>\n<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; J.T. Havill, C. Bhatiya, S.M. Johnson, J.D. Sheets, and J.S. Thompson.  <a href=\"http://bioinformatics.oxfordjournals.org/content/30/21/3012\">A new approach for detecting riboswitches in DNA sequences</a>.  <i>Bioinformatics</i> 30(21):3012-3019, 2014.\n</small>\n</td>\n<td style=\"text-align:right\">\n<p><a href='http://denison.edu/academics/computer-science'>Department of Mathematics and Computer Science</a>\n<br />Denison University</a>\n<br />Granville, Ohio, USA\n</td>\n</tr>\n</table>\n<!-- Footer Table End -->\n" << endl;

			if (debug)
				cout << "   WRITE RESULTS COMPLETE" << endl;
			
			//Remove temporary output directory 
			tercmd = "rm -rf " + path + "/output";
			system(tercmd.c_str());
			results.close();
			
			//delete orfS, orfF, ribo, and args
			delete [] orfS;
			delete [] orfF;    
			delete ribo;	
			for (int i = 1; i < numThreads; i++)
			{
				delete args[i];
			}
			delete [] args;
			
			if (debug) 
				cout << "   NUMBER RIBOSWITCHES FOUND: " << numRibo << endl;
		}

		//Check if no riboswitches were found for Truncated results page
		if (!ribos_found)
		{
			overall_trunc << "&nbsp;&nbsp;&nbsp;<tt><font color=\'red\'>No potential riboswitches were found</font></tt></td></tr>";
		}

		overall << "</table>" << endl;
		overall_trunc << "</table>" << endl;
		
		ribos_found = false;
	}
	

	//Add overall results page footer and close 
	overall << "<p></p>" << endl;
	overall << "<!-- Footer Table Begin -->\n<table width=\"100%\" border=\"0\" class=\"footer\" cellpadding=\"4px\">\n<tr style=\"vertical-align:bottom\">\n<td style=\"text-align:left\">\nIf you find this software useful, please cite:\n<small>\n<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; J.T. Havill, C. Bhatiya, S.M. Johnson, J.D. Sheets, and J.S. Thompson.  <a href=\"http://bioinformatics.oxfordjournals.org/content/30/21/3012\">A new approach for detecting riboswitches in DNA sequences</a>.  <i>Bioinformatics</i> 30(21):3012-3019, 2014.\n</small>\n</td>\n<td style=\"text-align:right\">\n<p><a href='http://denison.edu/academics/computer-science'>Department of Mathematics and Computer Science</a>\n<br />Denison University</a>\n<br />Granville, Ohio, USA\n</td>\n</tr>\n</table>\n<!-- Footer Table End -->\n" << endl;
	overall.close();

	//Add Truncated results page footer and close
	overall_trunc << "<p></p>" << endl;
	overall_trunc << "<!-- Footer Table Begin -->\n<table width=\"100%\" border=\"0\" class=\"footer\" cellpadding=\"4px\">\n<tr style=\"vertical-align:bottom\">\n<td style=\"text-align:left\">\nIf you find this software useful, please cite:\n<small>\n<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; J.T. Havill, C. Bhatiya, S.M. Johnson, J.D. Sheets, and J.S. Thompson.  <a href=\"http://bioinformatics.oxfordjournals.org/content/30/21/3012\">A new approach for detecting riboswitches in DNA sequences</a>.  <i>Bioinformatics</i> 30(21):3012-3019, 2014.\n</small>\n</td>\n<td style=\"text-align:right\">\n<p><a href='http://denison.edu/academics/computer-science'>Department of Mathematics and Computer Science</a>\n<br />Denison University</a>\n<br />Granville, Ohio, USA\n</td>\n</tr>\n</table>\n<!-- Footer Table End -->\n" << endl;
	overall_trunc.close();

	//Close the sequence file
	sequences.close();

	return 0;
}
