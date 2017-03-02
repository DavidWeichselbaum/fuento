/*
FUENTO - the FUnctional ENrichment TOol
Copyright 2016 David Weichselbaum

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Fuento was developed in the course of a Master's thesis in Bojan Zagrovic's group at the MFPL.
Exact Fisher's Test adapted from https://github.com/RabadanLab/SAVI/blob/master/bin/pileup2multiallele_vcf.cpp
Using UniProt and EBI APIs.
Thanks to Anton Polyanski and Matea Hajnic for beta testing.
*/

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <string.h>			// strerror(), strstr()
#include <fstream>			// filestream
#include <sstream>			// stringstream
#include <cmath>			// log(), exp()
#include <unistd.h>			// getcwd()
#include <sys/stat.h>			// mkdir(), stat()
#include <sys/types.h>			// stat()
#include <errno.h>			// error handling
#include <limits.h>			// realpath()
#include <time.h>			// datestamp
#include <curl/curl.h>			// downloading
#include <getopt.h>			// arguments
#include <boost/format.hpp>		// output
#include <boost/regex.hpp>		// regex_replace(), find
#include <boost/lexical_cast.hpp>	// string formating of numbers
#include <boost/algorithm/string/trim.hpp> // triming
#include <boost/algorithm/string.hpp>	// split()
#include <boost/range/algorithm_ext/erase.hpp>

using namespace std;

const string home 		= getenv("HOME");
const string directoryName 	= home + "/.fuento/";
const string database		= "http://www.geneontology.org/ontology/go.obo";
const string oboName 		= "go.obo";
const string backListName 	= "background.lst";
const string slimListName	= "GOslim.lst";
const string defaultName	= "defaults.dat";
      string oboFileName 	= directoryName + oboName;		// has to be non-constant, is changed if GO slim is used in analyis
const string backListFileName 	= directoryName + backListName;
const string slimListFileName	= directoryName + slimListName;
const string defaultFileName	= directoryName + defaultName;
string red, green, yellow, cyan, magenta, reset, redErr, resetErr; 	// for colors

void add_file_label(const string fileName, const string label);

void grow_log_facs(vector<double> &logFacs, const int n){					// establish list of logarithmic factorials 
	if(logFacs.size() > n){ return; }
	logFacs.reserve(n);
	for(int i=logFacs.size(); i<=n; ++i){
		logFacs.push_back( logFacs[i-1] + log(double(i)) );
	}
}

double log_hypergeometric_prob(const vector<double> &logFacs, const int a, const int b, const int c, const int d){	// log of hypergeometric probability P(a,b,c,d) = ((a+b)!*(c+d)!*(a+c)!*(b+d)!) / (a!*b!*c!*d!*(a+b+c+d)!)
	return logFacs[a+b] + logFacs[c+d] + logFacs[a+c] + logFacs[b+d] - logFacs[a] - logFacs[b] - logFacs[c] - logFacs[d] - logFacs[a+b+c+d];
}

double hypergeometric_prob(vector<double> &logFacs, const int a, const int b, const int c, const int d){	// hypergeometric probability P(a,b,c,d) = ((a+b)!*(c+d)!*(a+c)!*(b+d)!) / (a!*b!*c!*d!*(a+b+c+d)!)
	grow_log_facs(logFacs, a+b+c+d);								// ensure factorials list is big enough
	return exp(log_hypergeometric_prob(logFacs, a, b, c, d));
}

double fisher_test_less(vector<double> &logFacs, map<int,map<int,map<int,map<int,double> > > > &pBuffFisher, const int a, const int b, const int c, const int d){
	if(pBuffFisher[a][b][c][d] != 0){ return pBuffFisher[a][b][c][d]; }				// if table allready known, return its value
	grow_log_facs(logFacs, a+b+c+d);								// ensure factorials list is big enough
	double pSum = 0;
	for(int x=0; x <= a; ++x){								
		if( a+b-x >= 0 && a+c-x >= 0 && d-a+x >=0 ){					// guard against zero fields in contingency table
			pSum+= exp( log_hypergeometric_prob(logFacs, x, a+b-x, a+c-x, d-a+x) );	// add up probabilities more extreme than given
		}
	}
	pBuffFisher[a][b][c][d] = pSum;								// remember value and table
	return pSum;										
}

double fisher_test_greater(vector<double> &logFacs, map<int,map<int,map<int,map<int,double> > > > pBuffFisher, const int a, const int b, const int c, const int d){
	if(pBuffFisher[a][b][c][d] != 0){ return pBuffFisher[a][b][c][d]; }				// if table allready known, return its value
	grow_log_facs(logFacs, a+b+c+d);								// ensure factorials list is big enough
	double pSum = 0;
	for(int x=a; x <= (a+b+c+d); ++x){
		if( a+b-x >= 0 && a+c-x >= 0 && d-a+x >=0 ){                                    // guard against zero fields in contingency table
			pSum+= exp( log_hypergeometric_prob(logFacs, x, a+b-x, a+c-x, d-a+x) );   // add up probabilities more extreme than given
		}
	}
	pBuffFisher[a][b][c][d] = pSum;								// remember value and table
	return pSum;
}

double log_n_over_k(vector<double> &logFacs, const int n, const int k){	// (n,k) = n!/(k!(n-k)!) = exp( log(k!) - log(n!) - log((n-k)!) )
	if(k==0){ return 1; }
	grow_log_facs(logFacs, n);	// 0 <= k <= n	
	return logFacs[n] - logFacs[k] - logFacs[n-k];
}

double n_over_k(vector<double> &logFacs, const int n, const int k){	// (n,k) = n!/(k!(n-k)!) = exp( log(k!) - log(n!) - log((n-k)!) )
	return exp(log_n_over_k(logFacs, n, k));
}

double binomial_pmf(vector<double> &logFacs, const int n, const int k, const double p){		// pmf(n,k,p) = (n,k) p^k (1-p)^(n-k)
	return exp( log_n_over_k(logFacs, n, k) + log(p)*k + log(1-p)*(n-k) );
}

double binomial_test_less(vector<double> &logFacs, map<int,map<int,map<double,double> > > &pBuff, const int n, const int k, const double p){
	if(pBuff[n][k][p] != 0){ return pBuff[n][k][p]; }				// if table allready known, return its value
	double P = 0;
	for(int i=k; i<=n; i++){ P += binomial_pmf(logFacs, n, i, p); }
	pBuff[n][k][p] = P; 
	return P;
}

double binomial_test_greater(vector<double> &logFacs, map<int,map<int,map<double,double> > > &pBuff, const int n, const int k, const double p){
	if(pBuff[n][k][p] != 0){ return pBuff[n][k][p]; }				// if table allready known, return its value
	double P = 0;
	for(int i=0; i<=k; i++){ P += binomial_pmf(logFacs, n, i, p); }
	pBuff[n][k][p] = P; 
	return P;
}

// sort functions, implemented separately for performance
bool fisherCompare(const std::pair<int, vector <double> >& firstElem, const std::pair<int, vector <double> >& secondElem) {
	return firstElem.second[0] < secondElem.second[0];}

bool fisherCorrCompare(const std::pair<int, vector <double> >& firstElem, const std::pair<int, vector <double> >& secondElem) {
	if(firstElem.second[1] == secondElem.second[1]){ return firstElem.second[0] < secondElem.second[0]; }	// sort by non-corrected values if equal
	return firstElem.second[1] < secondElem.second[1];}

bool binomialCompare(const std::pair<int, vector <double> >& firstElem, const std::pair<int, vector <double> >& secondElem) {
	return firstElem.second[2] < secondElem.second[2];}

bool binomialCorrCompare(const std::pair<int, vector <double> >& firstElem, const std::pair<int, vector <double> >& secondElem) {
	if(firstElem.second[3] == secondElem.second[3]){ return firstElem.second[2] < secondElem.second[2]; }	// sort by non-corrected values if equal
	return firstElem.second[3] < secondElem.second[3];}

bool hyperCompare(const std::pair<int, vector <double> >& firstElem, const std::pair<int, vector <double> >& secondElem) {
	return firstElem.second[4] < secondElem.second[4];}

bool hyperCorrCompare(const std::pair<int, vector <double> >& firstElem, const std::pair<int, vector <double> >& secondElem) {
	if(firstElem.second[5] == secondElem.second[5]){ return firstElem.second[4] < secondElem.second[4]; }	// sort by non-corrected values if equal
	return firstElem.second[5] < secondElem.second[5];}

bool foldCompare(const std::pair<int, vector <double> >& firstElem, const std::pair<int, vector <double> >& secondElem) {
	return firstElem.second[6] > secondElem.second[6];}					// a higher fold-enrichent is more significant

bool backNumCompare(const std::pair<int, vector <double> >& firstElem, const std::pair<int, vector <double> >& secondElem) {
	return firstElem.second[7] > secondElem.second[7];}					// sorted for more functions

bool setNumCompare(const std::pair<int, vector <double> >& firstElem, const std::pair<int, vector <double> >& secondElem) {
	return firstElem.second[8] > secondElem.second[8];}					// sorted for more functions

bool funNameCompare(const std::pair<int, vector <double> >& firstElem, const std::pair<int, vector <double> >& secondElem, vector<int>& functions, vector<int>& GOIDs, vector<string>& GOnames){
	int firstGOindex = find(GOIDs.begin(), GOIDs.end(), functions[firstElem.first]) - GOIDs.begin();
		if(firstGOindex > GOIDs.size()){ cerr << "  Error: sorting by function names" << endl; exit(1); }
	string firstName = GOnames[firstGOindex];
	int secondGOindex = find(GOIDs.begin(), GOIDs.end(), functions[secondElem.first]) - GOIDs.begin();
		if(secondGOindex > GOIDs.size()){ cerr << "  Error: sorting by function names" << endl; exit(1);}
	string secondName = GOnames[secondGOindex];
	return firstName < secondName;
}
bool expNameCompare(const std::pair<int, vector <double> >& firstElem, const std::pair<int, vector <double> >& secondElem, vector<int>& functions, vector<int>& GOIDs, vector<string>& GOexpls){
	int firstGOindex = find(GOIDs.begin(), GOIDs.end(), functions[firstElem.first]) - GOIDs.begin(); 
		if(firstGOindex > GOIDs.size()){ cerr << "  Error: sorting by function explanations" << endl; exit(1); } 
	string firstName = GOexpls[firstGOindex];
	int secondGOindex = find(GOIDs.begin(), GOIDs.end(), functions[secondElem.first]) - GOIDs.begin();
		if(secondGOindex > GOIDs.size()){ cerr << "  Error: sorting by function explanations" << endl; exit(1); }
	string secondName = GOexpls[secondGOindex];
	return firstName < secondName;
}

void check_permission(const FILE* file, const string message){
	if(file == NULL){
		if(errno == EACCES){ cerr << redErr << message << endl << resetErr; }
		else{ cerr << redErr << "  Something went wrong: " << strerror(errno) << endl << resetErr; }
		exit(2);
	}
}

const string current_date(){
	time_t	now = time(0);
	struct 	tm tstruct;
	char	buf[80];
	tstruct = *localtime(&now);
	strftime(buf, sizeof(buf), "%Y-%m-%d", &tstruct);
	return buf;
}

void create_background(const string nameBackIDs, string nameBack, const string backListFileName="", const string column=""){
	string ID, function, url, line;
	string currentID = "";
	string IDlist = "";
	string tmpFileName = nameBack + ".tmp";
	string outFileName = nameBack;	
	vector <string> IDs;
	char fullOutFilePath[99999];
	FILE *tmpFile = fopen(tmpFileName.c_str(), "w");
		check_permission(tmpFile, "  Permission denied for saving .tmp in folder of ID file.\n");
	realpath(outFileName.c_str(), fullOutFilePath);
	CURL* c = curl_easy_init();
	CURLcode err;
	ifstream file(nameBackIDs.c_str());
		if(! file){ cerr << redErr << "  Error reading background-ID file: " << nameBackIDs << endl << resetErr; exit(5); }
	ofstream outStream(outFileName.c_str());
		if(! outStream){ cerr << redErr << "  Error writing background file: " << outFileName << endl << resetErr; exit(5); }
	while(file >> ID){
		if(find (IDs.begin(), IDs.end(), ID) != IDs.end()){ continue; }
		IDs.push_back(ID);
	}
	if(column == "geneontology" || column == ""){	// geneontology by default
		cout << "Downloading geneontology annotation..." << endl;
		outStream << "# TYPE:\tgeneontology\n";	// add type explanation
		for(int i=0; i<IDs.size(); i++){	// get GOIDs
			IDlist = IDlist + IDs[i];
			IDlist = IDlist + ',';
			if((i+1)%100==0 || i==IDs.size()-1){
				url = "https://www.ebi.ac.uk/QuickGO/GAnnotation?protein=" + IDlist + "&format=tsv&col=proteinID,goID";  // format: tab separated column of protein ID and genontology ID
				curl_easy_setopt(c, CURLOPT_URL, url.c_str() );
				curl_easy_setopt(c, CURLOPT_WRITEDATA, tmpFile);
				err = curl_easy_perform(c);
				if(err){ cerr << redErr << "  Error: online service not reachable. Try later." << endl << resetErr; exit(1); }
				IDlist = "";
				cerr << boost::format("%5d/%d\t%3.2f%%\r") % (i+1) % IDs.size() % ((i+1)*double(100)/IDs.size());
			}
		}
	}else{						// uniprot function
		cout << "Downloading UniProt annotation: " << column << "..." << endl;
		outStream << "# TYPE:\t" << column << endl;// add type explanation
		for(int i=0; i<IDs.size(); i++){
			IDlist = IDlist + "accession:" + IDs[i];
			IDlist = IDlist + "%20OR%20";			// '20%' as space between OR
			if((i+1)%100==0 || i==IDs.size()-1){
				url = "http://www.uniprot.org/uniprot/?query=" + IDlist + "accession:xxx&columns=id," + column + "&format=tab";  // format: tab separated colums of protein ID and function. Accession:xxx to end 'OR' string
				curl_easy_setopt(c, CURLOPT_URL, url.c_str() );
				curl_easy_setopt(c, CURLOPT_WRITEDATA, tmpFile);
				err = curl_easy_perform(c);
				if(err){ cerr << redErr << "  Error: online service not reachable. Try later." << endl << resetErr; exit(1); }
				IDlist = "";
				cerr << boost::format("%5d/%d\t%3.2f%%\r") % (i+1) % IDs.size() % ((i+1)*double(100)/IDs.size());
			}
		}
	}
	cerr << "\r" << endl;
	curl_easy_cleanup(c);
	fclose(tmpFile);	
	ifstream tmpStream(tmpFileName.c_str());
		if(! tmpStream){ cerr << redErr << "  Error reading temporary background file: " << tmpFileName << endl << resetErr; exit(5); }
	while(getline(tmpStream, line)){
		istringstream lineStream(line);
		lineStream >> ID;
		if(ID != currentID){								// if first annotation for this ID
  			if (find (IDs.begin(), IDs.end(), ID) != IDs.end()){			// ignore every line not headed by input ID
				if(currentID != ""){ outStream << '\n';}
				outStream << ID;
				currentID = ID;
				IDs.erase(remove(IDs.begin(), IDs.end(), ID), IDs.end());
			}
		}
		if(ID == currentID){			// intentionally new condition
			while(getline(lineStream, function, ',')){
			boost::trim_left(function);
			outStream << '\t' << function; }
		}
	}
	outStream << endl;
	if(IDs.size() != 0){
		cout << IDs.size() << " element(s) have zero functions annotated. Not UniProt IDs? \n";
		for(int i=0; i<IDs.size(); i++){
			outStream << IDs[i] << '\n'; 	// IDs have to be added without functions  
			cout << IDs[i] << '\t';
		}
		cout << endl;
	}
	remove(tmpFileName.c_str());
	outStream.close();
	if(backListFileName != ""){		// add to background list file
		ofstream backListFile(backListFileName.c_str(), std::fstream::out | fstream::app);
		backListFile << fullOutFilePath << endl;
		backListFile.close();
	}
}

void get_db(){
	FILE *oboFile = fopen(oboFileName.c_str(), "w");
		check_permission(oboFile, "  Permission denied for saving geneontology database.\nTry \"sudo fuento -g\"");
	CURL* c = curl_easy_init();
	CURLcode err;
	cout << "Downloading from:\t" << database << " ..." << endl; // save go.obo file
	curl_easy_setopt( c, CURLOPT_URL, database.c_str() );
	curl_easy_setopt( c, CURLOPT_WRITEDATA, oboFile);
	err = curl_easy_perform(c);
	if(err){ cerr << redErr << "Error: online service not reachable. Try later." << endl << resetErr; exit(1); }
	curl_easy_cleanup( c );
	cout  << "Database file created:\t" << oboFileName << endl;
	fclose(oboFile);	
	ifstream oboStream(oboFileName.c_str());

	string line;
	vector<string> words;
	vector<string> slimDefs;
	vector<string> slimNames;
	vector<string> slimFileNames;
	while(getline(oboStream, line)){	// extract go.slim definitions from go.obo
		if(line == ""){ break; }	// end at end of header
		if(line.substr(0, 17) == "subsetdef: goslim"){
			boost::split(words, line, boost::is_any_of(" "));
			slimDefs.push_back(words[1]);
			slimFileNames.push_back(directoryName + words[1] + ".obo");
			boost::split(words, line, boost::is_any_of("\""));
			slimNames.push_back(words[1]);
		}
	}
	for(int i=0; i<slimDefs.size(); i++){	// extract go.slims from go.obo
		ofstream slimStream(slimFileNames[i].c_str());
			if(! slimStream){ cerr << redErr << "Could not save gene ontology slim file: " << slimFileNames[i] << endl << resetErr; exit(5); }
		slimStream << "# DEF:\t" << slimNames[i] << endl;
		boost::split(words, slimNames[i], boost::is_any_of(" "));
		if(words.size() > 0){ slimStream << "# LAB:\t" << words[0] << endl; }	// add automatic label
		string entry = "";
		bool found = false;
		oboStream.clear();
		oboStream.seekg(0);
		while(getline(oboStream, line)){
			if(line == "[Term]"){
				if(found){
					slimStream << entry;
				}
				entry = "";
				found = false;
			}
			if(line == "subset: " + slimDefs[i]){ found = true; }
			entry += line + "\n";
		}
		slimStream.close();
		cout << "Gene ontology slim file " << green << "\"" << slimNames[i] << "\"" << reset << " created: " << slimFileNames[i] << endl;
	}
	oboStream.close();
	fstream listFileStream(slimListFileName.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
		if(! listFileStream){ cerr << redErr << "Could not access slim list file " << slimListFileName << resetErr << endl; exit(5); }
	vector<bool> slimFileFound (false, slimFileNames.size());
	while(getline(listFileStream, line)){		// check if slims allready in listfile, else add them
		for(int i=0; i<slimFileNames.size(); i++){
			if(line == slimFileNames[i]){ slimFileFound[i] = true; }
		}
	}
	listFileStream.clear();
	listFileStream.seekg(0);
	for(int i=0; i<slimFileNames.size(); i++){
		if( ! slimFileFound[i]){
			listFileStream << slimFileNames[i] << endl;
			cout << "Added gene ontology slim file " << slimFileNames[i] << " to slim list file." << endl;
		}
	}
	listFileStream.close();
}

string get_go_date(ifstream & file){
	string line;
	int i=0;
	while(getline(file, line) && i++ < 2){}
	return line;
}

string common_prefix(vector<string> names){
	int minLength = names[0].length();
	char current;
	for(int i=1; i<names.size(); i++){
		if(names[i].length() < minLength){ minLength = names[i].length(); }
	}
	for(int i=0; i<minLength; i++){
		current = names[0][i];
		for(int j=1; j<names.size(); j++){
			if(current != names[j][i]){
				return names[0].substr(0, i);
			}
		}
	}
	return names[0].substr(0, minLength);
}

void create_ID_file(const string nameBackground, const string nameIDfile){	// creates a uniprot ID file from a background, to make a new (updated) background
	string line;
	string ID;
	ifstream backFile(nameBackground.c_str());
	ofstream IDfile(nameIDfile.c_str());
		if(! IDfile){ cerr << redErr << "  Can't write to file: " << nameIDfile << endl << resetErr; exit(2); }
	while(getline(backFile, line)){
		if(line[0] == '#'){ continue; }		// ignore comments
		istringstream lineStream(line);
		lineStream >> ID;	
		IDfile << ID << endl;
	}
	backFile.close();
	IDfile.close();
}

void merge_backgrounds(const string backgr1name, const string backgr2name, const char mergeFlag, string backgrOutName){
	string line, token;
	vector<string> pros1;
	vector<string> pros2;
	vector<string> intersection;
	vector< vector<string> > funcs1;
	vector< vector<string> > funcs2;
	ifstream backgr1(backgr1name.c_str());						// check files and open
		if(! backgr1){ cerr << redErr << "  Error reading background file: " << backgr1name << endl << resetErr; exit(2); }
	ifstream backgr2(backgr2name.c_str());
		if(! backgr2){ cerr << redErr << "  Error reading background file: " << backgr2name << endl << resetErr; exit(2); }
	while(getline(backgr1, line)){							// read in set1 and set2
		vector<string> funcsColumn;
		boost::split(funcsColumn, line, boost::is_any_of("\t"));
		pros1.push_back(funcsColumn[0]);					// add gene ID
		funcsColumn.erase(funcsColumn.begin());
		funcs1.push_back(funcsColumn);						// add functions
	}
	while(getline(backgr2, line)){
		vector<string> funcsColumn;
		boost::split(funcsColumn, line, boost::is_any_of("\t"));
		pros2.push_back(funcsColumn[0]);					// add gene ID
		funcsColumn.erase(funcsColumn.begin());
		funcs2.push_back(funcsColumn);						// add functions
	}
	for(int i=0; i<pros1.size(); i++){			// standard behavior is extension of set1 by set2
		int j = find(pros2.begin(), pros2.end(), pros1[i]) - pros2.begin();
		if(j < pros2.size()){
			for(int k=0; k<funcs2[j].size(); k++){
				int l = find(funcs1[i].begin(), funcs1[i].end(), funcs2[j][k]) - funcs1[i].begin(); 	// check if func already annotated
				if(l >= funcs1[i].size()){ funcs1[i].push_back(funcs2[j][k]); }
			}
			pros2.erase(pros2.begin() + j);
			funcs2.erase(funcs2.begin() + j);
		}
		else if(mergeFlag=='I'){			// intersection: remove not found set1 proteins
			pros1.erase(pros1.begin() + i);
			funcs1.erase(funcs1.begin() + i);
		}
	}
	if(mergeFlag=='U'){ 					// union: add not found set2 proteins
		pros1.insert(pros1.end(), pros2.begin(), pros2.end());
		funcs1.insert(funcs1.end(), funcs2.begin(), funcs2.end());
	}
	backgr1.close();
	backgr2.close();
	ofstream backgrOut(backgrOutName.c_str());
		if(! backgrOut){ cerr << redErr << "  Error writing background file: " << backgrOutName << endl << resetErr; exit(2); }	
	for(int i=0; i<pros1.size(); i++){			// write to file
		backgrOut << pros1[i] << "\t";
		for(int j=0; j<funcs1[i].size(); j++){ backgrOut << funcs1[i][j] << "\t"; }
		backgrOut << endl;
	}
	backgrOut.close();
	add_file_label(backgrOutName, "");
}

vector<string> get_file_types(const string fileName){
	string line;
	vector<string> out;
	ifstream file(fileName.c_str());
		if(! file){ cerr << redErr << "  Error reading background-ID file: " << fileName << endl << resetErr; exit(5); }
	while(getline(file, line)){
		if(line.substr(0,7) == "# TYPE:"){
			string rest = line.substr(8);
			boost::split(out, rest, boost::is_any_of("\t"));
			if(out[out.size()-1] == ""){ out.pop_back(); }		// prevent empty last type 
			file.close();
			return out;
		}
	}
}

int get_gofile_enties(const string fileName){
	ifstream file(fileName.c_str());
		if(! file){ cerr << redErr << "  Error reading GO slim file: " << fileName << endl << resetErr; exit(5); }
	int n=0;
	string line;
	while(getline(file, line)){ if(line == "[Term]"){ n++; }} 
	file.close();
	return n;
}

string get_file_definition(const string fileName){
	ifstream file(fileName.c_str());
		if(! file){ cerr << redErr << "  Error reading file: " << fileName << endl << resetErr; exit(5); }
	string line;
	string def = "Custom GO slim";
	vector<string> words;
	while(getline(file, line)){
		if(line.substr(0,6) == "# DEF:"){
			boost::split(words, line, boost::is_any_of("\t"));
			if(words.size() >= 2){		// saveguard against missing definition
				def = words[1];
				break;
			}
		}
	} 
	file.close();
	return def;
}

string get_file_label(const string fileName){
	ifstream file(fileName.c_str());
		if(! file){ cerr << redErr << "  Error reading file: " << fileName << endl << resetErr; exit(5); }
	string line;
	string label = "";
	vector<string> words;
	while(getline(file, line)){
		if(line.substr(0,6) == "# LAB:"){
			boost::split(words, line, boost::is_any_of("\t"));
			if(words.size() >= 2){		// saveguard against missing label
				label = words[1];
				break;
			}
		}
	} 
	file.close();
	return label;
}

int get_file_lines(const string fileName){
	ifstream file(fileName.c_str());
		if(! file){ cerr << redErr << "  Error reading file: " << fileName << endl << resetErr; exit(5); }
	string line;
	int count = 0;
	while(getline(file, line)){
		if(line[0] != '#'){
			count++;
		}
	} 
	file.close();
	return count;
}

void add_file_label(const string fileName, const string label){
	string tmpFileName = fileName + ".tmp";
	ifstream file(fileName.c_str());
		if(! file){ cerr << redErr << "  Error reading file: " << fileName << endl << resetErr; exit(5); }
	ofstream tmpFile(tmpFileName.c_str());
		if(! file){ cerr << redErr << "  Error reading file: " << tmpFileName << endl << resetErr; exit(5); }
	tmpFile << "# LAB:\t" << label << endl;
	string line;
	while(getline(file, line)){
		if(line.substr(0,6) == "# LAB:"){ continue; }		// leave out old label
		tmpFile << line << endl;
	}
	file.close();
	tmpFile.close();
	if( rename(tmpFileName.c_str(), fileName.c_str()) != 0){ cerr << redErr << "  Error renaming file: " << tmpFileName << " to: " << fileName << endl << resetErr; exit(5); }
}

size_t write_data(void *buffer, size_t size, size_t nmemb, void *userp) 		// used with curl to prevent output, when only redirect is used
{ return size * nmemb; }

string map_ids(const string fileName, const string mapID){
	string tmpFileName = fileName + ".tmp";
	string outFileName = fileName + ".acc";
	bool filterMessage = false;
	char *redUrl;
	string list, ID, ID1, ID2, url, line;
	string IDstring = "";
	string filters = "+";
	string currentID = "";
	CURL* c;
	CURLcode err;
	boost::regex exp("http://www.uniprot.org/mapping/([A-Z0-9]*?).tab");
	boost::smatch match;
	vector <string> IDs;
	vector <string> IDsBuff;
	vector <string> currentMaps;

	ifstream file(fileName.c_str());
		if(! file){ cout << redErr << "  Error reading ID file: " << fileName << endl << resetErr; exit(5); }
	ofstream outFile(outFileName.c_str());
		if(! outFile){ cerr << redErr << "  Error opening mapped background-ID file: " << outFileName << endl << resetErr; exit(5); }
	while(getline(file, line)){
		istringstream lineStream(line);
		lineStream >> ID;
		if(ID == "#"){
			filters = line.substr(1);				// remove comment symbol
			cout << "Applying filters:" << filters << endl;
			replace(filters.begin(), filters.end(), ' ', '+'); 	// format filter string
			filters = filters + "+";
			continue;
		}
		if(find (IDs.begin(), IDs.end(), ID) != IDs.end()){	// prevent double entries
			cout << "  Found " << mapID << " entry \"" << ID << "\" multiple times. It will only be mapped to once." << endl;
			continue;
		}
		IDs.push_back(ID);
	}

	cout << "Mapping " << mapID << " entries from " << fileName << " to UniProt Accession (ACC)" << endl;
	for(int i=0; i<IDs.size(); i++){
		IDstring = IDstring + IDs[i] + ",";
		IDsBuff.push_back(IDs[i]);
		if((i+1)%50==0 || i==IDs.size()-1){
			FILE *tmpFile = fopen(tmpFileName.c_str(), "w");
				check_permission(tmpFile, "  Permission denied for saving .tmp in folder of ID file.\n");
			c = curl_easy_init();
			url = "http://www.uniprot.org/mapping/?from=" + mapID + "&to=ACC&format=tab&query=" + IDstring;	// gene id mapping
			curl_easy_setopt(c, CURLOPT_URL, url.c_str() );
			curl_easy_setopt(c, CURLOPT_FOLLOWLOCATION, 1L);				// only interested in redirect
			curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, write_data);				// prevent unused output
			err = curl_easy_perform(c);
				if(err){ cerr << redErr << "  Error: online service not reachable. Try later." << endl << resetErr; exit(1); }
			err = curl_easy_getinfo(c, CURLINFO_EFFECTIVE_URL, &redUrl);
				if(err){ cerr << redErr << "  Error: online service not reachable. Try later." << endl << resetErr; exit(1); }
			if(boost::regex_search(string(redUrl), match, exp)){
				list = string(match[1].first, match[1].second);							// get uniprot list string
			} else { cerr << "Error: could not resolve UniProt mapping request. Try again later." << endl; exit(0);} 

			c = curl_easy_init();
			url = "http://www.uniprot.org/uniprot/?query=yourlist%3A" + list + filters + "&format=tab&columns=yourlist(" + list + "),id&sort=yourlist:" + list; // use uniprot list as reference for filtering, output with column of mapped-from ID and sort for it
			curl_easy_setopt(c, CURLOPT_URL, url.c_str() );
			curl_easy_setopt(c, CURLOPT_FOLLOWLOCATION, 1L);
			curl_easy_setopt(c, CURLOPT_WRITEDATA, tmpFile);
			err = curl_easy_perform(c);
				if(err){ cerr << redErr << "  Error: online service not reachable. Try later." << endl << resetErr; exit(1); }
			fclose(tmpFile);

			ifstream tmpStream(tmpFileName.c_str());
				if(! tmpStream){ cerr << redErr << "  Error reading temporary background file: " << tmpFileName << endl << resetErr; exit(5); }
			while(getline(tmpStream, line)){
				istringstream lineStream(line);
				if(line.substr(0,8) == "yourlist"){ continue; } 						// ignore header
				lineStream >> ID1 >> ID2;
				IDsBuff.erase(remove(IDsBuff.begin(), IDsBuff.end(), ID1), IDsBuff.end());
				if(ID1 != currentID){
					if(currentMaps.size() > 1){
						cout << "  Supplied ID \"" << currentID << "\" of type \"" << mapID << "\" is ambiguous and has multiple entries (first will be picked): " << boost::algorithm::join(currentMaps, "; ") << endl;
						filterMessage = true;
					}
					if(currentMaps.size() > 0){ outFile << currentMaps[0] << endl; }						// write mapped ID file
					currentID = ID1;
					currentMaps.clear();
				}
				currentMaps.push_back(ID2);
			}
			if(currentMaps.size() > 1){
				cout << "  Supplied ID \"" << currentID << "\" of type \"" << mapID << "\" is ambiguous and has multiple entries: " << boost::algorithm::join(currentMaps, "; ") << endl;
				filterMessage = true;
			}
			if(currentMaps.size() > 0){ outFile << currentMaps[0] << endl; }						// write mapped ID file
			currentMaps.clear();

			if(IDsBuff.size()>0){
				cout << "  " << IDsBuff.size() << " identifiers could not be mapped to UniProt ACC: " << boost::algorithm::join(IDsBuff, "; ") << endl;;
			}
			
			cerr << boost::format("%5d/%d\t%3.2f%%\r") % (i+1) % IDs.size() % ((i+1)*double(100)/IDs.size()); 	// display progress
			IDstring = "";
			IDsBuff.clear();
		}
	}
	curl_easy_cleanup(c);
	remove(tmpFileName.c_str());
	outFile.close();
	if(filterMessage){ cout << "Ambiguous entries occured. Maybe add filters in header of ID file? See --help." << endl;}
	cout << "Created file: " << outFileName << endl;
	return outFileName;
}

void copy_file(const string srcName, const string dstName){
	ifstream  src(srcName.c_str(), ios::binary);
		if(! src){ cerr << redErr << "  Could not open file: " << srcName << resetErr; exit(1); }
	ofstream  dst(dstName.c_str(), ios::binary);
		if(! dst){ cerr << redErr << "  Could not open file: " << dstName << resetErr; exit(1); }
	dst << src.rdbuf();
	src.close();
	dst.close();
}

int main(int argc, char *argv[]){
	bool pipeTrue 		= ! isatty(STDIN_FILENO);		// checking for pipe input
	if(isatty(STDOUT_FILENO)){					// only use colors if output to terminal
		green 	= "\033[1;32m";
		cyan	= "\033[0;36m";
		magenta	= "\033[0;35m";
		yellow	= "\033[1;33m";
		red 	= "\033[0;31m";
		reset	= "\033[0m";
	}else{								// if output to pipe, use no colors
		green 	= "";
		cyan	= "";
		magenta	= "";
		yellow	= "";
		red 	= "";
		reset	= "";
	}
	if(isatty(2)){							// if sterr is written to terminal, use red errors
		redErr	= "\033[0;31m";
		resetErr= "\033[0m";
	}
	else{
		redErr	= "";						// if sterr to pipe, use no color
		resetErr= "";
	}
	
	bool namespaceArg	= false;
	bool columnsArg		= false;
	bool filterArg		= false;
	bool reverseChoice 	= false;
	bool fileChoice 	= false;
	bool obsoleteChoice 	= false;
	bool backRefChoice 	= false;
	bool printProt 		= false;
	bool printGOIDs 	= false;
	bool printExpl 		= false;
	bool printReadable 	= true;
	bool printMin		= false;
	bool listChoice		= false;
	bool mergeChoice	= false;
	bool mergeAddChoice	= true;
	bool getDbChoice	= false;
	bool updateChoice	= false;
	bool filterIDtrue	= false;
	bool filterStringTrue	= false;
	bool modusF		= true;
	bool modusC		= false;
	bool modusP		= false;
	bool modusA		= true;
	bool fisherTrue		= true;
	bool binomialTrue	= false;
	bool hyperTrue		= false;
	bool foldTrue		= true;
	bool backNumTrue	= true; 
	bool setNumTrue		= true;
	bool headerTrue		= true;
	bool fisherCorr		= true;
	bool binomialCorr	= false;
	bool hyperCorr		= false;
	bool sortReverse	= false;
	bool allTrue		= false;
	bool defaultChoice	= false;
	bool defaultFileTrue 	= false;
	char option;
	char modus;
	char mergeFlag;
	char* fullPath;
	char sortChar		= 'F';
	char corr		= 'B';
	int status;
	int backNum;
	int sortIndex		= 0;
	int option_index 	= 0;
	int slimNum		= 0;
	int maxNum 		= 0;
	int setNum 		= 0;
	int index		= 0;
	int randNum 		= 100;
	double definedCutoff 	= -1.0;
	double upperBbound 	= 0.0;
	double fdr		= 0.05;
	string nameBackgr;
	string nameSlim;
	string modName;
	string line;
	string filterString;
	string nameSlimShort	= "";
	string outFileName 	= "";
	string newBackName	= "";
	string backFunction	= "";
	string addListFileName	= "";
	string addListName	= "";
	string mergeBackgr1name	= "";
	string mergeBackgr2name	= "";
	string mergeBackgrOutName = "";
	string mapID	 	= "";
	string sep		= "\t";
	string modi		= "FA";
	string columns		= "1FfNnE2x";
	string displayColumns	= "FfNnEx";
	string goDate		= "no date";
	string namespaceIn	= modi;
	string columnsIn	= columns;
	string filterIn		= "";
	vector<bool> filterStringsFound;
	vector<bool> filtersFound;
	vector<int> filters;
	vector<string> backgroundList;
	vector<string> backgroundLabList;
	vector<string> slimList;
	vector<string> slimLabList;
	vector<string> filterList;
	vector<string> fileNames;
	vector<string> filterStrings;

	// check for saved defaults
	ifstream defaultFile(defaultFileName.c_str());
	if(defaultFile){
		defaultFileTrue = true;
		string rest;
		while(getline(defaultFile, line)){
			rest = line.substr(2);
			istringstream lineStream(rest);
			switch(line[0]){
				case 'n':
					if(rest.size() == 0){ break; } 
					namespaceArg = true;
					namespaceIn = rest;
					break;
				case 'c':
					if(rest.size() == 0){ break; } 
					columnsArg = true;
					columnsIn = rest;
					break;
				case 'C':
					lineStream >> corr;
					lineStream >> fdr;
					break;
				case 'f':
					if(rest.size() == 0){ break; } 
					filterArg = true;
					filterIn = rest;
					break;
				case 's':
					nameSlim = rest;
					break;
				case 'x':
					definedCutoff = atof(rest.c_str());
					break;
				case 'a':
					allTrue = atoi(rest.c_str());
					break;
				case 'r':
					reverseChoice = atoi(rest.c_str());
					break;
				case 'm':
					maxNum = atoi(rest.c_str());
					break;
				case 't':
					randNum = atoi(rest.c_str());
					break;
				case 'G':
					printMin = atoi(rest.c_str());
					break;
				case 'H':
					headerTrue = atoi(rest.c_str());
					break;
				case 'e':
					obsoleteChoice = atoi(rest.c_str());
					break;
				case 'A':
					backRefChoice = atoi(rest.c_str());
					break;
				case 'S':
					sep = rest;
					break;
				default:
					cerr << redErr << "Error: default file " << defaultFileName << " corrupted. Deleting it." << endl << resetErr;
					remove(defaultFileName.c_str());
					exit(1);
					break;
			}

		}
	}
	defaultFile.close();

	ostringstream usageStream;
	usageStream << 		"FUENTO - the FUnctional ENrichment TOol       GNU Public Licence, David Weichselbaum 2016\n"
				"\n"
				"Tests sets of genes for significant enrichment of gene functions against a given background.\n"
				"The user can generate and archive backgrounds automatically. In- and output is subject to a variety of filters.\n"
				"Genes are anotated by the gene ontology by default, but any annotation can be used. Most gene IDs can be mapped\n"
				"to the internal standard 'UniProt accession ID'. Sets and backgrounds are not restricted to enrichment of genes\n"
				"and can contain any string identifier instead of both genes and functions. By default, fuento contains output\n"
				"to results which are significant in an permutation test using 100 random protein sets.\n"
				"\n"
				"usage:  fuento [OPTION...] [ BACKGROUND [ SET... ] ]\n"
				"  options:\n"
				"\n";
	
	usageStream <<	boost::format(
				"   -n --namespace <STRING>    Gene Ontology namespace(s) displayed.\n"
				"      default: '--namespace %s'. Components of STRING:\n"
				"        F:  molecular function\n"
				"        C:  cellular component\n"
				"        P:  biological process\n"
				"        A:  aberrant function\n") % namespaceIn;
	usageStream <<	boost::format(
				"   -c --columns <STRING>      columns to display and statistical test(s) to perform. Results will be sorted by first test column.\n"
				"      default: '--columns %s'. Components of STRING:\n" 
				"        F/f:  Fisher's exact test / with multiple-hypothesis correction\n"			// sort index: 0/1
				"        B/b:  binomial test       / with multiple-hypothesis correction\n"			//             2/3
				"        H/h:  hypergeometric test / with multiple-hypothesis correction\n"			//             4/5
				"         E :  fold-enrichment\n"								//              6
				"        N/n:  number of functions in Background / number of functions in set\n"		//             7/8
				"        x/X:  short/long explanation of function\n"						//	       9/10
				"         P :  protein IDs\n"
				"         G :  Gene Ontology IDs\n"
				"        1-5:  highlight next column: 1:green 2:cyan 3:magenta 4:yellow 5:red\n") % columnsIn;
	string corrStr;
	if(corr == 'B' || corr == 'b'){ corrStr = (boost::format("%c") % corr).str(); }
	else{ 				corrStr = (boost::format("%c %.1e") % corr % fdr).str(); }
	usageStream << boost::format(
				"   -C --correction <STRING>   multiple-hypothesis correction \n"
				"      default: '--correction %s'. Possible values for STRING:\n" 
				"         B:           Bonferroni correction, returns corrected p-values\n"
				"         C <NUMBER>:  Benjamini-Hochberg FDR (False Discovery rate) correction, returns corrected p-values for a given FDR=NUMBER\n"
				"         A <NUMBER>:  Benjamini-Hochberg-Yekutieli FDR adjustment, returns adjusted p-values for a given FDR=NUMBER\n"
				"         T <NUMBER>:  Benjamini-Hochberg test, returns only 0=significant, 1=insignificant for a given FDR=NUMBER.\n"
				"                      If a Benjamini-Hochberg test corrected column leads, it disables permutation tests, and displays all significant numbers.\n") % corrStr;
	string filterStr = "";
	if(filterArg){ filterStr = (boost::format("default: \'%s\'. ") % filterIn).str(); }
	usageStream << boost::format(
				"   -f --filter <STRING>       filter for GOIDs or plaintext functions separated by semicolon.\n"
				"      %sIf FILTER(s) contains whitespace, enclose in quotes. Set '--all' to see all possible functions.\n"
				"      Plaintext functions can be filtered by regex.\n") % filterStr;
	string slimStr = "u"; // for the word 'use' in the below usage string
	if(nameSlim != ""){ slimStr = (boost::format("default: \'%s\'. U") % nameSlim).str(); }
	usageStream << boost::format(
				"   -s --slim <FILE/NUMBER/LABEL>\n"
				"      %sse GO slim instead of full go.obo. Slims can be viewed with '-l'\n" ) % slimStr;
	string cutoffStr = "";
	if(definedCutoff >= 0){ cutoffStr = (boost::format("Default: \'%f\'. ") % definedCutoff).str(); }
	usageStream << boost::format(
				"   -x --cutoff <NUMBER>       cutoff for first column (p-value, fold-enrichment, function number)\n"
				"      %sCan be supplied in scientific or standard notation. Applied only to first test column\n") % cutoffStr;
	string reverseStr = "r";
	if(reverseChoice){ reverseStr = "default: set. R"; }
	usageStream << boost::format(
				"   -r --reverse               %severse enrichment, display functional depletion for all columns.\n") % reverseStr;
	string allStr = "d";
	if(allTrue){ allStr = "default: set. D"; }
	usageStream << boost::format(
				"   -a --all                   %siplay all enrichments.\n") % allStr;
	string maxStr = "m";
	if(maxNum){ maxStr = (boost::format("default: %d. M") % maxNum).str(); }
	usageStream << boost::format(
				"   -m --max <NUMBER>          %sax number of functions displayed (overrides -t).\n") % maxStr;
	string trialStr = "r";
	if(randNum){ trialStr = (boost::format("default: %d. R") % randNum).str(); }
	usageStream << boost::format(
				"   -t --trial_number <NUMBER> %sandom trial number.\n") % trialStr;
	string globalStr = "p";
	if(printMin){ globalStr = "default: set. P"; }
	usageStream << boost::format(
				"   -G --global_minimum        %srint only minimum p-value of each set\n") % globalStr;
	string headStr = "d";
	if(! headerTrue){ headStr = "default: set. D"; }
	usageStream << boost::format(
				"   -H --no_header             %so not display header explaining columns. Will print only filename of set instead.\n") % headStr;
	string obsoleteStr = "d";
	if(obsoleteChoice){ obsoleteStr = "default: set. D"; }
	usageStream << boost::format(
				"   -e --obsolete              %so not ignore obsolete go entries\n") % obsoleteStr;
	string backStr = "b";
	if(backRefChoice){ backStr = "default: set. B"; }
	usageStream << boost::format(
				"   -A --back_reference        %sack-reference functions ('is_a:' marker in go.obo)\n"
				"      Only recommended if background was made by hand, since \'-b\' accounts for that.\n") % backStr;
	if(sep == "\t"){ usageStream <<
				"   -S --sep <STRING>          specify separator (default: \'-sep \"\\t\"\')\n";
	}else{            usageStream << boost::format(
				"   -S --sep <STRING>          specify separator (default: \'-sep \"\%s\"\')\n") % sep; }
	usageStream <<
				"   -d --defaults              set current options as custom defaults. Use without options to restore defaults.\n"
				"        Defaults are setable for all arguments listed above, but not for the ones listed below.\n"
				"        Currently set defaults can be viewed with \'-h\'.\n"
				"\n"
				"   -b --background <FILE>     create background file from uniprot-accid file.\n"
				"      Downloads newest annotation from ebi server. This may take a while.\n"
				"   -B --background_function <FILE> <STRING>\n"
				"      Create background file from uniprot-accid file with function type STRING.\n"
				"      STRING should be a uniprotkb column name explained here: 'www.uniprot.org/help/uniprotkb_column_names'.\n"
				"   -M --map <STRING> [<FILE(S)>]\n"
				"      maps gene IDs of type STRING to uniprot-accid when creating backgrounds with '-b' or '-B', or maps FILE(S) to\n" 
				"      FILES_STRING. Available gene ID types are found here: 'http://www.uniprot.org/help/programmatic_access#id_mapping_examples'\n"
				"      If several IDs can be mapped to, the first in the list is used.\n"
				"   -l --list                  lists available background files and slims. columns: number, entries, explanation, [label], file\n"
				"   -L --add <FILE> <LABEL>    adds FILE to list of background/slim files with custom LABEL.\n"
				"      Slim FILEs must have the ending '.obo'. to remove entry, delete or rename file.\n"
				"   -R --merge <FILE_A> <FILE_B> <OPERATOR> <STRING>\n"
				"      merge two backgrounds FILEs A/B and save as file STRING. Possible values for OPERATOR:\n"
				"          'u': union of a and b, keep all IDs and functions.\n"
				"          'i': intersection of a and b, keep IDs that occur in both sets.\n"
				"          'e': extend a by b, keep only IDs of a, extend functions from b.\n"
				"      If OPERATOR is lowercase, background is not added to list (useful for scripting).\n"
				"\n"
				"   -g --get                   get newest Gene Ontology database (go.obo) from uniprot server. Saves to \"~/.fuento\"\n"
				"   -u --update                update everything, including Gene Ontology database and backgrounds.\n"
				"        Generates new background files in the same folder with a new timestamp.\n"
				"   -o --output <FILE>         specify output FILE instead of stdout.\n"
				"   -O --auto_output           generate output file of the longest common prefix of set names.\n"
				"      Will be saved as \".enr\".\n"
				"   -v --version               display version and exit\n"
				"   -h --help                  display this and exit\n"
				"\n"
				" background: <FILE/NUMBER/LABEL>\n"
				"   Background files can be specified by a filename, number of the background in the list or label (see with -l).\n"
				"   Labels can be added to backgrounds and go slims with the '-l' option or by heading the file with: '# lab:<tab>label'\n"
				" set: <FILES/FILES*>\n"
				"   Using wildcards and globbing is possible since all strings following the background will be assumed sets.\n"
				"   Piping whitespace separated uniprot IDs into fuento is an alternative to specifying protein sets.\n"
				"   If piping is used, sets can be divided by the word \'END\'.\n";
	const string usage = usageStream.str();
	const string version = "fuento 0.17\nDavid Weichselbaum 31.1.2017\n";
	
	static struct option long_options[] =
		{
	        {"namespace",		required_argument,		0, 'n'},
	        {"columns",		required_argument,		0, 'c'},
	        {"correction",		no_argument,			0, 'C'},
	        {"filter",		required_argument,		0, 'f'},
	        {"slim",		required_argument,		0, 's'},
	        {"cutoff",		required_argument,		0, 'x'},
	        {"all", 		no_argument,       		0, 'a'},
	        {"reverse", 		no_argument,       		0, 'r'},
	        {"max",			required_argument,		0, 'm'},
	        {"trial_number",	required_argument,		0, 't'},
	        {"global_minimum", 	no_argument,       		0, 'G'},
	        {"no_header", 		no_argument,       		0, 'H'},
	        {"obsolete", 		no_argument,       		0, 'e'},
	        {"back_reference",	no_argument,       		0, 'A'},
	        {"sep",			required_argument,		0, 'S'},
	        {"defaults",		no_argument,       		0, 'd'},
	        {"background",		required_argument,		0, 'b'},
	        {"background_function",	no_argument,			0, 'B'},
	        {"map",			required_argument,		0, 'M'},
	        {"list",		no_argument,       		0, 'l'},
	        {"add",			no_argument,			0, 'L'},
	        {"merge",		no_argument,			0, 'R'},
	        {"get",			no_argument,       		0, 'g'},
	        {"update",		no_argument,       		0, 'u'},
	        {"output",		required_argument,		0, 'o'},
	        {"auto_output",		no_argument,			0, 'O'},
	        {"version",		no_argument,			0, 'v'},
	        {"help",		no_argument,       		0, 'h'},
		{0,         		0,                 		0,  0 }
		};

	////////////////// get arguments
	if(argc > 1 && strstr(argv[1], "--help")){ cout << usage; return 0; }
	while((option = getopt_long(argc, argv, "n:c:Cf:s:x:arm:t:GHeAS:b:BM:lLRgudo:Ovh", long_options, &option_index)) != -1){
		switch (option){
		case 'n': // --namespace <STRING>
			namespaceArg = true; 
			namespaceIn = optarg;
			break;
		case 'c': // --columns <STRING>
			columnsArg = true;
			columnsIn = optarg;
			break;
		case 'C': // --correction
			if(optind >= argc){ cerr << redErr << "  Error correction argument\n" << resetErr; return 1; }
			corr = argv[optind++][0];
			if( ! (corr == 'B' || corr == 'b' || corr == 'C' || corr == 'c' || corr == 'A' || corr == 'a' || corr == 'T' || corr == 't')){ cerr << redErr << "  Error correction argument\n" << resetErr; return 1; }
			if(corr == 'C' || corr == 'c' || corr == 'A' || corr == 'a' || corr == 'T' || corr == 't'){
				if(optind >= argc){ cerr << redErr << "  Error correction argument\n" << resetErr; return 1; }
				fdr = atof(argv[optind++]);
				if(fdr <= 0){ cerr << redErr << "  Error correction argument: p-value has to be positive\n" << resetErr; return 1; }
			}
			break; 
		case 'f': // --filter <STRING>
			filterArg = true;
			filterIn = optarg;
			break;
		case 's': // --slim <FILE/NUMBER/LABEL>
			nameSlim = optarg;
			break;
		case 'x': // --cutoff <NUMBER>
			definedCutoff = atof(optarg);
			if(definedCutoff < 0){ cerr << redErr << "  Error cutoff argument\n" << resetErr; return 1; } 
			randNum = 0;		// no random trial
			break;
		case 'a': // --all
			allTrue = true;
			randNum = 0;		// no random trial
			break;
		case 'r': // --reverse
			reverseChoice = true;
			break;
		case 'm': // --max <NUMBER>
			maxNum = atoi(optarg);
			if(maxNum <= 0){ cerr << redErr << "  Error max display argument\n" << resetErr; return 1; }
			randNum = 0;		// no random trial
			break;
		case 't': // --trial_number <NUMBER> 
			randNum = atoi(optarg);
			if(randNum < 0){ cerr << redErr << "  Error random trial number argument\n" << resetErr; return 1; }
			break;
		case 'G': // --global_minimum
			printMin = true;
			randNum = 0;		// no random trial
			break;
		case 'H': // --no_header
			headerTrue = false;
			break;
		case 'e': // --obsolete
			obsoleteChoice = true;
			break;
		case 'A': // --back_reference
			backRefChoice = true;
			break;
		case 'S': // --sep <STRING>
			sep = optarg;
			break;
		case 'd': // --defaults
			defaultChoice = true;
			break;
		case 'b': // --background <FILE>
			newBackName = optarg;
			break;
		case 'B': // --background_function <FILE> <STRING>
			if(optind < argc && *argv[optind] != '-'){ newBackName = argv[optind++];  }else{ cerr << redErr << "  Error background argument\n" << resetErr; return 1; }
			if(optind < argc && *argv[optind] != '-'){ backFunction = argv[optind++]; }else{ cerr << redErr << "  Error background argument\n" << resetErr; return 1; }
			break;
		case 'M': // --map <STRING>
			mapID = optarg;
			break;
		case 'l': // --list
			listChoice = true;
			break;
		case 'L': // --add <FILE> <LABEL>
			if(optind < argc && *argv[optind] != '-'){ addListFileName = argv[optind++]; }else{ cerr << redErr << "  Error \"add file to list\" argument A\n" << resetErr; return 1;}
			if(optind < argc && *argv[optind] != '-'){ addListName = argv[optind++];     }else{ cerr << redErr << "  Error \"add file to list\" argument A\n" << resetErr; return 1;}
			break;
		case 'R': // --merge <FILE_A> <FILE_B> <OPERATOR> <STRING>
			mergeChoice = true;
			if(optind < argc && *argv[optind] != '-'){ mergeBackgr1name   = argv[optind++];  }else{ cerr << redErr << "  Error merge argument A\n" << resetErr; return 1;}
			if(optind < argc && *argv[optind] != '-'){ mergeBackgr2name   = argv[optind++];  }else{ cerr << redErr << "  Error merge argument B\n" << resetErr; return 1;}
			if(optind < argc && *argv[optind] != '-'){ mergeFlag   	      = *argv[optind++]; }else{ cerr << redErr << "  Error merge argument flag\n" << resetErr; return 1;}
			if(mergeFlag=='u' || mergeFlag=='i' || mergeFlag=='e'){ mergeAddChoice = false; }
			if( ! (mergeFlag=='U' || mergeFlag=='I' || mergeFlag=='E' || mergeFlag=='u' || mergeFlag=='i' || mergeFlag=='e') ){ cerr << redErr << "  Error merge argument flag. Allowed charters: U,I,E\n" << resetErr; return 1;}
			if(optind < argc && *argv[optind] != '-'){ mergeBackgrOutName = argv[optind++];  }else{ cerr << redErr << "  Error merge argument string\n" << resetErr; return 1;}
			break;
		case 'g': // --get
			getDbChoice = true;
			break;
		case 'u': // --update
			updateChoice = true;
			break;
		case 'o': // --output <FILE>
			if(optarg[0] == '-' || ((argc - optind > 0) && strstr(optarg, ".bkg") != NULL)){ cerr << redErr << "  Error output file argument\n" << resetErr; return 1; }
			fileChoice = true;
			outFileName = optarg;
			red = ""; green = ""; yellow = ""; cyan = ""; magenta = ""; reset = "";	// reset color codes
			break;
		case 'O': // --auto_output
			fileChoice = true;
			red = ""; green = ""; yellow = ""; cyan = ""; magenta = ""; reset = "";	// reset color codes
			break;
		case 'v': // --help
			cout << version;
			return 0;
		case 'h': // --help
			cout << usage;
			return 0;
		default: 
			cerr << "Error argument" << endl;
			return 0;
		}
	}
	if(namespaceArg){
		modusF = false; modusC = false; modusP = false; modusA = false; modi = "";
		if(strstr(namespaceIn.c_str(), "F") != NULL || strstr(namespaceIn.c_str(), "f") != NULL){ modusF = true; modi += "F"; }
		if(strstr(namespaceIn.c_str(), "C") != NULL || strstr(namespaceIn.c_str(), "c") != NULL){ modusC = true; modi += "C"; }
		if(strstr(namespaceIn.c_str(), "P") != NULL || strstr(namespaceIn.c_str(), "p") != NULL){ modusP = true; modi += "P"; }
		if(strstr(namespaceIn.c_str(), "A") != NULL || strstr(namespaceIn.c_str(), "a") != NULL){ modusA = true; modi += "A"; }
		if(! (modusF || modusC || modusP || modusA)){ cerr << redErr << "  Error modus argument\n" << resetErr; return 1; }
	}
	if(columnsArg){
		columns = columnsIn;
		fisherTrue = false; binomialTrue = false; hyperTrue = false; foldTrue = false;
		displayColumns = columns;
		boost::range::remove_erase_if(displayColumns, boost::algorithm::is_any_of("12345")); // remove color markers from displayed columns string
		if (columns.find_first_not_of("FfBbHhENnxXPG12345") != std::string::npos){ cerr << redErr << "  Error column argument\n" << resetErr; return 1;}
		if(strstr(columnsIn.c_str(), "F") != NULL){ fisherTrue = true; }		// see if test needs to be performed
		if(strstr(columnsIn.c_str(), "f") != NULL){ fisherTrue = true; fisherCorr = true; }
		if(strstr(columnsIn.c_str(), "B") != NULL){ binomialTrue = true; }
		if(strstr(columnsIn.c_str(), "b") != NULL){ binomialTrue = true; binomialCorr = true;}
		if(strstr(columnsIn.c_str(), "H") != NULL){ hyperTrue = true; }
		if(strstr(columnsIn.c_str(), "h") != NULL){ hyperTrue = true; hyperCorr = true;}
		if(strstr(columnsIn.c_str(), "E") != NULL){ foldTrue = true; }
		if(strstr(columnsIn.c_str(), "X") != NULL){ printExpl = true; }		// function explanations are only saved when protein is printed
		if(strstr(columnsIn.c_str(), "P") != NULL){ printProt = true; }		// protein affiliations are only saved when protein is printed
		sortChar = columns[ columns.find_first_of("FfBbHhENnxX") ];	// get index of function to sort by
		if     (sortChar == 'F'){ sortIndex = 0; }
		else if(sortChar == 'f'){ sortIndex = 1; }
		else if(sortChar == 'B'){ sortIndex = 2; }
		else if(sortChar == 'b'){ sortIndex = 3; }
		else if(sortChar == 'H'){ sortIndex = 4; }
		else if(sortChar == 'h'){ sortIndex = 5; }
		else if(sortChar == 'E'){ sortIndex = 6; sortReverse = true;}	// higher protein numbers and enrichments are more significant
		else if(sortChar == 'N'){ sortIndex = 7; sortReverse = true;}
		else if(sortChar == 'n'){ sortIndex = 8; sortReverse = true;}
		else if(sortChar == 'x'){ sortIndex = 9; }
		else if(sortChar == 'X'){ sortIndex = 10; }
		if(sortIndex == 0){ fisherTrue = true; }			// prevent situations in which no proper sort column test is specified
	}
	if(filterArg){
		boost::split(filterList, filterIn, boost::is_any_of(";"));
		for(int i=0; i<filterList.size(); i++){
			if(boost::regex_search(filterList[i], boost::regex("GO:[0-9]{7}"))){		// filter as goid
				filters.push_back( atoi(filterList[i].substr(3, 7).c_str()) );
				filtersFound.push_back(false);
				filterIDtrue = true;
			}else{										// filter as string
				filterStrings.push_back(filterList[i]);
				filterStringsFound.push_back(false);
				filterStringTrue = true;
			}
		}
	}
	////////////////// checking directory integrity
	status = mkdir(directoryName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);						// checking or creating "~/.fuento/" directory
	if(! status){ cerr << "Directory created: " << directoryName << endl; }
	else{
		if(! errno == EEXIST){ cerr << redErr << status << "  Error creating directory: " << directoryName << endl << resetErr; exit(2); }
	}		

	if(getDbChoice){													// argument: get database file (must be possible before checking for db integrity) 
		get_db();
		return 0;
	}

	ifstream mappingFile(oboFileName.c_str());
	if(! mappingFile){ cerr << redErr << "  Error reading GO mapping file: " << oboFileName <<  "\n  Try \"fuengo -g\"\n"  << resetErr; return 2;}// check mapping file
	goDate = get_go_date(mappingFile);
	mappingFile.close();


	fstream checkFile1(backListFileName.c_str(), std::fstream::out | std::fstream::app);							// create list file if not there
	if(! checkFile1){ cerr << redErr << "  Error reading/creating list file: " << backListFileName << endl  << resetErr; return 2;}
	checkFile1.close();
	fstream backListFile(backListFileName.c_str(), std::fstream::in);
	int nBacks = 0;
	while(backListFile >> line){
		nBacks++;
		if(find (backgroundList.begin(), backgroundList.end(), line) != backgroundList.end()){ continue; }
	  	ifstream background(line.c_str());
		if(background){
			backgroundList.push_back(line);
			backgroundLabList.push_back(get_file_label(line));
		}
		background.close();
	}
	backListFile.close();
	if(nBacks != backgroundList.size()){											// update listfile to remove non-existing backgrounds
		backListFile.open(backListFileName.c_str(), std::fstream::out | std::fstream::trunc);				// create new background file
		for(int i=0; i<backgroundList.size(); i++){
			backListFile << backgroundList[i] << endl;
		}
		backListFile.close();
	}

	fstream checkFile2(slimListFileName.c_str(), std::fstream::out | std::fstream::app);							// create list file if not there
	if(! checkFile2){ cerr << redErr << "  Error reading/creating slim file: " << slimListFileName << endl  << resetErr; return 2;}
	checkFile2.close();
	fstream slimListFile(slimListFileName.c_str(), std::fstream::in);		// check slim file
	int nSlims = 0;
	while(slimListFile >> line){
		nSlims++;
	  	ifstream slim(line.c_str());
		if(find (slimList.begin(), slimList.end(), line) != slimList.end()){ continue; }
		if(slim){
			slimList.push_back(line);
			slimLabList.push_back(get_file_label(line));
		}
		slim.close();
	}
	slimListFile.close();
	if(nSlims != slimList.size()){												// remove unnecessary entries if existing
		slimListFile.open(slimListFileName.c_str(), std::fstream::out | std::fstream::trunc);				// create new slim file
		for(int i=0; i<slimList.size(); i++){
			slimListFile << slimList[i] << endl;
		}
		slimListFile.close();
	}

	////////////////// working arguments
	if(defaultChoice && optind == 2){ 	// if --default was only argument, delete defaults
		cout << "Deleting custom defaults." << endl;
		remove(defaultFileName.c_str());
		return 0;
	}
	if(defaultChoice){
		cout << "Saving custom defaults." << endl;
		ofstream defaultFile(defaultFileName.c_str());
			if(! defaultFile){ cerr << redErr << "Error: Not able to wirte to default file " << defaultFileName << endl << resetErr; exit(1); }
		defaultFile << boost::format("n %s\n") % namespaceIn;
		defaultFile << boost::format("c %s\n") % columnsIn;
		defaultFile << boost::format("C %c %f\n") % corr % fdr;
		defaultFile << boost::format("f %s\n") % filterIn;
		defaultFile << boost::format("s %s\n") % nameSlim;
		defaultFile << boost::format("x %f\n") % definedCutoff;
		defaultFile << boost::format("a %d\n") % allTrue;
		defaultFile << boost::format("r %d\n") % reverseChoice;
		defaultFile << boost::format("m %d\n") % maxNum;
		defaultFile << boost::format("t %d\n") % randNum;
		defaultFile << boost::format("G %d\n") % printMin;
		defaultFile << boost::format("H %d\n") % headerTrue;
		defaultFile << boost::format("e %d\n") % obsoleteChoice;
		defaultFile << boost::format("A %d\n") % backRefChoice;
		defaultFile << boost::format("S %s\n") % sep;
		defaultFile.close();
		return 0;
	}
	if(mergeChoice){													// merge backgound files
		if(mergeAddChoice){
			string date = current_date();
			mergeBackgrOutName = mergeBackgrOutName + "_" + date + ".bkg";
			addListFileName = mergeBackgrOutName;
		}
		merge_backgrounds(mergeBackgr1name, mergeBackgr2name, mergeFlag, mergeBackgrOutName);
		if(! mergeAddChoice){ return 0; }										// if new background should not be added, return
	}
	if(addListFileName != ""){												// for adding files
	  	ifstream background(addListFileName.c_str());
		if(background){
			fullPath = realpath(addListFileName.c_str(), NULL);
			if(boost::regex_search(addListFileName, boost::regex(".obo"))){						// if '.obo', add to slim list
				if(find (slimList.begin(), slimList.end(), fullPath) == slimList.end()){
					slimList.push_back(fullPath);
					slimListFile.open(slimListFileName.c_str(), std::fstream::out | std::fstream::trunc);				// create new background file
					for(int i=0; i<slimList.size(); i++){
						slimListFile << slimList[i] << endl;
					}
					slimListFile.close();
					cout << "Added file to GO slim archive: " << fullPath << endl;
				}
			}
			else{
				if(find (backgroundList.begin(), backgroundList.end(), fullPath) == backgroundList.end()){
					backgroundList.push_back(fullPath);
					backListFile.open(backListFileName.c_str(), std::fstream::out | std::fstream::trunc);				// create new background file
					for(int i=0; i<backgroundList.size(); i++){
						backListFile << backgroundList[i] << endl;
					}
					backListFile.close();
					cout << "Added file to background archive: " << fullPath << endl;
				}
			}
			if(addListName != ""){
				add_file_label(addListFileName, addListName);
			}
		}
		else{ cerr << redErr << "  Background/slim file not readable: " << fullPath << endl << resetErr; }
		background.close();
		return 0;
	}
	if(newBackName != ""){													// create background
		string backFunctionBuff = backFunction;
		string idFile = newBackName;
		if(mapID != ""){												// maping argument with background leads to mapping and creation of background
			idFile = map_ids(newBackName, mapID);
		}
		if(backFunctionBuff == ""){	// create background file name
			newBackName = newBackName + "_" + current_date() + ".bkg";
		}else{			// if column is specified, remove special characters from column string before making a new background file name
	        	backFunctionBuff.erase(remove_if(backFunctionBuff.begin(), backFunctionBuff.end(), [](char c) { return !isalpha(c); } ), backFunctionBuff.end());
			newBackName = newBackName + "_" + backFunctionBuff + "_" + current_date() + ".bkg";
		}
		create_background(idFile, newBackName, backListFileName, backFunction);
		cout << "Created file: " << newBackName << endl;
		return 0;
	}
	if(mapID != ""){													// map file(s)   -  mapping argument alone leads to maping of following sets
		if(argc - optind < 1)	{ cerr << redErr <<  "  mandatory argument missing\n" << resetErr; return 1; }		// check for files to map
		for(int i=optind; i<argc; i++){											// check maping files
			struct stat statBuff;
			if(stat(argv[i], &statBuff) == -1){ cerr << redErr << "  Can't stat file: " << argv[i] << resetErr << endl; exit(1); }
			if(statBuff.st_mode & S_IFMT != S_IFREG){ cerr << redErr << "  Error reading set file: " << argv[i] << endl << resetErr; exit(1); } // only map reall files
			map_ids(argv[i], mapID);
		}
		return 0;
	}
	if(listChoice){														// display listFile
		cout << "Gene ontology version " << goDate << endl;
		if(backgroundList.size() != 0){
			cout << green << "Available backgrounds:" << endl << reset;
			for(int i=0; i<backgroundList.size(); i++){
				int entries = get_file_lines(backgroundList[i]);
				vector<string> types = get_file_types(backgroundList[i]);
				string typeString = boost::algorithm::join(types, ", ");
				string label = get_file_label(backgroundList[i]);
				cout << boost::format("#%s%2d%s:%10d entries%35s%s%20s%s\t%s\n") % cyan % (i+1) % reset % entries % typeString % cyan % label % reset % backgroundList[i] ;
			}
		}
		if(slimList.size() != 0){
			cout << green << "Available GO slims:" << endl << reset;
			for(int i=0; i<slimList.size(); i++){
				int entries = get_gofile_enties(slimList[i]);
				string definition = get_file_definition(slimList[i]);
				string label = get_file_label(slimList[i]);
				cout << boost::format("#%s%2d%s:%10d entries%35s%s%20s%s\t%s\n") % cyan % (i+1) % reset % entries % definition % cyan % label % reset % slimList[i]; 
			}
		}
		return 0;
	}
	if(updateChoice){													// update data
		string tmpBack1name, tmpBack2name, tmpIDname, outdatedBackName;
		for(int i=0; i<backgroundList.size(); i++){
			cout << "Updating file:\t" << backgroundList[i] << endl;
			tmpBack1name = 		backgroundList[i] + ".1.tmp";
			tmpBack2name = 		backgroundList[i] + ".2.tmp";
			tmpIDname = 		backgroundList[i] + ".IDs.tmp";
			outdatedBackName = 	backgroundList[i] + ".obs";
			if(boost::regex_search(backgroundList[i], boost::regex("[0-9]{4}-[0-9]{2}-[0-9]{2}"))){
				newBackName = boost::regex_replace(backgroundList[i], boost::regex("[0-9]{4}-[0-9]{2}-[0-9]{2}"), current_date()); // update date in filename
			}
			else{
				newBackName = backgroundList[i] + "_" + current_date() + ".bkg";	// if no date previously, add it
			}

			vector<string> types = get_file_types(backgroundList[i]);
			if(types.size() == 0){ cout << "Ommitting custom background file " << backgroundList[i] << endl; continue; }
			string label = get_file_label(backgroundList[i]);
			create_ID_file(backgroundList[i], tmpIDname);	// make ID file for use in new background
			for(int j=0; j<types.size(); j++){
				create_background(tmpIDname, tmpBack2name, "", types[j]);
				if(j==0){ copy_file(tmpBack2name, tmpBack1name); continue; }
				merge_backgrounds(tmpBack1name, tmpBack2name, 'U', tmpBack1name);
			}

			if(rename(backgroundList[i].c_str(), outdatedBackName.c_str()) == 0){		// if a file is renamed it is removed from the listfile at next execution
				cout << "Renamed file "	<< backgroundList[i] << " to " << outdatedBackName << endl;
			} else {
				cerr << redErr << "  Could not write background file " << outdatedBackName << endl;
				exit(1);
			}
			if(rename(tmpBack1name.c_str(), newBackName.c_str()) == 0){
				cout << "Updated file "	<< backgroundList[i] << " to " << newBackName << endl;
			} else {
				cerr << redErr << "  Could not write background file " << newBackName << endl;
				exit(1);
			}

			add_file_label(newBackName, label);
			backListFile.open(backListFileName.c_str(), std::fstream::out | std::fstream::app);		// add to background list file
			fullPath = realpath(newBackName.c_str(), NULL);
			backListFile << newBackName << endl;
			backListFile.close();
			remove(tmpIDname.c_str());
			remove(tmpBack1name.c_str());
			remove(tmpBack2name.c_str());
		}
		get_db();
		return 0;
	}
	if(corr == 'T' && (sortIndex == 1 || sortIndex == 3 || sortIndex == 5)){ // if sorted for Benjamini-Hochberg test, unset permutation test and set cutoff so all significant values are displayed
		randNum = 0;
		definedCutoff = 0;
	}
	////////////////// mandatory arguments
	if(!pipeTrue && argc - optind < 2)	{ cerr << redErr <<  "  Mandatory argument missing\n" << resetErr; return 1; }				// check for missing arguments
	if(pipeTrue && argc - optind < 1)	{ cerr << redErr <<  "  Mandatory argument missing\n" << resetErr; return 1; }				// check for missing arguments
	nameBackgr = argv[optind++];											// get background string
	for(int i=0; i<backgroundLabList.size(); i++){									// check for background specified as label
		if(backgroundLabList[i] == nameBackgr){
			nameBackgr = backgroundList[i];
			break;
		}
	}
	if(backNum = atoi(nameBackgr.c_str())){										// check for background specified as number
		if(backNum < 0 || backNum > backgroundList.size()){ cerr << redErr << "  No background file No. " << backNum << endl << resetErr; return 2; }
		nameBackgr = backgroundList[backNum-1];
		cerr << redErr << "Using background file: " << nameBackgr << endl << resetErr;
	}
	ifstream backgroundFile(nameBackgr.c_str());
	if(! backgroundFile){ cerr << redErr << "  Error reading background file: " << nameBackgr << endl << resetErr; return 2; }			// check background file
	nameBackgr = nameBackgr.substr(nameBackgr.find_last_of("/") + 1);							// make short background name string for output
	if(nameSlim != ""){													// check slim if specified
		for(int i=0; i<slimLabList.size(); i++){									// check for background specified as label
			if(slimLabList[i] == nameSlim){
				nameSlim = slimList[i];
				cerr << redErr << "Using slim file: " << nameSlim << endl << resetErr;
				break;
			}
		}
		if(slimNum = atoi(nameSlim.c_str())){										// check for slim specified as number
			if(slimNum < 0 || slimNum > slimList.size()){ cerr << redErr << "  No GO slim file No. " << slimNum << endl << resetErr; return 2; }
			nameSlim = slimList[slimNum-1];
			cerr << redErr << "Using slim file: " << nameSlim << endl << resetErr;
		}
		ifstream slimFile(nameSlim.c_str());
		if(! slimFile){ cerr << redErr << "  Error reading slim file: " << nameSlim << endl << resetErr; return 2; }				// check slim file
		slimFile.close();
		nameSlimShort = nameSlim.substr(nameSlim.find_last_of("/") + 1);						// make short slim name string for output
		oboFileName = nameSlim; 
	}	
	for(int i=optind; i<argc; i++){												// check set files
		struct stat statBuff;
		if(stat(argv[i], &statBuff) == -1){ cerr << redErr << "  Can't stat file: " << argv[i] << resetErr << endl; exit(1); }
		if(statBuff.st_mode & S_IFMT != S_IFREG){ cerr << redErr << "  Error reading set file: " << argv[i] << endl << resetErr; exit(1); } // only add reall files
		fileNames.push_back(argv[i]);
	}
	setNum = fileNames.size();
	if(pipeTrue){
		fileNames.push_back("STANDARD_INPUT");
		setNum = 1;
	}
	if(!pipeTrue && fileChoice && outFileName == ""){											// make output name if not specified
		outFileName = common_prefix(fileNames);
		if(outFileName[outFileName.size()-1] == '.'){ 	outFileName = outFileName + "enr"; }		// make shure there is no double points
		else{				outFileName = outFileName + ".enr"; }
	}
	if(pipeTrue && fileChoice && outFileName == ""){											// make output name if not specified from date if pipe
		string date = current_date();
		outFileName = date + ".enr";
	}
	if(defaultFileTrue){ cerr << redErr << "Using custom defaults." << endl << resetErr; }							// informing user aubout default setting
	////////////////// defining working variables
	bool found, obsolete;
	int h,i,j,k,l,m, ID, t00, t01, t10, t11, geneCountSet, unmatchedSet, function, functionIndex, functionBuffer, GOindex, aberIndex, displayNum, filterIDindex, filterStringIndex;
	int totalFunctionCount = 0;
	int geneCountAll = 0;
	int unmatchedAll = 0;
	int maxFun = 0;
	int countAbbIDs = -1; 			// aberrant 'GOIDs' start at -1, normal GOIDS at 1 
	double out, cutoff, cutoffSum, p;
	double minTot = 1.0;
	string identifier, IDstring, name, modString, explString, token, nameProt, functionString, nameSet, aberFun;
	vector<int> backRefIndices;
	vector<int> countAll;
	vector<int> countSet;
	vector<int> indicesTrial; 
	vector<int> backRefGOIDs;
	vector<int> GOIDs;
	vector<int> functionsBuffer;
	vector<int> functions;
	vector<double> logFacs (1, 0);		// logFacs[0]=0
	vector<char> GOnamespace;
	vector<string> GOnames;
	vector<string> GOexpls;
	vector<string> tokens;
	vector<string> proteinNames;
	vector<string> wrongIDs;
	vector<string> aberFunBuffer;
	vector< vector<int> > functionsAll;
	vector< vector<int> > affiliationSet;
	vector< pair< int, vector<double> > > enrichments;
	vector< pair< int, vector<double> > > enrichmentsDisplay;
	map<int,map<int,map<int,map<int,double> > > > pBuffFisher;		// container for allready generated p-values
	map<int,map<int,map<double,double> > > pBuffBinomial;			// container for allready generated binomial pValues
	ostringstream outStream;
	srand ( unsigned ( std::time(0) ) );
	////////////////// reading in ontology database
	i=-1;															// each ID adds 1, beginning should be 0
	mappingFile.open(oboFileName.c_str(), std::fstream::in);
	while(getline(mappingFile, line)){											// read in the GO mapping file
		istringstream lineStream(line);
		if(lineStream.str() == ""){ continue; }
		lineStream >> identifier;
		if("id:" == identifier){
			lineStream >> IDstring;
			ID = atoi(IDstring.substr(3, 7).c_str());								// if GOID is GO:0098501 we want 98501
			found = false;
		}
		else if("name:" == identifier){
			obsolete = (lineStream.str().substr(6, 8) == "obsolete");
			name = lineStream.str().substr(6);
		}
		else if("namespace:" == identifier){
			lineStream >> modString;
			if(obsolete && ! obsoleteChoice){ continue; }
			i++;
			GOIDs.push_back(ID);
			GOnames.push_back(name);
			GOnamespace.push_back(modString[0]);
			found = true;
		}
		else if(backRefChoice && found && "is_a:" == identifier){
			lineStream >> IDstring;
			ID = atoi(IDstring.substr(3, 7).c_str());								// if GOID is GO:0098501 we want 98501
			backRefGOIDs.push_back(ID);
			backRefIndices.push_back(i);										// at "namespace" iterated already
		}
		else if(printExpl && found && "def:" == identifier){
			explString = lineStream.str().substr(5);
			GOexpls.push_back(explString);
		}
		else if("[Typedef]" == identifier){
			break;
		}
	}
	mappingFile.close();
	////////////////// reading in background file
	while(getline(backgroundFile, line)){
		if(line[0] == '#'){ continue; }						// ignore comments
		indicesTrial.push_back(geneCountAll++);					// fill with gene indices
		tokens.clear(); functionsBuffer.clear(); aberFunBuffer.clear();
		istringstream lineStream(line);
		while(getline(lineStream, token, '\t')){ tokens.push_back(token); }	// read all tab-delimited into vector 
		nameProt = tokens[0];
		proteinNames.push_back(nameProt);
		vector<int> functionColumn; 
		for(int h=1; h<tokens.size(); h++){
			functionString = tokens[h];
			if(boost::regex_match(functionString, boost::regex("(GO:[0-9]{7})"))){						// if function is no GOID, save as aberrant function
				function = atoi(functionString.substr(3, 7).c_str());							// if GOID is GO:0098501 we want 98501
 				if(find(functionsBuffer.begin(), functionsBuffer.end(), function) != functionsBuffer.end()){ continue; }// guard against multiple functions for same protein
				else{ functionsBuffer.push_back(function); }
			}
			else{
 				if(find(aberFunBuffer.begin(), aberFunBuffer.end(), functionString) != aberFunBuffer.end()){ continue; }	// guard against multiple aberrant functions for same protein
				else{ aberFunBuffer.push_back(functionString); }
			}
		}
		for(int i=0; i<aberFunBuffer.size(); i++){				// handle aberrant functions => functions
			aberFun = aberFunBuffer[i];	
			aberIndex = find(GOnames.begin(), GOnames.end(), aberFun) - GOnames.begin();
			if(aberIndex >= GOnames.size()){				// if not found yet, append
				function = countAbbIDs;
				GOIDs.push_back( countAbbIDs-- );			// aberrant 'GOIDs' count backwards from -1
				GOnames.push_back(aberFun);
				GOnamespace.push_back('a');
				GOexpls.push_back("aberrant function");
			}
			else{	function = GOIDs[aberIndex]; }
			functionsBuffer.push_back(function);				// append rectified function
		}
		for(int i=0; i<functionsBuffer.size(); i++){
			function = functionsBuffer[i];	
			if(backRefChoice){				// backreference handling
				for(int j=0; j<GOIDs.size(); j++){
					if(function == GOIDs[j]){
						for(int k=0; k<backRefIndices.size(); k++){
							if(backRefIndices[k] == j){
								functionBuffer = backRefGOIDs[k];
								if(find(functionsBuffer.begin(), functionsBuffer.end(), functionBuffer) == functionsBuffer.end()){
									functionsBuffer.push_back( functionBuffer ); 
								}
							}
						}
						break;
					}
				}
			}
			functionIndex = -1;
			for(GOindex=0; GOindex < functions.size(); GOindex++){		// if function id is allready found
				if(function == functions[GOindex]){
					countAll[GOindex] ++;				// count it
					functionIndex = GOindex;
				}
			}
			if(functionIndex < 0){						// if funciton id is not found allredy
				functions.push_back(function);
				countAll.push_back(1);					// make new counter
				countSet.push_back(0);
				functionIndex = functions.size() - 1;
				if(printProt){ affiliationSet.push_back( vector<int>() ); }	// for each function prepare vector of indices for proteins linked to that function
			}
			functionColumn.push_back(functionIndex);
			totalFunctionCount++;
		}
		functionsAll.push_back(functionColumn);
		if(functionsBuffer.size()==0){unmatchedAll++;}
	}
	if(totalFunctionCount == 0){
		cerr << redErr << "  Empty background file: " << nameBackgr << endl << resetErr;
		return 2;
	}
	if(countAbbIDs == -1){ modusA = false; }		// no aberrant output needed if none is found

	int ignoreFunCount = 0;
	for(int i=0; i<functions.size(); i++){ if(countAll[i] <= 1){ ignoreFunCount++; }}	// count functions with only 1 occurrence, ignore in bonferroni correction of p-values 

	////////////////// reading all possible sets
	istream *in;
	for(int h=0; h<setNum; h++){
		ifstream comparisonFile(fileNames[h].c_str());
		if(! pipeTrue){
			in = &comparisonFile;
			nameSet = fileNames[h];
			nameSet = nameSet.substr(nameSet.find_last_of("/") + 1);							// truncate full filename, get rid of path
		}
		else{
			in = &cin;
			nameSet = "STANDARD_IN_SET_" + boost::lexical_cast<std::string>(h+1);
		}
		wrongIDs.clear(); enrichmentsDisplay.clear(); 
		if(printProt){ for(int i=0; i<affiliationSet.size(); i++){ affiliationSet[i].clear(); } }				// delete protein indices previously linked to function i
		
		fill(countSet.begin(), countSet.end(), 0);
		geneCountSet = 0; unmatchedSet = 0;
		while(*in >> nameProt){
			found = false;
			if(nameProt == "END"){
				setNum++;
				break;
			}
			for(int i=0; i<proteinNames.size(); i++){
				if(proteinNames[i] == nameProt){
					if(functionsAll[i].size() == 0){ unmatchedSet++;}
					for(int j=0; j<functionsAll[i].size(); j++){
						countSet[functionsAll[i][j]]++;
						if(printProt){ affiliationSet[ functionsAll[i][j] ].push_back(i); }			// fill affiliationSet with number (i) corresponding to protein
					}
					geneCountSet++;
					found = true;
					break;
				}
			}
			if(! found){ wrongIDs.push_back(nameProt); }
		}
		if(wrongIDs.size()){
			cerr << redErr <<  endl << "  Set file " << nameSet << " contains " << wrongIDs.size() << " IDs without matching background IDs:" << endl << resetErr;
			for(int i=0; i<wrongIDs.size(); i++){
				cerr << redErr << wrongIDs[i] << " " << resetErr;
			}
			cerr << redErr << endl << endl << resetErr;
		}
		if(! geneCountSet){
			cerr << redErr << "  Ignoring file without IDs matching background IDs: " << nameSet << endl << resetErr;
			continue;
		}
		////////// calculate enrichments for all columns
		cutoffSum = 0;
		for(int i=0; i<randNum+1; i++){ 						// calculate enrichments at least once, more often for permutation test
			bool testFisherTrue = fisherTrue; bool testBinomialTrue = binomialTrue; bool  testHyperTrue = hyperTrue; bool testFoldTrue = foldTrue; char testCorr = corr;	// lokal changable variables 
			enrichments.clear();
			if(i == 1){								// when permutation test
				switch(sortIndex){						// only calculate what you need 
					case 0:		testFisherTrue = true; testBinomialTrue = false; testHyperTrue = false; testFoldTrue = false; testCorr = ' '; 	break;
					case 1:		testFisherTrue = true; testBinomialTrue = false; testHyperTrue = false; testFoldTrue = false;		 	break;
					case 2:		testFisherTrue = false; testBinomialTrue = true; testHyperTrue = false; testFoldTrue = false; testCorr = ' '; 	break;
					case 3:		testFisherTrue = false; testBinomialTrue = true; testHyperTrue = false; testFoldTrue = false;		 	break;
					case 4:		testFisherTrue = false; testBinomialTrue = false; testHyperTrue = true; testFoldTrue = false; testCorr = ' '; 	break;
					case 5:		testFisherTrue = false; testBinomialTrue = false; testHyperTrue = true; testFoldTrue = false;		 	break;
					case 6:		testFisherTrue = false; testBinomialTrue = false; testHyperTrue = false; testFoldTrue = true; testCorr = ' '; 	break;
					case 7:		testFisherTrue = false; testBinomialTrue = false; testHyperTrue = false; testFoldTrue = false; testCorr = ' '; 	break;
					case 8:		testFisherTrue = false; testBinomialTrue = false; testHyperTrue = false; testFoldTrue = false; testCorr = ' '; 	break;

				}
				random_shuffle(indicesTrial.begin(), indicesTrial.end());	// shuffle all protein indices
				fill(countSet.begin(), countSet.end(), 0);			// reset all function counters
				if(sortReverse){ cutoff=0; }
				if(! sortReverse){ cutoff=999999999; }
				for(int j=0; j<geneCountSet; j++){				// get as many protein indices as in current set
					k = indicesTrial[j];
					for(int l=0; l<functionsAll[k].size(); l++){
						countSet[functionsAll[k][l]]++;			// fill function counter according to function of gene index
					}
				}
			}
			for(int j=0; j<functions.size(); j++){
				if(countSet[j]+countAll[j] <= 1){ continue; }			// enrichment for functions with only one entry is a logic impossibility
				double fisherOut = -1;
				double binomialOut = -1;
				double hyperOut = -1;
				double foldOut = -1;
	 			if(testFisherTrue){
					t00 = countAll[j];
					t10 = countSet[j];
					t01 = geneCountAll - countAll[j];
					t11 = geneCountSet - countSet[j];
		      			if(!reverseChoice){ fisherOut = fisher_test_less(logFacs, pBuffFisher, t00, t01, t10, t11); }
		      			else              { fisherOut = fisher_test_greater(logFacs, pBuffFisher, t00, t01, t10, t11); }
	 			}
	 			if(testBinomialTrue){
					double p = countAll[j] / double(geneCountAll);
					int n = geneCountSet;
					int k = countSet[j];
	 				if(!reverseChoice){ binomialOut = binomial_test_less(logFacs, pBuffBinomial, n, k, p); }
	 				else		  { binomialOut = binomial_test_greater(logFacs, pBuffBinomial, n, k, p); }
	 			}
				if(testHyperTrue){
					t00 = countAll[j];
					t10 = countSet[j];
					t01 = geneCountAll - countAll[j];
					t11 = geneCountSet - countSet[j];
	 				if(!reverseChoice){ hyperOut = hypergeometric_prob(logFacs, t00, t01, t10, t11); }
	 				else		  { hyperOut = 1 - hypergeometric_prob(logFacs, t00, t01, t10, t11); }
				}
	 			if(testFoldTrue){
	 				if(!reverseChoice){ foldOut = (countSet[j]/double(geneCountSet)) / (countAll[j]/double(geneCountAll)); }
	 				else 		  { foldOut = (countAll[j]/double(geneCountAll)) / (countSet[j]/double(geneCountSet)); }
	 			}
	 			vector <double> enrichment (9, -1.0);		// initalize all values with -1
	 			enrichment[0] = fisherOut;
	 			enrichment[1] = fisherOut;	// for testCorr
	 			enrichment[2] = binomialOut;
	 			enrichment[3] = binomialOut;	// for testCorr
	 			enrichment[4] = hyperOut;
	 			enrichment[5] = hyperOut;	// for testCorr
	 			enrichment[6] = foldOut;
	 			enrichment[7] = countAll[j];
	 			enrichment[8] = countSet[j];
				enrichments.push_back(make_pair(j, enrichment));
			}
			///////////////// multiple hypothesis testCorrection
			/// Bonferroni (standard)
			if(testCorr == 'B' || testCorr == 'b' || testCorr == 'C' || testCorr == 'c' || testCorr == 'A' || testCorr == 'a'){
				for(int j=0; j<enrichments.size(); j++){
		 			enrichments[j].second[1] *= (functions.size()-ignoreFunCount);
		 			enrichments[j].second[3] *= (functions.size()-ignoreFunCount);
		 			enrichments[j].second[5] *= (functions.size()-ignoreFunCount);
				}
			}
			if(testCorr == 'C' || testCorr == 'c' || testCorr == 'A' || testCorr == 'a'){
				if(fisherCorr){
					/// FDR testCorrection
					sort(enrichments.begin(), enrichments.end(), fisherCompare);
					for(j=0; j<enrichments.size(); j++){ enrichments[j].second[1] /= (j+1); }
					/// FDR adjustment
					if(testCorr == 'A' || testCorr == 'a'){
						vector <double> qbuff (enrichments.size(), -1);
						for(j=0; j<enrichments.size(); j++){
							double minQ = enrichments[j].second[1];
							for(k=j+1; k<enrichments.size(); k++){ if(minQ > enrichments[k].second[1]){ minQ = enrichments[k].second[1]; }}
							qbuff[j] = minQ;
						}
						for(j=0; j<enrichments.size(); j++){ enrichments[j].second[1] = qbuff[j]; }
					}
				}
				if(binomialCorr){
					/// FDR testCorrection
					sort(enrichments.begin(), enrichments.end(), binomialCompare);
					for(j=0; j<enrichments.size(); j++){ enrichments[j].second[3] /= (j+1); }
					/// FDR adjustment
					if(testCorr == 'A' || testCorr == 'a'){
						vector <double> qbuff (enrichments.size(), -1);
						for(j=0; j<enrichments.size(); j++){
							double minQ = enrichments[j].second[3];
							for(k=j+1; k<enrichments.size(); k++){ if(minQ > enrichments[k].second[3]){ minQ = enrichments[k].second[3]; }}
							qbuff[j] = minQ;
						}
						for(j=0; j<enrichments.size(); j++){ enrichments[j].second[3] = qbuff[j]; }
					}
				}
				if(hyperCorr){
					/// FDR testCorrection
					sort(enrichments.begin(), enrichments.end(), hyperCompare);
					for(j=0; j<enrichments.size(); j++){ enrichments[j].second[5] /= (j+1); }
					/// FDR adjustment
					if(testCorr == 'A' || testCorr == 'a'){
						vector <double> qbuff (enrichments.size(), -1);
						for(j=0; j<enrichments.size(); j++){
							double minQ = enrichments[j].second[5];
							for(k=j+1; k<enrichments.size(); k++){ if(minQ > enrichments[k].second[5]){ minQ = enrichments[k].second[5]; }}
							qbuff[j] = minQ;
						}
						for(j=0; j<enrichments.size(); j++){ enrichments[j].second[5] = qbuff[j]; }
					}
				}
			}
			/// Benjamini–Hochberg test
			if(testCorr == 'T' || testCorr == 't'){
				if(fisherCorr){
					sort(enrichments.begin(), enrichments.end(), fisherCompare);
					for(k=0; k<enrichments.size(); k++){ if(enrichments[k].second[0] > fdr*(k+1)/double(enrichments.size())){ break; } } // k+1 to start at 1 	enrichment.size() allready is reduced by ignoreFunCount
					for(j=0; j<=k; j++){ enrichments[j].second[1] = 0; }
					for(j=k+1; j<enrichments.size(); j++){ enrichments[j].second[1] = 1; }
				}
				if(binomialCorr){
					sort(enrichments.begin(), enrichments.end(), binomialCompare);
					for(k=0; k<enrichments.size(); k++){ if(enrichments[k].second[2] > fdr*(k+1)/double(enrichments.size())){ break; } }
					for(j=0; j<=k; j++){ enrichments[j].second[3] = 0; }
					for(j=k+1; j<enrichments.size(); j++){ enrichments[j].second[3] = 1; }
				}
				if(hyperCorr){
					sort(enrichments.begin(), enrichments.end(), hyperCompare);
					for(k=0; k<enrichments.size(); k++){ if(enrichments[k].second[4] > fdr*(k+1)/double(enrichments.size())){ break; } }
					for(j=0; j<=k; j++){ enrichments[j].second[5] = 0; }
					for(j=k+1; j<enrichments.size(); j++){ enrichments[j].second[5] = 1; }
				}
			}
			///////////////// sort enrichments by first test column
			switch(sortIndex){
				case 0: sort(enrichments.begin(), enrichments.end(), fisherCompare); 	 	break;
				case 1: sort(enrichments.begin(), enrichments.end(), fisherCorrCompare); 	break;
				case 2: sort(enrichments.begin(), enrichments.end(), binomialCompare); 		break;
				case 3: sort(enrichments.begin(), enrichments.end(), binomialCorrCompare); 	break;
				case 4: sort(enrichments.begin(), enrichments.end(), hyperCompare); 		break;
				case 5: sort(enrichments.begin(), enrichments.end(), hyperCorrCompare);		break;
				case 6: sort(enrichments.begin(), enrichments.end(), foldCompare);	 	break;
				case 7: sort(enrichments.begin(), enrichments.end(), backNumCompare);	 	break;
				case 8: sort(enrichments.begin(), enrichments.end(), setNumCompare);	 	break;
				case 9: sort(enrichments.begin(), enrichments.end(), bind(funNameCompare, placeholders::_1, placeholders::_2, functions, GOIDs, GOnames));	break;
				case 10: sort(enrichments.begin(), enrichments.end(), bind(expNameCompare, placeholders::_1, placeholders::_2, functions, GOIDs, GOexpls));	break;
			}
			///////////////// differentiate between calculting enrichments for display and permutation test
			if(i == 0){							// enrichments for display
				for(j=0; j<enrichments.size(); j++){			// copy enrichments
	 				vector <double> enrichment (9, -1.0);		// initalize all values with -1
					for(int k=0; k<enrichments[j].second.size(); k++){ enrichment[k] = enrichments[j].second[k]; }		// copy p-values
					enrichmentsDisplay.push_back(make_pair(enrichments[j].first, enrichment));				// copy function number
				}
			}else{
				cutoff = enrichments[0].second[sortIndex];
				cutoffSum += cutoff;
			}
		}
		if(randNum)			{ cutoff = cutoffSum / randNum; }
		else if(definedCutoff >= 0)	{ cutoff = definedCutoff; }		// definedCutoff is initialized as -1 so a cutoff of 0 can be used 
		else				{ cutoff = 1; }				// if no p-value bound is set, ignore it 	
		///////////////// output
		////// print maximum enrichment of first test column
		if(printMin){
			cout << enrichmentsDisplay[0].second[sortIndex] << sep.c_str() << cyan.c_str() << nameSet.c_str() << reset.c_str() << endl;
			continue;
		}
		////// header
		if(headerTrue){
			outStream << boost::format("# File: %s%s%s  Background: %s  Modus: %s  Columns: %s  Correction: %c  ") % cyan.c_str() % nameSet.c_str() % reset.c_str() % nameBackgr.c_str() % modi % displayColumns % corr;
			if(nameSlimShort != ""){ outStream << "GO-Slim: " << nameSlimShort << "  "; }
			if(reverseChoice){ outStream << "Reverse  "; }
			if(backRefChoice){ outStream << "BackRef  "; }
			outStream << endl;
			outStream << boost::format("# Proteins bkg/set: %d / %d  Functions used/ignored: %d / %d  Functionless Proteins bkg/set: %d / %d") % geneCountAll % geneCountSet % functions.size() % ignoreFunCount % unmatchedAll % unmatchedSet;
			if(randNum)			{ outStream << boost::format("  %s-cutoff: %.2e (%d Tries)\n") % sortChar % cutoff % randNum; } 
			else if(definedCutoff >= 0)	{ outStream << boost::format("  user-defined cutoff: %.2e\n") % definedCutoff; }
			else if(allTrue)		{ outStream << "  no cutoff\n"; }
			else				{ outStream << endl; }
			bool color = false;
			char colorChar;
			int printCol = 1;
			for(int k=0; k<columns.length(); k++){
				char column = columns[k];
				switch(column){
					case 'F': outStream << printCol << ": Fisher Exact Test";  break;
					case 'f': outStream << printCol << ": Fisher Exact Test Corrected";  break;
					case 'B': outStream << printCol << ": Binomial Test";  break;
					case 'b': outStream << printCol << ": Binomial Test Corrected";  break;
					case 'H': outStream << printCol << ": Hypergeometric Test";  break;
					case 'h': outStream << printCol << ": Hypergeometric Test Corrected";  break;
					case 'E': outStream << printCol << ": Fold-Enrichment";  break;
					case 'N': outStream << printCol << ": Function Number in Background";  break;
					case 'n': outStream << printCol << ": Function Number in Set";  break;
					case 'x': outStream << printCol << ": Name of Function";  break;
					case 'X': outStream << printCol << ": Explanation of Function";  break;
					case 'P': outStream << printCol << ": Protein ID(s)";  break;
					case 'G': outStream << printCol << ": Gene Ontology ID"; break;
					case '1':case'2':case'3':case'4':case'5': color = true; colorChar = column; printCol--; break;
					default: break;
				}
				if(! color){ outStream << boost::format("%s%.2e") % reset % sep.c_str(); }
				else{	switch(colorChar){	case '1': outStream << green; break;
								case '2': outStream << cyan; break;
								case '3': outStream << magenta; break;
								case '4': outStream << yellow; break;
								case '5': outStream << red; break; } }
				color = false;
				printCol++;
			}
			outStream << endl;
		} else {
			outStream << boost::format("# %s%s%s\n") % cyan.c_str() % nameSet.c_str() % reset.c_str();
		}
		////// body
		outStream << reset;
		if(filterIDtrue){ for(int k=0; k<filters.size(); k++){ filtersFound[k] = false; } }		// reset filter recognition
		if(filterStringTrue){ for(int k=0; k<filterStrings.size(); k++){ filterStringsFound[k] = false; } }		// reset filter recognition
		for(int i=0; i<4; i++){										// go through all modi
			if(headerTrue && modusF+modusC+modusP+modusA > 1){							// check if namespace headers needed
				outStream << magenta;									
				if(i==0 && modusF){ modus = 'm'; outStream << "# molecular function:\n"; }	
				if(i==1 && modusC){ modus = 'c'; outStream << "# cellular component:\n"; }
				if(i==2 && modusP){ modus = 'b'; outStream << "# biological process:\n"; }	// modus=='b' because 'biological process' ~ P
				if(i==3 && modusA){ modus = 'a'; outStream << "# aberrant functions:\n"; }
				outStream << reset;
			}
			modus = 'x';										// dummy modus
			if(i==0 && modusF){ modus = 'm'; }	
			if(i==1 && modusC){ modus = 'c'; }
			if(i==2 && modusP){ modus = 'b'; }							// modus=='b' because 'biological process' ~ P
			if(i==3 && modusA){ modus = 'a'; }
			displayNum = 0;
			for(int j=0; j<enrichmentsDisplay.size(); j++){
				///// preparation
				index 	= enrichmentsDisplay[j].first;
				p 	= enrichmentsDisplay[j].second[sortIndex];
				if(!allTrue &&  sortReverse && p < cutoff){ break; }
				if(!allTrue && !sortReverse && p > cutoff){ break; }
				GOindex = find(GOIDs.begin(), GOIDs.end(), functions[index]) - GOIDs.begin(); if(GOindex >= GOIDs.size()){ continue; }		// look up GOindex, if not found, continue
				if(modus != GOnamespace[GOindex]){ continue; }											// continue if wrong modus
				///// filters
				if(filterIDtrue){
					filterIDindex = find(filters.begin(), filters.end(), functions[index]) - filters.begin();
					if(filterIDindex == filters.size()){ continue; }									// if filter is present, only output filter GOIDs
					else{filtersFound[filterIDindex] = true; }										// remove found filters
				}
				if(filterStringTrue){
					found = false;
					for(int k=0; k<filterStrings.size(); k++){ if(boost::regex_search(GOnames[GOindex], boost::regex(filterStrings[k]))){ found = true; } }
					if(! found){ continue; }												// if filter is present, only output filter functions
					else{filterStringsFound[filterStringIndex] = true; }									// remove found filterStrings
				}
				///// columns
				bool color = false;
				char colorChar;
				int printCol = 1;
				for(int k=0; k<columns.length(); k++){
					char column = columns[k];
					switch(column){
						case 'F': outStream << boost::format("%.2e") % enrichmentsDisplay[j].second[0]; break;
						case 'f': outStream << boost::format("%.2e") % enrichmentsDisplay[j].second[1]; break;
						case 'B': outStream << boost::format("%.2e") % enrichmentsDisplay[j].second[2]; break;
						case 'b': outStream << boost::format("%.2e") % enrichmentsDisplay[j].second[3]; break;
						case 'H': outStream << boost::format("%.2e") % enrichmentsDisplay[j].second[4]; break;
						case 'h': outStream << boost::format("%.2e") % enrichmentsDisplay[j].second[5]; break;
						case 'E': outStream << boost::format("%.2f") % enrichmentsDisplay[j].second[6]; break;
						case 'N': outStream << boost::format("%5d")  % enrichmentsDisplay[j].second[7]; break;
						case 'n': outStream << boost::format("%5d")  % enrichmentsDisplay[j].second[8]; break;
						case 'x': outStream << boost::format("%s")   % GOnames[GOindex].c_str(); break;
						case 'X': if(GOexpls.size() > GOindex){ 						// guard against displaying non-existant explanations
							  outStream << boost::format("%s")   % GOexpls[GOindex].c_str();
							  }else{ outStream << "missing" << endl; } break;
						case 'P': for(int l=0; l<affiliationSet[index].size(); l++){ 
							  outStream << proteinNames[affiliationSet[index][l]] << ";"; } break;
						case 'G': if(GOIDs[GOindex] > 0){ 							// guard against displaying abberant IDs
						          outStream << boost::format("%sGO:%07d") % sep.c_str() % GOIDs[GOindex];
							  }else{ outStream << "missing" << endl; } break;
						case '1':case'2':case'3':case'4':case'5': color = true; colorChar = column; break;
						default: break;
					}
					if(! color){ outStream << boost::format("%s%.2e") % reset % sep.c_str(); }
					else{	switch(colorChar){	case '1': outStream << green; 	break;
									case '2': outStream << cyan; 	break;
									case '3': outStream << magenta; break;
									case '4': outStream << yellow; 	break;
									case '5': outStream << red; 	break; } }
					color = false;
				}
				///// column end
				outStream << endl;
				if(maxNum && ++displayNum >= maxNum){ break; }					// increments displayed functions and checks maxNum after every other chance to continue or break
			}
			if(! fileChoice && ! printMin){								// output immediately if no file specified
				cout << outStream.str();
				outStream.str("");
			}
			if(fileChoice && ! printMin){
				cerr << boost::format("%5d/%d\t%3.2f%%\r") % (h+1) % setNum % ((h+1)*double(100)/setNum);
			}
		}
		if(filterIDtrue){	for(int j=0; j<filters.size(); j++){ if( ! filtersFound[j]){ outStream << red << boost::format("# filter 'GO:%07d' not found\n") % filters[j] << reset; } } } // note empty filters
		if(filterStringTrue){	for(int j=0; j<filterStrings.size(); j++){ if( ! filterStringsFound[j]){ outStream << red << boost::format("# filter '%s' not found\n") % filterStrings[j] << reset; } } } // note empty filters
	}
	////////////////// output in file
	if(fileChoice){
		cout << "Generating File " << outFileName << endl;
		ofstream outFile(outFileName.c_str());
			if(! outFile){ cerr << redErr << "  Error writing output file: " << outFileName << endl << resetErr; exit(5); }
		outFile << outStream.str();
	}
}
