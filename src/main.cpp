/*
 * main.cpp
 * 
 * 
 * Modified from:
 * 
 * http://www.pcf.leeds.ac.uk/research/highlight/view/4
 *
 * Polymer & Complex Fluids Group
 * School of Physics & Astronomy
 * University of Leeds
 *
 * This program converts compliance measurements into storage and loss
 * moduli, using the method described by 
 * R M L Evans, Manlio Tassieri, Dietmar Auhl and Thomas A Waigh in 
 * Phys. Rev. E 80, 012501 (2009).
 *
 * Coded by Chirag Kalelkar, Complex Fluids and Polymer Engineering Group,
 * National Chemical Laboratory, Pune, India,
 * Modified from code by David Pearce, University of Manchester, U.K.
 * Edited by R M L Evans, University of Leeds, U.K.
 * This code by Jonathan Bramble, University of Leeds, UK
 * 
 * Moduli is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * frap-tool is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "../include/moduli.h"
#include <iostream>     // cout, endl
#include <fstream>      // fstream
#include <boost/tokenizer.hpp>

using namespace std; 
using namespace boost;

int main(int argc, char* argv[]) {

	int c, index;
	char* datafile; 
	char* results;
	char* appname = argv[0];

	while((c = getopt(argc, argv, "d:r:")) != -1){			//taken from gnu example
		switch(c){
		case 'd':
			datafile = optarg;
			break;
		case 'r':
			results = optarg;
			break;
		default:
			fprintf(stderr, "usage: %s [-d csvfilename] [-r resultscsvfilename] \n", appname);
			exit(-1);
		}
	
	}
	
	char_separator<char> sep(",");				//csv files are separated with ,  can replace for otherfile types
	typedef tokenizer< char_separator<char> > Tokenizer;	//define the tokenizer, abit OTT here
	std::vector< string > vec;				//a standard vector of strings
	string line;						

	char * pEnd;						// for convertions
	int rc = 0;						// counters
	int ck = 0;

	ifstream infile;					// initilise a stream
	infile.open(datafile);					// here we read the file once to get the number of data points
    	if (!infile.is_open()) {
		cout << "File not opened" << endl;
		return 1;
	}
	while (getline(infile,line)){rc++;}
	infile.close();

	cout << "Number of data points: " << rc << endl;
	
	boost::numeric::ublas::vector<double> J(rc);		// initialise the data storage
	boost::numeric::ublas::vector<double> t(rc);
	
	infile.open(datafile);				// reopen the file

    	while (getline(infile,line) && ck < (rc-1))
    	{
        	Tokenizer tok(line, sep);
        	vec.assign(tok.begin(),tok.end());		
		t(ck) = strtod(vec[0].c_str(),&pEnd);		// convert the strings into data and add to arrays
		J(ck) = strtod(vec[1].c_str(),&pEnd);
		ck++;	
    	}

	infile.close();						// close the file

	Moduli *moduli; 					// initialise the moduli class
	moduli = new Moduli();

	int pts = 3000;
		
	moduli->J0(0.1);					// set params
	moduli->eta(1/0.00003);
	moduli->dpoints(pts);
	moduli->maxfreq(100);
	moduli->minfreq(1/t(ck-1));				// the -1 is there to stop it reading passed the end of the file
	moduli->N(rc);
	
	moduli->setJ(J);	
	moduli->sett(t);

	boost::numeric::ublas::vector<double> gp(pts);		//allocate data arrays for solution (dpoints)
	boost::numeric::ublas::vector<double> gpp(pts);
	boost::numeric::ublas::vector<double> ww(pts);

	moduli->run();						// run calculations

	moduli->getgp(gp);					// retrieve results
	moduli->getgpp(gpp);
	moduli->getgpp(ww);

	ofstream outfile(results);

	if (!outfile.is_open()) return 1;
	for(int i = 0; i < pts; i ++){
		outfile << ww[i] << "," << gp[i] << "," << gpp[i] << endl;
	}
	outfile.close();
	

}
