/*
 * moduli.h
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

#ifndef MODULI_H
#define MODULI_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;

class Moduli{
	
public:
	//getter and setters
	double J0(){return _J0;}
	void J0(double J0_){_J0=J0_;}
	
	double eta(){return _eta;}
	void eta(double eta_){_eta=eta_;}

	int dpoints(){return _dpoints;}
	void dpoints(int dpoints_){_dpoints=dpoints_;}

	double maxfreq(){return _maxfreq;}
	void maxfreq(double maxfreq_){_maxfreq=maxfreq_;}

	double minfreq(){return _minfreq;}
	void minfreq(double minfreq_){_minfreq=minfreq_;}

	int N(){return _N;}
	void N(int N_){_N=N_;}

	void run();

	void getgp(boost::numeric::ublas::vector<double>&);
	void getgpp(boost::numeric::ublas::vector<double>&);

	void setJ(boost::numeric::ublas::vector<double>&);
	void sett(boost::numeric::ublas::vector<double>&);

private:
	void frequencysweep();
	int _dpoints, _N;
	double _J0, _eta, _maxfreq, _minfreq;
	double smin, smax, step, w, ReSum, ImSum, x, y; 

	boost::numeric::ublas::vector<double> J, t;
	boost::numeric::ublas::vector<double> Gp, Gpp;
};


#endif


