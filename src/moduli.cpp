/*
 * moduli.cpp
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

#include "../include/moduli.h"

void Moduli::frequencysweep(){
	smin = log(_minfreq);  //get from end vector
	smax = log(_maxfreq);
	step = (smax-smin)/(_dpoints-1);
	ReSum = 0;
	ImSum = 0;

	//std::cout << _eta << std::endl;
	std::cout << smin << std::endl;
	std::cout << _minfreq << std::endl;
	std::cout << smax << std::endl;

	for(int i=0; i< _dpoints; ++ i)
	{
		w = pow(10,(smin + i*step));	// logspace
		ww(i) = w;

		//std::cout << w << std::endl;

		for(int k=2; k<(_N-1); ++k){	// messed with length here
			 ReSum += (cos(w*t(k-1))-cos(w*t(k))) * (J(k)-J(k-1))/(t(k)-t(k-1));
      			 ImSum += (-sin(w*t(k-1))+sin(w*t(k))) * (J(k)-J(k-1))/(t(k)-t(k-1));
		}

		x = ((1.-cos(w*t(1)))*(J(1)-_J0)/t(1)) + (cos(w*t(_N-1))/_eta) + ReSum;   // check indexes of J
    		y = (w*_J0) - (sin(w*t(1))*(J(1)-_J0)/t(1)) - (sin(w*t(_N-1))/_eta) + ImSum;  // check first sign of sin
		
		Gp(i) = w*y/(x*x+y*y); //complex conj multiply
		Gpp(i) = w*x/(x*x+y*y);

		//std::cout << Gp(i) << std::endl;
	}

}

void Moduli::run(){
	std::cout << "Running calculations..." << std::endl;   // add validations here

	Gp = boost::numeric::ublas::vector<double>(_dpoints);
	Gpp = boost::numeric::ublas::vector<double>(_dpoints);
	ww = boost::numeric::ublas::vector<double>(_dpoints);

	frequencysweep();
	std::cout << "...Complete" << std::endl;
}

void Moduli::getgp(boost::numeric::ublas::vector<double>& ret_data){
	ret_data=Gp;
}

void Moduli::getgpp(boost::numeric::ublas::vector<double>& ret_data){
	ret_data=Gpp;
}

void Moduli::getww(boost::numeric::ublas::vector<double>& ret_data){
	ret_data = ww;	
}

void Moduli::setJ(boost::numeric::ublas::vector<double>& _J){
	J=_J;
}

void Moduli::sett(boost::numeric::ublas::vector<double>& _t){
	t=_t;
}




