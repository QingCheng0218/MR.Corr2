#ifndef pdsoft_hpp
#define pdsoft_hpp

#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <string.h>

using namespace std;
using namespace arma;

struct pdsoftObj{
	mat omega;
	mat sigma;
	mat theta;
	mat theta_inv;
};

pdsoftObj  pdsoft(mat s, mat lam, double tau = 1e-8, string init = "soft",
	bool standard = 1, double tolin = 1e-8, double tolout = 1e-8, int maxitin = 1e4,
	int maxitout = 1e3, bool quiet = 1);

#endif 
