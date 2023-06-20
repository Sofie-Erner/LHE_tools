#include <iostream>
#include <fstream>
#include <string>
#include <list>

#include "LHEF/LHEF.hpp"

// ***** Functions *****
// ***** Add to Histogram
std::vector<int> into_hist( std::vector<int>& hist, int spin );

// ***** Write Vector to File
template <typename T1 > void vec_to_file(std::vector<T1>& vec, std::string ofile, std::string vec_name );
void vec_hist_to_file( std::vector<double>& vec_s, std::vector<int>& vec_h, std::string ofile, std::string vec_name_s, std::string vec_name_h, int n_events );

// ***** Vector stuff
std::vector<double> linspace(double start_in, double end_in, int num_in);
//template <typename T2> 
void print_vector(std::vector<double>& vec);
void printNestedVector(std::vector< std::vector<int> >& nested_vector);
void saveNestedVector(std::vector< std::vector<int> >& nested_vector, std::string ofile, std::string vec_name);
void saveNestedVector(std::vector< std::vector<double> >& nested_vector, std::string ofile, std::string vec_name);

// Sort cos numbers into vectors
std::vector< std::vector<int> > sort_cos(int cos_n, float cos_min, float cos_step, double cos_val, std::vector< std::vector<int> > cos_vec, int sAi );

// ***** Basic Kinematics (Based on Mathematica file)
// *** Functions for event.PUP[i] (Fourvector of the ith particle)
double ptOf(const std::vector<double>& vec);
double enOf(const std::vector<double>& vec);
double thetaOf(const std::vector<double>& vec);
double phiOf(const std::vector<double>& vec);
double phifn(const std::vector<double>& vec);
double rapOf(const std::vector<double>& vec);
double etaOf(const std::vector<double>& vec);
double FourLengthSq(const std::vector<double>& vec);

std::vector<double> ptVec(const std::vector<double>& vec);
std::vector<double> ThreeVectorFrom(const std::vector<double>& vec);
std::vector<double> addVecs(const std::vector<double>& vec1, const std::vector<double>& vec2);
