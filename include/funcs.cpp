#include <iostream>
#include <fstream>
#include <string>
#include <list>

#include "LHEF/LHEF.hpp"
#include "header.h"

// ***** Functions *****
// ***** Add to Histogram
std::vector<int> into_hist( std::vector<int>& hist, int spin ){
  if ( spin == 1 ){
    hist[1] += 1;
  }
  else if ( spin == -1 ){
    hist[0] += 1;
  }
  return hist;
}

// ***** Write Vector to File
template <typename T1 > void vec_to_file(std::vector<T1>& vec, std::string ofile, std::string vec_name){
  std::ofstream out_file;
  out_file.open(ofile, std::fstream::app); //append to file

  if ( vec_name.size() != 0 ){
    out_file << vec_name << ": ";
  }

  for ( auto s: vec ){ // loop over vecor entries
    out_file << s << " ";
  }
  out_file << "\n";
  out_file.close();
}

void vec_hist_to_file( std::vector<double>& vec_s, std::vector<int>& vec_h, std::string ofile, std::string vec_name_s, std::string vec_name_h, int n_events ){
  if ( n_events < 25 ){ //arbitrary so that the output file is managable
    vec_to_file(vec_s,ofile, vec_name_s);
    vec_to_file(vec_h,ofile, vec_name_h);
  } else{ // if many events, only do histograms
  std::cout << "Too many events, only outputting histograms to file" << std::endl;
    vec_to_file(vec_h,ofile, vec_name_h);
  }
}

// ***** Create vector of evenly spaced numbers
std::vector<double> linspace(double start_in, double end_in, int num_in){

  std::vector<double> linspaced;
  double num = static_cast<double>(num_in);

  if (num == 0){ return linspaced; }
  if (num == 1){
      linspaced.push_back(start_in);
      return linspaced;
    }

  double delta = (end_in - start_in)/(num - 1);

  for(int i=0; i < num-1; ++i){
      linspaced.push_back(start_in + delta * i);
    }
  linspaced.push_back(end_in); // That start and end are exactly the same as the input
  return linspaced;
}

// ***** Print vector & nested list
//template <typename T2> 
void print_vector(std::vector<double>& vec){
  //std::cout << "size: " << vec.size() << std::endl;
  std::cout << " [";
  for (auto d: vec)
    std::cout << d << " ";
  std::cout << "]\n";
}
void printNestedVector(std::vector< std::vector<int> >& nested_vector){
  std::cout << "[\n";
 
  // Print the nested_list
  for ( std::vector<int> s1: nested_vector ){
    std::cout << " [";
    for (int s2: s1 ){
      std::cout <<" " << s2 << " ";
    }
    std::cout << "]\n";
  }

  std::cout << "]\n";
}
void saveNestedVector(std::vector< std::vector<double> >& nested_vector, std::string ofile, std::string vec_name){
  std::ofstream out_file;
  out_file.open(ofile, std::fstream::app); //append to file

  if ( vec_name.size() != 0 ){
    out_file << vec_name << ": ";
  }
 
  // Save the nested_list
  for ( std::vector<double> s1: nested_vector ){
    for (double s2: s1 ){
      out_file << s2 << " ";
    }
    out_file << "\n";
  }

  out_file << "\n";
  out_file.close();
}
void saveNestedVector(std::vector< std::vector<int> >& nested_vector, std::string ofile, std::string vec_name){
  std::ofstream out_file;
  out_file.open(ofile, std::fstream::app); //append to file

  if ( vec_name.size() != 0 ){
    out_file << vec_name << ": ";
  }
 
  // Save the nested_list
  for ( std::vector<int> s1: nested_vector ){
    for (int s2: s1 ){
      out_file << s2 << " ";
    }
    out_file << "\n";
  }

  out_file << "\n";
  out_file.close();
}

// Sort cos numbers into vectors
std::vector< std::vector<int> > sort_cos(int cos_n, float cos_min, float cos_step, double cos_val, std::vector< std::vector<int> > cos_vec, int sAi ){
  for (int i =0; i < cos_n; ++i){
    if ( cos_min + i*cos_step  <= cos_val && cos_val < cos_min + (i+1)*cos_step ){
      cos_vec[i][0] += 1;
      if ( sAi == 1 ){
        cos_vec[i][1] += 1;
      }
    }
  }
  return cos_vec;
}

// ***** Basic Kinematics (Based on Mathematica file)
// *** Functions for event.PUP[i] (Fourvector of the ith particle)
// *** Structure is p_lab = (Px, Py, Pz, E and M in GeV)
double ptOf(const std::vector<double>& vec ){ return sqrt(vec[0]*vec[0] + vec[1]*vec[1]); }
double enOf(const std::vector<double>& vec){ return vec[3]; }
double thetaOf(const std::vector<double>& vec){ return acos(vec[2]/(sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]))); }
double phiOf(const std::vector<double>& vec){ return phifn(ptVec(vec)); }
double phifn(const std::vector<double>& vec){
  double pi = 2 * acos(0.0);
  double px = vec[0], py = vec[1];

  if ( px > 0 && py > 0){ return atan(py/px); }
  else if ( px < 0 && py > 0 ){return atan(py/px) + pi; }
  else if ( px < 0 && py < 0 ){return pi + atan(py/px); }
  else if (px > 0 && py < 0 ){return 2*pi + atan(py/px); }
  else {return 0;}
}
double rapOf(const std::vector<double>& vec){ return (1.0/2.0)*log((vec[4]+vec[2])/(vec[4]-vec[2])); }
double etaOf(const std::vector<double>& vec){ return -log(tan((1.0/2.0)*thetaOf(vec))); }
double FourLengthSq(const std::vector<double>& vec){ return (vec[3]*vec[3] - vec[0]*vec[0] - vec[1]*vec[1] - vec[2]*vec[2]); }

std::vector<double> ptVec(const std::vector<double>& vec){ std::vector<double> out_vec{vec[1], vec[2]}; return out_vec; }
std::vector<double> ThreeVectorFrom(const std::vector<double>& vec){ std::vector<double> out_vec{vec[1], vec[2], vec[3]}; return out_vec; }
std::vector<double> addVecs(const std::vector<double>& vec1, const std::vector<double>& vec2){ 
  std::vector<double> out_vec{vec1[0] + vec2[0], vec1[1] + vec2[1], vec1[2] + vec2[2], vec1[3] + vec2[3],}; return out_vec; }

/* ***** Pre-Defined LHE Functions
- HEPEUP (http://home.thep.lu.se/~leif/LHEF/classLHEF_1_1HEPEUP.html#ae2cd4d1f01947bf6bd7c404cdefe2d8c)
  int 	NUP: The number of particle entries in the current event.
 
  int 	IDPRUP: The subprocess code for this event (as given in LPRUP).
 
  double 	XWGTUP: The weight for this event.
 
  std::pair< double, double > 	XPDWUP: The PDF weights for the two incoming partons.
 
  double 	SCALUP: The scale in GeV used in the calculation of the PDF's in this event.
 
  double 	AQEDUP: The value of the QED coupling used in this event.
 
  double 	AQCDUP: The value of the QCD coupling used in this event.
 
  std::vector< long > 	IDUP: The PDG id's for the particle entries in this event.
 
  std::vector< int > 	ISTUP: The status codes for the particle entries in this event.
 
  std::vector< std::pair< int, int > > 	MOTHUP: Indices for the first and last mother for the particle entries in this event.
 
  std::vector< std::pair< int, int > > 	ICOLUP: The colour-line indices (first(second) is (anti)colour) for the particle entries in this event.
 
  std::vector< std::vector< double > > 	PUP: Lab frame momentum (Px, Py, Pz, E and M in GeV) for the particle entries in this event.
 
  std::vector< double > 	VTIMUP:	Invariant lifetime (c*tau, distance from production to decay in mm) for the particle entries in this event.
 
  std::vector< double > 	SPINUP: Spin info for the particle entries in this event given as the cosine of the angle between the spin vector of a particle and the 3-momentum of the decaying particle, specified in the lab frame. 

*/