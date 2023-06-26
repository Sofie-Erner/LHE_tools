#include <iostream>
#include <fstream>
#include <string>
#include <list>

#include "LHEF/LHEF.hpp"
#include "header.h"
/*  *************************************
    ***** Analysing LHE event files ***** 
    Based on https://github.com/Hitham2496/simple-lhe-reader
    Using http://home.thep.lu.se/~leif/LHEF/index.html
    Adapted to match Mathematica file
    *************************************
*/

/*  *************************************
    ***** Output ***** 
    2 files;
    1) angular distribution for b quark
    2) invariant mass for b pair
    *************************************
*/

// ***** Main *****
int main(int argn, char** argv){
  // ***** Setup
  if (argn != 2){
    std::cerr << "Usage: " << argv[0] << " event_file.lhe\n";
    return EXIT_FAILURE;
  }
  using namespace LHEF;

  Reader reader(argv[1]);
  Writer writer(std::cout);

  // ***** Input File
  std::string file = argv[1];

  // ***** Output Files
  size_t found = 0, found_temp = 0;
  std::string expr = ".lhe"; // ending to remove
  while ( found_temp != std::string::npos ){
    found = found_temp;
    found_temp = file.find("/",found+1);
  }
  int pos_end = file.find(expr)-found;
  std::cout << file << std::endl;

  // Angular Distribution
  std::string cos_file = file.substr(found,pos_end);
  cos_file = "cos.dat";
  std::ofstream cos_f;
  cos_f.open(cos_file); //Clear file
  cos_f.close();

  // Invariant Mass
  std::string inv_file = file.substr(found,pos_end);
  inv_file = "inv_mass.dat";
  std::ofstream inv_f;
  inv_f.open(inv_file); //Clear file
  inv_f.close();

  // ***** Variables ***** 
  int event_number = 0; // number of events

  // ** Process Values
  // Center-of-mass energy
  int beam2_type = reader.headerBlock.find("beam 2 type");
  int ebeam1_start = reader.headerBlock.find(" = ebeam1");
  int ebeam2_start = reader.headerBlock.find(" = ebeam2");

  std::string temp_string = "beam 2 type";
  std::string sub_1 = reader.headerBlock.substr(beam2_type + temp_string.size(), ebeam1_start-beam2_type-temp_string.size());

  temp_string = " = ebeam1  ! beam 1 total energy in GeV";
  std::string sub_2 = reader.headerBlock.substr(ebeam1_start + temp_string.size() ,ebeam2_start-ebeam1_start-temp_string.size());

  int E1 = std::stod(sub_1), E2 = std::stod(sub_2);
  double Ecms = sqrt(4*E1*E2);
  //std::cout << Ecms << std::endl;

  // ** Angular distribution
  float cos_min = -1.0, cos_max = 1.0;
  int cos_n = 41;
  float cos_step = (cos_max - cos_min)/float(cos_n);
  std::vector<double> cos_vals(cos_n,0);
  double cos_val;

  // ** Invariant Mass
  float inv_min = 0, inv_max = 250; // should be ECMS but that's way too high
  int inv_n = 41;
  float inv_step = (inv_max - inv_min)/float(inv_n);
  std::vector<double> inv_vals(inv_n,0);
  double inv_val;

  
  // ***** Loop over Events ***** 
  while (reader.readEvent()){
    ++event_number;
    //std::cout << "---------------------" << std::endl;
    //std::cout << "event no. " << event_number << std::endl;
    if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;
    //writer.eventComments() << reader.eventComments;
    writer.hepeup = reader.hepeup;
    writer.hepeup.heprup = &writer.heprup;

    HEPEUP event =  writer.hepeup; // stores the event information

    // *** Variables
    std::vector<double> p3, p4, comb_vec; // four-vectors for momenta

    // *** Loop over particles in event
    //std::cout << "--------" << std::endl;
    for (int i=0; i < event.NUP; i++){ // loop over particles in event
      //std::cout << "particle: " << event.IDUP[i] << std::endl;

      if ( event.IDUP[i] == 5 ){ //outgoing b quark
        p3 = event.PUP[i];
      }
      else if ( event.IDUP[i] == -5 ){ // outgoing b-bar quark
        p4 = event.PUP[i];
      } 
    }

    // ***** Analysis *****
    cos_val = cos(thetaOf(p3));
    inv_val = sqrt(FourLengthSq(addVecs(p3,p4)));

    //print_vector(p3);
    //print_vector(p4);
    //std::cout << inv_val << std::endl;

    // ***** Fill Histograms *****
    // ** Angular distribution
    for (int i =0; i < cos_n; ++i){
        if ( cos_min + i*cos_step  <= cos_val && cos_val < cos_min + (i+1)*cos_step ){
          cos_vals[i] += 1;
        }
        if ( inv_min + i*inv_step  <= inv_val && inv_val < inv_min + (i+1)*inv_step ){
          inv_vals[i] += 1;
        }
    }

  } // DONE W. LOOPING OVER EVENTS
  
  // Normalise histograms
  for (int i =0; i < cos_n; ++i){
    cos_vals[i] = cos_vals[i]/event_number;
    inv_vals[i] = inv_vals[i]/event_number;
  }

  // ***** Output ***** 
  vec_to_file(cos_vals,cos_file,"");
  vec_to_file(inv_vals,inv_file,"");

  std::cout << "processed " << event_number << " events" << std::endl;
  return 0.;
}
