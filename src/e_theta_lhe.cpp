#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <cstdlib>

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
    2D histgram for photon angle versus energy
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

  // ** Get Cross-section
  int cross_start = reader.headerBlock.find("Integrated weight (pb)");
  std::string sub_headerblock = reader.headerBlock.substr(cross_start);
  int cross_end = sub_headerblock.find("</MGGenerationInfo>");
  std::string sub_sub_headerblock = sub_headerblock.substr(0,cross_end);
  int cross_mid = sub_sub_headerblock.find(":");
  double cross_sec = std::stod(sub_sub_headerblock.substr(cross_mid+1));

  // ***** Input File
  std::string file = argv[1];

  // ***** Output File
  // ** E & theta file
  std::string E_theta_file = "E_theta.dat";
  std::ofstream E_theta_f; 
  E_theta_f.open(E_theta_file); //Clear file
  E_theta_f.close();
  
  // ***** Variables ***** 
  int event_number = 0; // number of events
  int n_accept = 0;

  // ** Energy and Theta
  float t_1 = 12.4, t_2 = 31.4, t_3 = 32.2, t_4 = 128.7, t_5 = 130.7, t_6 = 155.1;

  float Emin = 1.8, Emax = 5.8, tmin = 10, tmax = 160;

  float d_E = 0.1, E1, E2;
  float E_bins = (Emax-Emin+d_E)/d_E;
  std::vector<double> E_vals = linspace(Emin, Emax, E_bins);
  double E_val;

  float d_t = 1;
  float t_bins = (tmax-tmin+d_t)/d_t;
  std::vector<double> theta_vals = linspace(tmin, tmax, t_bins);
  double theta_val;
  double t_rad_conv = 180.0/M_PI;

  std::vector< std::vector<double> > E_theta;
  std::vector<double> temp_theta(t_bins,0);

  for (int i=0; i<E_bins; i++){
    E_theta.push_back(temp_theta);
  }

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
    int n_photons = event.NUP;
    std::vector<double> q1, q2, p3, p4; // four-vectors for momenta
    std::vector<std::vector<double>> phot_p; //list of vectors for photon momenta
    std::string proc_type = "p"; // photon
    int phot=0;
    double temp_theta, temp_E, temp_pz;

    // *** Loop over particles in event

    //std::cout << "--------" << std::endl;
    for (int i=0; i < event.NUP; i++){ // loop over particles in event
      //std::cout << "particle: " << event.IDUP[i];
      //print_vector(event.PUP[i]);
      //std::cout << std::endl;

      if ( event.IDUP[i] == 11 ){
        if ( i > 1 ){ // only for outgoing electron
          p3 = event.PUP[i];
          proc_type = "e"; // electron and positron in final state
        }
        else { // incoming
          E2 = enOf(event.PUP[i]);
        }
      }
      else if ( event.IDUP[i] == -11 ){
        if ( i > 1 ){ // outgoing positron
          p4 = event.PUP[i];
        }
        else { // incoming
          E1 = enOf(event.PUP[i]);
        }
      } 
      else if ( event.IDUP[i] == 22 ){ // Photon
        phot_p.push_back(event.PUP[i]);
        phot += 1;
      }
      else if(event.IDUP[i]==14 || event.IDUP[i]==12 || event.IDUP[i]==16){ 
        p3=event.PUP[i];
        proc_type = "n"; // neutrinos
      } 
      else if(event.IDUP[i]==-14 || event.IDUP[i]==-12 || event.IDUP[i]==-16){ p4=event.PUP[i]; }
      else if(event.IDUP[i]==9900032){ proc_type = "DM"; }  // Hidden Photon
      else if(event.IDUP[i]==9000005 || event.IDUP[i]==9000006 ){ proc_type = "DM"; } // Axion
      else if(event.IDUP[i]==23){ // Z boson, NTGC process
        q2 = event.PUP[i];
        proc_type = "NTGC";
      }    
    }

    //std::cout << proc_type << std::endl;

    // lab to CM frame
    double vel = (E1 - E2)/(E1 + E2);

    // *** Determine detected vs undetected photons photon  *** 
    int photon_n = 0, un_photons = 0;  // undected photons
    n_photons = phot;

    for ( int i=0; i< n_photons; i++ ){
      // Photon Detection Efficiency
      //int rand_undetect = rand() % 101; // [0,100]
      temp_theta = 180-t_rad_conv*thetaOf(phot_p[i]);

      temp_E = 1/sqrt(1-vel*vel) * ( enOf(phot_p[i]) - vel*phot_p[i][2] );
      //temp_pz = 1/sqrt(1-vel*vel) * ( phot_p[i][2] - vel*enOf(phot_p[i]) );
      //temp_theta = 180 - t_rad_conv*acos(temp_pz/temp_E);

      //if (rand_undetect <= 95){ // photon is detected
        if ( (temp_theta > t_1 && temp_theta < t_2 ) //&& temp_E > 4) 
              || (temp_theta > t_3 && temp_theta < t_4) 
              || (temp_theta > t_5 && temp_theta < t_6) )
              {photon_n = i; } // photon within detector limits
        else { un_photons += 1; }
      //}else { photon_n = i; }//un_photons += 1; } // photon is undetected
    }
    
    //std::cout << n_photons-1 << " " << un_photons << std::endl;
    // ***** Fill Histograms *****
    if ( (n_photons-1) == un_photons ){ // Only if one photon is detected
      int on_off = 0;

      // ***** Note photon properties
      q1 = phot_p[photon_n]; // photon momentum
      theta_val = 180-t_rad_conv*thetaOf(q1); //photon cos value

      // Change to CMS
      //std::cout << vel << " " << q1[2] << std::endl;
      E_val = 1/sqrt(1-vel*vel) * ( enOf(q1) - vel*q1[2] );
      //temp_pz = 1/sqrt(1-vel*vel) * ( q1[2] - vel*enOf(q1) );
      //theta_val = 180-t_rad_conv*acos(temp_pz/E_val);


      // *****  With detector cuts
      if ( ( (theta_val > t_1 && theta_val < t_2 ) //&& E_val > 4) 
          || (theta_val > t_3 && theta_val < t_4) 
          || (theta_val > t_5 && theta_val < t_6) 
          ) && E_val > 1.8 ){
        if ( proc_type == "p"){ on_off = 1; }// only one photon detected
        else if ( proc_type == "e" ){ 
          //std::cout << (180-thetaOf(p3)*t_rad_conv) << " " << (180-thetaOf(p4)*t_rad_conv) << std::endl;
          if ( ((180-thetaOf(p3)*t_rad_conv) >= t_6 || (180-thetaOf(p3)*t_rad_conv) <= t_1 ) 
                && ( (180-thetaOf(p4)*t_rad_conv) >= t_6 || (180-thetaOf(p4)*t_rad_conv) <= t_1 )){
            // one photon detected, electron/positron not
            on_off = 1;
          }
        }
        else if ( proc_type == "n" ){ on_off = 1; } // can't detect neutrinos
        else if ( proc_type == "DM" ){ on_off = 1; }
        else if ( proc_type == "NTGC" ){ on_off = 1; }

        //std::cout << on_off << std::endl;
        if ( on_off ==1 ){ // If requirements met
          n_accept += 1;

          //E_val = 1/sqrt(1-vel*vel) * ( E_val - vel*q1[2] );
          //double temp_val = 1/sqrt(1-vel*vel) * ( q1[2] - vel*enOf(q1) );
          //theta_val = 180 - t_rad_conv * acos( temp_val/E_val );
          //theta_val = 180-t_rad_conv*thetaOf(q1); // back using theta lab

          int E_val_i = 100;
          for (int i =0; i < E_bins; ++i){
            //std::cout << Emin + i*d_E << " " << enOf(q1) << " " << Emin + (i+1)*d_E << std::endl;
            if ( Emin + i*d_E  <= E_val && E_val < Emin + (i+1)*d_E ){
              E_val_i = i;
            }
          }
          if (E_val_i != 100){
            for (int i =0; i < t_bins; ++i){
              //std::cout << tmin + i*d_t << " " << theta_val << std::endl;
              if ( tmin + i*d_t  <= theta_val && theta_val < tmin + (i+1)*d_t ){
                E_theta[E_val_i][i] += 1;
              }
            }
          }
        }
      }
    }
  } // DONE W. LOOPING OVER EVENTS

  // ***** Output ***** 
  std::cout << "processed " << event_number << " events" << std::endl;
  std::cout << "accepted # events: " << n_accept  << std::endl;
  std::cout << "Cross-section: " << cross_sec << std::endl;

  // scale & save histogram
  //saveNestedVector(E_theta,E_theta_file,"");

  E_theta_f.open(E_theta_file);

  //10^3 [fb/pb]
  
  double luminosity = 20;
  double sf = ((cross_sec*pow(10,3))/(event_number));

  for ( std::vector<double> s1: E_theta ){
    for (double s2: s1 ){
      E_theta_f << s2*sf << " ";
    }
    E_theta_f << "\n";
  }
  
  E_theta_f.close();

  return 0.;
}
