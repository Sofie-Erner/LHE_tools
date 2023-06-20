#include <iostream>
#include <fstream>
#include <string>
#include <list>
//#include <bits/stdc++.h>

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
    // *** Total # photons and #s(gamma)=+1
    // Photon Helicity
    - "hsA_no_cuts":        #s(gamma)=-1     #s(gamma)=+1 (w. cuts)
    - "hsA_w_cuts":         -- || -- (no cuts)

    // Electron Helicity
    - "hsE_all":            #s(e-)=-1     #s(e-)=+1 (for all events, no cuts at all)
    - "hsE_no_cuts":        -- || -- (no cuts) 
    - "hsE_w_cuts":         -- || -- (w. cuts) 

    // s(e-)=+1
    - "hsE1_no_cuts":       #s(gamma)=+-1     #s(gamma)=+1 (no cuts)
    - "hsE1_w_cuts":        -- || -- (w. cuts) 

    // s(e-)=-1
    - "hsEm1_no_cuts":      #s(gamma)=+-1     #s(gamma)=+1 (no cuts)
    - "hsEm1_w_cuts":        -- || -- (w. cuts)      

    // *** Angular distribution
    - "cos":                cos values ([-1,-0.95,.....,+0.95,+1])
    - "cos %":              cos distribution of probability of s(gamma)=+1, #s(gamma)=+1 / # s(gamma)=+-1 (w. cuts)
    - "cos % no cuts":      -- || -- (no cuts) 

    - "cos nums":           Total # photons and # s(gamma)=+1 cos distribution (w. cuts)
    - "cos nums no cuts":   -- || -- (no cuts)
    - "cos all no cuts":    Cos distibution of s(gamma)=+1 for both s(e-)=+- 1

    - "cos h minus":        s(e-)=-1 cos distribution of s(gamma)=+1

    // *** Other
    - "s(e+) vs s(e-)":     Comparison of s(e+) and s(e-): #s(e-)=-s(e+)     #s(e-)=s(e+)
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

  // ***** Output File
  // ** General
  std::ofstream out_file;

  size_t found = 0, found_temp = 0;
  std::string expr = ".lhe"; // ending to remove

  while ( found_temp != std::string::npos ){
    found = found_temp;
    found_temp = file.find("/",found+1);
  }

  int pos_end = file.find(expr)-found;

  std::cout << file << std::endl;
  std::string f_out = file.substr(found,pos_end);
  f_out=f_out+".dat";
  std::cout << f_out << std::endl;

  out_file.open(f_out);
  out_file.close();

  // ** E & cos relationship
  std::string E_cos_file = file.substr(found,pos_end);
  E_cos_file = "E_" + E_cos_file + ".dat";
  std::ofstream E_cos_f;
  //Clear file
  E_cos_f.open(E_cos_file);
  E_cos_f.close();
  
  // ***** Variables ***** 
  //double pi = 2 * acos(0.0);
  int event_number = 0; // number of events
  // A: photon, E: Electron, E1: Electrons w. spin proj. 1, Em1: Electrons w. spin proj. -1 
  std::vector<double> sA, sE, sE1, sEm1; // lists for spin projection values

  // ** lists for number of -1 and 1 spin projections
  std::vector<int> hsA_no_cuts {0,0}, hsE_no_cuts {0,0}, hsE1_no_cuts {0,0}, hsEm1_no_cuts {0,0}; 
  std::vector<int> hsA_w_cuts {0,0}, hsE_w_cuts {0,0}, hsE1_w_cuts {0,0}, hsEm1_w_cuts {0,0};
  std::vector<int> hsE_all {0,0};
  std::vector<int> hsZ {0,0,0};

  std::vector<int> n_sE_sP {0,0}; //number of times { s(e+) = s(e-), s(e+) != s(e-)}

  float ang_cut_max = -log(tan((17.0*M_PI/180.0)*(1.0/2.0))); //-log(tan((17.0*M_PI/180.0)*(1.0/2.0))); // min and max angular cuts from detector
  float ang_cut_min = -log(tan((140.0*M_PI/180.0)*(1.0/2.0))); //-log(tan((140.0*M_PI/180.0)*(1.0/2.0)));
  //std::cout << "eta min: " << ang_cut_min << " eta max: " << ang_cut_max << std::endl;

  // ** angular distributions
  float cos_min = -1.0, cos_max = 1.0;
  int cos_n = 41;
  float cos_step = (cos_max - cos_min)/float(cos_n);
  std::vector<double> cos_vals = linspace(cos_min, cos_max, cos_n);
  double cos_val;

  // create nested vector with [#s(E)=1, #s(A)=1] for each cos value
  std::vector< std::vector<int> > cos_sE_sA, cos_sE_sA_no_cuts, cos_all, cos_h_minus;
  std::vector<double> cos_per, cos_per_no_cuts;
  std::vector<int> temp_2 = {0,0};

  for (int i = 0; i < cos_n; i++){
    cos_sE_sA.push_back(temp_2);
    cos_sE_sA_no_cuts.push_back(temp_2);
    cos_all.push_back(temp_2);
    cos_h_minus.push_back(temp_2);
  }

  // ** E & cos relationship
  float E_min = 0.0, E_max = 10.0;
  int E_n = 21;
  float E_step = (E_max - E_min)/float(E_n);
  std::vector<double> E_vals = linspace(E_min, E_max, E_n);
  vec_to_file(E_vals,E_cos_file,"Es");
  vec_to_file(cos_vals,E_cos_file,"Cos");

  std::vector< std::vector<int> > E_cos;
  std::vector<int> temp_cos(cos_n,0);

  for (int i=0; i<E_n; i++){
    E_cos.push_back(temp_cos);
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
    int iE = 0; //, iA = 0;
    int sEi = 0, sPi = 0, sAi, sZ; // spin of the electron and positron for the current event
    int n_photons = event.NUP;
    std::vector<double> q1, q2, p3, p4; // four-vectors for momenta
    std::vector<std::vector<double>> phot_p; //list of vectors for photon momenta
    std::vector<double> phot_s, phot_cos;
    std::string proc_type = "p"; // photon
    int phot=0;

    // *** Loop over particles in event
    //std::cout << " --------- " << std::endl;
    //for(int i = 0; i < event.MOTHUP.size(); i++){
     //std::cout << event.MOTHUP[i].first << ", " << event.MOTHUP[i].second << std::endl;
    //}

    //std::cout << "--------" << std::endl;
    for (int i=0; i < event.NUP; i++){ // loop over particles in event
      //std::cout << "particle: " << event.IDUP[i] << std::endl;

      if ( event.IDUP[i] == 11 ){
        if ( i <= 1){ // only incoming electron
          iE = i;
          sEi = event.SPINUP[i];
          hsE_all = into_hist(hsE_all,sEi); // all incoming electron spin balance
        }
        else if ( i > 1 ){ // only for outgoing electron
          p3 = event.PUP[i];
          proc_type = "e"; // electron/ positron
        }
      }
      else if ( event.IDUP[i] == -11 ){
        if ( i > 1 ){ // outgoing positron
          p4 = event.PUP[i];
        }
        else{ sPi = event.SPINUP[i]; } // Incoming
      } 
      else if ( event.IDUP[i] == 22 ){ // Photon
        std::vector<double> p_lab = event.PUP[i];

        phot_p.push_back(p_lab);
        phot_s.push_back(event.SPINUP[i]);

        double temp_cos = cos(thetaOf(p_lab));
        phot_cos.push_back(temp_cos);
        phot += 1;
      }
      else if(event.IDUP[i]==14 || event.IDUP[i]==12 || event.IDUP[i]==16){ 
        p3=event.PUP[i];
        proc_type = "n"; // neutrinos
      } 
      else if(event.IDUP[i]==-14 || event.IDUP[i]==-12 || event.IDUP[i]==-16){ p4=event.PUP[i]; }
      else if(event.IDUP[i]==9900032){ // Hidden Photon
        proc_type = "DM";
        //ang_cut_min = -1.681; //ang_cut_max = 1.681;
      }
      else if(event.IDUP[i]==9000005 || event.IDUP[i]==9000006 ){ // Axion
       //std::cout << "TEST!" << std::endl;
        proc_type = "DM";
        //ang_cut_min = -1.681; //ang_cut_max = 1.681;
      }
      else if(event.IDUP[i]==23){ // Z boson, NTGC process
        q2 = event.PUP[i];
        sZ = event.SPINUP[i];
        proc_type = "NTGC";
      }    
    }

    //std::cout << proc_type << std::endl;

    // *** Compare electron and positron spin ***
    if (sEi != sPi){ n_sE_sP[1] += 1; }
    else{ n_sE_sP[0] += 1; }

    // *** Determine detected vs undetected photons photon  *** 
    int photon_n = 0, un_photons = 0;  // undected photons
    n_photons = phot;

    //std::cout << ang_cut_min << " " << etaOf(phot_p[0]) << " " << ang_cut_max << std::endl;

    for ( int i=0; i< n_photons; i++ ){
      if ( etaOf(phot_p[i]) > ang_cut_min && etaOf(phot_p[i]) < ang_cut_max ){ photon_n = i; }
      else if ( etaOf(phot_p[i]) < ang_cut_min || etaOf(phot_p[i]) > ang_cut_max ){ un_photons += 1; }
    }
    
    //std::cout << n_photons-1 << " " << un_photons << std::endl;
    // ***** Fill Histograms *****
    if ( 1 == 1){ //( (n_photons-1) == un_photons ){ // Only if one photon is detected
      int on_off = 0;

      // ***** Note photon properties
      q1 = phot_p[photon_n]; // photon momentum
      cos_val = phot_cos[photon_n]; //photon cos value
      sAi = phot_s[photon_n];
      sA.push_back(sAi);
      hsA_no_cuts = into_hist(hsA_no_cuts,sAi);

      // ***** Note incoming electron properties
      sE.push_back(sEi);
      hsE_no_cuts = into_hist(hsE_no_cuts,sEi);

      // *****  With detector cuts
      if (1 == 1 ){ //( etaOf(q1) > ang_cut_min && etaOf(q1) < ang_cut_max ){
        if ( proc_type == "p"){ on_off = 1; }// only one photon detected
        else if ( proc_type == "e" ){ //&& ptOf(p3)*ptOf(p3)<0.15 && ptOf(p4)*ptOf(p4)<0.15 ){
          float temp_min = ang_cut_min; //-log(tan((150.0*M_PI/180.0)*(1.0/2.0)));
          float temp_max = -log(tan((40.0*M_PI/180.0)*(1.0/2.0))); // Update cut

          if ( (etaOf(p3) > temp_max || etaOf(p3) < temp_min ) && ( etaOf(p4) > temp_max || etaOf(p4) < temp_min )){
            // one photon detected, electron/positron not
            on_off = 1;
            /*
            if (thetaOf(p3)>(pi/180) && thetaOf(p4)>(pi/180) ){ // \theta > 1 degree
              on_off = 1;
            }
            else if (((thetaOf(p3)>(pi/180) && thetaOf(p4)<(pi/180)) or (thetaOf(p3)<(pi/180) && thetaOf(p4)>(pi/180)))){
              on_off = 1;
            }*/
          }
        }
        else if ( proc_type == "n" ){ on_off = 1; } // can't detect neutrinos
        else if ( proc_type == "DM" ){ on_off = 1; }
        else if ( proc_type == "NTGC" ){ 
          on_off = 1; 

          if ( sZ == -1 ){
            hsZ[0] += 1;
          }
          else if ( sZ == 0 ){
            hsZ[1] += 1;
          }
          else if ( sZ == 1 ){
            hsZ[2] += 1;
          }
        }

        if ( on_off ==1 ){ // If requirements met
          // *** Note photon spin projection for each electron spin projection
          hsE_w_cuts = into_hist(hsE_w_cuts,sEi);
          hsA_w_cuts = into_hist(hsA_w_cuts,sAi);

          if ( event.SPINUP[iE] == 1 ){ // if s(e-) = 1
            sE1.push_back(sAi);
            hsE1_w_cuts = into_hist(hsE1_w_cuts,sAi);
            cos_sE_sA = sort_cos(cos_n, cos_min, cos_step, cos_val, cos_sE_sA, sAi);
          } 
          else if ( event.SPINUP[iE] == -1 ){// if s(e-) = -1
            sEm1.push_back(sAi);
            hsEm1_w_cuts = into_hist(hsEm1_w_cuts,sAi);
          }
        }
      }

      // *** No cuts
      // not for SM ee->eea(a) and ee->aa(a) as we require all photons but one to be undetected
      if ( proc_type == "DM" || proc_type == "n" || proc_type == "p" || proc_type == "NTGC" || (proc_type == "e" && on_off ==1) ){
        // cos & spin
        cos_all = sort_cos(cos_n, cos_min, cos_step, cos_val, cos_all, sAi);

        if ( event.SPINUP[iE] == 1 ){
          hsE1_no_cuts = into_hist(hsE1_no_cuts,sAi);
          cos_sE_sA_no_cuts = sort_cos(cos_n, cos_min, cos_step, cos_val, cos_sE_sA_no_cuts, sAi);
        } 
        else if ( event.SPINUP[iE] == -1 ){
          hsEm1_no_cuts = into_hist(hsEm1_no_cuts,sAi);
          cos_h_minus = sort_cos(cos_n, cos_min, cos_step, cos_val, cos_h_minus, sAi);
        }


        // E & cos
        /*
        int E_val_i = 100;
        for (int i =0; i < E_n; ++i){
          if ( E_min + i*E_step  <= q1[3] && q1[3] < E_min + (i+1)*E_step ){
            E_val_i = i;
          }
        }
        if (E_val_i == 100){std::cout << "Error in photon energy" << std::endl;}

        for (int i =0; i < cos_n; ++i){
          if ( cos_min + i*cos_step  <= cos_val && cos_val < cos_min + (i+1)*cos_step ){
            E_cos[E_val_i][i] += 1;
          }
        }
        */
      }
    }
  } // DONE W. LOOPING OVER EVENTS

  for (int i =0; i < E_n; ++i){
    std::string E_name = "E" + std::to_string(i);
    vec_to_file(E_cos[i],E_cos_file,E_name);
  }

  // *** Calculate percentages for angular distribution *** 
  for (long unsigned int i = 0; i < cos_sE_sA.size(); ++i){
    if ( cos_sE_sA[i][0] != 0 ){
      cos_per.push_back(double(cos_sE_sA[i][1])/cos_sE_sA[i][0]); // calculate ratio of s(A)=1 and s(A)=-1
    }
    else{ cos_per.push_back(0.0); }
    if ( cos_sE_sA_no_cuts[i][0] != 0 ){
      cos_per_no_cuts.push_back(double(cos_sE_sA_no_cuts[i][1])/cos_sE_sA_no_cuts[i][0]);
    }
    else{ cos_per_no_cuts.push_back(0.0); }
  }

  // ***** Output ***** 
  std::cout << "sA: " << sA.size() << std::endl;
  std::cout << "sE: " << sE.size() << std::endl;
  std::cout << "sE1: " << sE1.size() << std::endl;
  std::cout << "sEm1: " << sEm1.size() << std::endl;

  // *** sA & hsA
  vec_hist_to_file(sA,hsA_no_cuts,f_out,"sA","hsA_no_cuts",event_number);
  vec_hist_to_file(sA,hsA_w_cuts,f_out,"sA","hsA_w_cuts",event_number);
  
  // *** sE
  vec_hist_to_file(sE,hsE_all,f_out,"sE","hsE_all",event_number);
  vec_hist_to_file(sE,hsE_no_cuts,f_out,"sE","hsE_no_cuts",event_number);
  vec_hist_to_file(sE,hsE_w_cuts,f_out,"sE","hsE_w_cuts",event_number);

  // *** sE1
  vec_hist_to_file(sE1,hsE1_no_cuts,f_out,"sE1","hsE1_no_cuts",event_number);
  vec_hist_to_file(sE1,hsE1_w_cuts,f_out,"sE1","hsE1_w_cuts",event_number);

  // *** sEm1
  vec_hist_to_file(sEm1,hsEm1_no_cuts,f_out,"sEm1","hsEm1_no_cuts",event_number);
  vec_hist_to_file(sEm1,hsEm1_w_cuts,f_out,"sEm1","hsEm1_w_cuts",event_number);

  // *** Angular distribution
  vec_to_file(cos_vals,f_out,"cos");
  vec_to_file(cos_per,f_out,"cos %");
  vec_to_file(cos_per_no_cuts,f_out,"cos % no cuts");

  // Total number and +- photons
  saveNestedVector(cos_sE_sA,f_out,"cos nums");
  saveNestedVector(cos_sE_sA_no_cuts,f_out,"cos nums no cuts");

  // All photons irregardless of electron helicity
  saveNestedVector(cos_all,f_out,"cos all no cuts");

  // Minus electron helicity
  saveNestedVector(cos_h_minus,f_out,"cos h minus");

  // *** Comparison of s(e+) and s(e-)
  vec_to_file(n_sE_sP,f_out,"s(e+) vs s(e-)");

  // *** Z boson spins
  vec_to_file(hsZ,f_out,"Z spin");

  std::cout << "processed " << event_number << " events\n";
  return 0.;
}
