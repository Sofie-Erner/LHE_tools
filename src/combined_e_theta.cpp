#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <string.h>
#include <stdio.h>

#include "LHEF/LHEF.hpp"
#include "header.h"
/*  *************************************
    ***** Analysing LHE event files ***** 
    Based on https://github.com/Hitham2496/simple-lhe-reader
    Using https://twiki.cern.ch/twiki/pub/CMSPublic/SWGuideSubgroupMC/mergeLheFiles.cpp
    Using http://home.thep.lu.se/~leif/LHEF/index.html
    Adapted to match Mathematica file
    *************************************
*/

/*  *************************************
    ***** Output ***** 
    2D histgram for photon angle versus energy
    *************************************
*/

int main(int argc, char** argv) {
    // ***** Setup ***** 
    using namespace LHEF;
    Writer writer(std::cout);

    if(argc < 2) {
    std::cout << ">>>combined_e_theta.cpp::Usage:   " << argv[0] << " fileToAdd1.lhe   fileToAdd2.lhe ..." << std::endl;
    return -1;
    }

    // ***** Input Files ***** 
    std::vector<char*> fileToAddNames;

    for(int fileIt = 0; fileIt < argc-1; fileIt++) {
        fileToAddNames.push_back( argv[1+fileIt] );
        std::cout << "fileToAddName = " << fileToAddNames.at(fileIt) << std::endl ;
    }

    std::cout << "Merging " << argc-1 << " LHE files" << std::endl;

    // ***** Output File *****   
    // ** E & theta file
    std::string E_theta_file = "E_theta.dat";
    std::ofstream E_theta_f;
    E_theta_f.open(E_theta_file); //Clear file
    E_theta_f.close();

    // ***** Variables ***** 
    // ** Belle II angular Coverage
    float t_1 = 12.4, t_2 = 31.4, t_3 = 32.2, t_4 = 128.7, t_5 = 130.7, t_6 = 155.1;
    float d_t = 1, tmin = 10, tmax = 160;
    float t_bins = (tmax-tmin+d_t)/d_t;
    std::vector<double> theta_vals = linspace(tmin, tmax, t_bins);
    double theta_val, t_rad_conv = 180.0/M_PI;

    // ** Belle II energy values
    float d_E = 0.1, Emin = 1.8, Emax = 5.8;
    float E_bins = (Emax-Emin+d_E)/d_E;
    std::vector<double> E_vals = linspace(Emin, Emax, E_bins);
    double E_val;

    // ** Energy - Theta Histogram
    std::vector< std::vector<int> > E_theta;
    std::vector<int> temp_theta(t_bins,0);
    for (int i=0; i<E_bins; i++){
        E_theta.push_back(temp_theta);
    }

    // ***** Loop over LHE Files *****
    for(int fileIt = 0; fileIt < argc-1; fileIt++) {
        // ***** Get LHE file
        std::cout << " Opening " << fileToAddNames.at(fileIt) << std::endl;
        std::ifstream fileToAdd(fileToAddNames.at(fileIt), std::ios::in);
        Reader reader(fileToAdd);

        if (fileIt == 0) {
            writer.heprup = reader.heprup;
            writer.init();
        }

        // ***** Variables
        int event_number = 0;
        int n_accept = 0;

        // ** Get Cross-section
        int cross_start = reader.headerBlock.find("Integrated weight (pb)");
        std::string sub_headerblock = reader.headerBlock.substr(cross_start);
        int cross_end = sub_headerblock.find("</MGGenerationInfo>");
        std::string sub_sub_headerblock = sub_headerblock.substr(0,cross_end);
        int cross_mid = sub_sub_headerblock.find(":");
        double cross_sec = std::stod(sub_sub_headerblock.substr(cross_mid+1));

        // ** Energy - Theta Histogram (Temporary for each file)
        std::vector< std::vector<int> > E_theta_temp;
        for (int i=0; i<E_bins; i++){
            E_theta_temp.push_back(temp_theta);
        }

        // ***** Loop over Events *****
        while (reader.readEvent()){
            ++event_number;
            if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;
            writer.hepeup = reader.hepeup;
            writer.hepeup.heprup = &writer.heprup;
            HEPEUP event =  writer.hepeup; // stores the event information

            // ** Variables
            int n_photons = event.NUP;
            std::vector<double> q1, q2, p3, p4; // four-vectors for momenta
            std::vector<std::vector<double>> phot_p; //list of vectors for photon momenta
            std::vector<double> phot_theta;
            std::string proc_type = "p"; // photon
            int phot=0;

            // ** Loop over particles in event
            //std::cout << "--------" << std::endl;
            for (int i=0; i < event.NUP; i++){ // loop over particles in event
                //std::cout << "particle: " << event.IDUP[i];
                //print_vector(event.PUP[i]);
                //std::cout << std::endl;

                if ( event.IDUP[i] == 11 ){
                    if ( i > 1 ){ // only for outgoing electron
                    p3 = event.PUP[i];
                    proc_type = "e"; // electron/ positron
                    }
                }
                else if ( event.IDUP[i] == -11 ){
                    if ( i > 1 ){ // outgoing positron
                    p4 = event.PUP[i];
                    }
                } 
                else if ( event.IDUP[i] == 22 ){ // Photon
                    std::vector<double> p_lab = event.PUP[i];

                    phot_p.push_back(p_lab);

                    double temp_theta = 180-t_rad_conv*thetaOf(p_lab) ;
                    phot_theta.push_back(temp_theta);
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

            // *** Determine detected vs undetected photons photon  *** 
            int photon_n = 0, un_photons = 0;  // undected photons
            n_photons = phot;

            for ( int i=0; i< n_photons; i++ ){
                if ( (phot_theta[i] > t_1 && phot_theta[i] < t_2) || (phot_theta[i] > t_3 && phot_theta[i] < t_4) || (phot_theta[i] > t_5 && phot_theta[i] < t_6) )
                    {photon_n = i; }
                else { un_photons += 1; }
            }
            
            //std::cout << n_photons-1 << " " << un_photons << std::endl;
            // ***** Fill Histograms *****
            if ( (n_photons-1) == un_photons ){ // Only if one photon is detected
                int on_off = 0;

                // ***** Note photon properties
                q1 = phot_p[photon_n]; // photon momentum
                theta_val = phot_theta[photon_n]; //photon cos value

                //std::cout << theta_val << " " << enOf(q1) << std::endl;

                // *****  With detector cuts
                if ( ( (theta_val > t_1 && theta_val < t_2) || (theta_val > t_3 && theta_val < t_4) || (theta_val > t_5 && theta_val < t_6) ) && enOf(q1) > 1.8 ){
                    if ( proc_type == "p"){ on_off = 1; }// only one photon detected
                    else if ( proc_type == "e" ){ 
                        if ( ((180-thetaOf(p3)*t_rad_conv) <= t_6 || (180-thetaOf(p3)*t_rad_conv) >= t_1 ) 
                                && ( (180-thetaOf(p4)*t_rad_conv) <= t_6 || (180-thetaOf(p4)*t_rad_conv) >= t_1 )){
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
                        int E_val_i = 100000;
                        for (int i =0; i < E_bins; ++i){
                            //std::cout << Emin + i*d_E << " " << enOf(q1) << " " << Emin + (i+1)*d_E << std::endl;
                            if ( Emin + i*d_E  <= enOf(q1) && enOf(q1) < Emin + (i+1)*d_E ){ E_val_i = i;}
                        }
                        if (E_val_i != 100000){
                            for (int i =0; i < t_bins; ++i){
                                //std::cout << tmin + i*d_t << " " << theta_val << std::endl;
                                if ( tmin + i*d_t  <= theta_val && theta_val < tmin + (i+1)*d_t ){ E_theta_temp[E_val_i][i] += 1;}
                            }
                        }
                    }
                }
            }
        } // DONE W. LOOPING OVER EVENTS

        std::cout << "processed " << event_number << " events" << std::endl;
        std::cout << "accepted # events: " << n_accept  << std::endl;
        std::cout << "Cross-section: " << cross_sec << std::endl;
    }

    // ***** Output ***** 
    saveNestedVector(E_theta,E_theta_file,"");

    return 0;
}


