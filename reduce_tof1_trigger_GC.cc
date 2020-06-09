/* This file is part of MAUS: http://micewww.pp.rl.ac.uk/projects/maus
 *
 * MAUS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MAUS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MAUS.	If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdlib.h>
#include <string>
#include <cstddef>
#include <vector>
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TTree.h"
// Reconstruction specific headers
#include "src/common_cpp/DataStructure/TOFEvent.hh"
#include "src/common_cpp/DataStructure/SciFiEvent.hh"
#include "src/common_cpp/DataStructure/KLEvent.hh"
#include "src/common_cpp/DataStructure/CkovEvent.hh"
#include "src/common_cpp/DataStructure/EMREvent.hh"
#include "src/common_cpp/DataStructure/GlobalEvent.hh"
 // MC specific headers
#include "src/common_cpp/DataStructure/Primary.hh"
#include "src/common_cpp/DataStructure/VirtualHit.hh"
#include "src/common_cpp/DataStructure/Track.hh"
#include "src/common_cpp/DataStructure/Hit.hh"
#include "src/common_cpp/DataStructure/MCEvent.hh"
 // General headers
#include "src/common_cpp/DataStructure/Spill.hh"
#include "src/common_cpp/DataStructure/Data.hh"
#include "src/common_cpp/JsonCppStreamer/IRStream.hh"
#include "Collection.h"



/** Perform a toy analysis using ROOT
 *
 *	Plot digits for TOF1 plane 0 and plane 1. See documentation in
 *	maus_user_guide "Accessing Data"
 */
int main(int argc, char * argv[]) { //val

    std::vector < std::string > filelist;
    if (argc != 0) {
        std::cout << "This is a MAUS tree skim that selects events that \n";
        std::cout << "produce triggers in TOF2 and writes the reconstructed event \n";
        std::cout << "objects and the MCevent objects to a tree. The resulting \n";
        std::cout << "tree may be used with the \"tree->Draw()\" command on the \n";
        std::cout << "root command line or with the \"run_sp_ensemble_analysis\"\n";
        std::cout << "executable.\n" << std::endl;
        std::cout << "To run the command, provide one or more input files as arguements.\n";
    }
    std::string type = argv[2];
    std::string Fnumerate=argv[1];
    std::string fname = "reduced_tree_"+Fnumerate;
    for (int ifile = 2; ifile < argc; ifile++) {
   		std::string infilename = argv[ifile];
        TFile * file = new TFile(infilename.c_str()); // valgrind
        if (file -> IsZombie()) continue;
        TTree * Spill;
        Spill = (TTree * ) file -> Get("Spill");
        if (!Spill) continue;
        file -> Close(); // valgrind 
        std::cout << "Adding " << argv[ifile] << " to file list." << std::endl;
        filelist.push_back(argv[ifile]);
    }

    MAUS::Data * data = new MAUS::Data();
    double tof1pos = 12929.56;
    double tof2pos = 21139.3;
   

    fname += ".root";
    TFile * outfile = new TFile(fname.c_str(), "RECREATE");//valgring

    std::cout << "Make tree; to fill later on" << std::endl;
    outfile->cd();
    TTree * tree = new TTree("recon_reduced_tree", "Reduced MAUS Spill");
    tree -> SetMaxTreeSize(200000000);
    int RunNumber, SpillNumber;
    MAUS::TOFEvent * tofevent = new MAUS::TOFEvent();
    MAUS::SciFiEvent * scifievent = new MAUS::SciFiEvent();
    MAUS::CkovEvent * ckovevent = new MAUS::CkovEvent();
    MAUS::KLEvent * klevent = new MAUS::KLEvent();
    MAUS::EMREvent * emrevent = new MAUS::EMREvent();
    MAUS::MCEvent * mcevent = new MAUS::MCEvent();
    //MAUS::Primary* primary = new MAUS::Primary();
    //MAUS::SpecialVirtualHitArray* sphitarray = new MAUS::SpecialVirtualHitArray();
    MAUS::GlobalEvent* globalevent = new MAUS::GlobalEvent();
    // MAUS::TOFHitArray* tofhitarray = new MAUS::TOFHitArray();
    tree -> Branch("RunNumber", & RunNumber, "RunNumber/I");
    tree -> Branch("SpillNumber", & SpillNumber, "SpillNumber/I");
    tree -> Branch("TOFBranch", "MAUS::TOFEvent", & tofevent, 64000, 10); //valgring
    tree -> Branch("SciFiBranch", "MAUS::SciFiEvent", & scifievent, 64000, 10);
    tree -> Branch("CkovBranch", "MAUS::CkovEvent", & ckovevent, 64000, 10);
    tree -> Branch("KLBranch", "MAUS::KLEvent", & klevent, 64000, 10);
    tree -> Branch("EMRBranch", "MAUS::EMREvent", & emrevent, 64000, 10);
    tree->Branch("GlobalBranch","MAUS::GlobalEvent",&globalevent, 64000, 10);

    std::vector < std::string > ::iterator filename;

    int TOF0count = 0, TOF1count = 0, TOF2count = 0;
    for (filename = filelist.begin(); filename < filelist.end(); filename++) {
        // irstream infile((*filename).c_str(), "Spill");
        std::cout << "Considering file " << ( * filename).c_str() << std::endl;

        TFile infile(( * filename).c_str());
        TTree * Spill = (TTree * ) infile.Get("Spill");
        Spill -> SetBranchAddress("data", & data);
        // Iterate over all events
        int Nentries = Spill -> GetEntries();

        for (int i = 0; i < Nentries; i++) {
            Spill -> GetEntry(i);
            MAUS::Spill * spill = data -> GetSpill();
            // Iterate over all events; top level event in the "Spill" tree
            // corresponds to a daq_event; e.g. "start_of_burst" (spill start
            //	signal) or "physics_event" (spill data)
            if (spill != NULL && spill -> GetDaqEventType() == "physics_event") {
                // Each recon event corresponds to a particle trigger; data in the
                // recon event should all have the same trigger
                if (type.find("recon") != std::string::npos) {
                    RunNumber = spill -> GetRunNumber();
                    SpillNumber = spill -> GetSpillNumber();
                }
                int Rentries = spill -> GetReconEvents() -> size();
                for (size_t ii = 0; ii < Rentries; ++ii) {
                    // TOF event holds TOF information for the particle trigger
                    MAUS::TOFEvent * tof_event = ( * spill -> GetReconEvents())[ii] -> GetTOFEvent();
                    MAUS::TOFEventDigit digits;
                    if (type.find("mc") != std::string::npos) {
                        mcevent = ( * spill -> GetMCEvents())[ii];
                        // if(fabs(mcevent->GetPrimary()->GetParticleId()) != 13) continue;
                    }
                    if (tof_event != NULL)
                        digits = tof_event -> GetTOFEventDigit();
                    else
                        continue;
                    // std::cout<<"Found TOF digits at event "<<i<<std::endl;
                    // Iterate over TOF digits
                    bool plane0digit = false;
                    bool plane1digit = false;
                    for (size_t j = 0; j < digits.GetTOF0DigitArray().size(); ++j) {
                        MAUS::TOFDigit tof0_digit = digits.GetTOF0DigitArray()[j];
                        // get the slab where the digit was registered
                        if (tof0_digit.GetPlane() == 0) {
                            plane0digit = true;
                        } else {
                            plane1digit = true;
                        }
                    }
                    if (plane0digit && plane1digit) TOF0count++;
                    plane0digit = false;
                    plane1digit = false;
                    for (size_t j = 0; j < digits.GetTOF1DigitArray().size(); ++j) {
                        MAUS::TOFDigit tof1_digit = digits.GetTOF1DigitArray()[j];
                        // get the slab where the digit was registered
                        if (tof1_digit.GetPlane() == 0) {
                            plane0digit = true;
                        } else {
                            plane1digit = true;
                        }
                    }
                    if (plane0digit && plane1digit) {
                        TOF1count++;
                        tofevent = tof_event;
                        scifievent = ( * spill -> GetReconEvents())[ii] -> GetSciFiEvent();
                        ckovevent = ( * spill -> GetReconEvents())[ii] -> GetCkovEvent();
                        klevent = ( * spill -> GetReconEvents())[ii] -> GetKLEvent();
                        emrevent = ( * spill -> GetReconEvents())[ii] -> GetEMREvent();
			globalevent= (*spill->GetReconEvents())[ii]->GetGlobalEvent();

                        if (type.find("sim") != std::string::npos) {
                            mcevent = ( * spill -> GetMCEvents())[ii];
                        }
                        tree -> Fill();

                    }

                    plane0digit = false;
                    plane1digit = false;
                    for (size_t j = 0; j < digits.GetTOF2DigitArray().size(); ++j) {
                        MAUS::TOFDigit tof2_digit = digits.GetTOF2DigitArray()[j];
                        // get the slab where the digit was registered
                        if (tof2_digit.GetPlane() == 0) {
                            plane0digit = true;
                        } else {
                            plane1digit = true;
                        }
                    }
                    if (plane0digit && plane1digit) TOF2count++;

                }

            }


        }


    }
    outfile = tree->GetCurrentFile();
    outfile -> cd();
    tree -> Print();
    tree -> Write();

    std::cout << "TOF0 Count:" << TOF0count << std::endl;
    std::cout << "TOF1 Count:" << TOF1count << std::endl;
    std::cout << "TOF2 Count:" << TOF2count << std::endl;
    outfile -> Close("R");
}
