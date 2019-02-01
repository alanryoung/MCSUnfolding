#include <stdlib.h>
#include <string>
#include <vector>

// This is what we are here for
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

// Standard ROOT stuff
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraph.h"

// Read directly from the MAUS data structure.
#include "src/common_cpp/DataStructure/TOFEvent.hh"
#include "src/common_cpp/DataStructure/SciFiEvent.hh"
#include "src/common_cpp/DataStructure/SciFiTrack.hh"
#include "src/common_cpp/DataStructure/SciFiTrackPoint.hh"
#include "src/common_cpp/DataStructure/SciFiSpacePoint.hh"
#include "src/common_cpp/DataStructure/KLEvent.hh"
#include "src/common_cpp/DataStructure/CkovEvent.hh"
#include "src/common_cpp/DataStructure/EMREvent.hh"
#include "src/common_cpp/DataStructure/MCEvent.hh"
#include "src/common_cpp/DataStructure/VirtualHit.hh"
#include "src/common_cpp/DataStructure/Primary.hh"
#include "src/common_cpp/DataStructure/Spill.hh"
#include "src/common_cpp/DataStructure/Data.hh"
#include "src/common_cpp/JsonCppStreamer/IRStream.hh"

// And my own home brew object collection.
#include "Collection.h"

  
class MCSAnalysis {
 public:
  MCSAnalysis(std::string tree, std::string mctree, std::string outname);
  ~MCSAnalysis() {};

  void Write();

  TChain* GetTree(){ return chain; }
  TChain* GetMCTree(){ return mcchain; }
  void Execute();
  void dataSelection();
  void generateMCSResponse();
  void DoUnfolding();

 private:
  
  int jUS, jDS, kUS, kDS;
  
  double TOF_lower_limit;
  double TOF_upper_limit;

  double meanp;
  double sigmap;

  double USrefplaneZ;
  double DSrefplaneZ;
  int USrefplaneI;
  int DSrefplaneI;

  Collection _USset;
  Collection _DSset;

  RooUnfoldResponse resp_thetaX;
  RooUnfoldResponse resp_thetaY;
  RooUnfoldResponse resp_thetaScatt;

  int USabsPlaneI;
  int DSabsPlaneI;

  // Chain containing only the data of interest (not necessarily MC).
  TChain* chain;
  // Training tree containing the response to the volume of interest
  TChain* mcchain; 
  
  
  MAUS::TOFEvent* tofevent;
  MAUS::SciFiEvent* scifievent;
  MAUS::KLEvent* klevent;
  MAUS::CkovEvent* ckovevent;
  MAUS::EMREvent* emrevent;
  MAUS::MCEvent* mcevent;
  MAUS::Primary* primary;
  MAUS::SpecialVirtualHitArray* sphitarray;

  TFile* outfile;
  TH1D* tof10;
  TH1D* tof10_sel;
  TH1D* tof21;
  TH1D* tof21_sel;
  TH1D* calc_mom;
  
  bool MatchUSDS();
  bool PIDSelection();
  double MomentumFromTOF();
  void findVirtualPlanes();
  void FillMCSResponse(bool event_ok, Vars& US, Vars& DS, Vars& USMC, Vars& DSMC);
  void FillVarsVirtual(Vars& tmpvar, int j);
  void FillCollectionSciFi(Collection& Set, int j, int k);
  void FillVarsSciFi(Vars& tmpvar, int j, int k);
  
};
    
