#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <typeinfo>
#include <iostream>
#include <sstream>
#include <math.h>

// This is what we are here for
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

// Standard ROOT stuff
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraphErrors.h"
#include "TSpectrum.h"
#include "Fit/Fitter.h"
#include "Math/Functor.h"
#include "Math/RootFinder.h"
#include "TLine.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TFileInfo.h"
#include "TSystem.h"
#include "TString.h"
#include "TFrame.h"

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
#include "src/common_cpp/DataStructure/Hit.hh"
#include "src/common_cpp/DataStructure/Primary.hh"
#include "src/common_cpp/DataStructure/Spill.hh"
#include "src/common_cpp/DataStructure/Data.hh"
#include "src/common_cpp/DataStructure/GlobalEvent.hh"

// Borrowed from Chris Rogers script
#include "src/legacy/BeamTools/BTFieldConstructor.hh"
#include "src/common_cpp/Utils/JsonWrapper.hh"
#include "src/common_cpp/Utils/Globals.hh"
#include "src/common_cpp/Globals/GlobalsManager.hh"
#include "src/common_cpp/Recon/Global/GlobalTools.hh"
// #include "src/common_cpp/JsonCppStreamer/IRStream.hh"


// And my own home brew object collection.
#include "Collection.h"

  
class MCSAnalysis {
 public:
  MCSAnalysis(std::string tree, std::string mctree, std::string reftree, std::string outname, std::map<std::string, double> histlimits);
  ~MCSAnalysis();

  void Write();

  TChain* GetTree(){ return chain; }
  TChain* GetRefTree(){ return refchain; }
  TChain* GetMCTree(){ return mcchain; }
  TChain* GetMCEmptyTree(){ return mcemptychain; }
  //TFileInfo* GetFileInfo(){ return fileinfo;}
  void Execute(int mode);
  void retrieve(int retrieve_mode);
  void dataSelection(int mode);
  void dataSelection2();
  void referenceSelection();
  void referenceSelection2();
  void generateMCSResponse();
  void TruthData(int mode);
  void ConvolveWithInputDistribution(std::string distname);
  void ConvolveWithVirtualInputDistribution(std::string distname);
  void DoUnfolding();
  void DoDeconvolution(std::string model, int n_sel);
  void DoFFTDeconvolution();
  void UpdateRunInfo();
  void PlotRunInfo();
  void FitGaussian(std::string outfilename);
  void CalculateChi2(std::string outfilename, std::string distname);
  void print();
  TH1D* trkreffix(TH1D* h1);
  TH1D* trkreffiy(TH1D* h1);
  TH1D* trkreffiscatt(TH1D* h1);
  TH1D* trkreffi2scatt(TH1D* h1);
  double myfunc(double pz, double s1, double s2, double E, double delta);
  double myfunc_deriv(double pz, double s1, double s2, double E, double delta);
  double myfunc1(double x);
  double myfunc_deriv1(double x);
  double operator()(double x) const { return x*x; } 
  double Eval(double x) const { return x+x; }
  double Derivative(double x) const { return 2*x; }

  void Setangdef(double a){ angdef=a; }
  void SetTOFUpperLimit(double a){ TOF_upper_limit=a; }
  void SetTOFLowerLimit(double a){ TOF_lower_limit=a; }
  void SetMomentumUpperLimit(double a){ momentum_upper_limit=a; }
  void SetMomentumLowerLimit(double a){ momentum_lower_limit=a; }
  void SetRadialLimit(double a){ fiducial_cut_radius=a; }
  void SetGradientLimit(double a){ sigmap=a; }
  void SetModelFileName(std::string a){ modelfile=a; }
  void SetModelFileName2(std::string a){ modelfile2=a; }
  void SetModelName1(std::string a) {modelname1=a; }
  void SetModelName2(std::string a) {modelname2=a; }
  void SetModelName3(std::string a) {modelname3=a; }
  void SetMaterial(std::string a) {material=a; }
  void SetTrkrEffiName(std::string a) {trkreffiname=a; }
  void SetTrkrEffiEmptyName(std::string a) {trkreffiemptyname=a; }
  void SetUSDataname(std::string a) {USDataname=a; }
  void SetDSDataname(std::string a) {DSDataname=a; }
  void SetUSRefname(std::string a) {USRefname=a; }
  void SetDSRefname(std::string a) {DSRefname=a; }
  void SetEndofDStrackerrefname(std::string a) {EndofDStrackerrefname=a; }
  void SetEndofDStrackerdataname(std::string a) {EndofDStrackerdataname=a; }
  void SetDStrackerStation22refname(std::string a) {DStrackerStation22refname=a; }
  void SetDStrackerStation32refname(std::string a) {DStrackerStation32refname=a;  }
  void SetDStrackerStation42refname(std::string a) {DStrackerStation42refname=a; }
  void SetDStrackerStation22dataname(std::string a) {DStrackerStation22dataname=a; }
  void SetDStrackerStation32dataname(std::string a) {DStrackerStation32dataname=a; }
  void SetDStrackerStation42dataname(std::string a) {DStrackerStation42dataname=a; }
  void SetStartofDStrackerrefname(std::string a) {StartofDStrackerrefname=a; }
  void SetStartofDStrackerdataname(std::string a) {StartofDStrackerdataname=a; }
  void SetCollect_USTOF0refname(std::string a) {Collect_USTOF0refname=a; }
  void SetCollect_USTOF0dataname(std::string a) {Collect_USTOF0dataname=a; }
  void SetCollect_USTOF1refname(std::string a) {Collect_USTOF1refname=a; }
  void SetCollect_USTOF1dataname(std::string a) {Collect_USTOF1dataname=a; }
  void SetChi2dataname(std::string a) {Chi2dataname=a;}
  void SetChi2refname(std::string a) {Chi2refname=a;}
  void SetCollect_DiffuserEnddataname(std::string a) {Collect_DiffuserEnddataname=a;}
  void SetCollect_DiffuserEndrefname(std::string a) {Collect_DiffuserEndrefname=a;}
  void SetCollect_DiffuserStartdataname(std::string a) {Collect_DiffuserStartdataname=a;}
  void SetCollect_DiffuserStartrefname(std::string a) {Collect_DiffuserStartrefname=a;}
  void SetCollect_EndofUSTdataname(std::string a) {Collect_EndofUSTdataname=a;}
  void SetCollect_EndofUSTrefname(std::string a) {Collect_EndofUSTrefname=a;}
  void SetCollect_StartofUSTdataname(std::string a) {Collect_StartofUSTdataname=a;}
  void SetCollect_StartofUSTrefname(std::string a) {Collect_StartofUSTrefname=a;}
  void Setdofdataname(std::string a) {dofdataname=a;}
  void Setdofrefname(std::string a) {dofrefname=a;}
//  void Setpiddataname(std::string a) {piddataname=a;}
//  void Setpidrefname(std::string a) {pidrefname=a;}
  void SetTOF0dataname(std::string a) {TOF0dataname=a;}
  void SetTOF0refname(std::string a) {TOF0refname=a;}
  void SetTOF1dataname(std::string a) {TOF1dataname=a;}
  void SetTOF1refname(std::string a) {TOF1refname=a;}
  void SetTOF2dataname(std::string a) {TOF2dataname=a;}
  void SetTOF2refname(std::string a) {TOF2refname=a;}

  void SetParentGeometryFile(std::string a) {geometryfile=a; }
  void SetFFTBinLimit(int a) { binlimit=a; }
  void SetFileName(std::string a) {outfilename=a; }

  void SetCutTOF1SP1(bool a){cut_TOF1_1sp = a;}
  void SetCutTOF0SP1(bool a){cut_TOF0_1sp = a;}
  void SetCutTKU1track(bool a){cut_TKU_1track = a;}
  void SetCutCHI2DOF(bool a){cut_CHi2_dof = a;}
  void SetCutTKURadial(bool a){cut_TKU_radial = a;}
  void SetCutDiffuserRadial(bool a){cut_diffuser_radial = a;}
  void SetCutFiducial(bool a){cut_fiducial = a;}
  void SetCutTOF01(bool a){cut_TOF01 = a;}
  void SetCutTOF01ext(bool a){cut_TOF01ext = a;}
  void SetCutAbsorberMomentum(bool a){cut_absorber_momentum = a;}
  void SetCutExtrudeTKUTOF0(bool a){cut_ext_TKU_TOF0 = a;}


  void SetSysOffset(std::string eff, double val) { _sys[eff] = val; }

 private:
  
  int jUS, jDS, kUS, kDS;
  
  double angdef;
  double TOF_lower_limit;
  double TOF_upper_limit;
  double momentum_lower_limit;
  double momentum_upper_limit;
  double semom;
  double rmsmom;

  double errrmsmom;
  double bwmom;
  double meanp;

  bool cut_TOF1_1sp;
  bool cut_TOF0_1sp;
  bool cut_TKU_1track;
  bool cut_CHi2_dof;
  bool cut_TKU_radial;
  bool cut_diffuser_radial;
  bool cut_fiducial;
  bool cut_TOF01;
  bool cut_TOF01ext;
  bool cut_absorber_momentum;
  bool cut_ext_TKU_TOF0;

  
  double fiducial_cut_radius;
  double sigmap;
  int binlimit;
  
  std::string modelfile;
  std::string modelfile2;
  std::string modelname1;
  std::string modelname2;
  std::string modelname3;
  std::string material;
  std::string trkreffiname;
  std::string trkreffiemptyname;
  std::string geometryfile;
  std::string USDataname;
  std::string USRefname;
  std::string DSDataname;
  std::string DSRefname;
  std::string StartofDStrackerrefname;
  std::string StartofDStrackerdataname;
  std::string DStrackerStation22refname;
  std::string DStrackerStation32refname;
  std::string DStrackerStation42refname;
  std::string DStrackerStation22dataname;
  std::string DStrackerStation32dataname;
  std::string DStrackerStation42dataname;
  std::string EndofDStrackerrefname;
  std::string EndofDStrackerdataname;
  std::string Collect_USTOF0dataname;
  std::string Collect_USTOF0refname;
  std::string Collect_USTOF1dataname;
  std::string Collect_USTOF1refname;
  std::string Collect_StartofUSTdataname;
  std::string Collect_StartofUSTrefname;
  std::string Collect_EndofUSTdataname;
  std::string Collect_EndofUSTrefname;
  std::string Collect_DiffuserStartdataname;
  std::string Collect_DiffuserStartrefname;
  std::string Collect_DiffuserEnddataname;
  std::string Collect_DiffuserEndrefname;

  std::string TOF0dataname;
  std::string TOF1dataname;
  std::string TOF2dataname;
  std::string TOF0refname;
  std::string TOF1refname;
  std::string TOF2refname;
  std::string Chi2dataname;
  std::string Chi2refname;
  std::string dofdataname;
  std::string dofrefname;
//  std::string piddataname;
//  std::string pidrefname;

  double USrefplaneZ;
  double DSrefplaneZ;
  unsigned int USrefplaneI;
  unsigned int DSrefplaneI;

  Collection _USset;
  Collection _DSset;
  Collection _USAbsset;
  Collection _DSAbsset;
  Collection _USMCset;
  Collection _DSMCset;
  Collection _UStmpset;
  Collection USTruthSet;
  Collection DSTruthSet;

  RooUnfoldResponse resp_thetaX;
  RooUnfoldResponse resp_thetaY;
  RooUnfoldResponse resp_thetaScatt;
  RooUnfoldResponse resp_theta2Scatt;

  RooUnfoldResponse tresp_thetaX;
  RooUnfoldResponse tresp_thetaY;
  RooUnfoldResponse tresp_thetaScatt;
  RooUnfoldResponse tresp_theta2Scatt;

  RooUnfoldResponse mresp_thetaX;
  RooUnfoldResponse mresp_thetaY;
  RooUnfoldResponse mresp_thetaScatt;
  RooUnfoldResponse mresp_theta2Scatt;
  
  RooUnfoldResponse fresp_thetaX;
  RooUnfoldResponse fresp_thetaY;
  RooUnfoldResponse fresp_thetaScatt;
  RooUnfoldResponse fresp_theta2Scatt;
  unsigned int USabsPlaneI;
  unsigned int DSabsPlaneI;


  
  std::map<std::string, double> _histlimits;
  std::map<std::string, double> _sys;
  // Chain containing only the data of interest (not necessarily MC).
  TChain* chain;
  // Training tree containing the response to the volume of interest
  TChain* refchain; 
  // Truth tree containing the response to the volume of interest
  TChain* mcchain; 
  // Truth tree containing the response to the empty channel
  TChain* mcemptychain; 
  
  int runnumber, LastRunNumber; 
  int SpillNumber;
  MAUS::TOFEvent* tofevent;
  MAUS::SciFiEvent* scifievent;
  MAUS::KLEvent* klevent;
  MAUS::CkovEvent* ckovevent;
  MAUS::EMREvent* emrevent;
  MAUS::MCEvent* mcevent;
  MAUS::Primary* primary;
  MAUS::SpecialVirtualHitArray* sphitarray;
  MAUS::GlobalEvent* globalevent;

  TFile* outfile;
  std::string outfilename;
  TH1D* diffradius;
  TH1D* projradius;
  TH2D* trackno;
  TH2D* tofhitno;
  TH1D* pathlengthabs;
  TH1D* tof10;
  TH1D* tof10_sel;
  TH1D* tof21;
  TH1D* tof21_sel;
  TH1D* calc_mom;
  TH2D* t_cor;
  TH1D* residual01;
  TH1D* residual12;
  TH1D* residual12p6;
  TH1D* residual01j;
  TH1D* residual01short;
  TH1D* residual12p2;
  TH1D* residual01p6;
  TH1D* residualcobb;
  TH1D* TOF0Energy;
  TH1D* difEloss;
  TH2D* energyloss;
  TH1D* energylossproj;
  TProfile* energylosspro;
  TH1D* TOF1Energy;
  TH1D* TOF2Energy;
  TH2D* pzdEdx;
  TH2D* residual_plot;
  TH2D* TOF01vsTOF12;
  TH2D* TOF01forvsTOF01absfor;
  TH2D* TOFcom;
  TH2D* TOF01vsMCTruth;
  TH2D* TOF12vsMCTruth;
  TH2D* TOF12Paul6thforvsMCTruth;
  TH2D* TOF12longPaul6thforvsMCTruth;
  TH2D* TOF12Paul2ndforvsMCTruth;
  TH2D* TOF01Paul6thforvsMCTruth;
  TH2D* TOF01shortPaul6thforvsMCTruth;
  TH2D* TOF01fordownvsMCTruth;
  TH2D* TOF12fordownvsMCTruth;
  TH2D* TOF12cobbvsMCTruth;
  TH2D* TOF01forupvsTOF01fordown;
  TH2D* USTMomvsTOF01Mom;
  TH2D* USTMomvsTOF01Momref;

  TH1D* cor_mom;
  TH1D* mctof10;
  TH1D* mctof10_sel;
  TH1D* mctof21;
  TH1D* mctof21_sel;
  TH1D* mccalc_mom;
  TH1D* mctrue_mom;
  TH1D* Bfield_DStrackerX;
  TH1D* Bfield_DStrackerY;
  TH1D* Bfield_DStrackerZ;
  TH1D* cuts_accept;
  TH1D* mccuts_accept;
  TH1D* pass_Tof1;
  TH1D* pass_Tof0;
  TH1D* pass_TKU;
  TH1D* pass_TKUrad;
  TH1D* pass_Diff;
  TH1D* pass_Fid;
  TH1D* pass_Fid_wide;
  TH1D* pass_Fid_old;
  TH1D* pass_ExtTKU0;
  TH1D* pass_TOF01;
  TH1D* pass_TOF01Ex;
  TH1D* pass_mom;
  TH1D* pass_chi2;

  TH1D* Hist_TOF0Require1SP;
  TH1D* Hist_TOF1Require1SP;
  TH1D* Histref_TOF0Require1SP;
  TH1D* Histref_TOF1Require1SP;
  TH1D* Histref_TOF01;
  TH1D* Hist_TOF01;
  TH1D* Hist_Chi2;
  TH1D* Histref_Chi2;
  TH1D* Hist_TKU_max_radius;
  TH1D* Histref_TKU_max_radius;
  TH1D* Hist_Diff_max_radius;
  TH1D* Histref_Diff_max_radius;
  TH1D* Histref_TOF01Ext;
  TH1D* Hist_TOF01Ext;
  TH1D* Hist_USTOF0;
  TH1D* Histref_USTOF0;
  TH1D* Hist_Fiducial;
  TH1D* Histref_Fiducial;
  TH1D* Hist_MOM;
  TH1D* Histref_MOM;

  TH1D* Hist_NC_TOF0Require1SP;
  TH1D* Hist_NC_TOF1Require1SP;
  TH1D* Histref_NC_TOF0Require1SP;
  TH1D* Histref_NC_TOF1Require1SP;
  TH1D* Histref_NC_TOF01;
  TH1D* Hist_NC_TOF01;
  TH1D* Hist_NC_Chi2;
  TH1D* Histref_NC_Chi2;
  TH1D* Hist_NC_TKU_max_radius;
  TH1D* Histref_NC_TKU_max_radius;
  TH1D* Hist_NC_Diff_max_radius;
  TH1D* Histref_NC_Diff_max_radius;
  TH1D* Histref_NC_TOF01Ext;
  TH1D* Hist_NC_TOF01Ext;
  TH1D* Hist_NC_USTOF0;
  TH1D* Histref_NC_USTOF0;
  TH1D* Hist_NC_Fiducial;
  TH1D* Histref_NC_Fiducial;
  TH1D* Hist_NC_MOM;
  TH1D* Histref_NC_MOM;

  TH1D* mcpass_Tof1;
  TH1D* mcpass_Tof0;
  TH1D* mcpass_TKU;
  TH1D* mcpass_TKUrad;
  TH1D* mcpass_Diff;
  TH1D* mcpass_Fid;
  TH1D* mcpass_Fid_wide;
  TH1D* mcpass_Fid_old;
  TH1D* mcpass_ExtTKU0;
  TH1D* mcpass_TOF01;
  TH1D* mcpass_TOF01Ex;
  TH1D* mcpass_mom;
  TH1D* mcpass_chi2;

  TH1D* HPt_Tof1;
  TH1D* HPt_Tof0;
  TH1D* HPt_TKU;
  TH1D* HPt_TKUrad;
  TH1D* HPt_Diff;
  TH1D* HPt_Fid;
  TH1D* HPt_ExtTKU0;
  TH1D* HPt_TOF01;
  TH1D* HPt_TOF01Ex;
  TH1D* HPt_mom;
  TH1D* HPt_Ch12;


  TH1D* Residual_x;
  TH1D* Residual_y;
  TH1D* Residual_z;
  TH1D* Residual_px;
  TH1D* Residual_py;
  TH1D* Residual_pz;
  TH1D* Residual_dxdz;
  TH1D* Residual_dydz;
  TH1D* Residual_px_corrected;
  TH1D* Residual_py_corrected;
  TH1D* Residual_pz_corrected;
  TH1D* Residual_dxdz_corrected;
  TH1D* Residual_dydz_corrected;


  TH1D* scattering_pos_x;
  TH1D* scattering_pos_y;
  TH1D* theta_true_x_graph;
  TH1D* theta_true_x_bin;
  TH1D* theta_true_y_graph;
  TH1D* theta_true_scat_graph;

  TH2D* scattering_proj_x_R;
  TH2D* scattering_proj_y_R;
  
  TH2D* scattering_proj_x_resp;
  TH2D* scattering_proj_y_resp;

  TH1D* theta_meas_y_graph;
  TH1D* theta_meas_x_graph;

  std::vector<int> _run;
  std::vector<double> _emit;
  std::vector<double> _emitError;
  std::vector<double> _meanX;
  std::vector<double> _meanXerror;
  std::vector<double> _meanY;
  std::vector<double> _meanYerror;
  std::vector<double> _meandXdz;
  std::vector<double> _meandXdzError;
  std::vector<double> _meandYdz;
  std::vector<double> _meandYdzError;
  std::vector<double> path_length;
  
  bool MatchUSDS();
  bool PIDSelection(bool isdata);
  bool RadialSelection(double pz, double pos, double radius);
  bool TruthMatchUSDS();
  bool TruthTime(bool isdata);
  bool TruthRadialSelection(double pz, double pos, double radius, int j);
  bool MomentumSelection();
  std::vector<double> DefineProjectionAngles(Vars US, Vars DS);
  std::vector<double> RotDefineProjectionAngles(Vars US, Vars DS, int l);
  TH1D *defineHist2(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup);
  double MomentumFromTOF(bool isdata);
  double BetheBloch(double pz, double Imat, double Z, double A, double hw);
  double MostProbBB(double pz, double Imat, double Z, double A, double hw, double R, double z);
  double TimeofFlight();
  double TimeofFlight12();
  std::vector<double> CalculatePathLength(double pz);
  std::vector<double> rCalculatePathLength(double pz);
  double PathLengthInLH2(double pz);
  double CorMomFromTOF(double pz, double mat, double diff);
  bool findVirtualPlanes();
  void FillMuScattResponse(bool event_ok, Vars& US, Vars& DS, Vars& USMC, Vars& DSMC);
  void FillMCSResponse(bool event_ok, Vars& US, Vars& DS, Vars& USMC, Vars& DSMC);
  void TruthGraph(Collection& USMC, Collection& DSMC);
  void FillVarsVirtual(Vars& tmpvar, int j);
  void FillCollectionSciFi(Collection& Set, int j, int k, double pz, int isDS, bool project=false);
  void FillCollectionSciFi2(Collection& Set, Vars tmpVars, bool project=false);
  void FillVarsSciFi(Vars& tmpvar, int j, int k, double pz, int isDS);
  void make_beam_histograms(Collection Set, std::string desc, std::string suffix);
  void make_acceptance_histograms(Collection USset, Collection DSset, 
				  std::string desc, std::string suffix);
  void make_scattering_acceptance_histograms(Collection USset,
					     Collection DSset,
					     Collection DSrefset,
					     std::string desc,
					     std::string suffix);
  Json::Value SetupConfig(int verbose_level, std::string geometry);
  Vars PropagateVarsMu(Vars event, double z0);
  Vars PropagateVarsMuVerbose(Vars event, double z0);
  Vars FillVars(TLorentzVector a_track_mom, TLorentzVector a_track_pos, int pid);
  Vars FillVars3(MAUS::ThreeVector b_track_mom, MAUS::ThreeVector b_track_pos, int pid);
  multiVars read_globals(double DS11, double DS22, double DS32, double DS42, double DS53, double US53, double US11, double diffuserlow, double diffuserhigh, double TOF0pos, double TOF1pos);
  double Check_Radius_US(double low_z, double high_z);
  Vars reset_Vars();
  void display_Vars(Vars tmpVars);
  void ReturnWords(std::string inputstring,std::vector<std::string> &v);
  double MaxRadius(double bz,double startx, double starty, double startz, double startpx, double startpy, double endx, double endy, double endz, double endpx, double endpy, double endpz);


  /*
  static TVectorD p_vec;
  static TVectorD res;
  static TVectorD pStart_vec;
  static TVectorD pStart_vec_y;
  static TVectorD theta_true_x;
  static TVectorD theta_true_y;
  */
  double pStart[19];
  // double p_vec[19];
  
  int counter;
  int nSize_true;
  int nSize_true_y;
  double scattering_proj_x_entry;
  double scattering_proj_x_new;
  double scattering_proj_x_centre;
  double scattering_proj_y_entry;
  double scattering_proj_y_new;
  double scattering_proj_y_centre;
  
  struct SumDistance2 {
    TVectorD D_x_func;
    TMatrixD e_x_func;
    TMatrixD R_x_func;

  SumDistance2(TVectorD D, TMatrixD e, TMatrixD R) :
    D_x_func(D), e_x_func(e), R_x_func(R) {}
    double operator() (const double * par) {
      double sum = 0;
      TVectorD p_vec(19);
      int z = 1;
      for (int j=0; j<19; j++)
	p_vec[j] = par[j];
      for (int i=10; i<19; i++){
	p_vec[i] = par[i-z];
	z += 2;
      }
      TVectorD res = R_x_func * e_x_func * p_vec;
      for (int i = 0; i<30; i++) {
	if (D_x_func(i) != 0) {
	  double d = pow(D_x_func(i)-res(i),2) / D_x_func(i);
	  sum += d;
	}
	/*
	  if (first) {
	  std::cout << "point " << i << "\t"
	  << D_x_func(i) << "\t"
	  << res(i) << "\t" << std::endl;
	  }
	  }
	  if (first) {
	  std::cout << "Total Initial distance square = " << sum << std::endl;
	  first = false;
	*/
      }
      return sum;
   }
  };




	     
};
    
class MyFunction1D { 
	public:
		double E;
		double delta;
		double s1;
		double s2;
		double operator()(double pz) const { return s1*sqrt(pow(pz+delta,2)+pow(105.65,2))/(pz+delta)+s2*sqrt(pow(pz-delta,2)+pow(105.65,2))/(pz-delta); } 
		double Derivative(double pz) const { return -(s1 - s2) * delta * pow(105.65,2)/ (pow(pz,2)*E); }
};
 
