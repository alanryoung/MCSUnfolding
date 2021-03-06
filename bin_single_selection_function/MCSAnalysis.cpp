#include "MCSAnalysis.h"

// #include "TVirtualFFT.h"

MCSAnalysis::MCSAnalysis(std::string tree, std::string reftree, std::string outname, std::map<std::string, double> histlimits)
//: p_vec(TVectorD(19)), res(TVectorD(30)), pStart_vec(TVectorD(19)),
//    pStart_vec_y(TVectorD(19)), theta_true_x(TVectorD(19)), theta_true_y(TVectorD(19))
{
  chain_tree= tree;
  chain_reftree = reftree;
  outfile = new TFile(outname.c_str(), "recreate");

  _sys["alXUS"] = 0.0; //-0.0356;
  _sys["alYUS"] = 0.0; //1.13;
  _sys["alXDS"] = 0.0; //-0.126;
  _sys["alYDS"] = 0.0; //0.118;
  _sys["thXUS"] = 0.0; //-0.0031;
  _sys["thYUS"] = 0.0; //-0.005;
  _sys["thXDS"] = 0.0; //-0.003;
  _sys["thYDS"] = 0.0; //-0.0151;
  _sys["resX"] = 1.0;
  _sys["resY"] = 1.0;
  _sys["psel_lcut"] = -1.;
  _sys["psel_ucut"] = -1.;
  _sys["TOF0_z"]  = 5287.24720607;
  _sys["TOF1_z"]  = 12929.2608098;
  _sys["TOF2_z"]  = 21138.306769;
  _sys["niter"]   = 10;
  _sys["abspos"]  = 16952.5;
  _sys["Nevents"]  = -1.;
  _sys["FracEvents"] = 1.;
  _sys["project"] = true;

  _sys["dataDS11"] = 18836.9;
  _sys["dataDS53"] = 19937.9;
  _sys["dataUS53"] = 13967.7;
  _sys["dataUS11"] = 15068.9;
  _sys["datadiffuserlow"] = 13726.2;
  _sys["datadiffuserhigh"] = 13746.2;
  _sys["dataTOF0pos"] = 5287.25;
  _sys["dataTOF1pos"] = 12929.6;

  _sys["refDS11"] = 18848.4;
  _sys["refDS53"] = 19949.4;
  _sys["refUS53"] = 13961.7;
  _sys["refUS11"] = 15062.9;
  _sys["refdiffuserlow"] = 13720.2;
  _sys["refdiffuserhigh"] = 13740.2;
  _sys["refTOF0pos"] = 5287.25;
  _sys["refTOF1pos"] = 12929.6;



  modelname1 = "Cobb";
  modelname2 = "GEANT";
  USDataname = "USData.txt";
  DSDataname = "DSData.txt";
  USRefname = "USRef.txt";
  DSRefname = "DSRef.txt";
 
  _histlimits["NbinsXY"] = histlimits.count("NbinsXY") != 0 ? histlimits["NbinsXY"]: 200;
  _histlimits["minXY"] = histlimits.count("minXY") != 0 ? histlimits["minXY"]: -0.2;
  _histlimits["maxXY"] = histlimits.count("maxXY") != 0 ? histlimits["maxXY"]:  0.2;
  _histlimits["NbinsTh"] = histlimits.count("NbinsTh") != 0 ? histlimits["NbinsTh"]: 100;
  _histlimits["minTh"] = histlimits.count("minTh") != 0 ? histlimits["minTh"]:  0.0;
  _histlimits["maxTh"] = histlimits.count("maxTh") != 0 ? histlimits["maxTh"]:  0.2;
  _histlimits["NbinsTh2"] = histlimits.count("NbinsTh2") != 0 ? histlimits["NbinsTh2"]: 400;
  _histlimits["minTh2"] = histlimits.count("minTh2") != 0 ? histlimits["minTh2"]:  0.0;
  _histlimits["maxTh2"] = histlimits.count("maxTh2") != 0 ? histlimits["maxTh2"]:  0.004;
  
  resp_thetaX.Setup(_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  resp_thetaY.Setup(_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  resp_thetaScatt.Setup(_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
  resp_theta2Scatt.Setup(_histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);

  tresp_thetaX.Setup(_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tresp_thetaY.Setup(_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tresp_thetaScatt.Setup(_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
  tresp_theta2Scatt.Setup(_histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);

  const int NBINS = 19;
  const int NPBINS = 30;
  double scat_bin_array[NBINS+1] =
    {-0.1151, -0.10298, -0.09087, -0.07875,  -0.06664,
     -0.05452, -0.04241, -0.03029, -0.01817, -0.00606,
     0.006058,  0.018174, 0.030289,  0.042405,   0.054521,
     0.066637,   0.078753,  0.090868,  0.102984,   0.1151};
  double pos_array[NPBINS+1];
  for(int i=0; i<NPBINS+1; i++)  pos_array[i] = -100 + (6.67 * i);

  runnumber  = 0;
  LastRunNumber = 0;
  tofevent   = new MAUS::TOFEvent();
  scifievent = new MAUS::SciFiEvent();
  klevent    = new MAUS::KLEvent();
  ckovevent  = new MAUS::CkovEvent();
  emrevent   = new MAUS::EMREvent();
  mcevent    = new MAUS::MCEvent();
  primary    = new MAUS::Primary();
  sphitarray = new MAUS::SpecialVirtualHitArray();
  globalevent= new MAUS::GlobalEvent();

  // set some default values that may be adjusted later.
  TOF_lower_limit = 27;
  TOF_upper_limit = 28;
  momentum_lower_limit = 140;
  momentum_upper_limit = 300;

  //Set default data cuts
  cut_TOF1_1sp = true;
  cut_TOF0_1sp = true;
  cut_TKU_1track = true;
  cut_CHi2_dof = true;
  cut_TKU_radial = true;
  cut_diffuser_radial = true;
  cut_fiducial = true;
  cut_TOF01 = true;
  cut_TOF01ext = true;
  cut_absorber_momentum = true;
  cut_ext_TKU_TOF0 = true;


  meanp = 300.0;
  sigmap = 0.035;
  binlimit = 22;
  // int counter = 0;
  // histograms of selection variables
  
  t_cor = new TH2D("t_cor","t_cor", 1000, 0, 1000,100, 8.05, 9.15);
  tof10 = new TH1D("tof10","TOF Between Stations 1 and 0; t_{TOF1} - t_{TOF0} (ns)", 150, 10, 40);
  tof10_sel = new TH1D("tof10_sel","TOF Selection Between Stations 1 and 0; t_{TOF1} - t_{TOF0} (ns)", 150, 10, 40);
  tof21 = new TH1D("tof21","TOF Between Stations 2 and 1; t_{TOF2} - t_{TOF1} (ns)", 150, 0, 50);
  tof21_sel = new TH1D("tof21_sel","TOF Selection Between Stations 2 and 1; t_{TOF2} - t_{TOF1} (ns)", 150, 0, 50);
  mctof10 = new TH1D("mctof10","TOF Between Stations 1 and 0; t_{TOF1} - t_{TOF0} (ns)", 150, 10, 40);
  mctof10_sel = new TH1D("mctof10_sel","TOF Selection Between Stations 1 and 0; t_{TOF1} - t_{TOF0} (ns)", 150, 10, 40);
  mctof21 = new TH1D("mctof21","TOF Between Stations 2 and 1; t_{TOF2} - t_{TOF1} (ns)", 150, 0, 50);
  mctof21_sel = new TH1D("mctof21_sel","TOF Selection Between Stations 2 and 1; t_{TOF2} - t_{TOF1} (ns)", 150, 0, 50);
  calc_mom = new TH1D("calc_mom","Momentum Calculated from TOF; Momentum (MeV/c)", 400, 0, 400);
  cor_mom = new TH1D("cor_mom","Cor Momentum Calculated from TOF; Momentum (MeV/c)", 400, 0, 400);
  mccalc_mom = new TH1D("mccalc_mom","Momentum Calculated from TOF; Momentum (MeV/c)", 400, 0, 400);
  mctrue_mom = new TH1D("mctrue_mom","Momentum from Virtual Planes; Momentum (MeV/c)", 400, 0, 400);

// Cuts and accepted muon histograms for absorber present data
  cuts_accept   = new TH1D("cuts_accept", 
			   ";Selection Criteria; Surviving Events", 14, 0, 14  );
  pass_Tof1 = new TH1D("datapass_TOF1",";Passed Cut; Events",3, 0, 3);
  pass_Tof0 = new TH1D("datapass_TOF0",";Passed Cut; Events",3, 0, 3);
  pass_TKU = new TH1D("datapass_TKU",";Passed Cut; Events",3, 0, 3);
  pass_TKUrad = new TH1D("datapass_TKUrad",";Passed Cut; Events",100, 0, 200);
  pass_Diff = new TH1D("datapass_Diff",";Passed Cut; Events",100, 0, 150);
  pass_Fid = new TH1D("datapass_Fid",";Passed Cut; Events",100, 0, 200);
  pass_Fid_wide = new TH1D("datapass_Fid_wide",";Passed Cut; Events",100, 0, 2000);
  pass_ExtTKU0 = new TH1D("datapass_ExtTKU0",";Passed Cut; Events",3, 0, 3);
  pass_TOF01 = new TH1D("datapass_TOF01",";Passed Cut; Events",100, 23, 35);
  pass_TOF01Ex = new TH1D("datapass_TOF01Ex",";Passed Cut; Events",100, 23, 35);
  pass_mom = new TH1D("datapass_mom",";Passed Cut; Events",100, 100, 300);
  pass_chi2= new TH1D("datapass_chi2",";Chi2/dof; Events",50, 0, 10);

// Cuts and accepted muon histograms for no absorber reference data
  mccuts_accept = new TH1D("mccuts_accept", 
			   ";Selection Criteria; Surviving Events", 14, 0, 14);
  mcpass_Tof1 = new TH1D("refpass_TOF1",";Passed Cut; Events",3, 0, 3);
  mcpass_Tof0 = new TH1D("refpass_TOF0",";Passed Cut; Events",3, 0, 3);
  mcpass_TKU = new TH1D("refpass_TKU",";Passed Cut; Events",3, 0, 3);
  mcpass_TKUrad = new TH1D("refpass_TKUrad",";Passed Cut; Events",100, 0, 200);
  mcpass_Diff = new TH1D("refpass_Diff",";Passed Cut; Events",100, 0, 150);
  mcpass_Fid = new TH1D("refpass_Fid",";Passed Cut; Events",100, 0, 200);
  mcpass_ExtTKU0 = new TH1D("refpass_ExtTKU0",";Passed Cut; Events",3, 0, 3);
  mcpass_TOF01 = new TH1D("refpass_TOF01",";Passed Cut; Events",100, 23, 35);
  mcpass_TOF01Ex = new TH1D("refpass_TOF01Ex",";Passed Cut; Events",100, 23, 35);
  mcpass_mom = new TH1D("refpass_mom",";Passed Cut; Events",100, 100, 300);
  mcpass_chi2= new TH1D("refpass_chi2",";Chi2/dof; Events",50, 0, 10);

//Residual histograms

  Residual_x = new TH1D("Residual_x","Difference in Position (X); #Delta X; Events",
                        100,-100,100);
  Residual_y = new TH1D("Residual_y","Difference in Position (Y); #Delta Y; Events",
                        100,-100,100);
  Residual_z = new TH1D("Residual_z","Difference in Position (Z); #Delta Z; Events",
                        100,-100,100);  
  Residual_px = new TH1D("Residual_px","Difference in Momentum (PX); #Delta PX; Events",
                        100,-25,25);
  Residual_py = new TH1D("Residual_py","Difference in Momentum (PY); #Delta PY; Events",
                        100,-25,25);
  Residual_pz = new TH1D("Residual_pz","Difference in Momentum (PZ); #Delta PZ; Events",
                        100,-25,25);
  Residual_dxdz = new TH1D("Residual_dxdz","Difference in dxdz; #Delta dxdz; Events",
                        100,-0.4,0.4);
  Residual_dydz = new TH1D("Residual_dydz","Difference in dydz; #Delta dydz; Events",
                        100,-0.4,0.4);
  Residual_px_corrected = new TH1D("Residual_px_corrected","Difference in Momentum (PX); #Delta PX; Events",
                        100,-25,25);
  Residual_py_corrected = new TH1D("Residual_py_corrected","Difference in Momentum (PY); #Delta PY; Events",
                        100,-25,25);
  Residual_pz_corrected = new TH1D("Residual_pz_corrected","Difference in Momentum (PZ); #Delta PZ; Events",
                        100,-25,25);
  Residual_dxdz_corrected = new TH1D("Residual_dxdz_corrected","Difference in dxdz; #Delta dxdz; Events",
                        100,-0.4,0.4);
  Residual_dydz_corrected = new TH1D("Residual_dydz_corrected","Difference in dydz; #Delta dydz; Events",
                        100,-0.4,0.4);



  scattering_pos_x =
    new TH1D("scattering_pos_x",
	     "Change in Position (X);#Delta X; Events per mm", 
	     NBINS, scat_bin_array);
  scattering_pos_y =
    new TH1D("scattering_pos_y",
	     "Change in Position (Y);#Delta Y; Events per mm", 
	     NBINS, scat_bin_array);

  theta_true_x_graph =
    new TH1D("theta_true_x_graph",
	     "Change in Projected Angle (X);#Delta#theta_{X}; Events per radian", 
	     NBINS, scat_bin_array);
  theta_true_y_graph =
    new TH1D("theta_true_y_graph",
	     "Change in Projected Angle (Y);#Delta#theta_{Y}; Events per radian", 
	     NBINS, scat_bin_array);
  
  scattering_proj_x_resp =
    new TH2D("scattering_proj_x_resp",
	     "Change in Projected Angle (X);#Delta#theta^{Measured}_{X}; #Delta#theta^{True}_{X}", 
	     NBINS, scat_bin_array, NBINS, scat_bin_array);
  scattering_proj_y_resp =
    new TH2D("scattering_proj_y_resp",
	     "Change in Projected Angle (Y);#Delta#theta^{Measured}_{Y}; #Delta#theta^{True}_{Y}", 
	     NBINS, scat_bin_array, NBINS, scat_bin_array);
  
  scattering_proj_x_R =
    new TH2D("scattering_proj_x_R",
	     "Change in Projected Angle (X);#Delta#theta^{Measured}_{X}; #Delta#theta^{True}_{Y}", 
	     NPBINS, pos_array, NBINS, scat_bin_array);
  scattering_proj_y_R =
    new TH2D("scattering_proj_y_R",
	     "Change in Projected Angle (Y);#Delta#theta^{Measured}_{Y}; #Delta#theta^{True}_{X}", 
	     NPBINS, pos_array, NBINS, scat_bin_array);
  theta_meas_x_graph =
    new TH1D("theta_meas_x_graph","Accepted Events;#Delta#theta^{True}_X",
	     NBINS, scat_bin_array);
  theta_meas_y_graph =
    new TH1D("theta_meas_y_graph","Accepted Events;#Delta#theta^{True}_Y",
	     NBINS, scat_bin_array);
  jUS=-1;
  jDS=-1;
  kUS=-1;
  kDS=-1;
  
  USrefplaneI = -1;
  DSrefplaneI = -1;
  USrefplaneZ = 0.;
  DSrefplaneZ = 0.;
  USabsPlaneI = -1;
  DSabsPlaneI = -1;
}

MCSAnalysis::~MCSAnalysis(){
//  delete chain;
//  delete refchain;
  // delete tof10;
  // delete tof21;
  // delete tof10_sel;
  // delete tof21_sel;
  // delete calc_mom;
}

Vars MCSAnalysis::reset_Vars(){
// Set Vars to a default values
  double split_factor=0.816496580927726; // Sqrt of 2/3
  double target_momentum=200;
  Vars tmpVars;
  tmpVars.X=0.0;
  tmpVars.Y=0.0;
  tmpVars.Z=0.0;
  tmpVars.dXdz=1./2.;
  tmpVars.dYdz=1./2.;
  tmpVars.px=target_momentum*split_factor/2;
  tmpVars.py=target_momentum*split_factor/2;
  tmpVars.pz=target_momentum*split_factor;
  tmpVars.TOF12=8.22475 / 0.299792458;;
  tmpVars.TOF01=7.64186 / 0.299792458;
  tmpVars.isgood=false;
  tmpVars.pid=-13;
  tmpVars.TOFX=0.0;
  return tmpVars;
}

void MCSAnalysis::display_Vars(Vars tmpVars){
   std::cerr<<" X: "<< tmpVars.X
            <<" Y: "<< tmpVars.Y
            <<" Z: "<< tmpVars.Z
            <<" dXdz: "<< tmpVars.dXdz
            <<" dYdz: "<< tmpVars.dYdz
            <<" px: "<< tmpVars.px
            <<" py: "<< tmpVars.py
            <<" pz: "<< tmpVars.pz
            <<" TOF12: "<< tmpVars.TOF12
            <<" TOF01: "<< tmpVars.TOF01
            <<" isgood: "<< tmpVars.isgood 
            <<" Time: "<< tmpVars.TOFX <<std::endl;
}

void MCSAnalysis::Write(){
  
  outfile->cd();
  tof10->Write();
  tof21->Write();
  tof10_sel->Write();
  tof21_sel->Write();
  calc_mom->Write();
  TCanvas *c1 = new TCanvas();
  calc_mom->GetXaxis()->SetRangeUser(120,280);
  calc_mom->Draw();
  cor_mom->SetLineColor(2);
  cor_mom->Draw("SAME");
  c1->SaveAs("calc_mom.pdf");
  cor_mom->Write();
  c1->Clear();
  calc_mom->GetXaxis()->SetRangeUser(120,280);
  cor_mom->Draw();
  c1->SaveAs("cor_mom.pdf");
  mctof10->Write();
  mctof21->Write();
  mctof10_sel->Write();
  mctof21_sel->Write();
  mccalc_mom->Write();
  mctrue_mom->Write();
  cuts_accept->Write();
  mccuts_accept->Write();
  pass_Tof1->Write();
  pass_Tof0->Write();
  pass_TKU->Write();
  pass_TKUrad->Write();
  pass_Diff->Write();
  pass_Fid->Write();
  pass_Fid_wide->Write();
  pass_ExtTKU0->Write();
  pass_TOF01->Write();
  pass_TOF01Ex->Write();
  pass_mom->Write();
  pass_chi2->Write();
  mcpass_Tof1->Write();
  mcpass_Tof0->Write();
  mcpass_TKU->Write();
  mcpass_TKUrad->Write();
  mcpass_Diff->Write();
  mcpass_Fid->Write();
  mcpass_ExtTKU0->Write();
  mcpass_TOF01->Write();
  mcpass_TOF01Ex->Write();
  mcpass_mom->Write();
  mcpass_chi2->Write();
  Residual_x->Write();
  Residual_y->Write();
  Residual_z->Write();
  Residual_px->Write();
  Residual_py->Write();
  Residual_pz->Write();
  Residual_dxdz->Write();
  Residual_dydz->Write();
  Residual_px_corrected->Write();
  Residual_py_corrected->Write();
  Residual_pz_corrected->Write();
  Residual_dxdz_corrected->Write();
  Residual_dydz_corrected->Write();


  //outfile->Close();
  delete c1;
}

void MCSAnalysis::Execute(int mode=0){
  
  int verbose_level = 1;
  // MAUS::GlobalsManager::InitialiseGlobals(JsonWrapper::JsonToString(SetupConfig(verbose_level)));
  std::cout<<"Calling Globals Initialisation"<<std::endl;
  std::cout<<"ResX: "<<_sys["resX"]<<std::endl;
  MAUS::GlobalsManager::InitialiseGlobals(
             JsonWrapper::JsonToString(SetupConfig(verbose_level, geometryfile)));
  std::cerr<<"Mode: "<<mode<<std::endl;
  std::cerr<<"Completed Globals Initialisation"<<std::endl;
  std::cerr<<"Require 1 TOF1 spacepoint "<<cut_TOF1_1sp<<std::endl;
  std::cerr<<"Require exactly 1 TOF0 Spacepoint: "<<cut_TOF0_1sp<<std::endl;
  std::cerr<<"Require exactly 1 track in Upstream tracker: "<<cut_TKU_1track<<std::endl;
  std::cerr<<"Chi-squared per degree of freedom: "<<cut_CHi2_dof<<std::endl;
  std::cerr<<"Maximum radius of muon in upstream tracker: "<<cut_TKU_radial<<std::endl;
  std::cerr<<"Maximum radius of muon in diffuser: "<<cut_diffuser_radial<<std::endl;
  std::cerr<<"Fiducial: "<<cut_fiducial<<std::endl;
  std::cerr<<"TOF01 time consistent with Muon: "<<cut_TOF01<<std::endl;
  std::cerr<<"Extruded TOF01 time consistent with Muon: "<<cut_TOF01ext<<std::endl;
  std::cerr<<"Momentum at Absorber: "<<cut_absorber_momentum<<std::endl;
  std::cerr<<"Successfully extrude from Upstream tracker to TOF0: "<<cut_ext_TKU_TOF0<<std::endl;
  if (mode == 0){
    std::cerr<<"Calling Dataselection"<<std::endl;
    dataSelection(chain_reftree, _sys["refDS11"], _sys["refDS53"], _sys["refUS53"], _sys["refUS11"], _sys["refdiffuserlow"], _sys["refdiffuserhigh"], _sys["refTOF0pos"], _sys["refTOF1pos"], "ref");
//    dataSelection(chain_tree, _sys["dataDS11"], _sys["dataDS53"], _sys["dataUS53"], _sys["dataUS11"], _sys["datadiffuserlow"], _sys["datadiffuserhigh"], _sys["dataTOF0pos"], _sys["dataTOF1pos"], "data");
  }
  else {
    dataretrieve();
    std::cerr<<"Mode = 1 - reference Retrieve"<<std::endl;
    refretrieve();   
    std::cerr<<"Mode = 1 - DoUnfolding"<<std::endl;
    DoUnfolding();
    std::cerr<<"Mode = 1 - DoFFTDeconvolution"<<std::endl;
    DoFFTDeconvolution();
    std::cerr<<"Mode = 1 - ConvolveWithInputDistribution"<<std::endl;
    ConvolveWithInputDistribution(modelname1.c_str());
   std::cerr<<"Mode = 1 - DoDeconvolution"<<std::endl;
    DoDeconvolution(modelname1.c_str(), 1);
    std::cerr<<"Mode = 1 - ConvolveWithInputDistribution"<<std::endl;
    ConvolveWithInputDistribution(modelname2.c_str());
    std::cerr<<"Mode = 1 - DoDeconvolution"<<std::endl;
    DoDeconvolution(modelname2.c_str(), 1);

  }
}

/////////////////////////////////////////////////////////////////////
// Function Refretrieve
/////////////////////////////////////////////////////////////////////

void MCSAnalysis::refretrieve(){
  Collection USretrieved, DSretrieved, USPostMom, DSPostMom, DSsetUncorrected;
  USretrieved.load_collection(USRefname);
  DSretrieved.load_collection(DSRefname);
  std::cerr<<"Entries in ref USretrieved: "<< USretrieved.N()<<std::endl;
  std::cerr<<"Entries in ref DSretrieved: "<< DSretrieved.N()<<std::endl;
  outfile->cd();
//---------------------------------------------------------------------------------------
// Momentum at absorber cut. 
//---------------------------------------------------------------------------------------
  for(int i=0;i<USretrieved.N();i++){
    double momentum_absorber= sqrt((USretrieved.E(i).px*USretrieved.E(i).px)+
                                 (USretrieved.E(i).py*USretrieved.E(i).py)+ 
                                 (USretrieved.E(i).pz*USretrieved.E(i).pz));
    mcpass_mom->Fill(momentum_absorber);
//  std::cerr<<"Momentum at centre of absorber: "<<momentum_absorber<<std::endl;
    if (momentum_absorber < momentum_lower_limit || momentum_absorber > momentum_upper_limit) {
//    std::cerr<<"**Failed momentum Cut"<<std::endl;
      if (cut_absorber_momentum) continue;
    }

//    std::cerr<<"Adding to AllTOF"<<std::endl;
//    display_Vars(globalVars.UScentre_absorber);
    USPostMom.append_instance(USretrieved.E(i));
//    std::cerr<<"Size of USAllTOF: "<<USAllTOF.N()<<std::endl;
//    display_Vars(globalVars.DScentre_absorber);
    DSPostMom.append_instance(DSretrieved.E(i));
//    std::cerr<<"Size of DSAllTOF: "<<DSAllTOF.N()<<std::endl;
//    cuts_accept->Fill("Momentum",1);    
  
//---------------------------------------------------------------------------------------
// Misalignment Correction at centre of Absorber 
//---------------------------------------------------------------------------------------
    _USMCset.append_instance(USretrieved.E(i));
    DSsetUncorrected.append_instance(DSretrieved.E(i));
    if (USretrieved.E(i).isgood && DSretrieved.E(i).isgood){
      Residual_x->Fill(USretrieved.E(i).X-DSretrieved.E(i).X);
      Residual_y->Fill(USretrieved.E(i).Y-DSretrieved.E(i).Y);
      Residual_z->Fill(USretrieved.E(i).Z-DSretrieved.E(i).Z);
      Residual_px->Fill(USretrieved.E(i).px-DSretrieved.E(i).px);
      Residual_py->Fill(USretrieved.E(i).py-DSretrieved.E(i).py);
      Residual_pz->Fill(USretrieved.E(i).pz-DSretrieved.E(i).pz);
      Residual_dxdz->Fill(USretrieved.E(i).dXdz-DSretrieved.E(i).dXdz);
      Residual_dydz->Fill(USretrieved.E(i).dYdz-DSretrieved.E(i).dYdz);
    }
  }
  std::cerr<<"Reference Momentum Cut complete and uncorrected residuals calculated. Number of entries: "<<_USMCset.N()<<std::endl;
  double mean_residual_px=Residual_px->GetMean();
  double mean_residual_py=Residual_py->GetMean();
  double mean_residual_pz=Residual_pz->GetMean();
  for (int i=0; i<_USMCset.N(); i++){
     Vars tempvarDS=DSsetUncorrected.E(i);
     tempvarDS.px=DSsetUncorrected.E(i).px+mean_residual_px;
     tempvarDS.py=DSsetUncorrected.E(i).py+mean_residual_py;
     tempvarDS.pz=DSsetUncorrected.E(i).pz+mean_residual_pz;
     tempvarDS.dXdz=tempvarDS.px/tempvarDS.pz;
     tempvarDS.dYdz=tempvarDS.py/tempvarDS.pz;
     _DSMCset.append_instance(tempvarDS);
     if (_USMCset.E(i).isgood && tempvarDS.isgood){
        Residual_px_corrected->Fill(_USMCset.E(i).px-tempvarDS.px);
        Residual_py_corrected->Fill(_USMCset.E(i).py-tempvarDS.py);
        Residual_pz_corrected->Fill(_USMCset.E(i).pz-tempvarDS.pz);
        Residual_dxdz_corrected->Fill(_USMCset.E(i).dXdz-tempvarDS.dXdz);
        Residual_dydz_corrected->Fill(_USMCset.E(i).dYdz-tempvarDS.dYdz);
     }
  }
  mcpass_mom->Write();
  std::cerr<<"Reference Residual Correction complete"<<std::endl;
  make_beam_histograms(USPostMom, "Upstream, Reference", "ref_10US_PostMom");
  make_beam_histograms(DSPostMom, "Downstream, Reference", "ref_10DS_PostMom");
}


/////////////////////////////////////////////////////////////////////
// Function Dataretrieve
/////////////////////////////////////////////////////////////////////

void MCSAnalysis::dataretrieve(){
  Collection USretrieved, DSretrieved, USPostMom, DSPostMom, DSsetUncorrected;
  USretrieved.load_collection(USDataname);
  DSretrieved.load_collection(DSDataname);
  std::cerr<<"Entries in data USretrieved: "<< USretrieved.N()<<std::endl;
  std::cerr<<"Entries in data DSretrieved: "<< DSretrieved.N()<<std::endl;
  outfile->cd();

//---------------------------------------------------------------------------------------
// Momentum at absorber cut. 
//---------------------------------------------------------------------------------------
  for(int i=0;i<USretrieved.N();i++){
    double momentum_absorber= sqrt((USretrieved.E(i).px*USretrieved.E(i).px)+
                                 (USretrieved.E(i).py*USretrieved.E(i).py)+ 
                                 (USretrieved.E(i).pz*USretrieved.E(i).pz));
    pass_mom->Fill(momentum_absorber);
//  std::cerr<<"Momentum at centre of absorber: "<<momentum_absorber<<std::endl;
    if (momentum_absorber < momentum_lower_limit || momentum_absorber > momentum_upper_limit) {
//    std::cerr<<"**Failed momentum Cut"<<std::endl;
      if (cut_absorber_momentum) continue;
    }

//    std::cerr<<"Adding to AllTOF"<<std::endl;
//    display_Vars(globalVars.UScentre_absorber);
    USPostMom.append_instance(USretrieved.E(i));
//    std::cerr<<"Size of USAllTOF: "<<USAllTOF.N()<<std::endl;
//    display_Vars(globalVars.DScentre_absorber);
    DSPostMom.append_instance(DSretrieved.E(i));
//    std::cerr<<"Size of DSAllTOF: "<<DSAllTOF.N()<<std::endl;
//    cuts_accept->Fill("Momentum",1);    
  
//---------------------------------------------------------------------------------------
// Misalignment Correction at centre of Absorber 
//---------------------------------------------------------------------------------------
    _USset.append_instance(USretrieved.E(i));
    DSsetUncorrected.append_instance(DSretrieved.E(i));
    if (USretrieved.E(i).isgood && DSretrieved.E(i).isgood){
      Residual_x->Fill(USretrieved.E(i).X-DSretrieved.E(i).X);
      Residual_y->Fill(USretrieved.E(i).Y-DSretrieved.E(i).Y);
      Residual_z->Fill(USretrieved.E(i).Z-DSretrieved.E(i).Z);
      Residual_px->Fill(USretrieved.E(i).px-DSretrieved.E(i).px);
      Residual_py->Fill(USretrieved.E(i).py-DSretrieved.E(i).py);
      Residual_pz->Fill(USretrieved.E(i).pz-DSretrieved.E(i).pz);
      Residual_dxdz->Fill(USretrieved.E(i).dXdz-DSretrieved.E(i).dXdz);
      Residual_dydz->Fill(USretrieved.E(i).dYdz-DSretrieved.E(i).dYdz);
    }
  }
  std::cerr<<"Data Momentum Cut complete and uncorrected residuals calculated. Number of entries: "<<_USset.N()<<std::endl;
  double mean_residual_px=Residual_px->GetMean();
  double mean_residual_py=Residual_py->GetMean();
  double mean_residual_pz=Residual_pz->GetMean();
  for (int i=0; i<_USset.N(); i++){
     Vars tempvarDS=DSsetUncorrected.E(i);
     tempvarDS.px=DSsetUncorrected.E(i).px+mean_residual_px;
     tempvarDS.py=DSsetUncorrected.E(i).py+mean_residual_py;
     tempvarDS.pz=DSsetUncorrected.E(i).pz+mean_residual_pz;
     tempvarDS.dXdz=tempvarDS.px/tempvarDS.pz;
     tempvarDS.dYdz=tempvarDS.py/tempvarDS.pz;
     _DSset.append_instance(tempvarDS);
     if (_USset.E(i).isgood && tempvarDS.isgood){
        Residual_px_corrected->Fill(_USset.E(i).px-tempvarDS.px);
        Residual_py_corrected->Fill(_USset.E(i).py-tempvarDS.py);
        Residual_pz_corrected->Fill(_USset.E(i).pz-tempvarDS.pz);
        Residual_dxdz_corrected->Fill(_USset.E(i).dXdz-tempvarDS.dXdz);
        Residual_dydz_corrected->Fill(_USset.E(i).dYdz-tempvarDS.dYdz);
     }
  }
  pass_mom->Write();
  std::cerr<<"Data Residual Correction complete"<<std::endl;
  make_beam_histograms(USPostMom, "Upstream, Data", "data_10US_PostMom");
  make_beam_histograms(DSPostMom, "Downstream, Data", "data_10DS_PostMom");
}

/////////////////////////////////////////////////////////////////////
// Function Dataselection
/////////////////////////////////////////////////////////////////////


void MCSAnalysis::dataSelection(std::string chain_t, double dataDS11, double dataDS53, double dataUS53,double dataUS11, double datadiffuserlow, double datadiffuserhigh, double dataTOF0pos, double dataTOF1pos, std::string dataref){

  std::cerr<<"Data selection start"<<std::endl;
  TChain* chain   = new TChain(chain_t.c_str());
  std::cerr<<"Chain created"<<std::endl;
  if(dataref=="data"){
//     chain->Add(datafile.c_str);
     std::cerr<<"Working on Data file"<<std::endl;
  }
  else
  {
//     chain->Add(trainfile.c_str);
     std::cerr<<"Working on Reference file"<<std::endl;
  }

  chain->SetBranchAddress("RunNumber", &runnumber);
  chain->SetBranchAddress("TOFBranch", &tofevent);
  chain->SetBranchAddress("SciFiBranch", &scifievent);
  chain->SetBranchAddress("CkovBranch", &ckovevent);
  chain->SetBranchAddress("KLBranch", &klevent);
  chain->SetBranchAddress("EMRBranch", &emrevent);
  chain->SetBranchAddress("GlobalBranch", &globalevent);

//  double dataDS11 = _sys["dataDS11"];
//  double dataDS53 = _sys["dataDS53"];
//  double dataUS53 = _sys["dataUS53"];
//  double dataUS11 = _sys["dataUS11"];
//  double datadiffuserlow = _sys["datadiffuserlow"];
//  double datadiffuserhigh = _sys["datadiffuserhigh"];
//  double dataTOF0pos = _sys["dataTOF0pos"];
//  double dataTOF1pos = _sys["dataTOF1pos"];

  // Cuts to apply to data 

  int Nentries = chain->GetEntries();

  Collection USAllTOF, DSAllTOF, USPostTOF0, DSPostTOF0, USPostTKU, DSPostTKU, USPostTKUrad, DSPostTKUrad;
  Collection USPostDiff, DSPostDiff, USPostChi2, DSPostChi2, USPostFid, DSPostFid, USPostExtTKU0, DSPostExtTKU0;
  Collection USPostTOF01, DSPostTOF01, USPostTOF01Ex, DSPostTOF01Ex,DSsetUncorrected;
//  double pz = 0.;
//  int Nevents = 0;
  int Fiducial_Pass = 0;
  int Fiducial_Fail = 0;

  for (int i=0; i<Nentries; i++){
//    std::cerr<<"**********************************************"<< std::endl;
//    std::cerr<<"Data Event number " << i << " of " << Nentries << std::endl;
    chain->GetEntry(i);
    if (i%100000==0) std::cout<<"Event "<<i<<"\n";
//    if (i>=100) break;
    cuts_accept->Fill("All Events",1);

    multiVars globalVars;
    globalVars.USTOF0=reset_Vars();
    globalVars.USTOF1=reset_Vars();
    globalVars.UScentre_absorber=reset_Vars();
    globalVars.UStrackerUS11=reset_Vars();
    globalVars.UStrackerUS53=reset_Vars();
    globalVars.USend_of_DStracker=reset_Vars();
    globalVars.DScentre_absorber=reset_Vars();
//    std::cerr<<"centre_absorber_US created with default values"<<std::endl;
//    display_Vars(globalVars.UScentre_absorber);
//    std::cerr<<"centre_absorber_DS created with default values"<<std::endl;
//    display_Vars(globalVars.DScentre_absorber);


//----------------------------------------------------------------------
// Require 1 TOF1 Space Point
//----------------------------------------------------------------------

    double rawTOF1HitTime;
    if ( tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArraySize() == 1 ){
      rawTOF1HitTime = tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime();
//      std::cerr<<"TOF1 time: "<<rawTOF1HitTime<<std::endl;
    }else if (tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArraySize() > 1 ){
//      std::cerr<<"**Multiple TOF1 hits - cutting"<<std::endl;
      pass_Tof1->Fill("Fail",1);
      if (cut_TOF1_1sp) continue;
    }else{
//      std::cerr<<"**No TOF1 hit"<<std::endl;
      pass_Tof1->Fill("Fail",1);
      if (cut_TOF1_1sp) continue;
    }
//    std::cerr<<"Adding to AllTOF"<<std::endl;
//    display_Vars(globalVars.UScentre_absorber);
    USAllTOF.append_instance(globalVars.UScentre_absorber);
//    std::cerr<<"Size of USAllTOF: "<<USAllTOF.N()<<std::endl;
//    display_Vars(globalVars.DScentre_absorber);
    DSAllTOF.append_instance(globalVars.DScentre_absorber);
//    std::cerr<<"Size of DSAllTOF: "<<DSAllTOF.N()<<std::endl;
    cuts_accept->Fill("1 TOF1",1);
    pass_Tof1->Fill("Pass",1);


//----------------------------------------------------------------------
// Require 1 TOF0 Space Point
//----------------------------------------------------------------------

    double rawTOF0HitTime;
    if ( tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArraySize() == 1 ){
      rawTOF0HitTime = tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[0].GetTime();
//      std::cerr<<"TOF0 time: "<<rawTOF0HitTime<<std::endl;
    }else if (tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArraySize() > 1 ){
//      std::cerr<<"**Multiple TOF0 hits - cutting"<<std::endl;
      pass_Tof0->Fill("Fail",1);
      if (cut_TOF0_1sp) continue;
    }else{
//      std::cerr<<"**No TOF0 hit"<<std::endl;
      pass_Tof0->Fill("Fail",1);
      if (cut_TOF0_1sp) continue;
    }
    
//    std::cerr<<"Adding to PostTOF0"<<std::endl;
//    display_Vars(globalVars.UScentre_absorber);
    USPostTOF0.append_instance(globalVars.UScentre_absorber);
//    std::cerr<<"Size of USPostTOF0: "<<USPostTOF0.N()<<std::endl;
//    display_Vars(globalVars.DScentre_absorber);
    DSPostTOF0.append_instance(globalVars.DScentre_absorber);
//    std::cerr<<"Size of DSPostTOF0: "<<DSPostTOF0.N()<<std::endl;
    cuts_accept->Fill("1 TOF0",1);
      pass_Tof0->Fill("Pass",1);

//---------------------------------------------------------------------------------------
// Require only 1 US track
//---------------------------------------------------------------------------------------


//  std::cerr<<" Number of tracks: "<<scifievent->scifitracks().size()<<std::endl;
  int number_upstream_tracks=0;
  for(size_t j=0; j<scifievent->scifitracks().size(); j++){
      int tracker = scifievent->scifitracks()[j]->tracker();
      if(tracker==0) number_upstream_tracks++;  
  }
//  std::cerr<<"Number of Upstream tracks: "<< number_upstream_tracks << std::endl;
  if (number_upstream_tracks==0){
//    std::cerr<<"**Fail 1 US track cut, no tracks"<<std::endl;
    pass_TKU->Fill("Fail",1);
    if (cut_TKU_1track) continue;    
  } else if(number_upstream_tracks>1){
//    std::cerr<<"**Fail 1 US track cut, more than 1 track"<<std::endl; 
    pass_TKU->Fill("Fail",1);
    if (cut_TKU_1track) continue;
  } 
  
// Output track info
  for ( size_t j=0; j<scifievent->scifitracks().size(); j++){
//    int npoints = scifievent->scifitracks()[j]->scifitrackpoints().size();
    int tracker = scifievent->scifitracks()[j]->tracker();
    if(tracker==0){
//     std::cerr<<"Upstream track scifi points"<<std::endl;
    } else {
//     std::cerr<<"Downstream track scifi points"<<std::endl;
    }
/*    for ( int k=0; k<npoints; k++){
      double zpos = scifievent->scifitracks()[j]->scifitrackpoints()[k]->pos().z();
      double xpos = scifievent->scifitracks()[j]->scifitrackpoints()[k]->pos().x();
      double ypos = scifievent->scifitracks()[j]->scifitrackpoints()[k]->pos().y();
//      std::cerr<<zpos<<" "<<xpos<<" "<<ypos<<std::endl;
    }*/
  }

    
/*
  int nUStracks= 0;
  if (scifievent->scifitracks().size()>=1){
    for ( size_t j=0; j<scifievent->scifitracks().size(); j++){
      if(scifievent->scifitracks()[j]->tracker()==0) nUStracks++ ;
    }
  }else{
    std::cerr<<"**Fail 1 US track cut, no tracks"<<std::endl;
    pass_TKU->Fill("Fail",1);
    continue;
  }
  if (nUStracks!=1){
    std::cerr<<"**Fail 1 US track cut, more than 1 track"<<std::endl; 
    pass_TKU->Fill("Fail",1);
    continue;
  }
*/
//    std::cerr<<"Adding to AllTOF"<<std::endl;
//    display_Vars(globalVars.UScentre_absorber);
    USPostTKU.append_instance(globalVars.UScentre_absorber);
//    std::cerr<<"Size of USAllTOF: "<<USAllTOF.N()<<std::endl;
//    display_Vars(globalVars.DScentre_absorber);
    DSPostTKU.append_instance(globalVars.DScentre_absorber);
//    std::cerr<<"Size of DSAllTOF: "<<DSAllTOF.N()<<std::endl;
    cuts_accept->Fill("1 US track",1);
    pass_TKU->Fill("Pass",1);

//---------------------------------------------------------------------------------------
// CHi2/dof cut
//---------------------------------------------------------------------------------------

    std::vector<MAUS::SciFiTrack*> tracks = scifievent->scifitracks();
    std::vector<MAUS::SciFiTrack*>::iterator track_iter;
    if(tracks.size() != 0){ // use if only want to take the first track
      track_iter = tracks.begin();
      double TKU_chiSquare = (*track_iter)->chi2();
      int TKU_ndof = (*track_iter)->ndf();
      double TKU_chiSquare_per_ndof = TKU_chiSquare / (1.0*TKU_ndof);
      pass_chi2->Fill(TKU_chiSquare_per_ndof);  
      if (TKU_chiSquare_per_ndof>5){
//        std::cerr<<"**Failed Chi2/dof cut"<<std::endl;
        if (cut_CHi2_dof) continue;
      }
    }else{ 
//      std::cerr<<"************************No tracks - why not cut earlier**************************"<<std::endl;
      if (cut_CHi2_dof) continue;
    }

//    std::cerr<<"Adding to AllTOF"<<std::endl;
//    display_Vars(globalVars.UScentre_absorber);
    USPostChi2.append_instance(globalVars.UScentre_absorber);
//    std::cerr<<"Size of USAllTOF: "<<USAllTOF.N()<<std::endl;
//    display_Vars(globalVars.DScentre_absorber);
    DSPostChi2.append_instance(globalVars.DScentre_absorber);
//    std::cerr<<"Size of DSAllTOF: "<<DSAllTOF.N()<<std::endl;
    cuts_accept->Fill("Chi2",1);    
    
//----------------------------------------------------------------------
// Retrieve Global Data - US track and DS track if present. Cut if no US track
//----------------------------------------------------------------------

    globalVars=read_globals(dataDS11, dataDS53, dataUS53, dataUS11, datadiffuserlow, datadiffuserhigh, dataTOF0pos, dataTOF1pos);
//    std::cerr<<"centre_absorber_US has updated globals data"<<std::endl;
//    display_Vars(globalVars.UScentre_absorber);
//    display_Vars(globalVars.DScentre_absorber);
    if (!globalVars.UScentre_absorber.isgood) {
//      std::cerr<<"**Failed US track Cut - Why not cut earlier ********************"<<std::endl;
      continue;
    }
    if (!globalVars.DScentre_absorber.isgood) {
//      std::cerr<<"**No DS ********************"<<std::endl;
    }
    cuts_accept->Fill("US Track",1);

//------------------------------------------------------------------------------------------------------
// US Tracker Radial Cut
//------------------------------------------------------------------------------------------------------

    double TKU_max_rad = Check_Radius_US(_sys["US53"],_sys["US11"]);
    pass_TKUrad->Fill(TKU_max_rad);
    if (TKU_max_rad>150.0){
//      std::cerr<<"**Failed US Tracker Radial Cut"<<std::endl;
      if (cut_TKU_radial) continue;
    }
//    std::cerr<<"Adding to AllTOF"<<std::endl;
//    display_Vars(globalVars.UScentre_absorber);
    USPostTKUrad.append_instance(globalVars.UScentre_absorber);
//    std::cerr<<"Size of USAllTOF: "<<USAllTOF.N()<<std::endl;
//    display_Vars(globalVars.DScentre_absorber);
    DSPostTKUrad.append_instance(globalVars.DScentre_absorber);
//    std::cerr<<"Size of DSAllTOF: "<<DSAllTOF.N()<<std::endl;
    cuts_accept->Fill("US tracker Radius",1);

//------------------------------------------------------------------------------------------------------
// Diffuser Cut
//------------------------------------------------------------------------------------------------------

//Check_Radius_US(double low_z, double high_z, double radius)
    double diff_max_rad = Check_Radius_US(_sys["diffuserlow"],_sys["diffuserhigh"]);
    pass_Diff->Fill(diff_max_rad);
    if (diff_max_rad >90.0){
//      std::cerr<<"**Failed Diffuser Cut"<<std::endl;
      if (cut_diffuser_radial) continue;
    }
//    std::cerr<<"Adding to AllTOF"<<std::endl;
//    display_Vars(globalVars.UScentre_absorber);
    USPostDiff.append_instance(globalVars.UScentre_absorber);
//    std::cerr<<"Size of USAllTOF: "<<USAllTOF.N()<<std::endl;
//    display_Vars(globalVars.DScentre_absorber);
    DSPostDiff.append_instance(globalVars.DScentre_absorber);
//    std::cerr<<"Size of DSAllTOF: "<<DSAllTOF.N()<<std::endl;
    cuts_accept->Fill("Diffuser",1);

//------------------------------------------------------------------------------------------------------
// Fiducial Cut
//------------------------------------------------------------------------------------------------------

    double fiducial_radius = sqrt(globalVars.USend_of_DStracker.X*globalVars.USend_of_DStracker.X + globalVars.USend_of_DStracker.Y*globalVars.USend_of_DStracker.Y);
    pass_Fid->Fill(fiducial_radius);
    pass_Fid_wide->Fill(fiducial_radius);
//    std::cerr<<"Event "<<i<<" Fiducial Radius: "<<fiducial_radius;
    if ( fiducial_radius > meanp){
      Fiducial_Fail++;
/*      std::cerr<<"  **Failed Fiducial Cut"<<std::endl;
      display_Vars(globalVars.USTOF0);
      display_Vars(globalVars.USTOF1);
      display_Vars(globalVars.UStrackerUS53);
      display_Vars(globalVars.UStrackerUS11);
      display_Vars(globalVars.UScentre_absorber);
      display_Vars(globalVars.DScentre_absorber);
      display_Vars(globalVars.USend_of_DStracker);
      PropagateVarsMuVerbose(globalVars.UStrackerUS11, globalVars.USend_of_DStracker.Z);*/
      if (cut_fiducial) continue;
    } else {
        Fiducial_Pass++;
//        std::cerr<<std::endl;
//      display_Vars(globalVars.USend_of_DStracker);
//      std::cerr<<sqrt(globalVars.USend_of_DStracker.X*globalVars.USend_of_DStracker.X + globalVars.USend_of_DStracker.Y*globalVars.USend_of_DStracker.Y)<<std::endl;
    }
//    std::cerr<<"Adding to AllTOF"<<std::endl;
//    display_Vars(globalVars.UScentre_absorber);
    USPostFid.append_instance(globalVars.UScentre_absorber);
//    std::cerr<<"Size of USAllTOF: "<<USAllTOF.N()<<std::endl;
//    display_Vars(globalVars.DScentre_absorber);
    DSPostFid.append_instance(globalVars.DScentre_absorber);
//    std::cerr<<"Size of DSAllTOF: "<<DSAllTOF.N()<<std::endl;
    cuts_accept->Fill("Fiducial",1);



//---------------------------------------------------------------------------------------
// TOF01 cut. Require TOF0 and TOF1 hit and TOF01 in bounds set by Upper and lower limits
//---------------------------------------------------------------------------------------

    double rawTOF2HitTime;     
    if ( tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArraySize() >= 1 ){
      rawTOF2HitTime = tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray()[0].GetTime();
//      std::cerr<<"TOF2 time: "<<rawTOF2HitTime<<std::endl;
    }else{
//      std::cerr<<"No TOF2 hit?"<<std::endl;
      rawTOF2HitTime = rawTOF1HitTime + 100.0 * 8.22475 / 0.299792458; // ns.
    }
    tof10->Fill(rawTOF1HitTime - rawTOF0HitTime);
    tof21->Fill(rawTOF2HitTime - rawTOF1HitTime);
    double rawTOF01HitTime = rawTOF1HitTime - rawTOF0HitTime;
    pass_TOF01->Fill(rawTOF01HitTime);
    if ( rawTOF01HitTime < TOF_lower_limit ||
         rawTOF01HitTime > TOF_upper_limit){
//      std::cerr<<"**Failed TOF01 Cut"<<std::endl;
      if (cut_TOF01) continue;
    }

//    std::cerr<<"Adding to AllTOF"<<std::endl;
//    display_Vars(globalVars.UScentre_absorber);
    USPostTOF01.append_instance(globalVars.UScentre_absorber);
//    std::cerr<<"Size of USAllTOF: "<<USAllTOF.N()<<std::endl;
//    display_Vars(globalVars.DScentre_absorber);
    DSPostTOF01.append_instance(globalVars.DScentre_absorber);
//    std::cerr<<"Size of DSAllTOF: "<<DSAllTOF.N()<<std::endl;

    tof10_sel->Fill(rawTOF1HitTime - rawTOF0HitTime);
    tof21_sel->Fill(rawTOF2HitTime - rawTOF1HitTime);
    cuts_accept->Fill("TOF01 Selection",1);


//---------------------------------------------------------------------------------------
// TOF01 Extrapolated cut. 
//---------------------------------------------------------------------------------------

  double TOF01Ext = globalVars.USTOF1.TOFX-globalVars.USTOF0.TOFX;
  pass_TOF01Ex->Fill(TOF01Ext);
  if (TOF01Ext < TOF_lower_limit || TOF01Ext > TOF_upper_limit) {
//    std::cerr<<"**Failed TOF01 Extrapolation Cut"<<std::endl;
    if(cut_TOF01ext) continue;
  }

//    std::cerr<<"Adding to AllTOF"<<std::endl;
//    display_Vars(globalVars.UScentre_absorber);
    USPostTOF01Ex.append_instance(globalVars.UScentre_absorber);
//    std::cerr<<"Size of USAllTOF: "<<USAllTOF.N()<<std::endl;
//    display_Vars(globalVars.DScentre_absorber);
    DSPostTOF01Ex.append_instance(globalVars.DScentre_absorber);
//    std::cerr<<"Size of DSAllTOF: "<<DSAllTOF.N()<<std::endl;
    cuts_accept->Fill("Extrapolated TOF01",1);


//---------------------------------------------------------------------------------------
// Successful extrapolation from  US tracker to TOF0
//---------------------------------------------------------------------------------------
       
    if (globalVars.USTOF0.isgood==false){
//      std::cerr<<"**Failed extrapolation from TKU to TOF0"<<std::endl;
      pass_ExtTKU0->Fill("Fail",1);
      if (cut_ext_TKU_TOF0) continue;
    }
//    std::cerr<<"Adding to AllTOF"<<std::endl;
//    display_Vars(globalVars.UScentre_absorber);
    USPostExtTKU0.append_instance(globalVars.UScentre_absorber);
//    std::cerr<<"Size of USAllTOF: "<<USAllTOF.N()<<std::endl;
//    display_Vars(globalVars.DScentre_absorber);
    DSPostExtTKU0.append_instance(globalVars.DScentre_absorber);
//    std::cerr<<"Size of DSAllTOF: "<<DSAllTOF.N()<<std::endl;
    cuts_accept->Fill("TKU_TOF0",1);
    pass_ExtTKU0->Fill("Pass",1);

//------------------------------------------------------------------------------------------------------
// Saving of events that pass all cuts
//------------------------------------------------------------------------------------------------------

//    std::cerr<<"**Event accepted"<<std::endl;
    _USset.append_instance(globalVars.UScentre_absorber);
    _DSset.append_instance(globalVars.DScentre_absorber);
//    _DSset.append_instance(globalVars.DScentre_absorber);
    cuts_accept->Fill("All Cuts",1);
  }
  _USset.save_collection("US"+dataref+".txt");
  _DSset.save_collection("DS"+dataref+".txt");
//  std::cerr<<"Size of USAllTOF: "<<USAllTOF.N()<<std::endl;
//  std::cerr<<"Size of DSAllTOF: "<<DSAllTOF.N()<<std::endl;
//  std::cerr<<"Size of USPreRadSel: "<<USPreRadSel.N()<<std::endl;
//  std::cerr<<"Size of DSPreRadSel: "<<DSPreRadSel.N()<<std::endl;



//------------------------------------------------------------------------------------------------------
// Make beam histograms
//------------------------------------------------------------------------------------------------------
  

  make_beam_histograms(USAllTOF, "Upstream, "+dataref, dataref+"_01US_alltof");
  make_beam_histograms(DSAllTOF, "Downstream, "+dataref, dataref+"_01DS_alltof");
  make_beam_histograms(USPostTOF0, "Upstream, "+dataref, dataref+"_02US_posttof0");
  make_beam_histograms(DSPostTOF0, "Downstream, "+dataref, dataref+"_02DS_posttof0");
  make_beam_histograms(USPostTKU, "Upstream, "+dataref, dataref+"_03US_PostTKU");
  make_beam_histograms(DSPostTKU, "Downstream, "+dataref, dataref+"_03DS_PostTKU"); 
  make_beam_histograms(USPostTKUrad, "Upstream, "+dataref, dataref+"_04US_PostTKUrad");
  make_beam_histograms(DSPostTKUrad, "Downstream, "+dataref, dataref+"_04DS_PostTKUrad");
  make_beam_histograms(USPostDiff, "Upstream, "+dataref, dataref+"_05US_PostDiff");
  make_beam_histograms(DSPostDiff, "Downstream, "+dataref, dataref+"_05DS_PostDiff");
  make_beam_histograms(USPostFid, "Upstream, "+dataref, dataref+"_06US_PostFid");
  make_beam_histograms(DSPostFid, "Downstream, "+dataref, dataref+"_06DS_PostFid");
  make_beam_histograms(USPostExtTKU0, "Upstream, "+dataref, dataref+"_07US_alltof");
  make_beam_histograms(DSPostExtTKU0, "Downstream, "+dataref, dataref+"_07DS_alltof");
  make_beam_histograms(USPostTOF01, "Upstream, "+dataref, dataref+"_08US_PostTOF01");
  make_beam_histograms(DSPostTOF01, "Downstream, "+dataref, dataref+"_08DS_PostTOF01");
  make_beam_histograms(USPostTOF01Ex, "Upstream, "+dataref, dataref+"_09US_PostTOF01Ex");
  make_beam_histograms(DSPostTOF01Ex, "Downstream, "+dataref, dataref+"_09DS_PostTOF01Ex"); 
  std::cerr<<"Number of events Pass fiducial cut: "<<Fiducial_Pass<<std::endl;
  std::cerr<<"Number of events Fail fiducial cut: "<<Fiducial_Fail<<std::endl;
  delete chain;
}


std::vector<double> MCSAnalysis::DefineProjectionAngles(Vars US, Vars DS){

  std::vector<double> projTheta;
  double USnorm = 1./sqrt(1 + US.dXdz*US.dXdz + US.dYdz*US.dYdz);
  double u[3] = {US.dXdz*USnorm, US.dYdz*USnorm, USnorm};
  double w[3] = {-u[0]*u[1], (u[0]*u[0] + u[2]*u[2]), -u[1]*u[2]};
  double Wnorm  = 1./sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
  w[0] *= Wnorm;
  w[1] *= Wnorm;
  w[2] *= Wnorm;
  double DSnorm = 1./sqrt(1 + DS.dXdz*DS.dXdz + DS.dYdz*DS.dYdz);
  double d[3] = {DS.dXdz*DSnorm, DS.dYdz*DSnorm, DSnorm};
  projTheta.push_back( atan( (d[0]*w[0] + d[1]*w[1] + d[2]*w[2])/
			     (d[0]*u[0] + d[1]*u[1] + d[2]*u[2]) ));
  projTheta.push_back( atan( (d[0]*u[2] - u[0]*d[2])\
			     /(d[0]*u[0] + d[1]*u[1] + d[2]*d[2]) * 
			     1./sqrt(u[2]*u[2] + u[0]*u[0])) );
  // projTheta.push_back( atan( (d[0]*u[2] - u[0]*d[2])			\
  ///(d[0]*u[0] + d[1]*u[1] + d[2]*u[2]) * 
  //			     1./sqrt(u[0]*u[0] + u[2]*u[2])) );

  projTheta.push_back( acos( ( (1 + US.dXdz * DS.dXdz + US.dYdz * DS.dYdz )/
			       sqrt(1 + US.dXdz*US.dXdz + US.dYdz*US.dYdz)/
			       sqrt(1 + DS.dXdz*DS.dXdz + DS.dYdz*DS.dYdz))) );
 
  return projTheta;
}

void MCSAnalysis::ConvolveWithInputDistribution(std::string distname){
  int failcount = 0;
  int passcount = 0;
  int isfirst = 0;
  if (distname.find(modelname1.c_str()) != std::string::npos)
    isfirst = 1;

  std::cerr<<"Model file: "<<modelfile<<std::endl;
  TFile* infile = new TFile(modelfile.c_str());
  
  std::string tmpname = "thetaX_refconv_";
  tmpname += distname;
  TH1D* thetaX_refconv = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

  tmpname = "thetaY_refconv_";
  tmpname += distname;
  TH1D* thetaY_refconv = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

  tmpname = "thetaScatt_refconv_";
  tmpname += distname;
  TH1D* thetaScatt_refconv = 
    new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta_{Scatt}; Events per mrad",
	     _histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);


  tmpname = "theta2Scatt_refconv_";
  tmpname += distname;
  TH1D* theta2Scatt_refconv = 
    new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta^{2}_{Scatt}; Events per mrad",
	     _histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);

  tmpname = "thetaScatt_refconv_vp_";
  tmpname += distname;
  TH2D* thetaScatt_refconv_vp = 
    new TH2D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;Momentum (MeV/c); #theta_{Scatt}", 
	     200, 100, 300, _histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
    
  TH2D* thetaXUS_thetaXDS = 
    new TH2D("thetaXUS_thetaXDS","Upstream vs. Downstream Angle;#theta_{X}^{US}; #theta_{X}^{DS}",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"], 
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  TH2D* thetaYUS_thetaYDS = 
    new TH2D("thetaYUS_thetaYDS","Upstream vs. Downstream Angle;#theta_{X}^{US}; #theta_{X}^{DS}",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"],
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

  tmpname = "thetaX_";
  tmpname += distname;

  std::cerr<<"Convolution with "<<tmpname<<" from "<<modelfile<<std::endl;
  TH1D* hx = (TH1D*)infile->Get(tmpname.c_str());
  TH1D* hy = (TH1D*)infile->Get(tmpname.c_str());
  std::cerr<<"hx,hy histograms created. Number of entries in _USMCset"<< _USMCset.N() << std::endl;
  // Collection DSConvSet;
  // for (int l=0; l<10; l++){
  for (int i=0; i<_USMCset.N(); i++){
    std::cerr<<"_USMCset entry: "<< i << std::endl;
    for (int j=0; j<20; j++){
      double dthetaX = hx->GetRandom() * _sys["resX"];
      double dthetaY = hy->GetRandom() * _sys["resY"];
      // First project the upstream track to the absorber 
      double zabspos = _sys["abspos"] + 0.0;
      if (abs(_USMCset.E(i).Z-zabspos)>.1){
//       std::cerr<<"ConvolveWithInputDistribution - _USMCset.E(i) not at absorber: "<<_USMCset.E(i).Z<<std::endl;
      }
      Vars projvarAbs;
      if (abs(_USMCset.E(i).Z-zabspos)>1.0){
         std::cerr<<"Reference data not a centre of absorber z="<<_USMCset.E(i).Z<<" Centre of Absorber="<<zabspos<<std::endl;
         projvarAbs = PropagateVarsMu(_USMCset.E(i), zabspos);
      }else{
         projvarAbs=_USMCset.E(i);
      }
      double xabs = projvarAbs.X;  /// _USMCset.E(i).X + _USMCset.E(i).dXdz * dzabsUS;
      double yabs = projvarAbs.Y;  /// _USMCset.E(i).Y + _USMCset.E(i).dYdz * dzabsUS;
      // Now add the angle from the model to the downstream measurement.
      double dXdz_abs = _DSMCset.E(i).dXdz + tan(dthetaY);
      double dYdz_abs = _DSMCset.E(i).dYdz + tan(dthetaX);
      // double d_thetaY    = atan(dXdz_abs) - atan(_DSMCset.E(i).dXdz);
      // double d_thetaX    = atan(dYdz_abs) - atan(_DSMCset.E(i).dYdz);
      
      // Project the track into the downstream reference plane
      // double xref = xabs + dXdz_abs * dzabsDS;
      // double yref = yabs + dYdz_abs * dzabsDS;
      Vars tmpvar = _DSMCset.E(i);
      tmpvar.X = xabs;
      tmpvar.Y = yabs;
      tmpvar.Z = _sys["abspos"];
      tmpvar.dXdz = dXdz_abs;
      tmpvar.dYdz = dYdz_abs;
      tmpvar.px   = dXdz_abs * _USMCset.E(i).pz;
      tmpvar.py   = dYdz_abs * _USMCset.E(i).pz;
      tmpvar.pz   = _USMCset.E(i).pz;
      tmpvar.TOF12= _USMCset.E(i).TOF12;
      tmpvar.TOF01= _USMCset.E(i).TOF01;
      std::vector<double> projDTheta = DefineProjectionAngles(tmpvar, _DSMCset.E(i));
      double d_thetaX = projDTheta[0];
      double d_thetaY = projDTheta[1];

      //double cosdthetaScatt = ( (1 + dXdz_abs * _DSMCset.E(i).dXdz + dYdz_abs*_DSMCset.E(i).dYdz) / 
      //			sqrt( 1 + dXdz_abs*dXdz_abs + dYdz_abs*dYdz_abs) / 
      //			sqrt(1 + _DSMCset.E(i).dXdz*_DSMCset.E(i).dXdz + _DSMCset.E(i).dYdz * _DSMCset.E(i).dYdz) );
      double dthetaScatt = projDTheta[2];  /// acos(cosdthetaScatt);
      Vars projvar = PropagateVarsMu(tmpvar, _sys["abspos"] + 2993.05);
      // Remove events that do not pass through the tracker at its centre
      // double xtracker = xabs + dXdz_abs * (_sys["abspos"] + 549.95);
      // double ytracker = yabs + dYdz_abs * (dzabsDS + 549.95);
      
      if( sqrt(projvar.X*projvar.X + projvar.Y*projvar.Y) > meanp ){
        failcount++;
//        std::cerr<<"ConvolveWithInputDistribution - event failed projvar check. Passed: "<< passcount<<" Failed: "<<failcount<<std::endl;
	tmpvar.X = 0.0;
	tmpvar.Y = 0.0;
	tmpvar.Z = 0.0;
	tmpvar.dXdz = 1./2.;
	tmpvar.dYdz = 1./2.;
      } else {passcount++;}
      std::vector<double> projTheta = DefineProjectionAngles(_USMCset.E(i), tmpvar);
      double thetaY = projTheta[1];   /// atan(tmpvar.dXdz) - atan(_USMCset.E(i).dXdz);
      double thetaX = projTheta[0];   /// atan(tmpvar.dYdz) - atan(_USMCset.E(i).dYdz);
      // double cosScatt = ( (1 + _USMCset.E(i).dXdz * tmpvar.dXdz +
      //		   _USMCset.E(i).dYdz * tmpvar.dYdz )/
      //		  sqrt(1 + _USMCset.E(i).dXdz*_USMCset.E(i).dXdz +
      //		       _USMCset.E(i).dYdz*_USMCset.E(i).dYdz)/
      //		  sqrt(1 + tmpvar.dXdz*tmpvar.dXdz +
      //		       tmpvar.dYdz*tmpvar.dYdz));
      double thetaScatt = projTheta[2];  /// acos(cosScatt);
      thetaXUS_thetaXDS->Fill(atan(_USMCset.E(i).dXdz), atan(tmpvar.dXdz));
      thetaYUS_thetaYDS->Fill(atan(_USMCset.E(i).dYdz), atan(tmpvar.dYdz));
      thetaX_refconv->Fill(thetaX);
      thetaY_refconv->Fill(thetaY);
      thetaScatt_refconv->Fill(thetaScatt);
      theta2Scatt_refconv->Fill(thetaScatt*thetaScatt);
      thetaScatt_refconv_vp->Fill(_USMCset.E(i).pz, thetaScatt);
      
      isfirst == 1 ? tresp_thetaX.Fill(thetaX, d_thetaX) : resp_thetaX.Fill(thetaX, d_thetaX);
      isfirst == 1 ? tresp_thetaY.Fill(thetaY, d_thetaY) : resp_thetaY.Fill(thetaY, d_thetaY);
      isfirst == 1 ? tresp_thetaScatt.Fill(thetaScatt, dthetaScatt) : resp_thetaScatt.Fill(thetaScatt, dthetaScatt); 
      isfirst == 1 ? tresp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt) : resp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt); 
    }
  }

  

  outfile->cd();
  thetaXUS_thetaXDS->Write();
  thetaYUS_thetaYDS->Write();
  thetaX_refconv->Write();
  thetaY_refconv->Write();
  thetaScatt_refconv->Write();
  theta2Scatt_refconv->Write();
  thetaScatt_refconv_vp->Write();

}

void MCSAnalysis::DoUnfolding(){

  std::string  tmpname = "thetaX_recoGold";
  TH1D* thetaX_recoGold = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmpname = "thetaY_recoGold";
  TH1D* thetaY_recoGold = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmpname = "thetaScatt_recoGold";
  TH1D* thetaScatt_recoGold = 
    new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta_{Scatt}; Events per mrad",_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
  tmpname = "theta2Scatt_recoGold";
  TH1D* theta2Scatt_recoGold = 
    new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta^{2}_{Scatt}; Events per mrad",_histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);

  tmpname = "thetaX_measured";
  TH1D* thetaX_measured = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmpname = "thetaY_measured";
  TH1D* thetaY_measured = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmpname = "thetaScatt_measured";
  TH1D* thetaScatt_measured = 
    new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta_{Scatt}; Events per mrad",_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
  tmpname = "theta2Scatt_measured";
  TH1D* theta2Scatt_measured = 
    new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta^{2}_{Scatt}; Events per mrad",_histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);


  tmpname = "thetaX_response";
  TH1D* thetaX_response = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmpname = "thetaY_response";
  TH1D* thetaY_response = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmpname = "thetaScatt_response";
  TH1D* thetaScatt_response = 
    new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta_{Scatt}; Events per mrad",_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
  tmpname = "theta2Scatt_response";
  TH1D* theta2Scatt_response = 
    new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta^{2}_{Scatt}; Events per mrad",_histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);


  for(int i=0; i<_DSset.N(); i++){
    
    double thetaY = atan(_DSset.E(i).dXdz) - atan(_USset.E(i).dXdz);
    double thetaX = atan(_DSset.E(i).dYdz) - atan(_USset.E(i).dYdz);
    double cosScatt = ( (1 + _USset.E(i).dXdz * _DSset.E(i).dXdz +
			 _USset.E(i).dYdz * _DSset.E(i).dYdz )/
			sqrt(1 + _USset.E(i).dXdz*_USset.E(i).dXdz +
			     _USset.E(i).dYdz*_USset.E(i).dYdz)/
			sqrt(1 + _DSset.E(i).dXdz*_DSset.E(i).dXdz +
			     _DSset.E(i).dYdz*_DSset.E(i).dYdz));
    double thetaScatt = acos(cosScatt);

    thetaX_recoGold->Fill(thetaX);
    thetaY_recoGold->Fill(thetaY);
    thetaScatt_recoGold->Fill(thetaScatt);
    theta2Scatt_recoGold->Fill(thetaScatt*thetaScatt);

    thetaX_measured->Fill(thetaX);
    thetaY_measured->Fill(thetaY);
    thetaScatt_measured->Fill(thetaScatt);
    theta2Scatt_measured->Fill(thetaScatt*thetaScatt);

  }
  float* source_thetaX = new float[int(_histlimits["NbinsXY"])];
  float* source_thetaY = new float[int(_histlimits["NbinsXY"])];
  float* source_thetaScatt = new float[int(_histlimits["NbinsTh"])];
  float* source_theta2Scatt = new float[int(_histlimits["NbinsTh2"])];

  for(int i=0;i<_DSMCset.N(); i++){
    std::vector<double> projTheta = DefineProjectionAngles(_USMCset.E(i), _DSMCset.E(i));
    double thetaMCY = projTheta[1]; /// atan(_DSMCset.E(i).dXdz) - atan(_USMCset.E(i).dXdz);
    double thetaMCX = projTheta[0]; /// atan(_DSMCset.E(i).dYdz) - atan(_USMCset.E(i).dYdz);
    // double cosScattMC = ( (1 + _USMCset.E(i).dXdz * _DSMCset.E(i).dXdz +
    //_USMCset.E(i).dYdz * _DSMCset.E(i).dYdz )/
    //sqrt(1 + _USMCset.E(i).dXdz*_USMCset.E(i).dXdz +
    //_USMCset.E(i).dYdz*_USMCset.E(i).dYdz)/
    // sqrt(1 + _DSMCset.E(i).dXdz*_DSMCset.E(i).dXdz +
    // _DSMCset.E(i).dYdz*_DSMCset.E(i).dYdz));
    double thetaScattMC = projTheta[2];  /// acos(cosScattMC);
    thetaX_response->Fill(thetaMCX);
    thetaY_response->Fill(thetaMCY);
    thetaScatt_response->Fill(thetaScattMC);
    theta2Scatt_response->Fill(thetaScattMC*thetaScattMC);

  }
  
  float* response_thetaX = new float[int(_histlimits["NbinsXY"])];
  float* response_thetaY = new float[int(_histlimits["NbinsXY"])];
  float* response_thetaScatt = new float[int(_histlimits["NbinsTh"])];
  float* response_theta2Scatt = new float[int(_histlimits["NbinsTh2"])];

  for(int i=0; i<_histlimits["NbinsXY"]; i++){
    source_thetaX[i] = thetaX_recoGold->GetBinContent(i+1);
    source_thetaY[i] = thetaY_recoGold->GetBinContent(i+1);
    response_thetaX[i] = thetaX_response->GetBinContent(i+1);
    response_thetaY[i] = thetaY_response->GetBinContent(i+1);
  }
  for(int i=0; i<_histlimits["NbinsTh"]; i++){
    source_thetaScatt[i] = thetaScatt_recoGold->GetBinContent(i+1);
    response_thetaScatt[i] = thetaScatt_response->GetBinContent(i+1);
  }
  for(int i=0; i<_histlimits["NbinsTh2"]; i++){
    source_theta2Scatt[i] = theta2Scatt_recoGold->GetBinContent(i+1);
    response_theta2Scatt[i] = theta2Scatt_response->GetBinContent(i+1);
  }

  TSpectrum *sX = new TSpectrum();
  sX->Deconvolution(source_thetaX, response_thetaX, _histlimits["NbinsXY"], 20, 10, 0.0);

  TSpectrum *sY = new TSpectrum();
  sY->Deconvolution(source_thetaY, response_thetaY, _histlimits["NbinsXY"], 20, 10, 0.0);
  for(int i=int(_histlimits["NbinsXY"]*3./4.); 
      i<int(_histlimits["NbinsXY"]*(3./4. + 1./5.)); i++){
    thetaX_recoGold->SetBinContent(i-148, source_thetaX[i]);
    thetaY_recoGold->SetBinContent(i-148, source_thetaY[i]);
  }
  double maxX = thetaX_recoGold->GetMaximum();
  double maxY = thetaY_recoGold->GetMaximum();
  for(int i = 1; i <= _histlimits["NbinsXY"]; i++){
    thetaX_recoGold->SetBinContent(i, thetaX_recoGold->GetBinContent(i)/maxX);
    thetaY_recoGold->SetBinContent(i, thetaY_recoGold->GetBinContent(i)/maxY);
  }
  TSpectrum *sScatt = new TSpectrum();
  sScatt->Deconvolution(source_thetaScatt, response_thetaScatt, _histlimits["NbinsTh"], 20, 10, 0.0);
  for(int i=0; i<int(_histlimits["NbinsTh"]); i++){
    thetaScatt_recoGold->SetBinContent(i, source_thetaScatt[i]);
  }
  TSpectrum *sScatt2 = new TSpectrum();
  sScatt2->Deconvolution(source_theta2Scatt, response_theta2Scatt, _histlimits["NbinsTh2"], 20, 10, 0.0);
  
  for(int i=0; i<int(_histlimits["NbinsTh2"]*(3./4.  + 1./5.)); i++){
    theta2Scatt_recoGold->SetBinContent(i, source_theta2Scatt[i]);
  }
  // double maxX = thetaScatt_recoGold->GetMaximum();
  //for(int i = 1; i <= _histlimits["NbinsTh"]; i++){
  //  thetaScatt_recoGold->SetBinContent(i, thetaScatt_recoGold->GetBinContent(i)/maxX);
  // }

  outfile->cd();
  thetaX_recoGold->Write();
  thetaY_recoGold->Write();
  thetaScatt_recoGold->Write();
  theta2Scatt_recoGold->Write();

  thetaX_measured->Write();
  thetaY_measured->Write();
  thetaScatt_measured->Write();
  theta2Scatt_measured->Write();

  thetaX_response->Write();
  thetaY_response->Write();
  thetaScatt_response->Write();
  theta2Scatt_response->Write();
  

  delete thetaX_recoGold;
  delete thetaY_recoGold;
  delete thetaScatt_recoGold;
  delete theta2Scatt_recoGold;

  delete thetaX_measured;
  delete thetaY_measured;
  delete thetaScatt_measured;
  delete theta2Scatt_measured;

  delete thetaX_response;
  delete thetaY_response;
  delete thetaScatt_response;
  delete theta2Scatt_response;

}


void MCSAnalysis::DoDeconvolution(std::string model, int n_sel=1){
  // The basic methods associated with the RooUnfolding package
  // First generate a histogram of the measured data.
  
  int isfirst = 0;
  if (model.find(modelname1.c_str()) != std::string::npos)
    isfirst = 1;

  // int n_base = int(2 * _DSset.N() / (n_sel * (n_sel + 1)));
  // int k=0;
  // for (int j=0; j<n_sel; j++){
  // define the number of events under the current selection
  int curr_sel = int(_DSset.N());  // n_base * j;
  int curr_sel_empty = int(_DSMCset.N());
  int curr_k = 0;

  
    
  std::string  tmpname = "thetaX_measdata";
  tmpname += model;
  TH1D* thetaX_measdata = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmpname = "thetaY_measdata";
  tmpname += model;
  TH1D* thetaY_measdata = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#Delta#theta_{Y}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmpname = "thetaScatt_measdata";
  tmpname += model;
  TH1D* thetaScatt_measdata = 
    new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta_{Scatt}; Events per mrad",_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
  tmpname = "theta2Scatt_measdata";
  tmpname += model;
  TH1D* theta2Scatt_measdata = 
    new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta^{2}_{Scatt}; Events per mrad",_histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);
  tmpname = "thetaScatt_measdata_vp";
  tmpname += model;
  TH2D* thetaScatt_measdata_vp = 
    new TH2D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;Momentum (MeV/c); #theta_{Scatt}", 
	     400, 100, 300, _histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
  tmpname = "thetaXUS_thetaXDS";
  tmpname += model;
  
  TH2D* thetaXUS_thetaXDS = 
    new TH2D(tmpname.c_str(),"Upstream vs. Downstream Angle;#theta_{X}^{US}; #theta_{X}^{DS}",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"], _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmpname = "thetaYUS_thetaYDS";
  tmpname += model;
  TH2D* thetaYUS_thetaYDS = 
    new TH2D(tmpname.c_str(),"Upstream vs. Downstream Angle;#theta_{X}^{US}; #theta_{X}^{DS}",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"],_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  
  tmpname = "projposUSDSdiff_";
  tmpname += model;
  TH2D* projposUSDSdiff = 
    new TH2D(tmpname.c_str(),
	     "Difference of US and DS projections at absorber;#Delta x_{DS-US} (mm); #Delta y_{DS-US}",
	     180, -400, 400, 180, -400, 400);
  
  
  const Int_t NBINS = 19;
  Double_t scat_bin_array[NBINS + 1] = {-0.1151,-0.0938,-0.0754,-0.0597,-0.0463,-0.0347,-0.0248,-0.0162,-0.00895,-0.00269,0.00269,0.00895,0.0162,0.0248,0.0347,0.0463,0.0597,0.0754,0.0938,0.1151};
  TH1D* scattering_proj_x = new TH1D("scattering_proj_x_DC","Change in Projected Angle (X);#Delta#theta_{X}; Events per radian", 
				     NBINS, scat_bin_array);
  TH1D* scattering_proj_y = new TH1D("scattering_proj_y_DC","Change in Projected Angle (Y);#Delta#theta_{Y}; Events per radian", 
				     NBINS, scat_bin_array);

  TH1D* scattering_proj_x_empty = new TH1D("scattering_proj_x_empty","Change in Projected Angle (X);#Delta#theta_{X}; Events per radian", 
				     NBINS, scat_bin_array);
  TH1D* scattering_proj_y_empty = new TH1D("scattering_proj_y_empty","Change in Projected Angle (Y);#Delta#theta_{Y}; Events per radian", 
				     NBINS, scat_bin_array);


  
  for(int i=curr_k; i<curr_sel_empty; i++){
      std::vector<double> projTheta = DefineProjectionAngles(_USMCset.E(i), _DSMCset.E(i));
    double thetaY = projTheta[1];  /// atan(_DSset.E(i).dXdz) - atan(_USset.E(i).dXdz);
    double thetaX = projTheta[0];  /// atan(_DSset.E(i).dYdz) - atan(_USset.E(i).dYdz);
    scattering_proj_x_empty->Fill(thetaX);
    scattering_proj_y_empty->Fill(thetaY);
  }

  for(int i=curr_k; i<curr_sel; i++){
    std::vector<double> projTheta = DefineProjectionAngles(_USset.E(i), _DSset.E(i));
    double thetaY = projTheta[1];  /// atan(_DSset.E(i).dXdz) - atan(_USset.E(i).dXdz);
    double thetaX = projTheta[0];  /// atan(_DSset.E(i).dYdz) - atan(_USset.E(i).dYdz);
    // double cosScatt 
    double thetaScatt = projTheta[2];  /// acos(cosScatt);
    double xUSabs = _USset.E(i).X - _USset.E(i).dXdz * (_USset.E(i).Z - _sys["abspos"]);
    double yUSabs = _USset.E(i).Y - _USset.E(i).dYdz * (_USset.E(i).Z - _sys["abspos"]);
    double xDSabs = _DSset.E(i).X - _DSset.E(i).dXdz * (_DSset.E(i).Z - _sys["abspos"]);
    double yDSabs = _DSset.E(i).Y - _DSset.E(i).dYdz * (_DSset.E(i).Z - _sys["abspos"]);
    projposUSDSdiff->Fill(xDSabs - xUSabs, yDSabs - yUSabs);
    thetaXUS_thetaXDS->Fill(atan(_USset.E(i).dXdz), atan(_DSset.E(i).dXdz));
    thetaYUS_thetaYDS->Fill(atan(_USset.E(i).dYdz), atan(_DSset.E(i).dYdz));
    thetaX_measdata->Fill(thetaX);
    thetaY_measdata->Fill(thetaY);
    thetaScatt_measdata->Fill(thetaScatt);
    theta2Scatt_measdata->Fill(thetaScatt*thetaScatt);
    // double dt0 = 7.64186 / 0.299792458; // ns. 
    // double dt = _USset.E(i).TOF01;
    // if (_USset.E(i).TOF12 < 99.0 * 8.22475 / 0.299792458){
    //  dt0 = 8.22475 / 0.299792458; // ns.
    //  dt = _USset.E(i).TOF12;
    //}
    double pz  = _USset.E(i).pz;
    // 105.65 / sqrt(dt*dt/dt0/dt0 - 1.0);
    thetaScatt_measdata_vp->Fill(pz, thetaScatt);
    scattering_proj_x->Fill(thetaX);
    scattering_proj_y->Fill(thetaY);
    // k++;
  }
  /*
    for ( int itr=1; itr<25; itr++){    
    RooUnfoldBayes unfold_thetaX(isfirst==1 ? &tresp_thetaX : &resp_thetaX, thetaX_measdata, itr);
    unfold_thetaX.Print();
    char res[16];
    sprintf( res, "%d", itr);
    TH1D* thetaX_reco = (TH1D*)unfold_thetaX.Hreco();
    tmpname = "thetaX_reco";
    tmpname += model;
    tmpname += "_";
    tmpname += res;
    thetaX_reco->SetName(tmpname.c_str());
    thetaX_reco->SetTitle(";#Delta #theta_{X}");
    outfile->cd();
    thetaX_reco->Write();
    
    delete thetaX_reco;
    }
  */
  RooUnfoldBayes unfold_thetaX(isfirst==1 ? &tresp_thetaX : &resp_thetaX, thetaX_measdata, int(_sys["niter"]));
  unfold_thetaX.Print();
  RooUnfoldBayes unfold_thetaY(isfirst==1 ? &tresp_thetaY : &resp_thetaY, thetaY_measdata, int(_sys["niter"]));
  unfold_thetaY.Print();
  RooUnfoldBayes 
    unfold_thetaScatt(isfirst==1 ? &tresp_thetaScatt : &resp_thetaScatt, thetaScatt_measdata, int(_sys["niter"]));
  unfold_thetaScatt.Print();
  RooUnfoldBayes 
    unfold_theta2Scatt(isfirst==1 ? &tresp_theta2Scatt : &resp_theta2Scatt, theta2Scatt_measdata, int(_sys["niter"]));
  unfold_theta2Scatt.Print();
  
  for (int i=1; i<=19; i++){
    double x_content = scattering_proj_x->GetBinContent(i);
    scattering_proj_x->Fill(i, x_content / scat_bin_array[i-1]);
    double y_content = scattering_proj_y->GetBinContent(i);
    scattering_proj_y->Fill(i, y_content / scat_bin_array[i-1]);
    double x_content_empty = scattering_proj_x_empty->GetBinContent(i);
    scattering_proj_x_empty->Fill(i, x_content_empty / scat_bin_array[i-1]);
    double y_content_empty = scattering_proj_y_empty->GetBinContent(i);
    scattering_proj_y_empty->Fill(i, y_content_empty / scat_bin_array[i-1]);

  }
  
  TH1D* thetaX_reco = (TH1D*)unfold_thetaX.Hreco();
  tmpname = "thetaX_reco";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaX_reco->SetName(tmpname.c_str());
  thetaX_reco->SetTitle(";#Delta #theta_{X}");
  TH1D* thetaY_reco = (TH1D*)unfold_thetaY.Hreco();
  tmpname = "thetaY_reco";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaY_reco->SetName(tmpname.c_str());
  thetaY_reco->SetTitle(";#Delta #theta_{Y}");
  TH1D* thetaScatt_reco = (TH1D*)unfold_thetaScatt.Hreco();
  tmpname = "thetaScatt_reco";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaScatt_reco->SetName(tmpname.c_str());
  thetaScatt_reco->SetTitle(";#theta_{Scatt}");
  TH1D* theta2Scatt_reco = (TH1D*)unfold_theta2Scatt.Hreco();
  tmpname = "theta2Scatt_reco";
  tmpname += model;
  // if(j>0) tmpname += j;
  theta2Scatt_reco->SetName(tmpname.c_str());
  theta2Scatt_reco->SetTitle(";#theta^{2}_{Scatt}");
  
  TH2D* thetaX_response = (TH2D*)unfold_thetaX.response()->Hresponse();
  tmpname = "thetaX_response";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaX_response->SetName(tmpname.c_str());
  thetaX_response->SetTitle(";#Delta #theta^{rec}_{X}; #Delta #theta^{true}_{X}");
  TH2D* thetaY_response = (TH2D*)unfold_thetaY.response()->Hresponse();
  tmpname = "thetaY_response";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaY_response->SetName(tmpname.c_str());
  thetaY_response->SetTitle(";#Delta #theta^{rec}_{Y}; #Delta #theta^{true}_{Y}");
  TH2D* thetaScatt_response = (TH2D*)unfold_thetaScatt.response()->Hresponse();
  tmpname = "thetaScatt_response";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaScatt_response->SetName(tmpname.c_str());
  thetaScatt_response->SetTitle(";#theta^{tracker}_{Scatt};#theta^{abs}_{Scatt}");
  TH2D* theta2Scatt_response = (TH2D*)unfold_theta2Scatt.response()->Hresponse();
  tmpname = "theta2Scatt_response";
  tmpname += model;
  // if(j>0) tmpname += j;
  theta2Scatt_response->SetName(tmpname.c_str());
  theta2Scatt_response->SetTitle(";(#theta^{tracker}_{Scatt})^{2}; (#theta^{abs}_{Scatt})^{2}");
  
  TH1D* thetaX_truth = (TH1D*)unfold_thetaX.response()->Htruth();
  tmpname = "thetaX_";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaX_truth->SetName(tmpname.c_str());
  thetaX_truth->SetTitle(";#Delta #theta_{X}");
  TH1D* thetaY_truth = (TH1D*)unfold_thetaY.response()->Htruth();
  tmpname = "thetaY_";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaY_truth->SetName(tmpname.c_str());
  thetaX_truth->SetTitle(";#Delta #theta_{Y}");
  TH1D* thetaScatt_truth = (TH1D*)unfold_thetaScatt.response()->Htruth();
  tmpname = "thetaScatt_";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaScatt_truth->SetName(tmpname.c_str());
  thetaScatt_truth->SetTitle(";#theta_{Scatt}");
  TH1D* theta2Scatt_truth = (TH1D*)unfold_theta2Scatt.response()->Htruth();
  tmpname = "theta2Scatt_";
  tmpname += model;
  // if(j>0) tmpname += j;
  theta2Scatt_truth->SetName(tmpname.c_str());
  theta2Scatt_truth->SetTitle(";#theta^{2}_{Scatt}");
  
  TH1D* thetaX_measured = (TH1D*)unfold_thetaX.response()->Hmeasured();
  tmpname = "thetaX_measured";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaX_measured->SetName(tmpname.c_str());
  thetaX_measured->SetTitle(";#Delta #theta_{X}");
  TH1D* thetaY_measured = (TH1D*)unfold_thetaY.response()->Hmeasured();
  tmpname = "thetaY_measured";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaY_measured->SetName(tmpname.c_str());
  thetaY_measured->SetTitle(";#Delta #theta_{Y}");
  TH1D* thetaScatt_measured = (TH1D*)unfold_thetaScatt.response()->Hmeasured();
  tmpname = "thetaScatt_measured";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaScatt_measured->SetName(tmpname.c_str());
  thetaScatt_measured->SetTitle(";#theta_{Scatt}");
  TH1D* theta2Scatt_measured = (TH1D*)unfold_theta2Scatt.response()->Hmeasured();
  tmpname = "theta2Scatt_measured";
  tmpname += model;
  // if(j>0) tmpname += j;
  theta2Scatt_measured->SetName(tmpname.c_str());
  theta2Scatt_measured->SetTitle(";#theta^{2}_{Scatt}");
  
  outfile->cd();
  thetaXUS_thetaXDS->Write();
  thetaYUS_thetaYDS->Write();
  thetaX_measdata->Write();
  thetaY_measdata->Write();
  thetaScatt_measdata->Write();
  theta2Scatt_measdata->Write();
  thetaScatt_measdata_vp->Write();
  thetaX_reco->Write();
  thetaY_reco->Write();
  thetaScatt_reco->Write();
  theta2Scatt_reco->Write();
  thetaX_response->Write();
  thetaY_response->Write();
  thetaScatt_response->Write();
  theta2Scatt_response->Write();
  thetaX_truth->Write();
  thetaY_truth->Write();
  thetaScatt_truth->Write();
  theta2Scatt_truth->Write();
  thetaX_measured->Write();
  thetaY_measured->Write();
  thetaScatt_measured->Write();
  theta2Scatt_measured->Write();
  scattering_proj_x->Write();
  scattering_proj_y->Write();
  scattering_proj_x_empty->Write();
  scattering_proj_y_empty->Write();
  projposUSDSdiff->Write();
  
  TCanvas* c1 = new TCanvas();
  thetaX_measdata->Draw();
  c1->Print("thetaX_measdata.pdf");
  c1->Clear();
  thetaY_measdata->Draw();
  c1->Print("thetaY_measdata.pdf");
  c1->Clear();
  thetaScatt_measdata->Draw();
  c1->Print("thetaScatt_measdata.pdf");
  c1->Clear();
  thetaScatt_measdata_vp->Draw("colz");
  c1->Print("thetaScatt_measdata_vp.pdf");
  c1->Clear();
  thetaX_reco->Draw();
  c1->Print("thetaX_reco.pdf");
  c1->Clear();
  thetaY_reco->Draw();
  c1->Print("thetaY_reco.pdf");
  c1->Clear();
  thetaScatt_reco->Draw();
  c1->Print("thetaScatt_reco.pdf");
  c1->Clear();
  thetaX_response->Draw();
  c1->Print("thetaX_response.pdf");
  c1->Clear();
  thetaY_response->Draw();
  c1->Print("thetaY_response.pdf");
  c1->Clear();
  thetaScatt_response->Draw();
  c1->Print("thetaScatt_response.pdf");
  c1->Clear();
  thetaX_measured->Draw();
  c1->Print("thetaX_measured.pdf");
  c1->Clear();
  thetaY_measured->Draw();
  c1->Print("thetaY_measured.pdf");
  c1->Clear();
  thetaScatt_measured->Draw();
  c1->Print("thetaScatt_measured.pdf");
  c1->Clear();
  thetaX_truth->Draw();
  c1->Print("thetaX_truth.pdf");
  c1->Clear();
  thetaY_truth->Draw();
  c1->Print("thetaY_truth.pdf");
  c1->Clear();
  thetaScatt_truth->Draw();
  c1->Print("thetaScatt_truth.pdf");
  c1->Clear();
  c1->SetLogy();
  scattering_proj_x->Draw();
  c1->Print("scattering_proj_x.pdf");
  c1->Clear();
  scattering_proj_y->Draw();
  c1->Print("scattering_proj_y.pdf");
  c1->Clear();
  
  delete c1;
  delete thetaX_measdata;
  delete thetaY_measdata;
  delete thetaScatt_measdata;
  delete thetaScatt_measdata_vp;
  delete thetaX_reco;
  delete thetaY_reco;
  delete thetaScatt_reco;
  delete thetaX_response;
  delete thetaY_response;
  delete thetaScatt_response;
  delete thetaX_truth;
  delete thetaY_truth;
  delete thetaScatt_truth;
  delete thetaX_measured;
  delete thetaY_measured;
  delete thetaScatt_measured;
  delete scattering_proj_x;
  delete scattering_proj_y;
  delete projposUSDSdiff;

}

void MCSAnalysis::FitGaussian(std::string outfilename){

	TFile *MyFile = new TFile(outfilename.c_str());
	TH1F * h1 = (TH1F*)MyFile->Get("thetaX_measdataCobb");
	h1->SetDirectory(0);
	h1->Draw();
	std::cout << "Fit result for thetaX_measdata " << std::endl;
	h1->Fit("gaus","","",-0.040,0.040);

	TH1F * h2 = (TH1F*)MyFile->Get("thetaY_measdataCobb");
	h2->SetDirectory(0);
	h2->Draw();
	std::cout << "Fit result for thetaY_measdata " << std::endl;
	h2->Fit("gaus","","",-0.040,0.040);

	TH1F * h3 = (TH1F*)MyFile->Get("thetaX_recoCobb");
	h3->SetDirectory(0);
	h3->Draw();
	std::cout << "Fit result for thetaX_reco " << std::endl;
	h3->Fit("gaus","","",-0.040,0.040);
	
	TH1F * h4 = (TH1F*)MyFile->Get("thetaY_recoCobb");
	h4->SetDirectory(0);
	h4->Draw();
	std::cout << "Fit result for thetaY_reco " << std::endl;
	h4->Fit("gaus","","",-0.040,0.040);
	
	TH1F * h5 = (TH1F*)MyFile->Get("thetaX_refconv_Cobb");
	h5->SetDirectory(0);
	h5->Draw();
	std::cout << "Fit result for thetaX_refconv_Cobb " << std::endl;
	h5->Fit("gaus","","",-0.040,0.040);
	
	TH1F * h6 = (TH1F*)MyFile->Get("thetaY_refconv_Cobb");
	h6->SetDirectory(0);
	h6->Draw();
	std::cout << "Fit result for thetaY_refconv_Cobb " << std::endl;
	h6->Fit("gaus","","",-0.040,0.040);
	
	MyFile->Close();
	delete h1;
	delete h2;
	delete h3;
	delete h4;
	delete h5;
	delete h6;
}

void MCSAnalysis::CalculateChi2(std::string outfilename, std::string distname){

	std::string tmpname = "thetaX_refconv_";
        tmpname += distname;
	TFile *MyFile = new TFile(outfilename.c_str());
        TH1F * hx = (TH1F*)MyFile->Get(tmpname.c_str());
	tmpname = "thetaX_measdata";
        tmpname += distname;
	TH1F * h1 = (TH1F*)MyFile->Get(tmpname.c_str());
	h1->SetDirectory(0);
	hx->SetDirectory(0);
	//hx->Rebin(GetNbinsX/40);
	h1->Draw();
	hx->Draw();
	std::cout << "Chi^2 for thetaX_reco and " << distname << std::endl;
	h1->Chi2Test(hx,"P");
	
        tmpname = "thetaY_refconv_";
        tmpname += distname;
	TH1F * hy = (TH1F*)MyFile->Get(tmpname.c_str());
        tmpname = "thetaY_measdata";
        tmpname += distname;
	TH1F * h2 = (TH1F*)MyFile->Get(tmpname.c_str());
	h2->SetDirectory(0);
	hy->SetDirectory(0);
	//hy->Rebin(GetNbinsX/40);
	h2->Draw();
	hy->Draw();
	std::cout << "Chi^2 for thetaY_reco and " << distname << std::endl;
	h2->Chi2Test(hy,"P");
	
	delete hx;
	delete h1;
	delete hy;
	delete h2;
        MyFile->Close();

}

void MCSAnalysis::DoFFTDeconvolution(){
  
  // The basic methods associated with the RooUnfolding package
  // First generate a histogram of the measured data.
  TH1D* thetaX_data = 
    new TH1D("thetaX_data","Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  TH1D* thetaY_data = 
    new TH1D("thetaY_data","Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  TH1D* thetaScatt_data = 
    new TH1D("thetaScatt_data","Scattering Angle between Momentum Vectors;#theta_{Scatt}; Events per mrad",_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
  TH1D* theta2Scatt_data = 
    new TH1D("theta2Scatt_data","Scattering Angle between Momentum Vectors;#theta_{Scatt}^{2}; Events per mrad",_histlimits["NbinsTh"], _histlimits["minTh2"], _histlimits["maxTh2"]);
  
  TH1D* thetaX_ref = 
    new TH1D("thetaX_ref","Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  TH1D* thetaY_ref = 
    new TH1D("thetaY_ref","Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  TH1D* thetaScatt_ref = 
    new TH1D("thetaScatt_ref","Scattering Angle between Momentum Vectors;#theta_{Scatt}; Events per mrad",_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
  TH1D* theta2Scatt_ref = 
    new TH1D("theta2Scatt_ref","Scattering Angle between Momentum Vectors;#theta_{Scatt}^{2}; Events per mrad",_histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);

  /*
  TH1D* thetaX_fft = 
    new TH1D("thetaX_fft","Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  TH1D* thetaY_fft = 
    new TH1D("thetaY_fft","Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  TH1D* thetaScatt_fft = 
    new TH1D("thetaScatt_fft","Scattering Angle between Momentum Vectors;#theta_{Scatt}; Events per mrad",_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);

  
  TH1D* fftX_data = 
    new TH1D("fftX_data","Change in Projected Angle (X);mode_{X}; Events per 4 mrad",400,0, 400);
  TH1D* fftY_data = 
    new TH1D("fftY_data","Change in Projected Angle (Y);mode_{Y}; Events per 4 mrad",400,0, 400);
  TH1D* fftScatt_dataRe = 
    new TH1D("fftScatt_dataRe","Scattering Angle between Momentum Vectors;mode_{Scatt}; Events per mrad",400,0,400);
  TH1D* fftScatt_dataIm = 
    new TH1D("fftScatt_dataIm","Scattering Angle between Momentum Vectors;mode_{Scatt}; Events per mrad",400,0,400);
  
  TH1D* fftX_ref = 
    new TH1D("fftX_ref","Change in Projected Angle (X);mode_{X}; Events per 4 mrad",400,0, 400);
  TH1D* fftY_ref = 
    new TH1D("fftY_ref","Change in Projected Angle (Y);mode_{Y}; Events per 4 mrad",400,0, 400);
  TH1D* fftScatt_refRe = 
    new TH1D("fftScatt_refRe","Scattering Angle between Momentum Vectors;mode_{Scatt}; Events per mrad",400,0, 400);
  TH1D* fftScatt_refIm = 
    new TH1D("fftScatt_refIm","Scattering Angle between Momentum Vectors;mode_{Scatt}; Events per mrad",400,0, 400);

  TH1D* fftX_div = 
    new TH1D("fftX_div","Change in Projected Angle (X);mode_{X}; Events per 4 mrad",400, 0, 400);
  TH1D* fftY_div = 
    new TH1D("fftY_div","Change in Projected Angle (Y);mode_{Y}; Events per 4 mrad",400, 0, 400);
  TH1D* fftScatt_divRe = 
    new TH1D("fftScatt_divRe","Scattering Angle between Momentum Vectors;mode_{Scatt}; Events per mrad",400, 0, 400);
  TH1D* fftScatt_divIm = 
    new TH1D("fftScatt_divIm","Scattering Angle between Momentum Vectors;mode_{Scatt}; Events per mrad",400, 0, 400);
  */
  TH2D* thetaScatt_measdata_vp = 
    new TH2D("thetaScatt_measdata_vp","Scattering Angle between Momentum Vectors;Momentum (MeV/c); #theta_{Scatt}", 
	     400, 100, 300, _histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
  
  const Int_t NBINS = 19;
  Double_t scat_bin_array[NBINS + 1] = {-0.1151,-0.0938,-0.0754,-0.0597,-0.0463,-0.0347,-0.0248,-0.0162,-0.00895,-0.00269,0.00269,0.00895,0.0162,0.0248,0.0347,0.0463,0.0597,0.0754,0.0938,0.1151};
  TH1D* scattering_proj_x_FFT= new TH1D("scattering_proj_x_FFT_DC","Change in Projected Angle (X);#Delta#theta_{X}; Events per radian", 
				     NBINS, scat_bin_array);
  TH1D* scattering_proj_y_FFT = new TH1D("scattering_proj_y_FFT_DC","Change in Projected Angle (Y);#Delta#theta_{Y}; Events per radian", 
				     NBINS, scat_bin_array);
  
  TH2D* projposUSDSdiff = 
    new TH2D("projposUSDSdiff_ref",
	     "Difference of US and DS projections at absorber;#Delta x_{DS-US} (mm); #Delta y_{DS-US}",
	     180, -400, 400, 180, -400, 400);

  for(int i=0; i<_DSset.N(); i++){
    std::vector<double> projTheta = DefineProjectionAngles(_USset.E(i), _DSset.E(i));
    double thetaY = projTheta[1];  // atan(_DSset.E(i).dXdz) - atan(_USset.E(i).dXdz);
    double thetaX = projTheta[0];  // atan(_DSset.E(i).dYdz) - atan(_USset.E(i).dYdz);
    // double cosScatt = ( (1 + _USset.E(i).dXdz * _DSset.E(i).dXdz +
    //_USset.E(i).dYdz * _DSset.E(i).dYdz )/
    //			sqrt(1 + _USset.E(i).dXdz*_USset.E(i).dXdz +
    //			      _USset.E(i).dYdz*_USset.E(i).dYdz)/
    //			sqrt(1 + _DSset.E(i).dXdz*_DSset.E(i).dXdz +
    //			     _DSset.E(i).dYdz*_DSset.E(i).dYdz));
    double thetaScatt = projTheta[2];  /// acos(cosScatt);
    
    thetaX_data->Fill(thetaX);
    thetaY_data->Fill(thetaY);
    thetaScatt_data->Fill(thetaScatt);
    theta2Scatt_data->Fill(thetaScatt*thetaScatt);
    /*
      double dt0 = 7.64186 / 0.299792458; // ns. 
    double dt = _USset.E(i).TOF01;
    if (_USset.E(i).TOF12 < 99.0 * 8.22475 / 0.299792458){
      dt0 = 8.22475 / 0.299792458; // ns.
      dt = _USset.E(i).TOF12;
    }
    */
    double pz  = _USset.E(i).pz;
    // 105.65 / sqrt(dt*dt/dt0/dt0 - 1.0);
    thetaScatt_measdata_vp->Fill(pz, thetaScatt);
    scattering_proj_x_FFT->Fill(thetaX);
    scattering_proj_y_FFT->Fill(thetaY);
  }

  for(int i=0; i<_DSMCset.N(); i++){
    
    double thetaY = atan(_DSMCset.E(i).dXdz) - atan(_USMCset.E(i).dXdz);
    double thetaX = atan(_DSMCset.E(i).dYdz) - atan(_USMCset.E(i).dYdz);
    double cosScatt = ( (1 + _USMCset.E(i).dXdz * _DSMCset.E(i).dXdz +
			 _USMCset.E(i).dYdz * _DSMCset.E(i).dYdz )/
			sqrt(1 + _USMCset.E(i).dXdz*_USMCset.E(i).dXdz +
			      _USMCset.E(i).dYdz*_USMCset.E(i).dYdz)/
			sqrt(1 + _DSMCset.E(i).dXdz*_DSMCset.E(i).dXdz +
			     _DSMCset.E(i).dYdz*_DSMCset.E(i).dYdz));
    double thetaScatt = acos(cosScatt);
    
    double xUSabs = _USMCset.E(i).X - _USMCset.E(i).dXdz * (_USMCset.E(i).Z - _sys["abspos"]);
    double yUSabs = _USMCset.E(i).Y - _USMCset.E(i).dYdz * (_USMCset.E(i).Z - _sys["abspos"]);
    double xDSabs = _DSMCset.E(i).X - _DSMCset.E(i).dXdz * (_DSMCset.E(i).Z - _sys["abspos"]);
    double yDSabs = _DSMCset.E(i).Y - _DSMCset.E(i).dYdz * (_DSMCset.E(i).Z - _sys["abspos"]);
    projposUSDSdiff->Fill(xDSabs - xUSabs, yDSabs - yUSabs);

    thetaX_ref->Fill(thetaX);
    thetaY_ref->Fill(thetaY);
    thetaScatt_ref->Fill(thetaScatt);
    theta2Scatt_ref->Fill(thetaScatt*thetaScatt);
  }
  
  /*
  int NXFFT=thetaX_data->Integral();
  thetaX_data->Sumw2();
  thetaX_data->Scale(1./NXFFT);
  thetaX_data->FFT(fftX_data,"Re R2R_3");
  int NYFFT=thetaY_data->Integral();
  thetaY_data->Sumw2();
  thetaY_data->Scale(1./NYFFT);
  thetaY_data->FFT(fftY_data,"Re R2R_3");
  int NScattFFT=thetaScatt_data->Integral();
  thetaScatt_data->Sumw2();
  thetaScatt_data->Scale(1./NScattFFT);
  thetaScatt_data->FFT(fftScatt_dataRe,"Re R2C P");
  thetaScatt_data->FFT(fftScatt_dataIm,"Im R2C P");
    
  for (int i=1; i<fftX_data->GetNbinsX() + 1; i++){
    double c = fftX_data->GetBinContent(i);
    if (i > 40 ) c = 0.;
    fftX_data->SetBinContent(i, fabs(c));
  }
  for (int i=1; i<fftY_data->GetNbinsX() + 1; i++){
    double c = fftY_data->GetBinContent(i);
    if (i > 40 ) c = 0.;
    fftY_data->SetBinContent(i, fabs(c));
  }
  // Supress modes higher than 10 if modes go to negative values.
  // bool iszero=false;
  //for (int i=1; i<fftScatt_data->GetNbinsX() + 1; i++){
  //  if(fftScatt_data->GetBinContent(i) < 1e-3)
      // iszero=true;
      // if(iszero)
  //    fftScatt_data->SetBinContent(i, 0.0);
  //}
  // Smooth the data
  
  TH1D* thetaScatt_smooth = (TH1D*)thetaScatt_data->Clone();
  thetaScatt_smooth->SetName("thetaScatt_smooth");
  fftScatt_data->FFT(thetaScatt_smooth,"Re R2R_7");
  thetaScatt_smooth->Scale(1./2./thetaScatt_data->GetNbinsX());
  thetaScatt_smooth->FFT(fftScatt_data,"Re R2R_7");
  
  NXFFT=thetaX_ref->Integral();
  thetaX_ref->Scale(1./NXFFT);
  thetaX_ref->FFT(fftX_ref,"Re R2R_3");
  NYFFT=thetaY_ref->Integral();
  thetaY_ref->Scale(1./NYFFT);
  thetaY_ref->FFT(fftY_ref,"Re R2R_3");
  NScattFFT=thetaScatt_ref->Integral();
  thetaScatt_ref->Scale(1./NScattFFT);
  thetaScatt_ref->FFT(fftScatt_refRe,"Re R2C P");
  thetaScatt_ref->FFT(fftScatt_refIm,"Im R2C P");
  

  // Correct the normalization after the FFT
  // fftX_ref->Scale(1./2./fftX_ref->GetNbinsX());
  // fftY_ref->Scale(1./2./fftY_ref->GetNbinsX());
  // fftScatt_ref->Scale(1./2./fftScatt_ref->GetNbinsX());
  
  for (int i=1; i<fftX_ref->GetNbinsX() + 1; i++){
    double c = fftX_ref->GetBinContent(i);
    if (i > 40 ) c = 0.;
    fftX_ref->SetBinContent(i, fabs(c));
  }
  for (int i=1; i<fftY_ref->GetNbinsX() + 1; i++){
    double c = fftY_ref->GetBinContent(i);
    if (i > 40 ) c = 0.;
    fftY_ref->SetBinContent(i, fabs(c));
  }
  
  // Supress modes higher than 10 if modes go to negative values.
  // iszero=false;
  
  double *re_full = new double[fftScatt_refRe->GetNbinsX()];
  double *im_full = new double[fftScatt_refRe->GetNbinsX()];
  for (int i=1; i<fftScatt_refRe->GetNbinsX() + 1; i++){
    // if(fftScatt_ref->GetBinContent(i) < 1e-3)
      // iszero=true;
      // if(iszero)
    //  fftScatt_data->SetBinContent(i, 0.0);
    double f0 = fftScatt_dataRe->GetBinContent(i);
    double g0 = fftScatt_dataIm->GetBinContent(i);
    double f1 = fftScatt_refRe->GetBinContent(i);
    double g1 = fftScatt_refIm->GetBinContent(i);
    double mag2 = f1*f1 + g1*g1;
    double rediv = (f1*f0 + g1*g0)/mag2;
    double imdiv = (g0*f1 - f0*g1)/mag2;
    if (mag2 < 1e-4 || fabs(f0) < 0.001 || fabs(g0) < 0.001 || fabs(f1) < 0.001 || fabs(g1) < 0.001){
      rediv = 0.0;
      imdiv = 0.0;
    }
    if ( i > binlimit && i < fftScatt_refRe->GetNbinsX()-binlimit ){	
      rediv = 0.0;
      imdiv = 0.0;
    }
    fftScatt_divRe->SetBinContent(i,rediv);
    fftScatt_divIm->SetBinContent(i,imdiv);
    re_full[i-1] = rediv;
    im_full[i-1] = imdiv;
  }
  int n = fftScatt_divRe->GetNbinsX();
  TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &n, "C2R M K");
  fft_back->SetPointsComplex(re_full,im_full);
  fft_back->Transform();
  TH1 *hb = 0;
  hb = TH1::TransformHisto(fft_back, hb, "Re");
  for (int i=1; i<fftScatt_refRe->GetNbinsX() + 1;i++){
    thetaScatt_fft->SetBinContent(i,hb->GetBinContent(i));
  }
  
  thetaScatt_smooth->SetName("thetaScatt_smooth");
  fftScatt_ref->FFT(thetaScatt_smooth,"Re R2R_7");
  thetaScatt_smooth->Scale(1./2./thetaScatt_ref->GetNbinsX());
  thetaScatt_smooth->FFT(fftScatt_ref,"Re R2R_7");
  
  // Divide the fft of the data by the fft of the reference
  fftX_div->Divide(fftX_data, fftX_ref);
  fftY_div->Divide(fftY_data, fftY_ref);

  // Undo the fft
  fftX_div->FFT(thetaX_fft,"Re R2R_3");
  fftY_div->FFT(thetaY_fft,"Re R2R_3");

  // Correct the normalization after the FFT
  thetaX_fft->Scale(1./2./fftX_data->GetNbinsX());
  thetaY_fft->Scale(1./2./fftY_data->GetNbinsX());
  thetaScatt_fft->Scale(1./2./thetaScatt_data->GetNbinsX());
  */
  outfile->cd();

  thetaX_data->Write();
  thetaY_data->Write();
  thetaScatt_data->Write();
  theta2Scatt_data->Write();
  thetaX_ref->Write();
  thetaY_ref->Write();
  thetaScatt_ref->Write();
  theta2Scatt_ref->Write();
  projposUSDSdiff->Write();
  scattering_proj_x_FFT->Write();
  scattering_proj_y_FFT->Write();
  /*
  thetaX_fft->Write();
  thetaY_fft->Write();
  thetaScatt_fft->Write();
 
  fftX_data->Write();
  fftY_data->Write();
  fftScatt_dataRe->Write();
  fftScatt_dataIm->Write();
  fftX_ref->Write();
  fftY_ref->Write();
  fftScatt_refRe->Write();
  fftScatt_refIm->Write();
  fftX_div->Write();
  fftY_div->Write();
  fftScatt_divRe->Write();
  fftScatt_divIm->Write();
  */
  delete projposUSDSdiff;
  delete scattering_proj_x_FFT;
  delete scattering_proj_y_FFT;
}

void MCSAnalysis::UpdateRunInfo(){
  _run.push_back(runnumber);
  _emit.push_back(_UStmpset.emittance(4, 1));
  _emitError.push_back(_emit.back() / sqrt(_UStmpset.N()));
  _meanX.push_back(_UStmpset.mean().X);  
  _meanXerror.push_back(sqrt(_UStmpset.cov()[0][0]));  
  _meanY.push_back(_UStmpset.mean().Y);  
  _meanYerror.push_back(sqrt(_UStmpset.cov()[2][2]));  
  _meandXdz.push_back(_UStmpset.mean().dXdz);  
  _meandXdzError.push_back(sqrt(_UStmpset.cov()[1][1]));  
  _meandYdz.push_back(_UStmpset.mean().dYdz);  
  _meandYdzError.push_back(sqrt(_UStmpset.cov()[3][3]));  
  
  _UStmpset.Set().clear();
  LastRunNumber = runnumber;
}

void MCSAnalysis::PlotRunInfo(){
  TGraphErrors* gemitt = new TGraphErrors(_run.size());
  gemitt->SetName("gemitt");
  gemitt->SetTitle(";Run Number;Emittance (mm)");
  TGraphErrors* gmeanX = new TGraphErrors(_run.size());
  gmeanX->SetName("gmeanX");
  gmeanX->SetTitle(";Run Number;Mean X (mm)");
  TGraphErrors* gmeanY = new TGraphErrors(_run.size());
  gmeanY->SetName("gmeanY");
  gmeanY->SetTitle(";Run Number;Mean Y (mm)");
  TGraphErrors* gmeandXdz = new TGraphErrors(_run.size());
  gmeandXdz->SetName("gmeandXdz");
  gmeandXdz->SetTitle(";Run Number;Mean dXdz (mm)");
  TGraphErrors* gmeandYdz = new TGraphErrors(_run.size());
  gmeandYdz->SetName("gmeandYdz");
  gmeandYdz->SetTitle(";Run Number;Mean dYdz (mm)");

  for (size_t q=0; q<_run.size(); q++){
    gemitt->SetPoint(q, _run.at(q), _emit.at(q));
    gemitt->SetPointError(q, 0, _emitError.at(q));
    gmeanX->SetPoint(q, _run.at(q), _meanX.at(q));
    gmeanX->SetPointError(q, 0, _meanXerror.at(q));
    gmeanY->SetPoint(q, _run.at(q), _meanY.at(q));
    gmeanY->SetPointError(q, 0, _meanYerror.at(q));
    gmeandXdz->SetPoint(q, _run.at(q), _meandXdz.at(q));
    gmeandXdz->SetPointError(q, 0, _meandXdzError.at(q));
    gmeandYdz->SetPoint(q, _run.at(q), _meandYdz.at(q));
    gmeandYdz->SetPointError(q, 0, _meandYdzError.at(q));
  }
  outfile->cd();
  gemitt->Write();
  gmeanX->Write();
  gmeanY->Write();
  gmeandXdz->Write();
  gmeandYdz->Write();

  _run.clear();
  _emit.clear();
  _meanX.clear();
  _meanY.clear();
  _meandXdz.clear();
  _meandYdz.clear();
  
}

bool MCSAnalysis::findVirtualPlanes(){
  bool planesfound = false;
  // The virtual planes should be the same event to event so once they
  // are found once they should be consistent
  if (USrefplaneI > 0 && USabsPlaneI >0 &&
      DSrefplaneI > 0 && DSabsPlaneI >0 ){
    if (static_cast<unsigned int>(USrefplaneI) < mcevent->GetVirtualHits()->size() &&
	static_cast<unsigned int>(DSrefplaneI) < mcevent->GetVirtualHits()->size() &&
	static_cast<unsigned int>(DSabsPlaneI) < mcevent->GetVirtualHits()->size() &&
	static_cast<unsigned int>(USabsPlaneI) < mcevent->GetVirtualHits()->size()){
      planesfound = true;
      return planesfound;
    }
  }
  if (1) { 
    double USmindiff = 9999, DSmindiff=9999;
    int USminI=-1, DSminI=-1;
    double USrefpos = _sys["abspos"] - 149.;
    double DSrefpos = _sys["abspos"] + 149.;
    for ( size_t j=0; j < mcevent->GetVirtualHits()->size(); j++){
      double ztest = mcevent->GetVirtualHits()->at(j).GetPosition().z();
      double USzdiff = fabs(ztest - USrefpos);
      double DSzdiff = fabs(ztest - DSrefpos);
      // std::cout<<j<<"\t"<<ztest<<"\t"<<USzdiff<<"\t"<<DSzdiff<<std::endl;
      if(USzdiff < USmindiff){
	
	USmindiff = USzdiff;
	USminI = j;
      }
      if(DSzdiff < DSmindiff){
	//std::cout<<j<<"\t"<<ztest<<"\t"<<DSzdiff<<std::endl;
	DSmindiff = DSzdiff;
	DSminI = j;
      }
    }
    if (USmindiff < 150.){
      planesfound = true;
      USabsPlaneI = USminI;
      // std::cout<<"US absorber plane.\t"<<USrefpos<<"\t"<<USmindiff<<"\t"<<USabsPlaneI<<"\t"<<mcevent->GetVirtualHits()->at(USabsPlaneI).GetPosition().z()<<std::endl;
    }
    if (DSmindiff < 150. && planesfound==true){
      DSabsPlaneI = DSminI;
      // std::cout<<"DS absorber plane.\t"<<DSrefpos<<"\t"<<DSmindiff<<"\t"<<DSabsPlaneI<<"\t"<<mcevent->GetVirtualHits()->at(DSabsPlaneI).GetPosition().z()<<std::endl;
    }
    else 
      planesfound = false;
    
    // Now I will fill indicies corresponding to the tracker
    // reference planes
    USrefpos = USrefplaneZ;
    DSrefpos = DSrefplaneZ;
    
    USmindiff = 9999, DSmindiff=9999;
    USminI=-1, DSminI=-1;
    
    for ( size_t j=0; j < mcevent->GetVirtualHits()->size(); j++){
      double ztest = mcevent->GetVirtualHits()->at(j).GetPosition().z();
      double USzdiff = fabs(ztest - USrefpos);
      if(USzdiff < USmindiff){
	USmindiff = USzdiff;
	USminI = j;
    }
      double DSzdiff = fabs(ztest - DSrefpos);
      if(DSzdiff < DSmindiff){
	DSmindiff = DSzdiff;
	DSminI = j;
      }
    }
    if (planesfound == true){
      if (USmindiff < 50. ){
	USrefplaneI = USminI;
	// std::cout<<"US reference plane.\t"<<USrefplaneZ<<"\t"<<USmindiff<<"\t"<<USrefplaneI<<"\t"<<mcevent->GetVirtualHits()->at(USrefplaneI).GetPosition().z()<<std::endl;
      }
      else  
	planesfound = false;
      if (DSmindiff < 50. ){
	DSrefplaneI = DSminI;
	// std::cout<<"DS reference plane.\t"<<DSrefplaneZ<<"\t"<<DSmindiff<<"\t"<<DSrefplaneI<<"\t"<<mcevent->GetVirtualHits()->at(DSrefplaneI).GetPosition().z()<<std::endl;
      }
      else 
	planesfound = false;
    }
    return planesfound;
  }
}
void MCSAnalysis::FillMuScattResponse(bool event_ok, Vars& US, Vars& DS, Vars& USMC, Vars& DSMC){
  double thetaYMC = atan(DSMC.dXdz) - atan(USMC.dXdz);
  double thetaXMC = atan(DSMC.dYdz) - atan(USMC.dYdz);
  double XMC = DSMC.X;
  double YMC = DSMC.Y;
  theta_true_x_graph->Fill(thetaXMC);
  theta_true_y_graph->Fill(thetaYMC);
  if ( event_ok ){
    scattering_proj_y_R->Fill(XMC, thetaYMC);
    scattering_proj_x_R->Fill(YMC, thetaXMC);
    theta_meas_x_graph->Fill(thetaXMC);
    theta_meas_y_graph->Fill(thetaYMC);
  }
}

void MCSAnalysis::FillMCSResponse(bool event_ok, Vars& US, Vars& DS, Vars& USMC, Vars& DSMC){
  
  double thetaYMC = atan(DSMC.dXdz) - atan(USMC.dXdz);
  double thetaXMC = atan(DSMC.dYdz) - atan(USMC.dYdz);
  double costhetaScattMC = ( (1 + USMC.dXdz * DSMC.dXdz + USMC.dYdz * DSMC.dYdz )/
			      sqrt(1 + USMC.dXdz*USMC.dXdz + USMC.dYdz*USMC.dYdz)/
			      sqrt(1 + DSMC.dXdz*DSMC.dXdz + DSMC.dYdz*DSMC.dYdz));
  double thetaScattMC = acos(costhetaScattMC);
  if(event_ok){
    double thetaY = atan(DS.dXdz) - atan(US.dXdz);
    double thetaX = atan(DS.dYdz) - atan(US.dYdz);
    double cosScatt = ( (1 + US.dXdz * DS.dXdz + US.dYdz * DS.dYdz )/
			sqrt(1 + US.dXdz*US.dXdz + US.dYdz*US.dYdz)/
			sqrt(1 + DS.dXdz*DS.dXdz + DS.dYdz*DS.dYdz));
    double thetaScatt = acos(cosScatt);
    resp_thetaX.Fill(thetaX, thetaXMC);
    resp_thetaY.Fill(thetaY, thetaYMC);
    resp_thetaScatt.Fill(thetaScatt, thetaScattMC);
    resp_theta2Scatt.Fill(thetaScatt, thetaScattMC * thetaScattMC);

    scattering_proj_x_R->Fill(DS.X - US.X, thetaYMC);
    scattering_proj_y_R->Fill(DS.Y - US.Y, thetaXMC);
    
    theta_true_x_graph->Fill(thetaXMC);
    theta_true_y_graph->Fill(thetaYMC);
  }
  else{
    resp_thetaX.Miss(thetaXMC);
    resp_thetaY.Miss(thetaYMC);
    resp_thetaScatt.Miss(thetaScattMC);
    resp_theta2Scatt.Miss(thetaScattMC * thetaScattMC);

    theta_true_x_graph->Fill(thetaXMC);
    theta_true_y_graph->Fill(thetaYMC);
  }
  
}
 
void MCSAnalysis::FillVarsVirtual(Vars &tmpvar, int j){
  
  tmpvar.X  = mcevent->GetVirtualHits()->at(j).GetPosition().x();
  tmpvar.Y  = mcevent->GetVirtualHits()->at(j).GetPosition().y();
  tmpvar.Z  = mcevent->GetVirtualHits()->at(j).GetPosition().z();
  tmpvar.dXdz = mcevent->GetVirtualHits()->at(j).GetMomentum().x()/
    mcevent->GetVirtualHits()->at(j).GetMomentum().z();
  tmpvar.dYdz = mcevent->GetVirtualHits()->at(j).GetMomentum().y()/
    mcevent->GetVirtualHits()->at(j).GetMomentum().z();
 // double px = mcevent->GetVirtualHits()->at(j).GetMomentum().x();
 // double py = mcevent->GetVirtualHits()->at(j).GetMomentum().y();
  double pz = mcevent->GetVirtualHits()->at(j).GetMomentum().z();
  tmpvar.TOF12 = 27.5 * sqrt(1 + 105.65*105.65/pz/pz);
  tmpvar.TOF01 = 26.5 * sqrt(1 + 105.65*105.65/pz/pz);
  tmpvar.pz = pz;
  tmpvar.pid = mcevent->GetVirtualHits()->at(j).GetParticleId();
  tmpvar.isgood = true;
    // mcevent->GetVirtualHits()->at(j).GetMomentum().z();
  /*
  while(_SPsets.size() <= j){
    Ensemble tmpens;`
    _SPsets.push_back(tmpens);
  }
  _SPsets[j].append_instance(tmpvar); 
  */
}


void MCSAnalysis::FillVarsSciFi(Vars& tmpvar, int j, int k, double pz, int isDS){
  double alX=0.0, alY=0.0, thX=0.0, thY=0.0;
  if(isDS==0){
    alX = _sys["alXUS"];
    alY = _sys["alYUS"];
    thX = _sys["thYUS"];
    thY = _sys["thXUS"];
  } else {
    alX = _sys["alXDS"];
    alY = _sys["alYDS"];
    thX = _sys["thYDS"];
    thY = _sys["thXDS"];
  }

  if(j < int(scifievent->scifitracks().size()) && j != -1){
    if(k < int(scifievent->scifitracks()[j]->scifitrackpoints().size()) && k != -1){
      tmpvar.X  = scifievent->scifitracks()[j]->scifitrackpoints()[k]->pos().x() + alX;
      tmpvar.Y  = scifievent->scifitracks()[j]->scifitrackpoints()[k]->pos().y() + alY;
      tmpvar.Z  = scifievent->scifitracks()[j]->scifitrackpoints()[k]->pos().z();
      tmpvar.dXdz = scifievent->scifitracks()[j]->scifitrackpoints()[k]->gradient().x() + tan(thX*atan(1.)/45.0);
      // scifievent->scifitracks()[j]->scifitrackpoints()[k]->mom().z();
      tmpvar.dYdz = scifievent->scifitracks()[j]->scifitrackpoints()[k]->gradient().y() + tan(thY*atan(1.)/45.0);
      tmpvar.pz = pz * sqrt(1 + pow(tmpvar.dXdz, 2) + pow(tmpvar.dYdz, 2));
      tmpvar.px = scifievent->scifitracks()[j]->scifitrackpoints()[k]->mom().x() + tan(thX*atan(1.)/45.0)*tmpvar.pz;
      tmpvar.py = scifievent->scifitracks()[j]->scifitrackpoints()[k]->mom().y() + tan(thY*atan(1.)/45.0)*tmpvar.pz;
      tmpvar.pid = -13;
      tmpvar.isgood = true;
      // scifievent->scifitracks()[j]->scifitrackpoints()[k]->mom().z();
    }
  }
  if( j == -1 || k == -1){
    // Vars tmpvar;
    tmpvar.X = 0.0;
    tmpvar.Y = 0.0;
    tmpvar.Z = 0.0;
    tmpvar.dXdz = 1./2.;
    tmpvar.dYdz = 1./2.;
    tmpvar.px   = pz/2.;
    tmpvar.py   = pz/2.;
    tmpvar.pz   = pz;
    tmpvar.pid = -13;
    tmpvar.isgood = false;
  }
  if( tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size() > 0 &&
      tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray().size() > 0)
    tmpvar.TOF12 = tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray()[0].GetTime() 
      - tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime();
  else
    tmpvar.TOF12 = 100.0 * 8.22475 / 0.299792458;
  if( tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray().size() > 0 &&
      tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size() > 0)
    tmpvar.TOF01 = tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime() 
      - tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[0].GetTime();
  else
    tmpvar.TOF01 = 100.0 * 7.64186 / 0.299792458;
}

void MCSAnalysis::make_beam_histograms(Collection Set, std::string desc, std::string suffix){

//  std::cerr<<"Start of Make beam Histograms";    
  std::string tmptitle = desc + ";X (mm); Y (mm)";
  std::string tmpname  = suffix + "_XY";
  TH2D* XY = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, -225, 225);
  tmptitle = desc + ";X (mm); dXdz";
  tmpname  = suffix + "_XdXdz";
  TH2D* XdXdz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, -0.125, 0.125);
  tmptitle = desc + ";X (mm); dYdz";
  tmpname  = suffix + "_XdYdz";
  TH2D* XdYdz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, -0.125, 0.125);
  tmptitle = desc + ";Y (mm); dXdz";
  tmpname  = suffix + "_YdXdz";
  TH2D* YdXdz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, -0.125, 0.125);
  tmptitle = desc + ";Y (mm); dYdz";
  tmpname  = suffix + "_YdYdz";
  TH2D* YdYdz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, -0.125, 0.125);
  tmptitle = desc + ";dXdz; dYdz";
  tmpname  = suffix + "_dXdzdYdz";
  TH2D* dXdzdYdz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -1.25, 1.25, 90, -1.25, 1.25);
  tmptitle = desc + ";X (mm); TOF12 (ns)";
  tmpname  = suffix + "_XTOF12";
  TH2D* XTOF12 = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, 20.0, 45.0);
  tmptitle = desc + ";Y (mm); TOF12 (ns)";
  tmpname  = suffix + "_YTOF12";
  TH2D* YTOF12 = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, 20.0, 45.0);
  tmptitle = desc + ";dXdz; TOF12 (ns)";
  tmpname  = suffix + "_dXdzTOF12";
  TH2D* dXdzTOF12 = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -0.125, 0.125, 90, 20.0, 45.0);
  tmptitle = desc + ";dYdz; TOF12 (ns)";
  tmpname  = suffix + "_dYdzTOF12";
  TH2D* dYdzTOF12 = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -0.125, 0.125, 90, 20.0, 45.0);
  
//  std::cerr<<" - Part 1 Complete";

  for(int i=0; i<Set.N(); i++){
    if(Set.E(i).X==0 && Set.E(i).Y==0) continue;
    XY->Fill(Set.E(i).X,Set.E(i).Y);
    XdXdz->Fill(Set.E(i).X,Set.E(i).dXdz);
    XdYdz->Fill(Set.E(i).X,Set.E(i).dYdz);
    YdXdz->Fill(Set.E(i).Y,Set.E(i).dXdz);
    YdYdz->Fill(Set.E(i).Y,Set.E(i).dYdz);
    dXdzdYdz->Fill(Set.E(i).dXdz,Set.E(i).dYdz);
    XTOF12->Fill(Set.E(i).X,Set.E(i).TOF12);
    YTOF12->Fill(Set.E(i).Y,Set.E(i).TOF12);
    dXdzTOF12->Fill(Set.E(i).dXdz,Set.E(i).TOF12);
    dYdzTOF12->Fill(Set.E(i).dYdz,Set.E(i).TOF12);
  }
//  std::cerr<<" - Part 2 Complete";
  
  XY->GetXaxis()->SetLabelSize(0.05);
  XY->GetXaxis()->SetTitleSize(0.05);
  XY->GetYaxis()->SetLabelSize(0.05);
  XY->GetYaxis()->SetTitleSize(0.05);
  XdXdz->GetXaxis()->SetLabelSize(0.05);
  XdXdz->GetXaxis()->SetTitleSize(0.05);
  XdXdz->GetYaxis()->SetLabelSize(0.05);
  XdXdz->GetYaxis()->SetTitleSize(0.05);
  XdYdz->GetXaxis()->SetLabelSize(0.05);
  XdYdz->GetXaxis()->SetTitleSize(0.05);
  XdYdz->GetYaxis()->SetLabelSize(0.05);
  XdYdz->GetYaxis()->SetTitleSize(0.05);
  YdXdz->GetXaxis()->SetLabelSize(0.05);
  YdXdz->GetXaxis()->SetTitleSize(0.05);
  YdXdz->GetYaxis()->SetLabelSize(0.05);
  YdXdz->GetYaxis()->SetTitleSize(0.05);
  YdYdz->GetXaxis()->SetLabelSize(0.05);
  YdYdz->GetXaxis()->SetTitleSize(0.05);
  YdYdz->GetYaxis()->SetLabelSize(0.05);
  YdYdz->GetYaxis()->SetTitleSize(0.05);
  dXdzdYdz->GetXaxis()->SetLabelSize(0.05);
  dXdzdYdz->GetXaxis()->SetTitleSize(0.05);
  dXdzdYdz->GetYaxis()->SetLabelSize(0.05);
  dXdzdYdz->GetYaxis()->SetTitleSize(0.05);
  XTOF12->GetXaxis()->SetLabelSize(0.05);
  XTOF12->GetXaxis()->SetTitleSize(0.05);
  XTOF12->GetYaxis()->SetLabelSize(0.05);
  XTOF12->GetYaxis()->SetTitleSize(0.05);
  YTOF12->GetXaxis()->SetLabelSize(0.05);
  YTOF12->GetXaxis()->SetTitleSize(0.05);
  YTOF12->GetYaxis()->SetLabelSize(0.05);
  YTOF12->GetYaxis()->SetTitleSize(0.05);
  dXdzTOF12->GetXaxis()->SetLabelSize(0.05);
  dXdzTOF12->GetXaxis()->SetTitleSize(0.05);
  dXdzTOF12->GetYaxis()->SetLabelSize(0.05);
  dXdzTOF12->GetYaxis()->SetTitleSize(0.05);
  dYdzTOF12->GetXaxis()->SetLabelSize(0.05);
  dYdzTOF12->GetXaxis()->SetTitleSize(0.05);
  dYdzTOF12->GetYaxis()->SetLabelSize(0.05);
  dYdzTOF12->GetYaxis()->SetTitleSize(0.05);

//  std::cerr<<" - Part 3 Complete ";

  outfile->cd();
//  std::cerr<<" - i";
  XY->Write();
//  std::cerr<<" - ii";
  XdXdz->Write();
//  std::cerr<<" - iii";
  XdYdz->Write();
//  std::cerr<<" - iv";
  XTOF12->Write();
//  std::cerr<<" - v";
  YdXdz->Write();
//  std::cerr<<" - vi";
  YdYdz->Write();
//  std::cerr<<" - vii";
  YTOF12->Write();
//  std::cerr<<" - viii";
  dXdzdYdz->Write();
//  std::cerr<<" - ix";
  dXdzTOF12->Write();
//  std::cerr<<" - x";
  dYdzTOF12->Write();

//  std::cerr<<" - Part 3 a Complete"; 
 
  TCanvas* c1 = new TCanvas();
  XY->Draw("colz");
  std::string tmpfile = suffix + "_XY.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();

//  std::cerr<<" - i";

  XdXdz->Draw("colz");
  tmpfile = suffix + "_XdXdz.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();

//  std::cerr<<" - ii";

  XdYdz->Draw("colz");
  tmpfile = suffix + "_XdYdz.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();

//  std::cerr<<" - iii";

  YdXdz->Draw("colz");
  tmpfile = suffix + "_YdXdz.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();

//  std::cerr<<" - iv";

  YdYdz->Draw("colz");
  tmpfile = suffix + "_YdYdz.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();

//  std::cerr<<" - v";

  dXdzdYdz->Draw("colz");
  tmpfile = suffix + "_dXdzdYdz.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();

//  std::cerr<<" - vi";

  XTOF12->Draw("colz");
  tmpfile = suffix + "_XTOF12.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();

//  std::cerr<<" - vii";

  YTOF12->Draw("colz");
  tmpfile = suffix + "_YTOF12.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();

//  std::cerr<<" - viii";

  dXdzTOF12->Draw("colz");
  tmpfile = suffix + "_dXdzTOF12.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();

//  std::cerr<<" - ix";

  dYdzTOF12->Draw("colz");
  tmpfile = suffix + "_dYdzTOF12.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();

//  std::cerr<<" - Part 4 Complete"<<std::endl;  

  delete c1;
  delete XY;
  delete XdXdz;
  delete XdYdz;
  delete YdXdz;
  delete YdYdz;
  delete XTOF12;
  delete YTOF12;
  delete dXdzdYdz;
  delete dXdzTOF12;
  delete dYdzTOF12;
}

void MCSAnalysis::make_acceptance_histograms(Collection USset, Collection DSset, 
					     std::string desc, std::string suffix){
  std::string tmptitle = desc + ";X (mm); Y (mm)";
  std::string tmpname  = suffix + "_posaccXY";
  TH2D* posaccXY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 75, -300, 300, 75, -300, 300);
  tmptitle = desc + ";X (mm); Y (mm)";
  tmpname  = suffix + "_divaccXY";
  TH2D* divaccXY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 75, -300, 300, 75, -300, 300);
  tmptitle = desc + ";X (mm); Y (mm)";
  tmpname  = suffix + "_posresXY";
  TH2D* posresXY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 75, -300, 300, 75, -300, 300);
  tmptitle = desc + ";X (mm); Y (mm)";
  tmpname  = suffix + "_divresXY";
  TH2D* divresXY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 75, -300, 300, 75, -300, 300);
  tmptitle = desc + ";X (mm); Y (mm)";
  tmpname  = suffix + "_posres2XY";
  TH2D* posres2XY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 75, -300, 300, 75, -300, 300);
  tmptitle = desc + ";X (mm); Y (mm)";
  tmpname  = suffix + "_divres2XY";
  TH2D* divres2XY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 75, -300, 300, 75, -300, 300);
  tmptitle = desc + ";X (mm); Y (mm)";
  tmpname  = suffix + "_posXY";
  TH2D* posXY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 75, -300, 300, 75, -300, 300);
  tmptitle = desc + ";dX/dz ; dY/dz";
  tmpname  = suffix + "_posaccDXDY";
  TH2D* posaccDXDY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 90, -0.125, 0.125, 90, -0.125, 0.125);
  tmptitle = desc + ";dX/dz ; dY/dz";
  tmpname  = suffix + "_divaccDXDY";
  TH2D* divaccDXDY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 90, -0.125, 0.125, 90, -0.125, 0.125);
  tmptitle = desc + ";dX/dz ; dY/dz";
  tmpname  = suffix + "_posresDXDY";
  TH2D* posresDXDY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 90, -0.125, 0.125, 90, -0.125, 0.125);
  tmptitle = desc + ";dX/dz ; dY/dz";
  tmpname  = suffix + "_divresDXDY";
  TH2D* divresDXDY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 90, -0.125, 0.125, 90, -0.125, 0.125);
  tmptitle = desc + ";dX/dz ; dY/dz";
  tmpname  = suffix + "_posres2DXDY";
  TH2D* posres2DXDY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 90, -0.125, 0.125, 90, -0.125, 0.125);
  tmptitle = desc + ";dX/dz ; dY/dz";
  tmpname  = suffix + "_divres2DXDY";
  TH2D* divres2DXDY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 90, -0.125, 0.125, 90, -0.125, 0.125);
  tmptitle = desc + ";dX/dz ; dY/dz";
  tmpname  = suffix + "_posDXDY";
  TH2D* posDXDY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 90, -0.125, 0.125, 90, -0.125, 0.125);

  tmptitle = desc + ";|#vec{r}^{proj#to DS}_{US} - #vec{r}^{meas}_{DS}";
  tmpname  = suffix + "_appRef";
  TH1D* appRef = new TH1D(tmpname.c_str(), tmptitle.c_str(), 200, 0, 200);
  tmptitle = desc + ";X_{US}^{proj#to DS} - X_{DS}^{meas};Y_{US}^{proj#to DS} - Y_{DS}^{meas}";
  tmpname  = suffix + "_projRefDiff";
  TH2D* projRefDiff = new TH2D(tmpname.c_str(), tmptitle.c_str(), 200, -100, 100, 200, -100, 100);

  double latchZ = 0.0;
  Vars DSproj;
  for (int i=0; i<USset.N(); i++){
    if (DSset.E(i).Z != 0 && latchZ == 0.0){
      latchZ = DSset.E(i).Z;
    } 
    if (latchZ != 0.0 && fabs(DSset.E(i).Z - latchZ) < 5){
//      std::cerr<<"make_acceptance_histograms - Propagating to DSset.E(i).Z"<<std::endl;
      DSproj = PropagateVarsMu(USset.E(i), DSset.E(i).Z);
    } else if(latchZ != 0.0){
//      std::cerr<<"make_acceptance_histograms - Propagating to latchZ"<<std::endl;
      DSproj = PropagateVarsMu(USset.E(i), latchZ);
    } else {
      continue;
    }
    
    Vars predDiff = DSproj - DSset.E(i);
    double approach = sqrt(predDiff.X * predDiff.X + 
			   predDiff.Y * predDiff.Y +
			   predDiff.Z * predDiff.Z);
    appRef->Fill(approach); // can be very large if the ds track is lost
    projRefDiff->Fill(predDiff.X, predDiff.Y);
    posXY->Fill(DSproj.X, DSproj.Y);
    posDXDY->Fill(DSproj.dXdz, DSproj.dYdz);
    
    if ( DSset.E(i).Z > USset.E(i).Z ){ // lost ds tracks are set to z = 0
      posaccXY->Fill(DSproj.X, DSproj.Y);
      posaccDXDY->Fill(DSproj.dXdz, DSproj.dYdz);
      if ( approach < 36.) { // 36 mm is approximately the resolution of the projection
	posresXY->Fill(DSproj.X, DSproj.Y);
	posresDXDY->Fill(DSproj.dXdz, DSproj.dYdz);
      }
      if ( approach < 2.*36.) { // 36 mm is approximately the resolution of the projection
	posres2XY->Fill(DSproj.X, DSproj.Y);
	posres2DXDY->Fill(DSproj.dXdz, DSproj.dYdz);
      }
    }
  }
  posXY->Sumw2();
  posaccXY->Sumw2();
  divaccXY->Divide(posaccXY, posXY);
  posresXY->Sumw2();
  divresXY->Divide(posresXY, posXY);
  posres2XY->Sumw2();
  divres2XY->Divide(posres2XY, posXY);
  posaccDXDY->Sumw2();
  divaccDXDY->Divide(posaccDXDY, posDXDY);
  posresDXDY->Sumw2();
  divresDXDY->Divide(posresDXDY, posDXDY);
  posres2DXDY->Sumw2();
  divres2DXDY->Divide(posres2DXDY, posDXDY);

  posXY->GetXaxis()->SetLabelSize(0.05);
  posXY->GetXaxis()->SetTitleSize(0.05);
  posXY->GetYaxis()->SetLabelSize(0.05);
  posXY->GetYaxis()->SetTitleSize(0.05);
  posDXDY->GetXaxis()->SetLabelSize(0.05);
  posDXDY->GetXaxis()->SetTitleSize(0.05);
  posDXDY->GetYaxis()->SetLabelSize(0.05);
  posDXDY->GetYaxis()->SetTitleSize(0.05);
  posaccXY->GetXaxis()->SetLabelSize(0.05);
  posaccXY->GetXaxis()->SetTitleSize(0.05);
  posaccXY->GetYaxis()->SetLabelSize(0.05);
  posaccXY->GetYaxis()->SetTitleSize(0.05);
  divaccXY->GetXaxis()->SetLabelSize(0.05);
  divaccXY->GetXaxis()->SetTitleSize(0.05);
  divaccXY->GetYaxis()->SetLabelSize(0.05);
  divaccXY->GetYaxis()->SetTitleSize(0.05);
  posresXY->GetXaxis()->SetLabelSize(0.05);
  posresXY->GetXaxis()->SetTitleSize(0.05);
  posresXY->GetYaxis()->SetLabelSize(0.05);
  posresXY->GetYaxis()->SetTitleSize(0.05);
  divresXY->GetXaxis()->SetLabelSize(0.05);
  divresXY->GetXaxis()->SetTitleSize(0.05);
  divresXY->GetYaxis()->SetLabelSize(0.05);
  divresXY->GetYaxis()->SetTitleSize(0.05);
  posres2XY->GetXaxis()->SetLabelSize(0.05);
  posres2XY->GetXaxis()->SetTitleSize(0.05);
  posres2XY->GetYaxis()->SetLabelSize(0.05);
  posres2XY->GetYaxis()->SetTitleSize(0.05);
  divres2XY->GetXaxis()->SetLabelSize(0.05);
  divres2XY->GetXaxis()->SetTitleSize(0.05);
  divres2XY->GetYaxis()->SetLabelSize(0.05);
  divres2XY->GetYaxis()->SetTitleSize(0.05);
  posaccDXDY->GetXaxis()->SetLabelSize(0.05);
  posaccDXDY->GetXaxis()->SetTitleSize(0.05);
  posaccDXDY->GetYaxis()->SetLabelSize(0.05);
  posaccDXDY->GetYaxis()->SetTitleSize(0.05);
  divaccDXDY->GetXaxis()->SetLabelSize(0.05);
  divaccDXDY->GetXaxis()->SetTitleSize(0.05);
  divaccDXDY->GetYaxis()->SetLabelSize(0.05);
  divaccDXDY->GetYaxis()->SetTitleSize(0.05);
  posresDXDY->GetXaxis()->SetLabelSize(0.05);
  posresDXDY->GetXaxis()->SetTitleSize(0.05);
  posresDXDY->GetYaxis()->SetLabelSize(0.05);
  posresDXDY->GetYaxis()->SetTitleSize(0.05);
  divresDXDY->GetXaxis()->SetLabelSize(0.05);
  divresDXDY->GetXaxis()->SetTitleSize(0.05);
  divresDXDY->GetYaxis()->SetLabelSize(0.05);
  divresDXDY->GetYaxis()->SetTitleSize(0.05);
  posres2DXDY->GetXaxis()->SetLabelSize(0.05);
  posres2DXDY->GetXaxis()->SetTitleSize(0.05);
  posres2DXDY->GetYaxis()->SetLabelSize(0.05);
  posres2DXDY->GetYaxis()->SetTitleSize(0.05);
  divres2DXDY->GetXaxis()->SetLabelSize(0.05);
  divres2DXDY->GetXaxis()->SetTitleSize(0.05);
  divres2DXDY->GetYaxis()->SetLabelSize(0.05);
  divres2DXDY->GetYaxis()->SetTitleSize(0.05);

  outfile->cd();
  posXY->Write();
  posaccXY->Write();
  posresXY->Write();
  divres2XY->Write();
  divaccXY->Write();
  divresXY->Write();
  divres2XY->Write();
  posDXDY->Write();
  posaccDXDY->Write();
  posresDXDY->Write();
  divres2DXDY->Write();
  divaccDXDY->Write();
  divresDXDY->Write();
  divres2DXDY->Write();
  appRef->Write();
  projRefDiff->Write();

  delete posXY;
  delete posDXDY;
  delete posaccXY;
  delete posaccDXDY;
  delete posresXY;
  delete posresDXDY;
  delete posres2XY;
  delete posres2DXDY;
  delete divaccXY;
  delete divaccDXDY;
  delete divresXY;
  delete divresDXDY;
  delete divres2XY;
  delete divres2DXDY;
  delete appRef;
  delete projRefDiff;
}

void MCSAnalysis::make_scattering_acceptance_histograms(Collection USset,
							Collection DSset,
							Collection DSrefset,
							std::string desc,
							std::string suffix){
  
  std::string tmptitle = desc + ";#theta_Y (radians)";
  std::string tmpname = "thetaY_all_";
  tmpname += suffix;
  TH1D* thetaY_all = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmptitle = desc + ";#theta_X (radians)";
  tmpname = "thetaX_all_";
  tmpname += suffix;
  TH1D* thetaX_all = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

  tmptitle = desc + ";#theta_Y (radians)";
  tmpname = "thetaY_acc_";
  tmpname += suffix;
  TH1D* thetaY_acc = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmptitle = desc + ";#theta_X (radians)";
  tmpname = "thetaX_acc_";
  tmpname += suffix;
  TH1D* thetaX_acc = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

  tmptitle = desc + ";#theta_Y (radians)";
  tmpname = "thetaY_divacc_";
  tmpname += suffix;
  TH1D* thetaY_divacc = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmptitle = desc + ";#theta_X (radians)";
  tmpname = "thetaX_divacc_";
  tmpname += suffix;
  TH1D* thetaX_divacc = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

  tmptitle = desc + ";#theta_Y (radians)";
  tmpname = "thetaY_res_";
  tmpname += suffix;
  TH1D* thetaY_res = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmptitle = desc + ";#theta_X (radians)";
  tmpname = "thetaX_res_";
  tmpname += suffix;
  TH1D* thetaX_res = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

  tmptitle = desc + ";#theta_Y (radians)";
  tmpname = "thetaY_divres_";
  tmpname += suffix;
  TH1D* thetaY_divres = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmptitle = desc + ";#theta_X (radians)";
  tmpname = "thetaX_divres_";
  tmpname += suffix;
  TH1D* thetaX_divres = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

  tmptitle = desc + ";#theta_Y (radians)";
  tmpname = "thetaY_res2_";
  tmpname += suffix;
  TH1D* thetaY_res2 = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmptitle = desc + ";#theta_X (radians)";
  tmpname = "thetaX_res2_";
  tmpname += suffix;
  TH1D* thetaX_res2 = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

  tmptitle = desc + ";#theta_Y (radians)";
  tmpname = "thetaY_divres2_";
  tmpname += suffix;
  TH1D* thetaY_divres2 = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmptitle = desc + ";#theta_X (radians)";
  tmpname = "thetaX_divres2_";
  tmpname += suffix;
  TH1D* thetaX_divres2 = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  
  for (int i=0; i<USset.N(); i++){
    if (i >= DSrefset.N()) 
      break;
    std::vector<double> projTheta = DefineProjectionAngles(USset.E(i), DSset.E(i));
    double thetaX = projTheta[0];
    double thetaY = projTheta[1];
    // double thetaScatt = projTheta[2];
    Vars DSdiff = DSset.E(i) - DSrefset.E(i);
    double approach = sqrt(DSdiff.X * DSdiff.X + 
			   DSdiff.Y * DSdiff.Y +
			   DSdiff.Z * DSdiff.Z);
    thetaX_all->Fill(thetaX);
    thetaY_all->Fill(thetaY);
    if (DSrefset.E(i).Z != 0){
      thetaX_acc->Fill(thetaX);
      thetaY_acc->Fill(thetaY);
    
      if (approach < 36.){
	thetaX_res->Fill(thetaX);
	thetaY_res->Fill(thetaY);
      }
      if (approach < 2* 36.){
	thetaX_res2->Fill(thetaX);
	thetaY_res2->Fill(thetaY);
      }
    }
  }
  thetaX_all->Sumw2();
  thetaY_all->Sumw2();
  thetaX_acc->Sumw2();
  thetaY_acc->Sumw2();
  thetaX_divacc->Divide(thetaX_acc, thetaX_all);
  thetaY_divacc->Divide(thetaY_acc, thetaY_all);
  thetaX_res->Sumw2();
  thetaY_res->Sumw2();
  thetaX_divres->Divide(thetaX_res, thetaX_all);
  thetaY_divres->Divide(thetaY_res, thetaY_all);
  thetaX_res2->Sumw2();
  thetaY_res2->Sumw2();
  thetaX_divres2->Divide(thetaX_res2, thetaX_all);
  thetaY_divres2->Divide(thetaY_res2, thetaY_all);

  thetaX_all->GetXaxis()->SetLabelSize(0.05);
  thetaX_all->GetXaxis()->SetTitleSize(0.05);
  thetaX_all->GetYaxis()->SetLabelSize(0.05);
  thetaX_all->GetYaxis()->SetTitleSize(0.05);
  thetaX_acc->GetXaxis()->SetLabelSize(0.05);
  thetaX_acc->GetXaxis()->SetTitleSize(0.05);
  thetaX_acc->GetYaxis()->SetLabelSize(0.05);
  thetaX_acc->GetYaxis()->SetTitleSize(0.05);
  thetaX_res->GetXaxis()->SetLabelSize(0.05);
  thetaX_res->GetXaxis()->SetTitleSize(0.05);
  thetaX_res->GetYaxis()->SetLabelSize(0.05);
  thetaX_res->GetYaxis()->SetTitleSize(0.05);
  thetaX_res2->GetXaxis()->SetLabelSize(0.05);
  thetaX_res2->GetXaxis()->SetTitleSize(0.05);
  thetaX_res2->GetYaxis()->SetLabelSize(0.05);
  thetaX_res2->GetYaxis()->SetTitleSize(0.05);

  thetaX_divacc->GetXaxis()->SetLabelSize(0.05);
  thetaX_divacc->GetXaxis()->SetTitleSize(0.05);
  thetaX_divacc->GetYaxis()->SetLabelSize(0.05);
  thetaX_divacc->GetYaxis()->SetTitleSize(0.05);
  thetaX_divres->GetXaxis()->SetLabelSize(0.05);
  thetaX_divres->GetXaxis()->SetTitleSize(0.05);
  thetaX_divres->GetYaxis()->SetLabelSize(0.05);
  thetaX_divres->GetYaxis()->SetTitleSize(0.05);
  thetaX_divres2->GetXaxis()->SetLabelSize(0.05);
  thetaX_divres2->GetXaxis()->SetTitleSize(0.05);
  thetaX_divres2->GetYaxis()->SetLabelSize(0.05);
  thetaX_divres2->GetYaxis()->SetTitleSize(0.05);


  thetaY_all->GetXaxis()->SetLabelSize(0.05);
  thetaY_all->GetXaxis()->SetTitleSize(0.05);
  thetaY_all->GetYaxis()->SetLabelSize(0.05);
  thetaY_all->GetYaxis()->SetTitleSize(0.05);
  thetaY_acc->GetXaxis()->SetLabelSize(0.05);
  thetaY_acc->GetXaxis()->SetTitleSize(0.05);
  thetaY_acc->GetYaxis()->SetLabelSize(0.05);
  thetaY_acc->GetYaxis()->SetTitleSize(0.05);
  thetaY_res->GetXaxis()->SetLabelSize(0.05);
  thetaY_res->GetXaxis()->SetTitleSize(0.05);
  thetaY_res->GetYaxis()->SetLabelSize(0.05);
  thetaY_res->GetYaxis()->SetTitleSize(0.05);
  thetaY_res2->GetXaxis()->SetLabelSize(0.05);
  thetaY_res2->GetXaxis()->SetTitleSize(0.05);
  thetaY_res2->GetYaxis()->SetLabelSize(0.05);
  thetaY_res2->GetYaxis()->SetTitleSize(0.05);

  thetaY_divacc->GetXaxis()->SetLabelSize(0.05);
  thetaY_divacc->GetXaxis()->SetTitleSize(0.05);
  thetaY_divacc->GetYaxis()->SetLabelSize(0.05);
  thetaY_divacc->GetYaxis()->SetTitleSize(0.05);
  thetaY_divres->GetXaxis()->SetLabelSize(0.05);
  thetaY_divres->GetXaxis()->SetTitleSize(0.05);
  thetaY_divres->GetYaxis()->SetLabelSize(0.05);
  thetaY_divres->GetYaxis()->SetTitleSize(0.05);
  thetaY_divres2->GetXaxis()->SetLabelSize(0.05);
  thetaY_divres2->GetXaxis()->SetTitleSize(0.05);
  thetaY_divres2->GetYaxis()->SetLabelSize(0.05);
  thetaY_divres2->GetYaxis()->SetTitleSize(0.05);

  outfile->cd();
  thetaX_all->Write();
  thetaY_all->Write();
  thetaX_acc->Write();
  thetaY_acc->Write();
  thetaX_res->Write();
  thetaY_res->Write();
  thetaX_res2->Write();
  thetaY_res2->Write();

  thetaX_divacc->Write();
  thetaY_divacc->Write();
  thetaX_divres->Write();
  thetaY_divres->Write();
  thetaX_divres2->Write();
  thetaY_divres2->Write();

  delete thetaX_all;
  delete thetaY_all;
  delete thetaX_acc;
  delete thetaY_acc;
  delete thetaX_res;
  delete thetaY_res;
  delete thetaX_res2;
  delete thetaY_res2;
  delete thetaX_divacc;
  delete thetaY_divacc;
  delete thetaX_divres;
  delete thetaY_divres;
  delete thetaX_divres2;
  delete thetaY_divres2;
}

Json::Value MCSAnalysis::SetupConfig(int verbose_level, std::string geometryfile) {
  std::cerr << "Running with verbose level " << verbose_level << std::endl;
  Json::Value config(Json::objectValue);
  config["maximum_module_depth"] = 50;
  config["check_volume_overlaps"] = true;
  config["reconstruction_geometry_filename"] = geometryfile;
  config["simulation_geometry_filename"] = geometryfile;
  config["simulation_geometry_debug"] = false;
  config["maximum_number_of_steps"] = 10000;
  config["will_do_stack_trace"] = true;
  config["keep_tracks"] = true;
  config["keep_steps"] = true;
  config["verbose_level"] = verbose_level;
  config["log_level"] = 1;
  config["geant4_visualisation"] = false;
  config["physics_model"] = "QGSP_BERT";
  config["reference_physics_processes"] = "mean_energy_loss";
  config["physics_processes"] = "standard";
  config["particle_decay"] = true;
  config["polarised_decay"] = true;
  config["charged_pion_half_life"] = -1.;
  config["muon_half_life"] = -1.;
  config["production_threshold"] = 0.5;
  config["fine_grained_production_threshold"] = Json::Value(Json::objectValue);
  config["default_keep_or_kill"] = true;
  config["spin_tracking"] = true;
  config["keep_or_kill_particles"] = "{\"neutron\":False}";
  config["kinetic_energy_threshold"] = 0.1;
  config["max_step_length"] = 100.;
  config["max_track_time"] = 1.e9;
  config["max_track_length"] = 1.e8;
  config["simulation_reference_particle"] = JsonWrapper::StringToJson(
			      std::string("{\"position\":{\"x\":0.0,\"y\":-0.0,\"z\":-5500.0},")+
			      std::string("\"spin\":{\"x\":0.0,\"y\":-0.0,\"z\":1.0},")+
			      std::string("\"momentum\":{\"x\":0.0,\"y\":0.0,\"z\":1.0},")+
			      std::string("\"particle_id\":-13,\"energy\":226.0,\"time\":0.0,")+
			      std::string("\"random_seed\":10}")
			      );
  config["stepping_algorithm"] = "ClassicalRK4";
  config["delta_one_step"] = -1.;
  config["delta_intersection"] = -1.;
  config["epsilon_min"] = -1.;
  config["epsilon_max"] = -1.;
  config["miss_distance"] = -1.;
  config["everything_special_virtual"] = false;
  config["field_tracker_absolute_error"] = 1.e-4;
  config["field_tracker_relative_error"] = 1.e-4;
  config["data_maximum_reference_count"] = 200;
  return config;
}

Vars MCSAnalysis::PropagateVarsMu(Vars event, double z0){
  
  Vars prop;
  double mass = 105.658; // muon mass [MeV/c^2]
  double time = event.TOFX; // time in lab frame [ns]
  double x = event.X; // horizontal position [mm]
  double y = event.Y; // vertical position [mm]
  double z = event.Z; // longitudinal position [mm]
  double px = event.px; // horizontal momentum component [MeV/c]
  double py = event.py; // vertical momentum component [MeV/c]
  double pz = event.pz;
  
  double energy = std::sqrt(pz*pz + px*px + py*py + mass*mass); // Total energy [MeV]
  
  double event_vector[8] = { time, x, y, z, energy, px, py, pz };
  
  BTField* field = dynamic_cast<BTField*>(MAUS::Globals::GetMCFieldConstructor());
  if (pz!=0){
    try {
//      std::cerr<<" - Calling GlobalTools propagate";
      MAUS::GlobalTools::propagate(event_vector, z0, field, 10., 
				 MAUS::DataStructure::Global::kMuPlus, 1, 1);
//      std::cerr<<" - Completed GlobalTools propagate";
    } catch (...){
//       std::cerr<<"**********Propagate failed**********"<<std::endl;
       event.TOF12 = -1;
    }
  }
  if ((event_vector[3]-z0)>.1){
//    std::cerr<<"Not Propagated to centre of Absorber"<<std::endl;
    prop = reset_Vars();
  }else{
    prop.X = event_vector[1];
    prop.Y = event_vector[2];
    prop.Z = event_vector[3];
    prop.px = event_vector[5];
    prop.py = event_vector[6];
    prop.pz = event_vector[7];
    prop.dXdz = event_vector[5]/event_vector[7];
    prop.dYdz = event_vector[6]/event_vector[7];
    prop.TOF12 = event.TOF12;
    prop.TOF01 = event.TOF01;
    prop.TOFX = event_vector[0];
  }
/*
  //straight track
  Vars proj;
  proj.X = px/pz * (z0 - z) + x;
  proj.Y = py/pz * (z0 - z) + y;
  proj.Z = z0;
  proj.px = px;
  proj.py = py;
  proj.pz = pz;
  proj.dXdz = px/pz;
  proj.dYdz = py/pz;
  proj.TOF12 = event.TOF12;
  proj.TOF01 = event.TOF01;

  Vars diffpp;
  diffpp.X = prop.X-proj.X;
  diffpp.Y = prop.Y-proj.Y;
  diffpp.Z = prop.Z-proj.Z;
  diffpp.px = prop.px-proj.px;
  diffpp.py = prop.py-proj.py;
  diffpp.pz = prop.pz-proj.pz;
  
//  if (diffpp.X > 1e-10){
//      std::cerr<<" X: "<<event.X<<" Y: "<<event.Y<<" Z: "<<event.Z<<" px: "<<event.px<<" py: "<<event.py<<" pz: "<<event.pz<<" z0: "<<z0<<std::endl;
//      std::cerr<<" X: "<<prop.X<<" Y: "<<prop.Y<<" Z: "<<prop.Z<<" px: "<<prop.px<<" py: "<<prop.py<<" pz: "<<prop.pz<<std::endl;
//      std::cerr<<" X: "<<proj.X<<" Y: "<<proj.Y<<" Z: "<<proj.Z<<" px: "<<proj.px<<" py: "<<proj.py<<" pz: "<<proj.pz<<std::endl;
//      std::cerr<<" X: "<<diffpp.X<<" Y: "<<diffpp.Y<<" Z: "<<diffpp.Z<<" px: "<<diffpp.px<<" py: "<<diffpp.py<<" pz: "<<diffpp.pz<<std::endl;

//    exit(0);}
*/
  return prop;
}

Vars MCSAnalysis::PropagateVarsMuVerbose(Vars event, double z0){
//  std::cerr<<"Verbose Propagate"<<std::endl;
//  display_Vars(event);
  Vars prop;
  double mass = 105.658; // muon mass [MeV/c^2]
  double time = event.TOFX; // time in lab frame [ns]
  double x = event.X; // horizontal position [mm]
  double y = event.Y; // vertical position [mm]
  double z = event.Z; // longitudinal position [mm]
  double px = event.px; // horizontal momentum component [MeV/c]
  double py = event.py; // vertical momentum component [MeV/c]
  double pz = event.pz;
  
  double energy = std::sqrt(pz*pz + px*px + py*py + mass*mass); // Total energy [MeV]
  
  double event_vector[8] = { time, x, y, z, energy, px, py, pz };
  double PropagationDirection;
  if (z<z0){
    PropagationDirection=1.0;
  } else {
    PropagationDirection=-1.0;
  }
  BTField* field = dynamic_cast<BTField*>(MAUS::Globals::GetMCFieldConstructor());
  if (pz!=0){
    std::cerr<<"Starting Z: "<<event_vector[3]<<" Target Z: "<<z0<<std::endl;
    while ((PropagationDirection*(z0-event_vector[3]))>.1){
      double z1=event_vector[3]+(PropagationDirection*50);
      if ((PropagationDirection==1.0&&z1>z0)||(PropagationDirection==-1.0&&z1<z0)){
        z1=z0;
      }
      try {
//      std::cerr<<" - Calling GlobalTools propagate.";
        MAUS::GlobalTools::propagate(event_vector, z1, field, 10., 
				 MAUS::DataStructure::Global::kMuPlus, 1, 1);
//      std::cerr<<" - Completed GlobalTools propagate ";
      } catch (...){
         std::cerr<<"**********Propagate failed**********"<<std::endl;
         event.TOF12 = -1;
      }
      if (abs(event_vector[3]-z1)>.1){
        std::cerr<<"Error in propagation"<<std::endl;
        break;
      } else {
        double total_momentum=sqrt(event_vector[5]*event_vector[5]+event_vector[6]*event_vector[6]+event_vector[7]*event_vector[7]);
        double radius=sqrt(event_vector[1]*event_vector[1]+event_vector[2]*event_vector[2]);
        std::cerr<<event_vector[3]<<" "<<event_vector[1]<<" "<<event_vector[2]<<" "<<radius<<" "<<event_vector[5]<<" "<<event_vector[6]<<" "<<event_vector[7]<<" "<<total_momentum<<std::endl;        
      }
    }
  }
  if (abs(event_vector[3]-z0)>.1){
    std::cerr<<"Not Propagated to Desired axial position: "<<z0<<std::endl;
    prop = reset_Vars();
  }else{
    prop.X = event_vector[1];
    prop.Y = event_vector[2];
    prop.Z = event_vector[3];
    prop.px = event_vector[5];
    prop.py = event_vector[6];
    prop.pz = event_vector[7];
    prop.dXdz = event_vector[5]/event_vector[7];
    prop.dYdz = event_vector[6]/event_vector[7];
    prop.TOF12 = event.TOF12;
    prop.TOF01 = event.TOF01;
    prop.TOFX = event_vector[0];
  }
/*
  //straight track
  Vars proj;
  proj.X = px/pz * (z0 - z) + x;
  proj.Y = py/pz * (z0 - z) + y;
  proj.Z = z0;
  proj.px = px;
  proj.py = py;
  proj.pz = pz;
  proj.dXdz = px/pz;
  proj.dYdz = py/pz;
  proj.TOF12 = event.TOF12;
  proj.TOF01 = event.TOF01;

  Vars diffpp;
  diffpp.X = prop.X-proj.X;
  diffpp.Y = prop.Y-proj.Y;
  diffpp.Z = prop.Z-proj.Z;
  diffpp.px = prop.px-proj.px;
  diffpp.py = prop.py-proj.py;
  diffpp.pz = prop.pz-proj.pz;
  
//  if (diffpp.X > 1e-10){
//      std::cerr<<" X: "<<event.X<<" Y: "<<event.Y<<" Z: "<<event.Z<<" px: "<<event.px<<" py: "<<event.py<<" pz: "<<event.pz<<" z0: "<<z0<<std::endl;
//      std::cerr<<" X: "<<prop.X<<" Y: "<<prop.Y<<" Z: "<<prop.Z<<" px: "<<prop.px<<" py: "<<prop.py<<" pz: "<<prop.pz<<std::endl;
//      std::cerr<<" X: "<<proj.X<<" Y: "<<proj.Y<<" Z: "<<proj.Z<<" px: "<<proj.px<<" py: "<<proj.py<<" pz: "<<proj.pz<<std::endl;
//      std::cerr<<" X: "<<diffpp.X<<" Y: "<<diffpp.Y<<" Z: "<<diffpp.Z<<" px: "<<diffpp.px<<" py: "<<diffpp.py<<" pz: "<<diffpp.pz<<std::endl;

//    exit(0);}
*/
//  std::cerr<<"End of Verbose Propagation"<<std::endl;
  return prop;
}


Vars MCSAnalysis::FillVars(TLorentzVector a_track_mom, TLorentzVector a_track_pos, int pid){
  Vars tempVars=reset_Vars();
  tempVars.X = a_track_pos.X();
  tempVars.Y = a_track_pos.Y();
  tempVars.Z = a_track_pos.Z();
  tempVars.px=  a_track_mom.X();
  tempVars.py = a_track_mom.Y();
  tempVars.pz = a_track_mom.Z();
  tempVars.dXdz = tempVars.px/tempVars.pz;
  tempVars.dYdz = tempVars.py/tempVars.pz;
  tempVars.isgood=true;
  tempVars.pid=pid;
  if( int(tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()) > 0 &&
      int(tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray().size()) > 0)
      tempVars.TOF12 = tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray()[0].GetTime() 
      - tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime();
  else{
//                                  std::cerr<<"TOF1 Array size: "<<tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()
//                                           <<"TOF2 Array size: "<<tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray().size()<<std::endl;
    tempVars.TOF12 = 100*8.22475 / 0.299792458;
  }
  if( int(tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray().size()) > 0 &&
        int(tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()) > 0)
      tempVars.TOF01 = tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime() 
        - tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[0].GetTime();
  else{
//                                  std::cerr<<"TOF0 Array size: "<<tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray().size()
//                                           <<"TOF1 Array size: "<<tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()<<std::endl;

    tempVars.TOF01 = 100*7.64186 / 0.299792458;
  }
  tempVars.TOFX = a_track_pos.T();
  return tempVars;
}

Vars MCSAnalysis::FillVars3(MAUS::ThreeVector b_track_mom, MAUS::ThreeVector b_track_pos, int pid){
  Vars tempVars=reset_Vars();
  tempVars.X = b_track_pos.x();
  tempVars.Y = b_track_pos.y();
  tempVars.Z = b_track_pos.z();
  tempVars.px=  b_track_mom.x();
  tempVars.py = b_track_mom.y();
  tempVars.pz = b_track_mom.z();
  tempVars.dXdz = tempVars.px/tempVars.pz;
  tempVars.dYdz = tempVars.py/tempVars.pz;
  tempVars.isgood=true;
  tempVars.pid=pid;
  if( int(tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()) > 0 &&
      int(tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray().size()) > 0)
      tempVars.TOF12 = tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray()[0].GetTime() 
      - tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime();
  else{
//                                  std::cerr<<"TOF1 Array size: "<<tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()
//                                           <<"TOF2 Array size: "<<tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray().size()<<std::endl;
    tempVars.TOF12 = 100*8.22475 / 0.299792458;
  }
  if( int(tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray().size()) > 0 &&
        int(tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()) > 0)
      tempVars.TOF01 = tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime() 
        - tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[0].GetTime();
  else{
//                                  std::cerr<<"TOF0 Array size: "<<tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray().size()
//                                           <<"TOF1 Array size: "<<tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()<<std::endl;

    tempVars.TOF01 = 100*7.64186 / 0.299792458;
  }
  tempVars.TOFX = 0;
  return tempVars;
}


//Vars Find_Track_Point(TLorentzVector pos, TLorentzVector mom, MAUS::DataStructure::Global::DetectorPoint Detector){
//}


multiVars MCSAnalysis::read_globals(double DS11, double DS53, double US53, double US11, double diffuserlow, double diffuserhigh, double TOF0pos, double TOF1pos){
// This is the main function to read the data from the root files
  double dz = 0.2;
  multiVars tempMV;
  tempMV.UScentre_absorber=reset_Vars();
  tempMV.USend_of_DStracker=reset_Vars();
  tempMV.USTOF0=reset_Vars();
  tempMV.USTOF1=reset_Vars();
  tempMV.DScentre_absorber=reset_Vars();
  tempMV.UStrackerUS11=reset_Vars();
  tempMV.UScentre_absorber.Z=_sys["abspos"];
  tempMV.USend_of_DStracker.Z= DS53;
  tempMV.USTOF0.Z= TOF0pos;
  tempMV.USTOF1.Z= TOF1pos;
  tempMV.DScentre_absorber.Z=_sys["abspos"];
  tempMV.UStrackerUS11.Z= US11;
  tempMV.UStrackerUS53.Z= US53;
     
    
//    std::cerr << "Globals..."<<std::endl;
  MAUS::DataStructure::Global::DetectorPoint virtual_plane = MAUS::DataStructure::Global::kVirtual;
  std::vector<MAUS::DataStructure::Global::PrimaryChain*>* pchains = globalevent->get_primary_chains();
//    std::vector<MAUS::DataStructure::Global::PrimaryChain*>::iterator pchains_iterator;


//    MAUS::DataStructure::Global::DetectorPoint diffuser_plane = MAUS::DataStructure::Global::kVirtual;
  int number_of_chains = 0;

  for (auto global_chain_iter = pchains->begin(); global_chain_iter != pchains->end(); ++global_chain_iter)
  {
     number_of_chains++;
  }
//---------------------------------------------------------------------------------------
// Read data from Globals. If available
//---------------------------------------------------------------------------------------
  if (number_of_chains>0){
//    std::cerr<<"Number of chains: " << number_of_chains << std::endl;

     for(auto pchains_iterator = pchains->begin(); pchains_iterator != pchains->end(); ++ pchains_iterator){
        MAUS::DataStructure::Global::PrimaryChain* chain = (*pchains_iterator);
//        std::cerr << "iterating? Chain: "<< chain->get_chain_type() <<" looking for "<< chain_type <<std::endl;
        std::cerr<<"Chain type: "<<chain->get_chain_type()<<std::endl;

// If available read for chain type 5 first

        if(chain->get_chain_type() == 5){
//            MAUS::DataStructure::Global::PrimaryChain* USchain=chain->GetUSDaughter();

//Chain 5 - Use Get data from Upstream track from primary chain

            std::vector<MAUS::DataStructure::Global::Track*> some_tracks = chain->GetMatchedTracks();
            std::vector<MAUS::DataStructure::Global::Track*>::iterator some_track_iterator;
//            std::cerr<<"Upstream Track - Through"<<std::endl;
            int track_number = 0;
            for(some_track_iterator = some_tracks.begin(); some_track_iterator < some_tracks.end(); ++some_track_iterator){
                track_number++;
                MAUS::DataStructure::Global::Track* a_track = (*some_track_iterator);
//                std::vector<MAUS::DataStructure::Global::DetectorPoint> some_detector_points = a_track->GetDetectorPoints();
//                std::cerr<<"Detector Points"<<std::endl;
//                for(int q=0;q<some_detector_points.size();q++){
//                  std::cerr<<"Detector Point: "<<some_detector_points.at(q)<<std::endl;;
                  if(a_track != NULL){
//                    std::cerr<<"Track Number: "<<track_number<<" Track Pid: "<<a_track->get_pid()<<std::endl;
                    if(a_track->get_pid() == -13 || a_track->get_pid() == 13){
                        std::vector< const MAUS::DataStructure::Global::TrackPoint * > track_points = a_track->GetTrackPoints(virtual_plane);
                        for(unsigned int p = 0; p < track_points.size(); p++){
                            TLorentzVector a_track_mom = track_points.at(p)->get_momentum();
                            TLorentzVector a_track_pos = track_points.at(p)->get_position();
//                            if ( track_points.at(p)->get_detector() == 1){
//                                std::cerr << "Detector type: " << track_points.at(p)->get_detector() << std::endl;
//                            }

//                            std::cerr<<"Time: "<< a_track_pos.T()<<std::endl;
                            std::cerr << " Looking for z=" << tempMV.UScentre_absorber.Z <<" and have z=" << a_track_pos.Z() << std::endl;
//                                      " or "<< tempMV.end_of_DStracker.Z <<" or "
//                                      << tempMV.Diffuser.Z <<" or "<< tempMV.TOF0.Z <<" and have z=" << a_track_pos.Z() << std::endl;
                            if(a_track_pos.Z() >= tempMV.UScentre_absorber.Z-dz && a_track_pos.Z() <= tempMV.UScentre_absorber.Z+dz){
//                                std::cerr << "-------------- FOUND Centre of Absorber!!!"<<std::endl;
//                                std::cerr << "Detector type: " << track_points.at(p)->get_detector() << std::endl;
                                tempMV.UScentre_absorber=FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                                std::cerr<<"Found global data from US track at centre of absorber"<<std::endl;
                                display_Vars(tempMV.UScentre_absorber);
                            } else if(a_track_pos.Z() >= tempMV.USend_of_DStracker.Z-dz && a_track_pos.Z() <= tempMV.USend_of_DStracker.Z+dz){
//                                std::cerr << "-------------- FOUND END OF DS Tracker!!!"<<std::endl;
//                                std::cerr << "Detector type: " << track_points.at(p)->get_detector() << std::endl;
                                tempMV.USend_of_DStracker = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.USTOF0.Z-dz && a_track_pos.Z() <= tempMV.USTOF0.Z+dz){
//                                std::cerr << "-------------- FOUND TOF0!!!"<<std::endl;
//                                std::cerr << "Detector type: " << track_points.at(p)->get_detector() << std::endl;
                                tempMV.USTOF0 = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.USTOF1.Z-dz && a_track_pos.Z() <= tempMV.USTOF1.Z+dz){
//                                std::cerr << "-------------- FOUND TOF1!!!"<<std::endl;
//                                std::cerr << "Detector type: " << track_points.at(p)->get_detector() << std::endl;
                                tempMV.USTOF1 = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.UStrackerUS11.Z-dz && a_track_pos.Z() <= tempMV.UStrackerUS11.Z+dz){
//                                std::cerr << "-------------- FOUND US 11!!!"<<std::endl;
//                                std::cerr << "Detector type: " << track_points.at(p)->get_detector() << std::endl;
                                tempMV.UStrackerUS11 = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.UStrackerUS53.Z-dz && a_track_pos.Z() <= tempMV.UStrackerUS53.Z+dz){
//                                std::cerr << "-------------- FOUND US 53!!!"<<std::endl;
//                                std::cerr << "Detector type: " << track_points.at(p)->get_detector() << std::endl;
                                tempMV.UStrackerUS53 = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            }
//                            std::cerr<<a_track_pos.Z()<<" "<<a_track_pos.X()<<" "<<a_track_pos.Y()<<std::endl;
//                            if (tempMV.UScentre_absorber.isgood==true && tempMV.USend_of_DStracker.isgood==true  
//                                 && tempMV.USTOF0.isgood==true && tempMV.UStrackerUS11.isgood==true && tempMV.UStrackerUS53.isgood==true) {
//                              std::cerr<<"Found all trackpoints"<<std::endl;
//                              break;
//                            }
                        }
                    }
                  }
//                }
            }

//Chain5 - Get Downstream data from DSDaughter

            MAUS::DataStructure::Global::PrimaryChain* DSchain=chain->GetDSDaughter();
            some_tracks = DSchain->GetMatchedTracks();
//            std::cerr<<"Downstream track - Through"<<std::endl;
            for(some_track_iterator = some_tracks.begin(); some_track_iterator < some_tracks.end(); ++some_track_iterator){
                MAUS::DataStructure::Global::Track* a_track = (*some_track_iterator);
                if(a_track != NULL){
                    if(a_track->get_pid() == -13 || a_track->get_pid() == 13){
                        std::vector< const MAUS::DataStructure::Global::TrackPoint * > track_points = a_track->GetTrackPoints(virtual_plane);
                        for(unsigned int p = 0; p < track_points.size(); p++){
                            TLorentzVector a_track_mom = track_points.at(p)->get_momentum();
                            TLorentzVector a_track_pos = track_points.at(p)->get_position();
//                            std::cerr << " Looking for z=" << tempMV.centre_absorber.Z << " or "<< tempMV.end_of_DStracker.Z <<" or "
//                                      << tempMV.Diffuser.Z <<" or "<< tempMV.TOF0.Z <<" and have z=" << a_track_pos.Z() << std::endl;
                            if((a_track_pos.Z() >= DS11-dz && a_track_pos.Z() <= DS11+dz)&&tempMV.DScentre_absorber.isgood==false){
//                                std::cerr << "-------------- FOUND Centre of Absorber from DS!!!"<<std::endl;
//                                std::cerr<<"Starting propagation from DS tracker to centre of absorber"<<std::endl;
//                                std::cerr<<a_track_pos.Z()<<" "<<a_track_pos.X()<<" "<<a_track_pos.Y()<<std::endl;
                                tempMV.DScentre_absorber=FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                                tempMV.DScentre_absorber=PropagateVarsMu(tempMV.DScentre_absorber,_sys["abspos"]);
                            }else if((a_track_pos.Z() >= DS11+0.65-dz && a_track_pos.Z() <= DS11+0.65+dz)&&tempMV.DScentre_absorber.isgood==false){
//                                std::cerr<<"Starting propagation from DS tracker to centre of absorber"<<std::endl;
//                                std::cerr<<a_track_pos.Z()<<" "<<a_track_pos.X()<<" "<<a_track_pos.Y()<<std::endl;
                                tempMV.DScentre_absorber=FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                                tempMV.DScentre_absorber=PropagateVarsMu(tempMV.DScentre_absorber,_sys["abspos"]);
                            }else if((a_track_pos.Z() >= DS11+1.3-dz && a_track_pos.Z() <= DS11+1.3+dz)&&tempMV.DScentre_absorber.isgood==false){
//                                std::cerr<<"Starting propagation from DS tracker to centre of absorber"<<std::endl;
//                                std::cerr<<a_track_pos.Z()<<" "<<a_track_pos.X()<<" "<<a_track_pos.Y()<<std::endl;
                                tempMV.DScentre_absorber=FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                                tempMV.DScentre_absorber=PropagateVarsMu(tempMV.DScentre_absorber,_sys["abspos"]);
                            }
//                            std::cerr<<a_track_pos.Z()<<" "<<a_track_pos.X()<<" "<<a_track_pos.Y()<<std::endl;
                        }
                    }
                }
           }

//If no chain 5 read data from chain 3, which is an upstream orphan

        }else if(chain->get_chain_type() == 3){
            std::vector<MAUS::DataStructure::Global::Track*> some_tracks = chain->GetMatchedTracks();
            std::vector<MAUS::DataStructure::Global::Track*>::iterator some_track_iterator;
            for(some_track_iterator = some_tracks.begin(); some_track_iterator < some_tracks.end(); ++some_track_iterator){
                MAUS::DataStructure::Global::Track* a_track = (*some_track_iterator);
                if(a_track != NULL){
                    if(a_track->get_pid() == -13 || a_track->get_pid() == 13){
                        std::vector< const MAUS::DataStructure::Global::TrackPoint * > track_points = a_track->GetTrackPoints(virtual_plane);
                        for(unsigned int p = 0; p < track_points.size(); p++){
                            TLorentzVector a_track_mom = track_points.at(p)->get_momentum();
                            TLorentzVector a_track_pos = track_points.at(p)->get_position();
//                            std::cerr << " Looking for z=" << tempMV.UStrackerUS11.Z <<" and have z=" << a_track_pos.Z() << std::endl;
//                            std::cerr << " Looking for z=" << tempMV.centre_absorber.Z << " or "<< tempMV.end_of_DStracker.Z <<" or "
//                                      << tempMV.Diffuser.Z <<" or "<< tempMV.TOF0.Z <<" and have z=" << a_track_pos.Z() << std::endl;
                            if(a_track_pos.Z() >= tempMV.UScentre_absorber.Z-dz && a_track_pos.Z() <= tempMV.UScentre_absorber.Z+dz){
//                                std::cerr << "-------------- FOUND Centre of Absorber!!!"<<std::endl;
//                                std::cerr << "Detector type: " << track_points.at(p)->get_detector() << std::endl;
                                tempMV.UScentre_absorber=FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.USend_of_DStracker.Z-dz && a_track_pos.Z() <= tempMV.USend_of_DStracker.Z+dz){
//                                std::cerr << "-------------- FOUND END OF DS Tracker!!!"<<std::endl;
//                                std::cerr << "Detector type: " << track_points.at(p)->get_detector() << std::endl;
                                tempMV.USend_of_DStracker = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.USTOF0.Z-dz && a_track_pos.Z() <= tempMV.USTOF0.Z+dz){
//                                std::cerr << "-------------- FOUND TOF0!!!"<<std::endl;
//                                std::cerr << "Detector type: " << track_points.at(p)->get_detector() << std::endl;
                                tempMV.USTOF0 = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.USTOF1.Z-dz && a_track_pos.Z() <= tempMV.USTOF1.Z+dz){
//                                std::cerr << "-------------- FOUND TOF1!!!"<<std::endl;
//                                std::cerr << "Detector type: " << track_points.at(p)->get_detector() << std::endl;
                                tempMV.USTOF1 = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.UStrackerUS11.Z-dz && a_track_pos.Z() <= tempMV.UStrackerUS11.Z+dz){
//                                std::cerr << "-------------- FOUND US 11!!!"<<std::endl;
//                                std::cerr << "Detector type: " << track_points.at(p)->get_detector() << std::endl;
                                tempMV.UStrackerUS11 = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.UStrackerUS53.Z-dz && a_track_pos.Z() <= tempMV.UStrackerUS53.Z+dz){
//                                std::cerr << "-------------- FOUND US 11!!!"<<std::endl;
//                                std::cerr << "Detector type: " << track_points.at(p)->get_detector() << std::endl;
                                tempMV.UStrackerUS53 = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            }
                            if (tempMV.UScentre_absorber.isgood==true && tempMV.USend_of_DStracker.isgood==true  
                                 && tempMV.USTOF0.isgood==true && tempMV.UStrackerUS11.isgood==true && tempMV.UStrackerUS53.isgood==true) {
//                              std::cerr<<"Found all trackpoints"<<std::endl;
                              break;
                            }
                        }
                    }
                }
            }
        }
     }
  }else{
//     std::cerr<<"No Globals Data present"<<std::endl;
//If no Globals execute this. First find the data from the US and DS trackers closest to the absorber

     if( scifievent->scifitracks().size() == 1 || 
        scifievent->scifitracks().size() == 2){
        for ( size_t j=0; j<scifievent->scifitracks().size(); j++){
           int npoints = scifievent->scifitracks()[j]->scifitrackpoints().size();
           double maxUS=0.0, minDS=44000;
           int tracker = scifievent->scifitracks()[j]->tracker();
           for ( int k=0; k<npoints; k++){
	      double zpos = scifievent->scifitracks()[j]->scifitrackpoints()[k]->pos().z();
	      if(tracker==0 && zpos > maxUS){
                 MAUS::ThreeVector b_track_pos=scifievent->scifitracks()[j]->scifitrackpoints()[k]->pos();
                 MAUS::ThreeVector b_track_mom=scifievent->scifitracks()[j]->scifitrackpoints()[k]->mom();
                 maxUS=b_track_pos.z();
                 tempMV.UStrackerUS11 = FillVars3(b_track_mom,b_track_pos,-13);
	      }
	      if(tracker==1 && zpos < minDS){      
                 MAUS::ThreeVector b_track_pos=scifievent->scifitracks()[j]->scifitrackpoints()[k]->pos();
                 MAUS::ThreeVector b_track_mom=scifievent->scifitracks()[j]->scifitrackpoints()[k]->mom();
                 minDS=b_track_pos.z();
                 tempMV.DScentre_absorber=FillVars3(b_track_mom,b_track_pos,-13);
	      }
           }
        }

// 
//        display_Vars(tempMV.UStrackerUS11);
        if ((tempMV.UStrackerUS11.isgood==true)&&abs(tempMV.UStrackerUS11.Z-US11)>dz){
           tempMV.UStrackerUS11=PropagateVarsMu(tempMV.UStrackerUS11,US11);
        }
        if ((tempMV.DScentre_absorber.isgood==true)&&abs(tempMV.DScentre_absorber.Z-_sys["abspos"])>dz){
              tempMV.DScentre_absorber=PropagateVarsMu(tempMV.DScentre_absorber,_sys["abspos"]);
        }
     }
  }
  if (tempMV.UStrackerUS11.isgood==true){
     if(tempMV.UScentre_absorber.isgood==false){
        tempMV.UScentre_absorber=PropagateVarsMu(tempMV.UStrackerUS11,_sys["abspos"]);
//        std::cerr<<"No Global Data - Propagated to centre of absorber"<<std::endl;
//        display_Vars(tempMV.UScentre_absorber);
     }
     if(tempMV.USend_of_DStracker.isgood==false){
        tempMV.USend_of_DStracker=PropagateVarsMu(tempMV.UStrackerUS11, DS53);
     }              
     if(tempMV.USTOF0.isgood==false){
        tempMV.USTOF0=PropagateVarsMu(tempMV.UStrackerUS11, TOF0pos);
     }              
     if(tempMV.USTOF1.isgood==false){
        tempMV.USTOF1=PropagateVarsMu(tempMV.UStrackerUS11, TOF1pos);
     }              
     if(tempMV.UStrackerUS53.isgood==false){
        tempMV.UStrackerUS53=PropagateVarsMu(tempMV.UStrackerUS11, US53);
     }              
  }
  return tempMV;
}

double MCSAnalysis::Check_Radius_US(double low_z, double high_z)
{
    double max_rad = 0.0;
//    std::cerr<<"Check Radius Called - ";
    double dz = 0.2;
//    bool found_tracks=false;
//    Vars tempVars=reset_Vars();
    std::vector<MAUS::DataStructure::Global::PrimaryChain*>* pchains = globalevent->get_primary_chains();
    for(auto pchains_iterator = pchains->begin(); pchains_iterator != pchains->end(); ++ pchains_iterator){
        MAUS::DataStructure::Global::PrimaryChain* chain = (*pchains_iterator);
        if((chain->get_chain_type() == 5) || (chain->get_chain_type() == 3) || (chain->get_chain_type() == 1)){
//            std::cerr<<"Found Chain "<< chain->get_chain_type()<<std::endl;
            MAUS::DataStructure::Global::PrimaryChain* USchain;
            if(chain->get_chain_type() == 5){
              USchain=chain->GetUSDaughter();
            }else{
              USchain=chain;
            }
            std::vector<MAUS::DataStructure::Global::Track*> some_tracks = USchain->GetMatchedTracks();
            std::vector<MAUS::DataStructure::Global::Track*>::iterator some_track_iterator;
            for(some_track_iterator = some_tracks.begin(); some_track_iterator < some_tracks.end(); ++some_track_iterator){
                MAUS::DataStructure::Global::Track* a_track = (*some_track_iterator);
                if(a_track != NULL){
//                    std::cerr<<"Tracks Found"<<std::endl;
                    if(a_track->get_pid() == -13 || a_track->get_pid() == 13){
//                        std::cerr<<"Correct pid"<<std::endl;
                        std::vector< const MAUS::DataStructure::Global::TrackPoint * > track_points = a_track->GetTrackPoints();
                        for(unsigned int p = 0; p < track_points.size(); p++){
 //                           std::cerr<<"Scanning track points"<<std::endl;
                            TLorentzVector a_track_pos = track_points.at(p)->get_position();
                            if(a_track_pos.Z() >= low_z-dz && a_track_pos.Z() <= high_z+dz){
//                              found_tracks=true;
                              double track_rad = sqrt((a_track_pos.X()*a_track_pos.X())+(a_track_pos.Y()*a_track_pos.Y()));
//                              std::cerr<<"track point found. Track radius is "<<track_rad<<std::endl;
                              if(track_rad >= max_rad){
                                max_rad = track_rad;
                              }
                            }
                        }
                    }
                }
             }
       }
    }
//    if (!found_tracks) std::cerr<<"Found no tracks for comparison"<<std::endl;
    return max_rad;
}



