#include "MCSAnalysis.h"

#include "TVirtualFFT.h"

MCSAnalysis::MCSAnalysis(std::string tree, std::string mctree, std::string reftree, std::string outname, std::map<std::string, double> histlimits)
//: p_vec(TVectorD(19)), res(TVectorD(30)), pStart_vec(TVectorD(19)),
//    pStart_vec_y(TVectorD(19)), theta_true_x(TVectorD(19)), theta_true_y(TVectorD(19))
{
  chain   = new TChain(tree.c_str());
  refchain   = new TChain(reftree.c_str());
  mcchain = new TChain(mctree.c_str());
  mcemptychain = new TChain(mctree.c_str());
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
  _sys["USref"]  = 15062;
  _sys["DSref"]  = 18849;
  _sys["niter"]   = 10;
  _sys["abspos"]  = 16952.5;
  _sys["diffpos"]  = 13620;
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
  modelname3 = "Moliere";
  USDataname = "USData.txt";
  DSDataname = "DSData.txt";
  USRefname = "USRef.txt";
  DSRefname = "DSRef.txt";
  EndofDStrackerrefname= "EndofDStrackerRef.txt";
  EndofDStrackerdataname= "EndofDStrackerData.txt";
  StartofDStrackerrefname= "StartofDStrackerRef.txt";
  StartofDStrackerdataname= "StartofDStrackerData.txt";
  Collect_USTOF0dataname= "Collect_USTOF0Data.txt";
  Collect_USTOF0refname= "Collect_USTOF0Ref.txt";
  Collect_USTOF1dataname= "Collect_USTOF1Data.txt";
  Collect_USTOF1refname= "Collect_USTOF1Ref.txt";
  TOF0dataname="TOF0data.txt";
  TOF1dataname="TOF1data.txt";
  TOF2dataname="TOF2data.txt";
  TOF0refname="TOF0ref.txt";
  TOF1refname="TOF1ref.txt";
  TOF2refname="TOF2ref.txt";
  Chi2dataname="Chi2data.txt";
  Chi2refname="Chi2ref.txt";
  dofdataname="dofdata.txt";
  dofrefname="dofref.txt";
  Collect_StartofUSTdataname="Collect_StartofUSTdata.txt";
  Collect_StartofUSTrefname="Collect_StartofUSTref.txt";
  Collect_EndofUSTdataname="Collect_EndofUSTdata.txt";
  Collect_EndofUSTrefname="Collect_EndofUSTref.txt"; 
  Collect_DiffuserStartdataname="Collect_DiffuserStartdata.txt"; 
  Collect_DiffuserEnddataname="Collect_DiffuserEnddata.txt"; 
  Collect_DiffuserStartrefname="Collect_DiffuserStartref.txt"; 
  Collect_DiffuserEndrefname="Collect_DiffuserEndref.txt"; 
 
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

  mresp_thetaX.Setup(_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  mresp_thetaY.Setup(_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  mresp_thetaScatt.Setup(_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
  mresp_theta2Scatt.Setup(_histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);
  
  fresp_thetaX.Setup(_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  fresp_thetaY.Setup(_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  fresp_thetaScatt.Setup(_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
  fresp_theta2Scatt.Setup(_histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);
  const int NBINS = 19;
  const int NPBINS = 30;
  double scat_bin_array[NBINS+1] =
    {-0.1151, -0.0938, -0.0754, -0.0597,  -0.0463,
     -0.0347, -0.0248, -0.0162, -0.00895, -0.00269,
     0.00269,  0.00895, 0.0162,  0.0248,   0.0347,
     0.0463,   0.0597,  0.0754,  0.0938,   0.1151};
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
  meanp = 300;
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


  fiducial_cut_radius = 300.0;
  sigmap = 0.035;
  binlimit = 22;
//  int counter = 0;
//  double semom = 0;
//  double rmsmom = 0;
//  double errrmsmom = 0;
//  double bwmom = 0;
  // histograms of selection variables
  
  diffradius = new TH1D("diffradius","diffradius", 500, 0, 500);
  projradius = new TH1D("projradius","projradius", 500, 0, 500);
  trackno = new TH2D("trackno","trackno",5,0,4,5,0,4);
  tofhitno = new TH2D("tofhitno","tofhitno",5,0,4,5,0,4);
  pathlengthabs = new TH1D("pathlengthabs","pathlengthabs", 500, 0, 500);
  t_cor = new TH2D("t_cor","t_cor", 1000, 0, 1000,100, 8.05, 9.15);
  TOF0Energy = new TH1D("TOF0Energy","TOF0Energy", 200, 150, 350);
  TOF1Energy = new TH1D("TOF1Energy","TOF1Energy", 200, 150, 350);
  TOF2Energy = new TH1D("TOF2Energy","TOF2Energy", 200, 150, 350);
  difEloss = new TH1D("difEloss","difEloss", 10, 0, 10);
  energylossproj = new TH1D("energylossproj","energylossproj", 40, 5000, 24000);
  //energylosspro = new TH1D("energylosspro","energylosspro", 40, 5000, 24000);
  residual01 = new TH1D("residual01","residual01", 36, -50, 100);
  residual12 = new TH1D("residual12","residual12", 36, -50, 100);
  residual12p6 = new TH1D("residual12","residual12", 36, -50, 100);
  residual01j = new TH1D("residual12","residual01", 36, -50, 100);
  residual12p2 = new TH1D("residual12","residual12", 36, -50, 100);
  residual01p6 = new TH1D("residual12","residual12", 36, -50, 100);
  residual01short = new TH1D("residual01short","residual01short", 36, -50, 100);
  residualcobb = new TH1D("residualcobb","residualcobb", 36, -50, 100);
  pzdEdx = new TH2D("pzdEdx","pzdEdz", 100, 1, 100,250, -500, 500);
  TOFcom = new TH2D("TOF01 vs TOF12","TOF01 vs TOF12", 33, 100, 350,33, 100, 350);
  energyloss = new TH2D("energyloss","energyloss", 40, 5000, 24000,350, 100, 450);
  TOF01vsMCTruth = new TH2D("TOF01vsMCTruth","TOF01vsMCTruth", 33, 100, 350,33, 100, 350);
  TOF12vsMCTruth = new TH2D("TOF12vsMCTruth","TOF12vsMCTruth", 33, 100, 350,33, 100, 350);
  TOF12Paul6thforvsMCTruth = new TH2D("TOF12forupvsMCTruth","TOF12forupvsMCTruth", 33, 100, 350, 33, 100, 350);
  TOF12longPaul6thforvsMCTruth = new TH2D("TOF12longforupvsMCTruth","TOF12longforupvsMCTruth", 33, 100, 350, 33, 100, 350);
  TOF12Paul2ndforvsMCTruth = new TH2D("TOF01forupvsMCTruth","TOF01forupvsMCTruth", 33, 100, 350, 33, 100, 350);
  TOF01Paul6thforvsMCTruth = new TH2D("TOF01forupvsMCTruth","TOF01forupvsMCTruth", 33, 100, 350, 33, 100, 350);
  TOF01shortPaul6thforvsMCTruth = new TH2D("TOF01shortforupvsMCTruth","TOF01shortforupvsMCTruth", 33, 100, 350, 33, 100, 350);
  TOF01fordownvsMCTruth = new TH2D("TOF01fordownvsMCTruth","TOF01fordownvsMCTruth", 33, 100, 350,33, 100, 350);
  TOF12fordownvsMCTruth = new TH2D("TOF12fordownvsMCTruth","TOF12fordownvsMCTruth", 33, 100, 350,33, 100, 350);
  TOF12cobbvsMCTruth = new TH2D("TOF12cobbvsMCTruth","TOF12cobbvsMCTruth", 33, 100, 350,33, 100, 350);
  TOF01forupvsTOF01fordown = new TH2D("TOF01forupvsTOF01fordown","TOF01forupvsTOF01fordown", 33, 100, 350,33, 100, 350);
  tof10 = new TH1D("tof10","TOF Between Stations 1 and 0; t_{TOF1} - t_{TOF0} (ns)", 150, 10, 40);
  tof10_sel = new TH1D("tof10_sel","TOF Selection Between Stations 1 and 0; t_{TOF1} - t_{TOF0} (ns)", 150, 10, 40);
  tof21 = new TH1D("tof21","TOF Between Stations 2 and 1; t_{TOF2} - t_{TOF1} (ns)", 150, 0, 50);
  tof21_sel = new TH1D("tof21_sel","TOF Selection Between Stations 2 and 1; t_{TOF2} - t_{TOF1} (ns)", 150, 0, 50);
  mctof10 = new TH1D("mctof10","TOF Between Stations 1 and 0; t_{TOF1} - t_{TOF0} (ns)", 150, 10, 40);
  mctof10_sel = new TH1D("mctof10_sel","TOF Selection Between Stations 1 and 0; t_{TOF1} - t_{TOF0} (ns)", 150, 10, 40);
  mctof21 = new TH1D("mctof21","TOF Between Stations 2 and 1; t_{TOF2} - t_{TOF1} (ns)", 150, 0, 50);
  mctof21_sel = new TH1D("mctof21_sel","TOF Selection Between Stations 2 and 1; t_{TOF2} - t_{TOF1} (ns)", 150, 0, 50);
  calc_mom = new TH1D("calc_mom","Momentum Calculated from TOF; Momentum (MeV/c)", 100, 0, 400);
  cor_mom = new TH1D("cor_mom","Cor Momentum Calculated from TOF; Momentum (MeV/c)", 100, 0, 400);
  mccalc_mom = new TH1D("mccalc_mom","Momentum Calculated from TOF; Momentum (MeV/c)", 400, 0, 400);
  mctrue_mom = new TH1D("mctrue_mom","Momentum from Virtual Planes; Momentum (MeV/c)", 400, 0, 400);
  Bfield_DStrackerX = new TH1D("BfieldDSTrackerX","BFieldDSTracker; Axial Position (m)",1000, _sys["dataDS11"], _sys["dataDS53"]);
  Bfield_DStrackerY = new TH1D("BfieldDSTrackerY","BFieldDSTracker; Axial Position (m)",1000, _sys["dataDS11"], _sys["dataDS53"]);
  Bfield_DStrackerZ = new TH1D("BfieldDSTrackerZ","BFieldDSTracker; Axial Position (m)",1000, _sys["dataDS11"], _sys["dataDS53"]);
  

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
  pass_Fid_old = new TH1D("datapass_Fid_old",";Passed Cut; Events",100, 0, 200);
  pass_ExtTKU0 = new TH1D("datapass_ExtTKU0",";Passed Cut; Events",3, 0, 3);
  pass_TOF01 = new TH1D("datapass_TOF01",";Passed Cut; Events",100, 23, 35);
  pass_TOF01Ex = new TH1D("datapass_TOF01Ex",";Passed Cut; Events",100, 23, 35);
  pass_mom = new TH1D("datapass_mom",";Passed Cut; Events",100, 100, 300);
  pass_chi2= new TH1D("datapass_chi2",";Chi2/dof; Events",50, 0, 10);

  Hist_TOF0Require1SP = new TH1D("Hist Data require 1 TOF0 SP",";Number of SP; Events", 10, 0, 10);
  Hist_TOF1Require1SP = new TH1D("Hist Data require 1 TOF1 SP",";Number of SP; Events", 10, 0, 10);
  Histref_TOF0Require1SP = new TH1D("Hist Ref require 1 TOF0 SP",";Number of SP; Events", 10, 0, 10);
  Histref_TOF1Require1SP = new TH1D("Hist Ref require 1 TOF1 SP",";Number of SP; Events", 10, 0, 10);
  Hist_TOF01 = new TH1D("Hist Data_TOF01",";TOF01 time; Events",100, 23, 40);
  Histref_TOF01 = new TH1D("Hist Ref TOF01",";TOF01 time; Events",100, 23, 40);
  Hist_Chi2 = new TH1D("Hist Data Chi2",";Chi2; Events",100, 0, 15);
  Histref_Chi2 = new TH1D("Hist Ref Chi2",";chi2; Events",100,0, 15);
  Hist_TKU_max_radius = new TH1D("Hist Data TKU_max_radius",";TKU_max_radius; Events",100, 0, 250);
  Histref_TKU_max_radius = new TH1D("Hist Ref TKU_max_radius",";TKU_max_radius; Events",100, 0, 250);
  Hist_Diff_max_radius = new TH1D("Hist Data Diff_max_radius",";Diff_max_radius; Events",100, 0, 200);
  Histref_Diff_max_radius = new TH1D("Hist Ref Diff_max_radius",";Diff_max_radius; Events",100, 0, 200);
  Hist_TOF01Ext = new TH1D("Hist Data_TOF01Ext",";TOF01 time; Events",100, 23, 35);
  Histref_TOF01Ext = new TH1D("Hist Ref_TOF01Ext",";TOF01 time; Events",100, 23, 35);
  Hist_USTOF0 = new TH1D("Hist Data USTOF0",";Successful Extrapolation; Events", 4, 0, 4);
  Histref_USTOF0 = new TH1D("Hist Ref USTOF0",";Successful Extrapolation; Events", 4, 0, 4);
  Hist_Fiducial = new TH1D("Hist Data Fiducial",";Fiducial radius; Events",100, 0, 250);
  Histref_Fiducial = new TH1D("Hist Ref Fiducial",";Fiducial radius; Events",100, 0, 250);
  Hist_MOM =  new TH1D("Hist Data Mom",";Absorber Momentum; Events",100, 0, 250);
  Histref_MOM =  new TH1D("Hist Ref Mom",";Absorber Momentum; Events",100, 0, 250);

// Cuts and accepted muon histograms for no absorber reference data
  mccuts_accept = new TH1D("mccuts_accept", 
			   ";Selection Criteria; Surviving Events", 14, 0, 14);
  mcpass_Tof1 = new TH1D("refpass_TOF1",";Passed Cut; Events",3, 0, 3);
  mcpass_Tof0 = new TH1D("refpass_TOF0",";Passed Cut; Events",3, 0, 3);
  mcpass_TKU = new TH1D("refpass_TKU",";Passed Cut; Events",3, 0, 3);
  mcpass_TKUrad = new TH1D("refpass_TKUrad",";Passed Cut; Events",100, 0, 200);
  mcpass_Diff = new TH1D("refpass_Diff",";Passed Cut; Events",100, 0, 150);
  mcpass_Fid = new TH1D("refpass_Fid",";Passed Cut; Events",100, 0, 200);
  mcpass_Fid_wide = new TH1D("refpass_Fid_wide",";Passed Cut; Events",100, 0, 2000);
  mcpass_Fid_old = new TH1D("refpass_Fid_old",";Passed Cut; Events",100, 0, 200);
  mcpass_ExtTKU0 = new TH1D("refpass_ExtTKU0",";Passed Cut; Events",3, 0, 3);
  mcpass_TOF01 = new TH1D("refpass_TOF01",";Passed Cut; Events",100, 23, 35);
  mcpass_TOF01Ex = new TH1D("refpass_TOF01Ex",";Passed Cut; Events",100, 23, 35);
  mcpass_mom = new TH1D("refpass_mom",";Passed Cut; Events",100, 100, 300);
  mcpass_chi2= new TH1D("refpass_chi2",";Chi2/dof; Events",50, 0, 10);

//Total Momentum Histogram prior to each cut
  HPt_Tof1 = new TH1D("01_HPt_tof1","; Total Momentum; Number", 100, 100, 300);
  HPt_Tof0 = new TH1D("02_HPt_Tof0","; Total Momentum; Number", 100, 100, 300);
  HPt_TKU = new TH1D("03_HPt_TKU","; Total Momentum; Number", 100, 100, 300);
  HPt_TKUrad = new TH1D("04_HPt_TKUrad","; Total Momentum; Number", 100, 100, 300);
  HPt_Diff = new TH1D("05_HPt_Diff","; Total Momentum; Number", 100, 100, 300);
  HPt_Fid = new TH1D("06_HPt_Fid","; Total Momentum; Number", 100, 100, 300);
  HPt_ExtTKU0 = new TH1D("07_HPt_ExtTKU0","; Total Momentum; Number", 100, 100, 300);
  HPt_TOF01 = new TH1D("08_HPt_TOF01","; Total Momentum; Number", 100, 100, 300);
  HPt_TOF01Ex = new TH1D("09_HPt_TOF01Ex","; Total Momentum; Number", 100, 100, 300);
  HPt_mom = new TH1D("10_HPt_mom","; Total Momentum; Number", 100, 100, 300);
  HPt_Ch12 = new TH1D("11_HPt_Ch12","; Total Momentum; Number", 100, 100, 300);


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
    new TH1D("theta_true_scat_graph","Scattering Angle between Momentum Vectors;#theta_{Scatt}; Events per mrad",
	     _histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
  
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
  delete chain;
  delete refchain;
  delete mcchain;
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
  
  TText t1 = TText(0.12,0.785,"MICE Internal");
  TText t3 = TText(0.12,0.75,"ISIS cycle");
  TText t2 = TText(0.12,0.715,"LiH, MAUS v3.1.0");
  t1.SetNDC(1);
  t1.SetTextSize(0.04);
  t1.SetTextFont(42);
  t2.SetNDC(1);
  t2.SetTextSize(0.03);
  t2.SetTextFont(42);
  t3.SetNDC(1);
  t3.SetTextSize(0.04);
  t3.SetTextFont(42);
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
  Bfield_DStrackerX->Write();
  Bfield_DStrackerY->Write();
  Bfield_DStrackerZ->Write();
  cuts_accept->Write();
  mccuts_accept->Write();
  trackno->Write();
  tofhitno->Write();
  diffradius->Write();
  projradius->Write();
//  TCanvas *c1 = new TCanvas();
  //calc_mom->GetXaxis()->SetRangeUser(120,280);
  //calc_mom->Draw();
  cor_mom->SetLineColor(2);
  cor_mom->Draw();
  semom = cor_mom->GetMeanError();
  rmsmom = cor_mom->GetRMS();
  errrmsmom = cor_mom->GetRMSError();
  cor_mom->Fit("gaus");
  c1->SaveAs("calc_mom.pdf");
  c1->SaveAs("calc_mom.root");
  cor_mom->Write();
  c1->Clear();
  calc_mom->GetXaxis()->SetRangeUser(120,280);
  cor_mom->Draw();
  c1->SaveAs("cor_mom.pdf");
  c1->Clear();
  TOFcom->GetXaxis()->SetTitle("pz analytic formula");
  TOFcom->GetYaxis()->SetTitle("pz Bethe-Bloch");
  TOFcom->Draw("colz");
  TLine *line2 = new TLine(100,100,350,350);
  line2->SetLineColor(kRed);
  line2->Draw();
  c1->SaveAs("TOFcom.pdf");
  c1->Clear();
  TOF01vsMCTruth->GetXaxis()->SetTitle("pz TOF01");
  TOF01vsMCTruth->GetYaxis()->SetTitle("MCTruth");
  TOF01vsMCTruth->Draw("colz");
  TLine *line3 = new TLine(100,100,350,350);
  line3->SetLineColor(kRed);
  line3->Draw();
  c1->SaveAs("TOF01vsMCTruth.pdf");
  c1->Clear();
  TOF12vsMCTruth->GetXaxis()->SetTitle("pz TOF12");
  TOF12vsMCTruth->GetYaxis()->SetTitle("MCTruth");
  TOF12vsMCTruth->Draw("colz");
  TLine *line4 = new TLine(100,100,350,350);
  line4->SetLineColor(kRed);
  line4->Draw();
  c1->SaveAs("TOF12vsMCTruth.pdf");
  c1->Clear();
  energyloss->GetXaxis()->SetTitle("position (mm)");
  energyloss->GetYaxis()->SetTitle("MCTruth energy (MeV)");
  energyloss->Draw("colz");
  c1->SaveAs("energyloss.pdf");
  c1->Clear();
  energylossproj=energyloss->ProjectionX();
  energylossproj->SetTitle("position (mm)");
  energylossproj->Draw();
  c1->SaveAs("energylossproj.pdf");
  c1->Clear();
  energylosspro=energyloss->ProfileX();
  //energylosspro->SetTitle("position (mm)");
  energylosspro->Draw();
  c1->SaveAs("energylosspro.pdf");
  c1->Clear();
  gStyle->SetOptStat(0);
  TOF12Paul6thforvsMCTruth->SetTitle("Corrected P downstream vs MC Truth");
  TOF12Paul6thforvsMCTruth->GetXaxis()->SetTitle("pz TOF12");
  TOF12Paul6thforvsMCTruth->GetYaxis()->SetTitle("pz MCTruth");
  TOF12Paul6thforvsMCTruth->Draw("colz");
  t1.Draw();
  t3.Draw();
  t2.Draw();
  TLine *line5 = new TLine(100,100,350,350);
  line5->SetLineColor(kRed);
  line5->Draw();
  c1->SaveAs("TOF12Paul6thforvsMCTruth.pdf");
  c1->Clear();
  TOF12longPaul6thforvsMCTruth->GetXaxis()->SetTitle("pz TOF12 long Paul 6th for");
  TOF12longPaul6thforvsMCTruth->GetYaxis()->SetTitle("pz MCTruth");
  TOF12longPaul6thforvsMCTruth->Draw("colz");
  TLine *line25 = new TLine(100,100,350,350);
  line25->SetLineColor(kRed);
  line25->Draw();
  c1->SaveAs("TOF12longPaul6thforvsMCTruth.pdf");
  c1->Clear();
  TOF12Paul2ndforvsMCTruth->GetXaxis()->SetTitle("pz TOF12 Paul 2nd for");
  TOF12Paul2ndforvsMCTruth->GetYaxis()->SetTitle("pz MCTruth");
  TOF12Paul2ndforvsMCTruth->Draw("colz");
  TLine *line15 = new TLine(100,100,350,350);
  line15->SetLineColor(kRed);
  line15->Draw();
  c1->SaveAs("TOF12Paul2ndforvsMCTruth.pdf");
  c1->Clear();
  TOF01Paul6thforvsMCTruth->GetXaxis()->SetTitle("pz TOF01 Paul 6th for");
  TOF01Paul6thforvsMCTruth->GetYaxis()->SetTitle("pz MCTruth");
  TOF01Paul6thforvsMCTruth->Draw("colz");
  TLine *line16 = new TLine(100,100,350,350);
  line16->SetLineColor(kRed);
  line16->Draw();
  c1->SaveAs("TOF01Paul6thforvsMCTruth.pdf");
  c1->Clear();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit();
  //TOF01shortPaul6thforvsMCTruth->Fit("pol1","","",170,270);
  //TF1 *f3 = TOF01shortPaul6thforvsMCTruth->GetFunction("pol1");
  TOF01shortPaul6thforvsMCTruth->SetTitle("Corrected P upstream vs MC Truth");
  TOF01shortPaul6thforvsMCTruth->GetXaxis()->SetTitle("pz TOF01");
  TOF01shortPaul6thforvsMCTruth->GetYaxis()->SetTitle("pz MCTruth");
  TOF01shortPaul6thforvsMCTruth->SetFillColor(kViolet+1);
  TOF01shortPaul6thforvsMCTruth->Draw("colz");
  c1->GetFrame()->SetFillColor(4);
  t1.Draw();
  t3.Draw();
  t2.Draw();
  TLine *line46 = new TLine(100,100,350,350);
  line46->SetLineColor(kRed);
  line46->Draw();
  c1->SaveAs("TOF01shortPaul6thforvsMCTruth.pdf");
  c1->Clear();
  TOF01fordownvsMCTruth->GetXaxis()->SetTitle("pz John's formula TOF01");
  TOF01fordownvsMCTruth->GetYaxis()->SetTitle("pz MCTruth");
  TOF01fordownvsMCTruth->Draw("colz");
  TLine *line6 = new TLine(100,100,350,350);
  line6->SetLineColor(kRed);
  line6->Draw();
  c1->SaveAs("TOF01fordownvsMCTruth.pdf");
  c1->Clear();
  TOF12fordownvsMCTruth->GetXaxis()->SetTitle("pz John's formula TOF12");
  TOF12fordownvsMCTruth->GetYaxis()->SetTitle("pz MCTruth");
  TOF12fordownvsMCTruth->Draw("colz");
  TLine *line26 = new TLine(100,100,350,350);
  line26->SetLineColor(kRed);
  line26->Draw();
  c1->SaveAs("TOF12fordownvsMCTruth.pdf");
  c1->Clear();
  TOF12cobbvsMCTruth->GetXaxis()->SetTitle("pz cobb's formula TOF12");
  TOF12cobbvsMCTruth->GetYaxis()->SetTitle("pz MCTruth");
  TOF12cobbvsMCTruth->Draw("colz");
  TLine *line36 = new TLine(100,100,350,350);
  line36->SetLineColor(kRed);
  line36->Draw();
  c1->SaveAs("TOF12cobbvsMCTruth.pdf");
  c1->Clear();
  TOF01forupvsTOF01fordown->GetXaxis()->SetTitle("pz Ryan's formula downstream TOF12");
  TOF01forupvsTOF01fordown->GetYaxis()->SetTitle("pz John's formula TOF01");
  TOF01forupvsTOF01fordown->Draw("colz");
  TLine *line7 = new TLine(100,100,350,350);
  line7->SetLineColor(kRed);
  line7->Draw();
  c1->SaveAs("TOF01forupvsTOF01fordown.pdf");
  c1->Clear();
  residual01p6->GetXaxis()->SetTitle("residual01p6");
  residual01p6->GetYaxis()->SetTitle("No. of events");
  residual01p6->Draw();
  c1->SaveAs("residual01p6.pdf");
  c1->Clear();
  residual01j->GetXaxis()->SetTitle("residual01j");
  residual01j->GetYaxis()->SetTitle("No. of events");
  residual01j->Draw();
  c1->SaveAs("residual01j.pdf");
  TText t11 = TText(0.63,0.785,"MICE Preliminary");
  TText t13 = TText(0.63,0.75,"ISIS cycle 2015/04");
  TText t12 = TText(0.63,0.715,"LiH, MAUS v3.1.2");
  t11.SetNDC(1);
  t11.SetTextSize(0.04);
  t11.SetTextFont(42);
  t12.SetNDC(1);
  t12.SetTextSize(0.03);
  t12.SetTextFont(42);
  t13.SetNDC(1);
  t13.SetTextSize(0.04);
  t13.SetTextFont(42);
  c1->Clear();
  residual12p2->GetXaxis()->SetTitle("residual12p2");
  residual12p2->GetYaxis()->SetTitle("No. of events");
  residual12p2->Draw();
  c1->SaveAs("residual12p2.pdf");
  c1->Clear();
  //gStyle->SetOptStat(01);
  residual12p6->SetTitle("");
  residual12p6->GetXaxis()->SetTitle("Momentum residual downstream");
  residual12p6->GetYaxis()->SetTitle("No. of events");
  residual12p6->GetYaxis()->SetTitleOffset(1.5);
  residual12p6->Draw();
  residual12p6->GetYaxis()->SetRangeUser(0,50000);
  t11.Draw();
  t13.Draw();  
  t12.Draw();
  residual12->Draw("SAMES");
  residual12->SetLineColor(kRed);
  gPad->Update();
  TLegend * legend2 = new TLegend(0.1,0.7,0.34,0.9);
  legend2->AddEntry(residual12p6,"Corrected momentum","l");
  legend2->AddEntry(residual12,"Naive momentum","l");
  legend2->Draw();
  /*
  TPaveStats *st = (TPaveStats*)residual12->FindObject("stats");
  st->SetY1NDC(0.4); 
  st->SetY2NDC(0.6); 
  */
  c1->SaveAs("residual12p6.pdf");
  c1->Clear();
  residualcobb->GetXaxis()->SetTitle("residualcobb");
  residualcobb->GetYaxis()->SetTitle("No. of events");
  residualcobb->Draw();
  c1->SaveAs("residualcobb.pdf");
  c1->Clear();
  residual01short->SetTitle("");
  residual01short->GetXaxis()->SetTitle("Momentum residual upstream");
  residual01short->GetYaxis()->SetTitle("No. of events");
  residual01short->GetYaxis()->SetTitleOffset(1.5);
  residual01short->Draw();
  //residual01short->GetYaxis()->SetRangeUser(0,35000);
  t11.Draw();
  t13.Draw();
  t12.Draw();
  residual01->Draw("SAMES");
  residual01->SetLineColor(kRed);
  TLegend * legend = new TLegend(0.1,0.7,0.34,0.9);
  legend->AddEntry(residual01short,"Corrected momentum","l");
  legend->AddEntry(residual01,"Naive Momentum","l");
  legend->Draw();
  gPad->Update();
  /*
  TPaveStats *st2 = (TPaveStats*)residual01->FindObject("stats");
  st2->SetY1NDC(0.4); 
  st2->SetY2NDC(0.6); 
  */
  c1->SaveAs("residual01short.pdf");
  c1->Clear();
  TOF0Energy->GetXaxis()->SetTitle("TOF0Energy");
  TOF0Energy->GetYaxis()->SetTitle("No. of events");
  TOF0Energy->Draw();
  c1->SaveAs("TOF0Energy.pdf");
  c1->Clear();
  TOF1Energy->GetXaxis()->SetTitle("TOF1Energy");
  TOF1Energy->GetYaxis()->SetTitle("No. of events");
  TOF1Energy->Draw();
  c1->SaveAs("TOF1Energy.pdf");
  c1->Clear();
  TOF2Energy->GetXaxis()->SetTitle("TOF2Energy");
  TOF2Energy->GetYaxis()->SetTitle("No. of events");
  TOF2Energy->Draw();
  c1->SaveAs("TOF2Energy.pdf");
  c1->Clear();
  difEloss->GetXaxis()->SetTitle("difEloss");
  difEloss->GetYaxis()->SetTitle("No. of events");
  difEloss->Draw();
  c1->SaveAs("difEloss.pdf");
  c1->Clear();
  gStyle->SetOptStat(1);
  pathlengthabs->GetXaxis()->SetTitle("Path Length in Absorber (mm)");
  pathlengthabs->GetYaxis()->SetTitle("No. of events");
  pathlengthabs->Draw();
  c1->SaveAs("pathlengthabs.pdf");
  c1->SaveAs("pathlengthabs.root");
  c1->Clear();
  gStyle->SetOptStat(1);
  gStyle->SetOptStat("nemR");
  theta_true_x_graph->GetXaxis()->SetTitle("#Delta#theta_{x} (mrad)");
  theta_true_x_graph->GetYaxis()->SetTitle("No. of events");
  //theta_true_x_graph->Sumw2();
  theta_true_x_graph->SetFillColor(kOrange-2);
  theta_true_x_graph->DrawNormalized();
  theta_true_x_graph->Write();
  TH1D* thetaX_reco = (TH1D*)outfile->Get("thetaX_recoGEANT");
  thetaX_reco->SetMarkerStyle(20);
  //thetaX_reco->Sumw2();
  thetaX_reco->DrawNormalized("e1 p SAMES");
  TLegend * legend20 = new TLegend(0.1,0.7,0.34,0.9);
  legend20->AddEntry(theta_true_x_graph,"MC","f");
  legend20->AddEntry(thetaX_reco,"Reco sim","p");
  legend20->Draw();
  gPad->Update();
  //TPaveStats *st = (TPaveStats*)thetaX_reco->FindObject("stats");
  //st->SetY1NDC(0.4); 
  //st->SetY2NDC(0.6); 
  double chi2 = theta_true_x_graph->Chi2Test(thetaX_reco,"WUP");
  TText t10 = TText(0.12,0.585,Form("Chi^2 = %f",chi2));
  t10.SetNDC(1);
  t10.SetTextSize(0.03);
  t10.SetTextFont(42);
  t10.Draw();
  double Kol = theta_true_x_graph->KolmogorovTest(thetaX_reco);
  std::cout << "Kolmogorov " << Kol << std::endl;
  TText t20 = TText(0.12,0.55,Form("Kolmogorov = %f",Kol));
  t20.SetNDC(1);
  t20.SetTextSize(0.03);
  t20.SetTextFont(42);
  t20.Draw();
  TText t100 = TText(0.652,0.585,"MICE Preliminary");
  TText t300 = TText(0.652,0.55,"ISIS cycle 2015/04");
  TText t200 = TText(0.652,0.515,"LiH, MAUS v3.1.2");
  t100.SetNDC(1);
  t100.SetTextSize(0.04);
  t100.SetTextFont(42);
  t200.SetNDC(1);
  t200.SetTextSize(0.03);
  t200.SetTextFont(42);
  t300.SetNDC(1);
  t300.SetTextSize(0.04);
  t300.SetTextFont(42);
  t100.Draw();
  t300.Draw();
  t200.Draw();
  c1->SaveAs("theta_true_x_graph.pdf");
  c1->SaveAs("theta_true_x_graph.root");
  c1->Clear();
  gStyle->SetOptStat(1);
  gStyle->SetOptStat("nemR");
  theta_true_y_graph->GetXaxis()->SetTitle("#Delta#theta_y (mrad)");
  theta_true_y_graph->GetYaxis()->SetTitle("No. of events");
  theta_true_y_graph->DrawNormalized();
  theta_true_y_graph->Write();
  theta_true_y_graph->SetFillColor(kOrange-2);
  //theta_true_y_graph->Sumw2();
  TH1D* thetaY_reco = (TH1D*)outfile->Get("thetaY_recoGEANT");
  thetaY_reco->SetMarkerStyle(20);
  //thetaY_reco->Sumw2();
  thetaY_reco->DrawNormalized("e1 p SAMES");
  TLegend * legend21 = new TLegend(0.1,0.7,0.34,0.9);
  legend21->AddEntry(theta_true_y_graph,"MC","f");
  legend21->AddEntry(thetaY_reco,"Reco sim","p");
  legend21->Draw();
  gPad->Update();
  //TPaveStats *st1 = (TPaveStats*)thetaY_reco->FindObject("stats");
  //st1->SetY1NDC(0.4); 
  //st1->SetY2NDC(0.6); 
  chi2 = theta_true_y_graph->Chi2Test(thetaY_reco,"WUP");
  TText t101 = TText(0.12,0.585,Form("Chi^2 = %.17g",chi2));
  t101.SetNDC(1);
  t101.SetTextSize(0.03);
  t101.SetTextFont(42);
  t101.Draw();
  Kol = theta_true_y_graph->KolmogorovTest(thetaY_reco);
  std::cout << "Kolmogorov " << Kol << std::endl;
  TText t21 = TText(0.12,0.55,Form("Kolmogorov = %.17g",Kol));
  t21.SetNDC(1);
  t21.SetTextSize(0.03);
  t21.SetTextFont(42);
  t21.Draw();
  c1->SaveAs("theta_true_y_graph.pdf");
  c1->SaveAs("theta_true_y_graph.root");

  pass_Tof1->Write();
  pass_Tof0->Write();
  pass_TKU->Write();
  pass_TKUrad->Write();
  pass_Diff->Write();
  pass_Fid->Write();
  pass_Fid_wide->Write();
  pass_Fid_old->Write();
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
  mcpass_Fid_wide->Write();
  mcpass_Fid_old->Write();
  mcpass_ExtTKU0->Write();
  mcpass_TOF01->Write();
  mcpass_TOF01Ex->Write();
  mcpass_mom->Write();
  mcpass_chi2->Write();
  HPt_Tof1->Write();
  HPt_Tof0->Write();
  HPt_TKU->Write();
  HPt_TKUrad->Write();
  HPt_Diff->Write();
  HPt_Fid->Write();
  HPt_ExtTKU0->Write();
  HPt_TOF01->Write();
  HPt_TOF01Ex->Write();
  HPt_mom->Write();
  HPt_Ch12->Write();

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

void MCSAnalysis::print(){
	std::cout << "standard error on mean momentum " << semom << std::endl;
	std::cout << "rms of momentum " << rmsmom << std::endl;
	std::cout << "Error on rms of momentum " << errrmsmom << std::endl;
	bwmom = (1.04e3/pow(TOF_upper_limit - 23.2,2)) - (1.04e3/pow(TOF_lower_limit - 23.2,2));
	std::cout << "width of bin mom " << bwmom << std::endl;
	std::cout << "bin width sqrt12 " << bwmom/sqrt(12) << std::endl;
}

void MCSAnalysis::Execute(int mode=0){
  
  int verbose_level = 1;
  // MAUS::GlobalsManager::InitialiseGlobals(JsonWrapper::JsonToString(SetupConfig(verbose_level)));
  std::cout<<"Calling Globals Initialisation"<<std::endl;
  std::cout<<"ResX: "<<_sys["resX"]<<std::endl;
  MAUS::GlobalsManager::InitialiseGlobals(
             JsonWrapper::JsonToString(SetupConfig(verbose_level, geometryfile)));
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
  std::cerr<<"End of DS tracker data file: "<<EndofDStrackerdataname<<std::endl;
  std::cerr<<"End of DS tracker ref file: "<<EndofDStrackerrefname<<std::endl;
  std::cerr<<"Start of DS tracker data file: "<<StartofDStrackerdataname<<std::endl;
  std::cerr<<"Start of DS tracker ref file: "<<StartofDStrackerrefname<<std::endl;
  std::cerr<<"US track to TOF0 data file: "<<Collect_USTOF0dataname<<std::endl;
  if (mode == 0){
    retrieve(0);
    generateMCSResponse();
    // DoUnfolding();
    DoDeconvolution(modelname2.c_str(),1);
  }
  else if (mode == 100){
    std::cerr<<"Running Pre-Data Analysis Phase"<<std::endl;
    dataSelection(mode);
    referenceSelection();
  }
  else if (mode == 101){
    std::cerr<<"Running Pre-Data Analysis Phase - Absorber Only"<<std::endl;
    dataSelection(mode);
  }
  else if (mode == 102){
    std::cerr<<"Running Pre-Data Analysis Phase - Empty Only"<<std::endl;
    referenceSelection();
  }
  else if (mode == 200){
     retrieve(0);
     retrieve(1);
//     Write();
  }
  else if (mode == 1){
    std::cerr<<"Running Main Analysis Phase"<<std::endl;
    retrieve(0);
    retrieve(1);   
    DoUnfolding();
    DoFFTDeconvolution();
    ConvolveWithInputDistribution(modelname1.c_str());
    DoDeconvolution(modelname1.c_str(), 1);
    ConvolveWithInputDistribution(modelname2.c_str());
    DoDeconvolution(modelname2.c_str(), 1);
    //FitGaussian(outfilename.c_str());
    //CalculateChi2(outfilename.c_str(), modelname2.c_str());
    Write();
  }
  else if (mode == 6){
    retrieve(1);
    Write();   
  }
  else if (mode == -1){
    retrieve(0);
    retrieve(1);   
    DoUnfolding();
    DoFFTDeconvolution();
    ConvolveWithInputDistribution(modelname1.c_str());
    DoDeconvolution(modelname1.c_str(), 1);
    ConvolveWithInputDistribution(modelname2.c_str());
    DoDeconvolution(modelname2.c_str(), 1);
    ConvolveWithInputDistribution(modelname3.c_str());
    DoDeconvolution(modelname3.c_str(), 1);
    TruthData(mode);
    //FitGaussian(outfilename.c_str());
    //CalculateChi2(outfilename.c_str(), modelname2.c_str());
    Write();
    print();
  }
  else if (mode == -2){
    retrieve(0);
    retrieve(1);   
    DoUnfolding();
    DoFFTDeconvolution();
    ConvolveWithInputDistribution(modelname1.c_str());
    DoDeconvolution(modelname1.c_str(), 1);
    ConvolveWithInputDistribution(modelname2.c_str());
    DoDeconvolution(modelname2.c_str(), 1);
    ConvolveWithInputDistribution(modelname3.c_str());
    DoDeconvolution(modelname3.c_str(), 1);
    TruthData(mode);
    //FitGaussian(outfilename.c_str());
    //CalculateChi2(outfilename.c_str(), modelname2.c_str());
    Write();
    print();
  }
  else if (mode == 2){
    retrieve(0);
    TruthData(mode);
    Write();
  } 
  else if (mode == 5){
    retrieve(0);
    retrieve(1);   
    DoUnfolding();
    DoFFTDeconvolution();
    TruthData(mode);
    ConvolveWithVirtualInputDistribution(modelname1.c_str());
    DoDeconvolution(modelname1.c_str(), 1);
    ConvolveWithVirtualInputDistribution(modelname2.c_str());
    DoDeconvolution(modelname2.c_str(), 1);
    ConvolveWithVirtualInputDistribution(modelname3.c_str());
    DoDeconvolution(modelname3.c_str(), 1);
    Write();
  }
  else {
    std::cout<<"Unknown operation mode"<<std::endl;
  } 

}

/////////////////////////////////////////////////////////////////////
// Function Retrieve
/////////////////////////////////////////////////////////////////////

void MCSAnalysis::retrieve(int retrieve_mode){
// retrieve_mode = 0 for absorber present and 1 for no absorber present.
  Collection USretrieved, DSretrieved, EndofDS, StartofDS, USPostFid, DSPostFid, USPostMom;
  Collection DSPostMom, DSsetUncorrected, USPostExtTKU0, DSPostExtTKU0, USPostTOF01Ex;
  Collection DSPostTOF01Ex, USPostDiff, DSPostDiff, USPostTKUrad, DSPostTKUrad, USPostChi2;
  Collection DSPostChi2, USPostTOF01, DSPostTOF01, USPostTOF0, DSPostTOF0, USAllTOF, DSAllTOF;
  Collection Collect_USTOF0, Collect_USTOF1, Collect_StartofUS, Collect_EndofUS, Collect_StartofUST;
  Collection Collect_EndofUST, Collect_DiffuserStart, Collect_DiffuserEnd;
  std::cerr<<"Retrieving Data from either data or ref"<<std::endl;
  if(retrieve_mode==0){
     std::cerr<<"Retrieving data from USretrieved: "<<USDataname<<std::endl;
     USretrieved.load_collection(USDataname);
     std::cerr<<"Retrieving data from DSretrieved: "<<DSDataname<<std::endl;
     DSretrieved.load_collection(DSDataname);
     std::cerr<<"Retrieving data from EndofDS: "<<EndofDStrackerdataname<<std::endl;
     EndofDS.load_collection(EndofDStrackerdataname);
     std::cerr<<"Retrieving data from StartofDS: "<<StartofDStrackerdataname<<std::endl;
     StartofDS.load_collection(StartofDStrackerdataname);
     std::cerr<<"Retrieving data from Collect_USTOF0: "<<Collect_USTOF0dataname<<std::endl;
     Collect_USTOF0.load_collection(Collect_USTOF0dataname);    
     std::cerr<<"Retrieving data from Collect_USTOF1: "<<Collect_USTOF1dataname<<std::endl;
     Collect_USTOF1.load_collection(Collect_USTOF1dataname);    
     std::cerr<<"Retrieving data from Collect_StartofUST: "<<Collect_StartofUSTdataname<<std::endl;
     Collect_StartofUST.load_collection(Collect_StartofUSTdataname);    
     std::cerr<<"Retrieving data from Collect_EndofUST: "<<Collect_EndofUSTdataname<<std::endl;
     Collect_EndofUST.load_collection(Collect_EndofUSTdataname);    
     std::cerr<<"Retrieving data from Collect_DiffuserStart: "<<Collect_DiffuserStartdataname<<std::endl;
     Collect_DiffuserStart.load_collection(Collect_DiffuserStartdataname);    
     std::cerr<<"Retrieving data from Collect_DiffuserEnd: "<<Collect_DiffuserEnddataname<<std::endl;
     Collect_DiffuserEnd.load_collection(Collect_DiffuserEnddataname);    
  }else{
     std::cerr<<"Retrieving ref from USretrieved: "<<USRefname<<std::endl;
     USretrieved.load_collection(USRefname);
     std::cerr<<"Retrieving ref from DSretrieved: "<<DSRefname<<std::endl;
     DSretrieved.load_collection(DSRefname);
     std::cerr<<"Retrieving ref from EndofDS: "<<EndofDStrackerrefname<<std::endl;
     EndofDS.load_collection(EndofDStrackerrefname);
     std::cerr<<"Retrieving ref from StartofDS: "<<StartofDStrackerrefname<<std::endl;
     StartofDS.load_collection(StartofDStrackerrefname);
     std::cerr<<"Retrieving ref from Collect_USTOF0: "<<Collect_USTOF0refname<<std::endl;
     Collect_USTOF0.load_collection(Collect_USTOF0refname);
     std::cerr<<"Retrieving ref from Collect_USTOF1: "<<Collect_USTOF1refname<<std::endl;
     Collect_USTOF1.load_collection(Collect_USTOF1refname);
     std::cerr<<"Retrieving data from Collect_StartofUST: "<<Collect_StartofUSTrefname<<std::endl;
     Collect_StartofUST.load_collection(Collect_StartofUSTrefname);    
     std::cerr<<"Retrieving data from Collect_EndofUST: "<<Collect_EndofUSTrefname<<std::endl;
     Collect_EndofUST.load_collection(Collect_EndofUSTrefname);    
     std::cerr<<"Retrieving data from Collect_DiffuserStart: "<<Collect_DiffuserStartrefname<<std::endl;
     Collect_DiffuserStart.load_collection(Collect_DiffuserStartrefname);    
     std::cerr<<"Retrieving data from Collect_DiffuserEnd: "<<Collect_DiffuserEndrefname<<std::endl;
     Collect_DiffuserEnd.load_collection(Collect_DiffuserEndrefname);    
  }
  std::cerr<<"Entries in USretrieved: "<< USretrieved.N()<<std::endl;
  std::cerr<<"Entries in DSretrieved: "<< DSretrieved.N()<<std::endl;
  std::cerr<<"Entries in EndofDS: "<< EndofDS.N()<<std::endl;
  std::cerr<<"Entries in startofDS: "<< StartofDS.N()<<std::endl;
  std::cerr<<"Entries in Collect_USTOF0: "<< Collect_USTOF0.N()<<std::endl; 
  std::cerr<<"Entries in Collect_USTOF1: "<< Collect_USTOF1.N()<<std::endl; 
  std::cerr<<"Entries in Collect_StartofUST: "<< Collect_StartofUST.N()<<std::endl; 
  std::cerr<<"Entries in Collect_EndofUST: "<< Collect_EndofUST.N()<<std::endl; 
  std::cerr<<"Entries in Collect_DiffuserStart: "<< Collect_DiffuserStart.N()<<std::endl; 
  std::cerr<<"Entries in Collect_DiffuserEnd: "<< Collect_DiffuserEnd.N()<<std::endl; 
  outfile->cd();
  BTField* field = dynamic_cast<BTField*>(MAUS::Globals::GetMCFieldConstructor());
  double pi = 3.14159265358979323846;
  double c_light=299.792458;
  double x = 0; double y = 0; double t = 0;

// Open files to load TOF data to.
   std::ifstream TOF0file;
   std::string TOF0line;
   std::vector<std::string> TOF0results;
   std::ifstream TOF1file;
   std::string TOF1line;
   std::vector<std::string> TOF1results;
   std::ifstream TOF2file;
   std::string TOF2line;
   std::vector<std::string> TOF2results;
   std::ifstream Chi2file;
   std::string Chi2line;
   std::vector<std::string> Chi2results;
   std::ifstream doffile;
   std::string dofline;
   std::vector<std::string> dofresults;

   if(retrieve_mode==0){   
      TOF0file.open(TOF0dataname);
      TOF1file.open(TOF1dataname);
      TOF2file.open(TOF2dataname);
      Chi2file.open(Chi2dataname);
      doffile.open(dofdataname);
   }else{
      TOF0file.open(TOF0refname);
      TOF1file.open(TOF1refname);
      TOF2file.open(TOF2refname);
      Chi2file.open(Chi2refname);
      doffile.open(dofrefname);
   }

//Plot axial magnetic field
  for(int i=0;i<1000;i++){
     double z=(_sys["dataDS11"]+((i+0.5)*(_sys["dataDS53"]-_sys["dataDS11"])/1000));
     std::vector<double> magfield=field->GetFieldValue2(x,y,z,t);
     Bfield_DStrackerX->Fill(z,magfield[0]);
     Bfield_DStrackerY->Fill(z,magfield[1]);
     Bfield_DStrackerZ->Fill(z,magfield[2]);
  }
  for(int i=0;i<USretrieved.N();i++){
     double AbsCenterpx=USretrieved.E(i).px;
     double AbsCenterpy=USretrieved.E(i).py;
     double AbsCenterpz=USretrieved.E(i).pz;
     double Ptot =sqrt(AbsCenterpx*AbsCenterpx+AbsCenterpy*AbsCenterpy+AbsCenterpz*AbsCenterpz);
//     if(i>100) break;

     getline(TOF0file,TOF0line);
     ReturnWords(TOF0line, TOF0results);
     getline(TOF1file,TOF1line);
     ReturnWords(TOF1line, TOF1results);
     getline(TOF2file,TOF2line);
     ReturnWords(TOF2line, TOF2results);
     getline(Chi2file,Chi2line);
     ReturnWords(Chi2line, Chi2results);
     getline(doffile,dofline);
     ReturnWords(dofline, dofresults);


//----------------------------------------------------------------------
// Require 1 TOF0 Space Point
//----------------------------------------------------------------------

     double rawTOF0HitTime=-100.;
     bool CutPassRequire1TOF0SP=false;
     int NumberTOF0SP=std::stoi(TOF0results[0]);
     if ( NumberTOF0SP == 1 ){
        rawTOF0HitTime = std::stod(TOF0results[1]);
        CutPassRequire1TOF0SP=true;
        if(retrieve_mode==0){pass_Tof1->Fill("Pass",1);}
        else{mcpass_Tof1->Fill("Pass",1);}
     }else{
        if(retrieve_mode==0){pass_Tof0->Fill("Fail",1);}
        else{mcpass_Tof0->Fill("Fail",1);}
     }
     if (CutPassRequire1TOF0SP){
        if(retrieve_mode==0){
           cuts_accept->Fill("1 TOF0",1);
           USPostTOF0.append_instance(USretrieved.E(i));
           DSPostTOF0.append_instance(DSretrieved.E(i));
        }else{
           mccuts_accept->Fill("1 TOF0",1);
           USPostTOF0.append_instance(USretrieved.E(i));
           DSPostTOF0.append_instance(DSretrieved.E(i));
        }  
     }
//    std::cerr<<"TOF0 Number of SP:"<<NumberTOF0SP<<" Pass or fail? "<<CutPassRequire1TOF0SP<<" TOF0 time: "<<rawTOF0HitTime<<std::endl;

//----------------------------------------------------------------------
// Require 1 TOF1 Space Point
//----------------------------------------------------------------------


     double rawTOF1HitTime=-100.;
     bool CutPassRequire1TOF1SP=false;
     int NumberTOF1SP=std::stoi(TOF1results[0]);
     if ( NumberTOF1SP == 1 ){
        rawTOF1HitTime = std::stod(TOF1results[1]);
        CutPassRequire1TOF1SP=true;
        if(retrieve_mode==0){pass_Tof1->Fill("Pass",1);}
        else{mcpass_Tof1->Fill("Pass",1);}
     }else{
        if(retrieve_mode==0){pass_Tof1->Fill("Fail",1);
        }else{mcpass_Tof1->Fill("Fail",1);}
     }
     if (CutPassRequire1TOF1SP && CutPassRequire1TOF0SP){
        if(retrieve_mode==0){
           cuts_accept->Fill("1 TOF1",1);
           USAllTOF.append_instance(USretrieved.E(i));
           DSAllTOF.append_instance(DSretrieved.E(i));
        }else{
           mccuts_accept->Fill("1 TOF1",1);
           USAllTOF.append_instance(USretrieved.E(i));
           DSAllTOF.append_instance(DSretrieved.E(i));
        }
     }
//     std::cerr<<"TOF1 Number of SP:"<<NumberTOF1SP<<" Pass or fail? "<<CutPassRequire1TOF1SP<<" TOF1 time: "<<rawTOF1HitTime<<std::endl;
//---------------------------------------------------------------------------------------
// TOF01 cut. Require TOF0 and TOF1 hit and TOF01 in bounds set by Upper and lower limits
//---------------------------------------------------------------------------------------

     double rawTOF2HitTime;
     if ( std::stoi(TOF2results[0]) >= 1 ){
        rawTOF2HitTime = std::stod(TOF2results[1]);
     }else{
        rawTOF2HitTime = rawTOF1HitTime + 100.0 * 8.22475 / 0.299792458; // ns.
     }

     tof10->Fill(rawTOF1HitTime - rawTOF0HitTime);
     tof21->Fill(rawTOF2HitTime - rawTOF1HitTime);
     double rawTOF01HitTime = rawTOF1HitTime - rawTOF0HitTime;
     bool CutPassTOF01=false;     
     if(retrieve_mode==0){
        pass_TOF01->Fill(rawTOF01HitTime);
     }else{
        mcpass_TOF01->Fill(rawTOF01HitTime);
     }
     if ( rawTOF01HitTime < TOF_lower_limit || rawTOF01HitTime > TOF_upper_limit){
        CutPassTOF01=false;
     }else{
        CutPassTOF01=true;
        tof10_sel->Fill(rawTOF1HitTime - rawTOF0HitTime);
        tof21_sel->Fill(rawTOF2HitTime - rawTOF1HitTime);
     }


     if (CutPassRequire1TOF1SP && CutPassRequire1TOF0SP && CutPassTOF01){
        if(retrieve_mode==0){
           cuts_accept->Fill("TOF01",1);
           USPostTOF01.append_instance(USretrieved.E(i));
           DSPostTOF01.append_instance(DSretrieved.E(i));
        }else{
           mccuts_accept->Fill("TOF01",1);
           USPostTOF01.append_instance(USretrieved.E(i));
           DSPostTOF01.append_instance(DSretrieved.E(i));
        }
     }
//     std::cerr<<"TOF01 :"<<rawTOF01HitTime<<" Pass or fail? "<<CutPassTOF01<<std::endl;
//---------------------------------------------------------------------------------------
// CHi2/dof cut
//---------------------------------------------------------------------------------------

     double TKU_chiSquare_per_ndof = std::stod(Chi2results[0]) / (1.0*std::stod(dofresults[0]));
     pass_chi2->Fill(TKU_chiSquare_per_ndof);  
     bool CutPassChi2=false;     
     if (TKU_chiSquare_per_ndof<=10){
        CutPassChi2=true;
     }
     if (CutPassRequire1TOF1SP && CutPassRequire1TOF0SP && CutPassTOF01 && CutPassChi2){    
        USPostChi2.append_instance(USretrieved.E(i));
        DSPostChi2.append_instance(DSretrieved.E(i));
        if(retrieve_mode==0){
           cuts_accept->Fill("Chi2",1);
        }else{
           mccuts_accept->Fill("Chi2",1); 
        } 
     }
    

//------------------------------------------------------------------------------------------------------
// US Tracker Radial Cut
//------------------------------------------------------------------------------------------------------
     double USTendx=Collect_EndofUST.E(i).X;
     double USTendy=Collect_EndofUST.E(i).Y;
     double USTendz=Collect_EndofUST.E(i).Z;
     double USTendpx=Collect_EndofUST.E(i).px;
     double USTendpy=Collect_EndofUST.E(i).py;
     double USTendpz=Collect_EndofUST.E(i).pz;
     double USTstartx=Collect_StartofUST.E(i).X;
     double USTstarty=Collect_StartofUST.E(i).Y;
     double USTstartz=Collect_StartofUST.E(i).Z;
     double USTstartpx=Collect_StartofUST.E(i).px;
     double USTstartpy=Collect_StartofUST.E(i).py;
     std::vector<double> USTmagfield=field->GetFieldValue2(x,y,USTstartz,t);
     double USTbz=USTmagfield[2];
     double TKU_max_rad=MaxRadius(USTbz, USTstartx, USTstarty, USTstartz, USTstartpx, USTstartpy, USTendx, USTendy, USTendz, USTendpx, USTendpy, USTendpz);
     bool CutPassTKU_max_rad=false;         
     if (retrieve_mode==0){
        pass_TKUrad->Fill(TKU_max_rad);
     }else{
        mcpass_TKUrad->Fill(TKU_max_rad);
     }
     if (TKU_max_rad<=150.0){
        CutPassTKU_max_rad=true;
     }

     if (CutPassRequire1TOF1SP && CutPassRequire1TOF0SP && CutPassTOF01 && CutPassChi2 && CutPassTKU_max_rad){    
        USPostTKUrad.append_instance(USretrieved.E(i));
        DSPostTKUrad.append_instance(DSretrieved.E(i));
        if(retrieve_mode==0){
           cuts_accept->Fill("US tracker Radius",1);
        }else{
           mccuts_accept->Fill("US tracker Radius",1); 
        } 
     }

//------------------------------------------------------------------------------------------------------
// Diffuser Cut
//------------------------------------------------------------------------------------------------------
     double diffuserendx=Collect_DiffuserEnd.E(i).X;
     double diffuserendy=Collect_DiffuserEnd.E(i).Y;
//     double diffuserendz=Collect_DiffuserEnd.E(i).Z;
//     double diffuserendpx=Collect_DiffuserEnd.E(i).px;
//     double diffuserendpy=Collect_DiffuserEnd.E(i).py;
//     double diffuserendpz=Collect_DiffuserEnd.E(i).pz;
     double diffuserstartx=Collect_DiffuserStart.E(i).X;
     double diffuserstarty=Collect_DiffuserStart.E(i).Y;
//     double diffuserstartz=Collect_DiffuserStart.E(i).Z;
//     double diffuserstartpx=Collect_DiffuserStart.E(i).px;
//     double diffuserstartpy=Collect_DiffuserStart.E(i).py;
//     std::vector<double> diffusermagfield=field->GetFieldValue2(x,y,diffuserstartz,t);
//     double diffuserbz=diffusermagfield[2];
//     double diff_max_rad=MaxRadius(diffuserbz, diffuserstartx, diffuserstarty, diffuserstartz, diffuserstartpx, diffuserstartpy, diffuserendx, diffuserendy, diffuserendz, diffuserendpx, diffuserendpy, diffuserendpz); 
     double diff_start_rad=sqrt(diffuserstartx*diffuserstartx+diffuserstarty*diffuserstarty);
     double diff_end_rad=sqrt(diffuserendx*diffuserendx+diffuserendy*diffuserendy);
     double diff_max_rad=max(diff_start_rad,diff_end_rad); 
//     std::cerr<<"Diff Start Rad: "<<sqrt(diffuserstartx*diffuserstartx+diffuserstarty*diffuserstartx)<<" Diff End Rad: "<<sqrt(diffuserendx*diffuserendx+diffuserendy*diffuserendx)<<" Diff Max Rad: "<<diff_max_rad<<" B-Field : "<<diffuserbz<<std::endl;
     bool CutPassDiff_max_rad=false;         
     if (retrieve_mode==0){
        pass_Diff->Fill(diff_max_rad);
     }else{
        mcpass_Diff->Fill(diff_max_rad);
     }
     if (diff_max_rad <=90.0){
        CutPassDiff_max_rad=true;         
     }

     if (CutPassRequire1TOF1SP && CutPassRequire1TOF0SP && CutPassTOF01 && CutPassChi2 && CutPassTKU_max_rad && CutPassTKU_max_rad && CutPassDiff_max_rad){    
        USPostDiff.append_instance(USretrieved.E(i));
        DSPostDiff.append_instance(DSretrieved.E(i));
        if(retrieve_mode==0){
           cuts_accept->Fill("Diffuser",1);
        }else{
           mccuts_accept->Fill("Diffuser",1);
        } 
     }

//---------------------------------------------------------------------------------------
// TOF01 Extrapolated cut. 
//---------------------------------------------------------------------------------------

  double TOF01Ext = Collect_USTOF1.E(i).TOFX-Collect_USTOF0.E(i).TOFX;
  bool CutPassTOF01Ext=true;
  if (TOF01Ext < TOF_lower_limit || TOF01Ext > TOF_upper_limit) {
     CutPassTOF01Ext=false;
//    std::cerr<<"**Failed TOF01 Extrapolation Cut"<<std::endl;
//    if(cut_TOF01ext) continue;
  }
  if(retrieve_mode==0){
     HPt_ExtTKU0->Fill(Ptot);
  }
  if (CutPassRequire1TOF1SP && CutPassRequire1TOF0SP && CutPassTOF01 && CutPassChi2 && CutPassTKU_max_rad && CutPassTKU_max_rad && CutPassDiff_max_rad && CutPassTOF01Ext) {    
     USPostTOF01Ex.append_instance(USretrieved.E(i));
     DSPostTOF01Ex.append_instance(DSretrieved.E(i));
     if(retrieve_mode==0){  
        pass_TOF01Ex->Fill(TOF01Ext);
        cuts_accept->Fill("Extrapolated TOF01",1);
     }else{
        mcpass_TOF01Ex->Fill(TOF01Ext); 
        mccuts_accept->Fill("Extrapolated TOF01",1);
     }
  }

//---------------------------------------------------------------------------------------
// Successful extrapolation from  US tracker to TOF0
//---------------------------------------------------------------------------------------
     
     bool CutPassUSTOF0=true;
     if (Collect_USTOF0.E(i).isgood==false){
        CutPassUSTOF0=false;
        if (retrieve_mode==0){
           pass_ExtTKU0->Fill("Fail",1);           
//           if (cut_ext_TKU_TOF0) continue;
	}else{
	   mcpass_ExtTKU0->Fill("Fail",1);
//	   if (cut_ext_TKU_TOF0) continue;
	}
     }

     if (CutPassRequire1TOF1SP && CutPassRequire1TOF0SP && CutPassTOF01 && CutPassChi2 && CutPassTKU_max_rad && CutPassTKU_max_rad && CutPassDiff_max_rad && CutPassTOF01Ext && CutPassUSTOF0 ) {    
        USPostExtTKU0.append_instance(USretrieved.E(i));
        DSPostExtTKU0.append_instance(DSretrieved.E(i));
        if (retrieve_mode==0){
          cuts_accept->Fill("TKU_TOF0",1);
          pass_ExtTKU0->Fill("Pass",1);
        }else{
           mccuts_accept->Fill("TKU_TOF0",1);
	   mcpass_ExtTKU0->Fill("Pass",1);
        }
     }

//------------------------------------------------------------------------------------------------------
// Fiducial Cut
//------------------------------------------------------------------------------------------------------
     double fiducialendx=EndofDS.E(i).X;
     double fiducialendy=EndofDS.E(i).Y;
     double fiducialendz=EndofDS.E(i).Z;
     double fiducialendpx=EndofDS.E(i).px;
     double fiducialendpy=EndofDS.E(i).py;
     double fiducialendpz=EndofDS.E(i).pz;
     double fiducialstartx=StartofDS.E(i).X;
     double fiducialstarty=StartofDS.E(i).Y;
     double fiducialstartz=StartofDS.E(i).Z;
     double fiducialstartpx=StartofDS.E(i).px;
     double fiducialstartpy=StartofDS.E(i).py;
     std::vector<double> fiducialmagfield=field->GetFieldValue2(x,y,fiducialstartz,t);
     double fiducialbz=fiducialmagfield[2];
     double fiducial_radius=MaxRadius(fiducialbz, fiducialstartx, fiducialstarty, fiducialstartz, fiducialstartpx, fiducialstartpy, fiducialendx, fiducialendy, fiducialendz, fiducialendpx, fiducialendpy, fiducialendpz);    
//------------------------------------------------------------------------------------------------------        
     double fiducial_radius2 = sqrt(EndofDS.E(i).X*EndofDS.E(i).X + EndofDS.E(i).Y*EndofDS.E(i).Y);
     bool CutPassFiducial=false;
     if (retrieve_mode==0){
        pass_Fid->Fill(fiducial_radius);
        pass_Fid_wide->Fill(fiducial_radius);
        pass_Fid_old->Fill(fiducial_radius2);
     }else{
        mcpass_Fid->Fill(fiducial_radius);
        mcpass_Fid_wide->Fill(fiducial_radius);
        mcpass_Fid_old->Fill(fiducial_radius2);     
     }
     if ( fiducial_radius <= fiducial_cut_radius){
        CutPassFiducial=true;
//       if (cut_fiducial) continue;

     }
     if (CutPassRequire1TOF1SP && CutPassRequire1TOF0SP && CutPassTOF01 && CutPassChi2 && CutPassTKU_max_rad && CutPassTKU_max_rad && CutPassDiff_max_rad && CutPassTOF01Ext && CutPassUSTOF0 && CutPassFiducial) {    
        if (retrieve_mode==0){
          cuts_accept->Fill("Fiducial",1);
        }else{
           mccuts_accept->Fill("Fiducial",1);
        }
     }
 

//---------------------------------------------------------------------------------------
// Momentum at absorber cut. 
//---------------------------------------------------------------------------------------
     double momentum_absorber= sqrt((USretrieved.E(i).px*USretrieved.E(i).px)+
                                 (USretrieved.E(i).py*USretrieved.E(i).py)+ 
                                 (USretrieved.E(i).pz*USretrieved.E(i).pz));
     bool CutPasMom=true;
     if (retrieve_mode==0){
        pass_mom->Fill(momentum_absorber);
     }else{
        mcpass_mom->Fill(momentum_absorber);
     }
     if (momentum_absorber < momentum_lower_limit || momentum_absorber > momentum_upper_limit) {
        CutPasMom=false;
//        if (cut_absorber_momentum) continue;
     }
//     std::cerr<<"Pass Momentum: "<<CutPasMom<<" momentum_absorber: "<<momentum_absorber<<" momentum_lower_limit: "<<momentum_lower_limit<<" momentum_upper_limit: "<<momentum_upper_limit<<std::endl;
     if (CutPassRequire1TOF1SP && CutPassRequire1TOF0SP && CutPassTOF01 && CutPassChi2 && CutPassTKU_max_rad && CutPassTKU_max_rad && CutPassDiff_max_rad && CutPassTOF01Ext && CutPassUSTOF0 && CutPassFiducial && CutPasMom) {    
        USPostMom.append_instance(USretrieved.E(i));
        DSPostMom.append_instance(DSretrieved.E(i));
        if (retrieve_mode==0){
          cuts_accept->Fill("Momentum",1);
        }else{
           mccuts_accept->Fill("Mmentum",1);
        }
 
     }
//---------------------------------------------------------------------------------------
// Apply Cuts 
//---------------------------------------------------------------------------------------

     if (retrieve_mode==0){                                 //Create Histograms for Data
        if (CutPassTOF01 && CutPassChi2 && CutPassTKU_max_rad && CutPassDiff_max_rad && CutPassTOF01Ext && CutPassUSTOF0 && CutPassFiducial && CutPasMom){
        std::cerr<<" CutPassRequire1TOF0SP: "<<CutPassRequire1TOF0SP<<" NumberTOF0SP: "<<NumberTOF0SP<<" CutPassRequire1TOF1SP: "<<CutPassRequire1TOF1SP<<" NumberTOF1SP: "<<NumberTOF1SP<<std::endl;
        }
//        if (CutPassRequire1TOF1SP && CutPassTOF01 && CutPassChi2 && CutPassTKU_max_rad && CutPassDiff_max_rad && CutPassTOF01Ext && CutPassUSTOF0 && CutPassFiducial && CutPasMom){                         
           Hist_TOF0Require1SP->Fill(NumberTOF0SP);         
//        }
//        if (CutPassRequire1TOF0SP && CutPassTOF01 && CutPassChi2 && CutPassTKU_max_rad && CutPassDiff_max_rad && CutPassTOF01Ext && CutPassUSTOF0 && CutPassFiducial && CutPasMom){                   
           Hist_TOF1Require1SP->Fill(NumberTOF1SP);         
//        }
//        if (CutPassRequire1TOF0SP && CutPassRequire1TOF1SP && CutPassChi2 && CutPassTKU_max_rad && CutPassDiff_max_rad && CutPassTOF01Ext && CutPassUSTOF0 && CutPassFiducial && CutPasMom){                  
           Hist_TOF01->Fill(rawTOF01HitTime);               
//        }
//        if (CutPassRequire1TOF0SP && CutPassRequire1TOF1SP && CutPassTOF01 && CutPassTKU_max_rad && CutPassDiff_max_rad && CutPassTOF01Ext && CutPassUSTOF0 && CutPassFiducial && CutPasMom){                  
           Hist_Chi2->Fill(TKU_chiSquare_per_ndof);               
//        }
//        if (CutPassRequire1TOF0SP && CutPassRequire1TOF1SP && CutPassTOF01 && CutPassTOF01 && CutPassDiff_max_rad && CutPassTOF01Ext && CutPassUSTOF0 && CutPassFiducial && CutPasMom){                  
           Hist_TKU_max_radius->Fill(TKU_max_rad);               
//        }
//        if (CutPassRequire1TOF0SP && CutPassRequire1TOF1SP && CutPassTOF01 && CutPassTOF01 && CutPassTKU_max_rad && CutPassTOF01Ext && CutPassUSTOF0 && CutPassFiducial && CutPasMom){                  
           Hist_Diff_max_radius->Fill(diff_max_rad);               
//        }
//        if (CutPassRequire1TOF0SP && CutPassRequire1TOF1SP && CutPassTOF01 && CutPassTOF01 && CutPassTKU_max_rad && CutPassDiff_max_rad && CutPassUSTOF0 && CutPassFiducial && CutPasMom){                  
           Hist_TOF01Ext->Fill(TOF01Ext);               
//        }
//        if (CutPassRequire1TOF0SP && CutPassRequire1TOF1SP && CutPassTOF01 && CutPassTOF01 && CutPassTKU_max_rad && CutPassDiff_max_rad && CutPassTOF01Ext && CutPassFiducial && CutPasMom){
           if (CutPassUSTOF0){                 
              Hist_USTOF0->Fill("true",1);               
           }else{
              Hist_USTOF0->Fill("false",1);               
//           }
        }
//        if (CutPassRequire1TOF0SP && CutPassRequire1TOF1SP && CutPassTOF01 && CutPassTOF01 && CutPassTKU_max_rad && CutPassDiff_max_rad && CutPassTOF01Ext && CutPassUSTOF0 && CutPasMom){                  
           Hist_Fiducial->Fill(fiducial_radius);               
//        }
//        if (CutPassRequire1TOF0SP && CutPassRequire1TOF1SP && CutPassTOF01 && CutPassTOF01 && CutPassTKU_max_rad && CutPassDiff_max_rad && CutPassTOF01Ext && CutPassUSTOF0 && CutPassFiducial){                  
           Hist_MOM->Fill(momentum_absorber);               
//        }
     }else{                                                 //Create Histograms for MC
//        if (CutPassTOF01 && CutPassChi2 && CutPassTKU_max_rad && CutPassDiff_max_rad && CutPassTOF01Ext && CutPassUSTOF0 && CutPassFiducial && CutPasMom){
        std::cerr<<" CutPassRequire1TOF0SP: "<<CutPassRequire1TOF0SP<<" NumberTOF0SP: "<<NumberTOF0SP<<" CutPassRequire1TOF1SP: "<<CutPassRequire1TOF1SP<<" NumberTOF1SP: "<<NumberTOF1SP<<std::endl;
//        }
//        if (CutPassRequire1TOF1SP && CutPassTOF01 && CutPassChi2 && CutPassTKU_max_rad && CutPassDiff_max_rad && CutPassTOF01Ext && CutPassUSTOF0 && CutPassFiducial && CutPasMom){                         
           Histref_TOF0Require1SP->Fill(NumberTOF0SP);       
//        }
//        if (CutPassRequire1TOF0SP && CutPassTOF01 && CutPassChi2 && CutPassTKU_max_rad && CutPassDiff_max_rad && CutPassTOF01Ext && CutPassUSTOF0 && CutPassFiducial && CutPasMom){                   
           Histref_TOF1Require1SP->Fill(NumberTOF1SP);       
//        }
//        if (CutPassRequire1TOF0SP && CutPassRequire1TOF1SP && CutPassChi2 && CutPassTKU_max_rad && CutPassDiff_max_rad && CutPassTOF01Ext && CutPassUSTOF0 && CutPassFiducial && CutPasMom){                   
           Histref_TOF01->Fill(rawTOF01HitTime);             
//        }
//        if (CutPassRequire1TOF0SP && CutPassRequire1TOF1SP && CutPassTOF01 && CutPassTKU_max_rad && CutPassDiff_max_rad && CutPassTOF01Ext && CutPassUSTOF0 && CutPassFiducial && CutPasMom){                   
           Histref_Chi2->Fill(TKU_chiSquare_per_ndof);             
//        }
//        if (CutPassRequire1TOF0SP && CutPassRequire1TOF1SP && CutPassTOF01 && CutPassChi2 && CutPassDiff_max_rad && CutPassTOF01Ext && CutPassUSTOF0 && CutPassFiducial && CutPasMom){                   
           Histref_TKU_max_radius->Fill(TKU_max_rad);             
//        }
//       if (CutPassRequire1TOF0SP && CutPassRequire1TOF1SP && CutPassTOF01 && CutPassChi2 && CutPassTKU_max_rad && CutPassTOF01Ext && CutPassUSTOF0 && CutPassFiducial && CutPasMom){                   
           Histref_Diff_max_radius->Fill(diff_max_rad);             
//        }
//        if (CutPassRequire1TOF0SP && CutPassRequire1TOF1SP && CutPassTOF01 && CutPassChi2 && CutPassTKU_max_rad && CutPassDiff_max_rad && CutPassUSTOF0 && CutPassFiducial && CutPasMom){                   
           Histref_TOF01Ext->Fill(TOF01Ext);             
//        }
//        if (CutPassRequire1TOF0SP && CutPassRequire1TOF1SP && CutPassTOF01 && CutPassChi2 && CutPassTKU_max_rad && CutPassDiff_max_rad && CutPassTOF01Ext && CutPassFiducial && CutPasMom){                   
           if (CutPassUSTOF0){                 
              Histref_USTOF0->Fill("true",1);               
           }else{
              Histref_USTOF0->Fill("false",1);               
           }
//        }
//        if (CutPassRequire1TOF0SP && CutPassRequire1TOF1SP && CutPassTOF01 && CutPassChi2 && CutPassTKU_max_rad && CutPassDiff_max_rad && CutPassTOF01Ext && CutPassUSTOF0 && CutPasMom){                   
           Histref_Fiducial->Fill(fiducial_radius);             
//        }
//        if (CutPassRequire1TOF0SP && CutPassRequire1TOF1SP && CutPassTOF01 && CutPassChi2 && CutPassTKU_max_rad && CutPassDiff_max_rad && CutPassTOF01Ext && CutPassUSTOF0 && CutPassFiducial){                   
           Histref_MOM->Fill(momentum_absorber);             
//       }
     }
  
//---------------------------------------------------------------------------------------
// Misalignment Correction at centre of Absorber 
//---------------------------------------------------------------------------------------
     if (retrieve_mode==0){
        _USset.append_instance(USretrieved.E(i));
     }else{
        _USMCset.append_instance(USretrieved.E(i));
     }
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
     if (runnumber != LastRunNumber){
        UpdateRunInfo();
     }
  }
  double mean_residual_px=Residual_px->GetMean();
  double mean_residual_py=Residual_py->GetMean();
  double mean_residual_pz=Residual_pz->GetMean();
  if (retrieve_mode==0){
     std::cerr<<"Data Momentum Cut complete and uncorrected residuals calculated. Number of entries: "<<_USset.N()<<std::endl;
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
     PlotRunInfo();
  }else{
     std::cerr<<"Reference Momentum Cut complete and uncorrected residuals calculated. Number of entries: "<<_USMCset.N()<<std::endl;
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
  if(retrieve_mode==0){
     Hist_TOF0Require1SP->Write();
     Hist_TOF1Require1SP->Write();
     Hist_TOF01->Write();
     Hist_Chi2->Write();
     Hist_TKU_max_radius->Write();
     Hist_Diff_max_radius->Write();
     Hist_TOF01Ext->Write();
     Hist_USTOF0->Write();
     Hist_Fiducial->Write();
     Hist_MOM->Write();
     cuts_accept->Write();
  }else{
     Histref_TOF0Require1SP->Write();
     Histref_TOF1Require1SP->Write();
     Histref_TOF01->Write();
     Histref_Chi2->Write();
     Histref_TKU_max_radius->Write();
     Histref_Diff_max_radius->Write();
     Histref_TOF01Ext->Write();
     Histref_USTOF0->Write();
     Histref_Fiducial->Write();
     Histref_MOM->Write();
     mccuts_accept->Write();
  }

  TOF0file.close();
  TOF1file.close();
  TOF2file.close();
  Chi2file.close();
  doffile.close();
}

/////////////////////////////////////////////////////////////////////
// Function Dataselection
/////////////////////////////////////////////////////////////////////


void MCSAnalysis::dataSelection(int mode){

  chain->SetBranchAddress("RunNumber", &runnumber);
  chain->SetBranchAddress("TOFBranch", &tofevent);
  chain->SetBranchAddress("SciFiBranch", &scifievent);
  chain->SetBranchAddress("CkovBranch", &ckovevent);
  chain->SetBranchAddress("KLBranch", &klevent);
  chain->SetBranchAddress("EMRBranch", &emrevent);
//  chain->SetBranchAddress("MCEvent", &mcevent);
  chain->SetBranchAddress("GlobalBranch", &globalevent);

  double dataDS11 = _sys["dataDS11"];
  double dataDS53 = _sys["dataDS53"];
  double dataUS53 = _sys["dataUS53"];
  double dataUS11 = _sys["dataUS11"];
  double datadiffuserlow = _sys["datadiffuserlow"];
  double datadiffuserhigh = _sys["datadiffuserhigh"];
  double dataTOF0pos = _sys["dataTOF0pos"];
  double dataTOF1pos = _sys["dataTOF1pos"];

  // Cuts to apply to data 
  // Restrict the number of entries to less than or equal to the mcchain entries
  int Nentries = chain->GetEntries();
  // < mcchain->GetEntries() ? 
  //  chain->GetEntries() : mcchain->GetEntries();
  // Loop over all tree entries.
  Collection USPostTKU, DSPostTKU;
  Collection Collect_EndOfDS, Collect_StartOfDS, Collect_USTOF0, Collect_USTOF1, Collect_StartofUST, Collect_EndofUST, Collect_DiffuserStart, Collect_DiffuserEnd;

// Open files to save TOF data to.
   std::ofstream TOF0file;
   TOF0file.open(TOF0dataname);
   std::ofstream TOF1file;
   TOF1file.open(TOF1dataname);
   std::ofstream TOF2file;
   TOF2file.open(TOF2dataname);
// Open files to save Chi2 and dof data to.
   std::ofstream Chi2file;
   Chi2file.open(Chi2dataname);
   std::ofstream doffile;
   doffile.open(dofdataname);

//  double pz = 0.;
  // Loop over all tree entries.
//  int Events_Pass = 0;
//  int Events_Fail = 0;
  std::cerr <<"Number of entries in data file: "<< Nentries << std::endl;
//  int Nevents = 0;
  // Nentries *= int(_sys["FracEvents"]) > 0 ? _sys["FracEvents"] : 1.0;
  for (int i=0; i<Nentries; i++){
    std::cerr<<"**********************************************"<< std::endl;
    std::cerr<<"Data Event number " << i << " of " << Nentries << std::endl;
    chain->GetEntry(i);
//     if (i%100!=0) {
//        continue;
//     }     
    std::cerr<<"Event Got"<<i<<"\n";
    if (i>=10000) break;
    // Set cuts based on the TOFs, ckov, kl, and EMR information
    // Locate the tracker reference planes. To be agnostic locate
    // the downstream most station of the upstream tracker and the
    // upstream most station of the downstream tracker.
     cuts_accept->Fill("All Events",1);
    // if (scifievent->scifitracks().size() != 2) continue;
    // if ( !SelectMomentum() ) continue;

     multiVars globalVars;
     globalVars.USTOF0=reset_Vars();
     globalVars.USTOF1=reset_Vars();
     globalVars.UScentre_absorber=reset_Vars();
     globalVars.UStrackerUS11=reset_Vars();
     globalVars.UStrackerUS53=reset_Vars();
     globalVars.DiffuserStart=reset_Vars();
     globalVars.DiffuserEnd=reset_Vars();
     globalVars.USend_of_DStracker=reset_Vars();
     globalVars.DScentre_absorber=reset_Vars();
     globalVars.USstart_of_DStracker=reset_Vars();
//    std::cerr<<"centre_absorber_US created with default values"<<std::endl;
//    display_Vars(globalVars.UScentre_absorber);
//    std::cerr<<"centre_absorber_DS created with default values"<<std::endl;
//    display_Vars(globalVars.DScentre_absorber);


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
  
//    std::cerr<<"Adding to AllTOF"<<std::endl;
//    display_Vars(globalVars.UScentre_absorber);
     USPostTKU.append_instance(globalVars.UScentre_absorber);
//    std::cerr<<"Size of USAllTOF: "<<USAllTOF.N()<<std::endl;
//    display_Vars(globalVars.DScentre_absorber);
     DSPostTKU.append_instance(globalVars.DScentre_absorber);
//    std::cerr<<"Size of DSAllTOF: "<<DSAllTOF.N()<<std::endl;
     cuts_accept->Fill("1 US track",1);
     pass_TKU->Fill("Pass",1);

   
//----------------------------------------------------------------------
// Retrieve Global Data - US track and DS track if present. Cut if no US track
//----------------------------------------------------------------------

    globalVars=read_globals(dataDS11, dataDS53, dataUS53, dataUS11, datadiffuserlow, datadiffuserhigh, dataTOF0pos, dataTOF1pos);

//    std::cerr<<"centre_absorber_US has updated globals data"<<std::endl;
//    display_Vars(globalVars.UStrackerUS53);
//    display_Vars(globalVars.DiffuserStart);
    if (!globalVars.UScentre_absorber.isgood) {
//      std::cerr<<"**Failed US track Cut - Why not cut earlier ********************"<<std::endl;
      continue;
    }
    if (!globalVars.DScentre_absorber.isgood) {
//      std::cerr<<"**No DS ********************"<<std::endl;
    }
    cuts_accept->Fill("US Track",1);

//----------------------------------------------------------------------
// Save TOF information for future processing
//----------------------------------------------------------------------

    int TOF0num = tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArraySize();
    int TOF1num = tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArraySize();
    int TOF2num = tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArraySize();
    TOF0file << TOF0num << " ";
    TOF1file << TOF1num << " ";
    TOF2file << TOF2num << " ";
    for (int TOFnum=0; TOFnum<TOF0num; TOFnum++){TOF0file << tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[TOFnum].GetTime()<< " ";}
    for (int TOFnum=0; TOFnum<TOF1num; TOFnum++){TOF1file << tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[TOFnum].GetTime()<< " ";}
    for (int TOFnum=0; TOFnum<TOF2num; TOFnum++){TOF2file << tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray()[TOFnum].GetTime()<< " ";}
    TOF0file << "\n";
    TOF1file << "\n";
    TOF2file << "\n";
    
//----------------------------------------------------------------------
// Save Chi2 and dof for future processing
//----------------------------------------------------------------------

    std::vector<MAUS::SciFiTrack*> tracks = scifievent->scifitracks();
    std::vector<MAUS::SciFiTrack*>::iterator track_iter;

    if(tracks.size() == 0){
        Chi2file << "10000.0" << "\n";
        doffile << "100" << "\n";
//      std::cerr<<"************************No tracks - why not cut earlier**************************"<<std::endl;
    }else{ 
      // use if only want to take the first track
        track_iter = tracks.begin();
        double TKU_chiSquare = (*track_iter)->chi2();
        int TKU_ndof = (*track_iter)->ndf();
        Chi2file << TKU_chiSquare << "\n";
        doffile << TKU_ndof << "\n";
    }



//------------------------------------------------------------------------------------------------------
// Saving of events that pass all cuts
//------------------------------------------------------------------------------------------------------

//    std::cerr<<"**Event accepted"<<std::endl;
    _USset.append_instance(globalVars.UScentre_absorber);
    _DSset.append_instance(globalVars.DScentre_absorber);
    Collect_EndOfDS.append_instance(globalVars.USend_of_DStracker);
    Collect_StartOfDS.append_instance(globalVars.USstart_of_DStracker);
    Collect_USTOF0.append_instance(globalVars.USTOF0);
    Collect_USTOF1.append_instance(globalVars.USTOF1);
    Collect_StartofUST.append_instance(globalVars.UStrackerUS53);
    Collect_EndofUST.append_instance(globalVars.UStrackerUS11);
    Collect_DiffuserStart.append_instance(globalVars.DiffuserStart);
    Collect_DiffuserEnd.append_instance(globalVars.DiffuserEnd);
    
//    _DSset.append_instance(globalVars.DScentre_absorber);
    cuts_accept->Fill("All Cuts",1);
 
  }
  _USset.save_collection(USDataname);
  _DSset.save_collection(DSDataname);
  Collect_EndOfDS.save_collection(EndofDStrackerdataname);
  Collect_StartOfDS.save_collection(StartofDStrackerdataname);
  Collect_USTOF0.save_collection(Collect_USTOF0dataname);
  Collect_USTOF1.save_collection(Collect_USTOF1dataname);
  Collect_StartofUST.save_collection(Collect_StartofUSTdataname);
  Collect_EndofUST.save_collection(Collect_EndofUSTdataname);
  Collect_DiffuserStart.save_collection(Collect_DiffuserStartdataname);
  Collect_DiffuserEnd.save_collection(Collect_DiffuserEnddataname);
// Close TOF, Chi2 and dof files
  TOF0file.close();
  TOF1file.close();
  TOF2file.close();
  Chi2file.close();
  doffile.close();

  
  std::cerr<<"Data: Size of USset: "<<_USset.N()<<std::endl;
  std::cerr<<"Data: Size of DSset: "<<_DSset.N()<<std::endl;
  std::cerr<<"Data: Size of Collect_EndOfDS: "<<Collect_EndOfDS.N()<<std::endl;
  std::cerr<<"Data: Collect_StartOfDS: "<<Collect_StartOfDS.N()<<std::endl;
  std::cerr<<"Data: Size of Collect_USTOF0: "<<Collect_USTOF0.N()<<std::endl;
  std::cerr<<"Data: Size of Collect_USTOF1: "<<Collect_USTOF1.N()<<std::endl;
  std::cerr<<"Data: Size of Collect_StartofUST: "<<Collect_StartofUST.N()<<std::endl;
  std::cerr<<"Data: Size of Collect_EndofUST: "<<Collect_EndofUST.N()<<std::endl;
  std::cerr<<"Data: Size of Collect_DiffuserStart: "<<Collect_DiffuserStart.N()<<std::endl;
  std::cerr<<"Data: Size of Collect_DiffuserEnd: "<<Collect_DiffuserEnd.N()<<std::endl;


//------------------------------------------------------------------------------------------------------
// Make beam histograms
//------------------------------------------------------------------------------------------------------
  

  make_beam_histograms(USPostTKU, "Upstream, Data", "data_03US_PostTKU");
  make_beam_histograms(DSPostTKU, "Downstream, Data", "data_03DS_PostTKU"); 
}

/////////////////////////////////////////////////////////////////////
// Function Reference Selection
/////////////////////////////////////////////////////////////////////


void MCSAnalysis::referenceSelection(){
  // Set addresses for tree selection

  std::cerr<<"Starting reference Selection"<<std::endl;
  refchain->SetBranchAddress("RunNumber", &runnumber);
  refchain->SetBranchAddress("TOFBranch", &tofevent);
  refchain->SetBranchAddress("SciFiBranch", &scifievent);
  refchain->SetBranchAddress("CkovBranch", &ckovevent);
  refchain->SetBranchAddress("KLBranch", &klevent);
  refchain->SetBranchAddress("EMRBranch", &emrevent);
  refchain->SetBranchAddress("GlobalBranch", &globalevent);

//  std::cerr<<"Branches set"<<std::endl;

  double refDS11 = _sys["refDS11"];
  double refDS53 = _sys["refDS53"];
  double refUS53 = _sys["refUS53"];
  double refUS11 = _sys["refUS11"];
  double refdiffuserlow = _sys["refdiffuserlow"];
  double refdiffuserhigh = _sys["refdiffuserhigh"];
  double refTOF0pos = _sys["refTOF0pos"];
  double refTOF1pos = _sys["refTOF1pos"];

  

  // Restrict the number of entries to less than or equal to the refchain entries
  
//  std::cerr<<"Getting Nentries"<<std::endl;
  int Nentries = refchain->GetEntries();
//  std::cerr<<"Nentries set"<<std::endl;
  Collection USPostTKU, DSPostTKU;
  Collection Collect_EndOfDS, Collect_StartOfDS, Collect_USTOF0, Collect_USTOF1, Collect_StartofUST, Collect_EndofUST, Collect_DiffuserStart, Collect_DiffuserEnd;

// Open files to save TOF data to.
   std::ofstream TOF0file;
   TOF0file.open(TOF0refname);
   std::ofstream TOF1file;
   TOF1file.open(TOF1refname);
   std::ofstream TOF2file;
   TOF2file.open(TOF2refname);
   std::ofstream Chi2file;
   Chi2file.open(Chi2refname);
   std::ofstream doffile;
   doffile.open(dofrefname);

//  double pz = 0.;
  std::cerr <<"Number of entries in data file: "<< Nentries << std::endl;
  // Loop over all tree entries.
  for (int i=0; i<Nentries; i++){
//    std::cerr<<"**********************************************"<< std::endl;
 //   std::cerr<<"Reference Event number " << i << " of " << Nentries << std::endl;
    refchain->GetEntry(i);
    if (i%100000==0) std::cout<<"Event "<<i<<"\n"; 
    if (i>=10000) break;
    mccuts_accept->Fill("All Events",1);

    multiVars globalVars;
    globalVars.USTOF0=reset_Vars();
    globalVars.USTOF1=reset_Vars();
    globalVars.UScentre_absorber=reset_Vars();
    globalVars.UStrackerUS11=reset_Vars();
    globalVars.UStrackerUS53=reset_Vars();
    globalVars.DiffuserStart=reset_Vars();
    globalVars.DiffuserEnd=reset_Vars();
    globalVars.USend_of_DStracker=reset_Vars();
    globalVars.DScentre_absorber=reset_Vars();
//    std::cerr<<"centre_absorber_US created with default values"<<std::endl;
//    display_Vars(globalVars.UScentre_absorber);
//    std::cerr<<"centre_absorber_DS created with default values"<<std::endl;
//    display_Vars(globalVars.DScentre_absorber);



//---------------------------------------------------------------------------------------
// Require only 1 US track
//---------------------------------------------------------------------------------------

//  std::cerr<<" Number of tracks: "<<scifievent->scifitracks().size()<<std::endl;
  int number_upstream_tracks=0;
  for(size_t j=0; j<scifievent->scifitracks().size(); j++){
      int tracker = scifievent->scifitracks()[j]->tracker();
      if(tracker==0) number_upstream_tracks++;  
  }
  if (number_upstream_tracks==0){
//    std::cerr<<"**Fail 1 US track cut, no tracks"<<std::endl;
    mcpass_TKU->Fill("Fail",1);
    if(cut_TKU_1track) continue;    
  } else if(number_upstream_tracks>1){
//    std::cerr<<"**Fail 1 US track cut, more than 1 track"<<std::endl; 
    mcpass_TKU->Fill("Fail",1);
    if(cut_TKU_1track) continue;   
  }

/*
  int nUStracks= 0;
  if (scifievent->scifitracks().size()>=1){
    for ( size_t j=0; j<scifievent->scifitracks().size(); j++){
      if(scifievent->scifitracks()[j]->tracker()==0) nUStracks++ ;
    }
  }else{
    std::cerr<<"**Fail 1 US track cut, no tracks"<<std::endl;
    mcpass_TKU->Fill("Fail",1);
    continue;
  }
  if (nUStracks!=1){
    std::cerr<<"**Fail 1 US track cut, more than 1 track"<<std::endl; 
    mcpass_TKU->Fill("Fail",1);
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
    mccuts_accept->Fill("1 US track",1);
    mcpass_TKU->Fill("Pass",1);
   
//----------------------------------------------------------------------
// Retrieve Global Data - US track and DS track if present. Cut if no US track
//----------------------------------------------------------------------

    globalVars=read_globals(refDS11, refDS53, refUS53, refUS11, refdiffuserlow, refdiffuserhigh, refTOF0pos, refTOF1pos);
//    std::cerr<<"centre_absorber_US has updated globals data"<<std::endl;
//    display_Vars(globalVars.UScentre_absorber);
//    display_Vars(globalVars.DScentre_absorber);
    if (!globalVars.UScentre_absorber.isgood) {
//    std::cerr<<"**Failed US track Cut - Why not cut earlier ********************"<<std::endl;
    continue;
    }
    double Ptot =sqrt(globalVars.UScentre_absorber.px*globalVars.UScentre_absorber.px+globalVars.UScentre_absorber.py*globalVars.UScentre_absorber.py+globalVars.UScentre_absorber.pz*globalVars.UScentre_absorber.pz);
    HPt_Ch12->Fill(Ptot);

    mccuts_accept->Fill("US Track",1);


//----------------------------------------------------------------------
// Save TOF information for future processing
//----------------------------------------------------------------------

    int TOF0num = tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArraySize();
    int TOF1num = tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArraySize();
    int TOF2num = tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArraySize();
    TOF0file << TOF0num << " ";
    TOF1file << TOF1num << " ";
    TOF2file << TOF2num << " ";
    for (int TOFnum=0; TOFnum<TOF0num; TOFnum++){TOF0file << tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[TOFnum].GetTime()<< " ";}
    for (int TOFnum=0; TOFnum<TOF1num; TOFnum++){TOF1file << tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[TOFnum].GetTime()<< " ";}
    for (int TOFnum=0; TOFnum<TOF2num; TOFnum++){TOF2file << tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray()[TOFnum].GetTime()<< " ";}
    TOF0file << "\n";
    TOF1file << "\n";
    TOF2file << "\n";

//----------------------------------------------------------------------
// Save Chi2 and dof for future processing
//----------------------------------------------------------------------

    std::vector<MAUS::SciFiTrack*> tracks = scifievent->scifitracks();
    std::vector<MAUS::SciFiTrack*>::iterator track_iter;

    if(tracks.size() == 0){
        Chi2file << "10000.0" << "\n";
        doffile << "100" << "\n";
//      std::cerr<<"************************No tracks - why not cut earlier**************************"<<std::endl;
    }else{ 
      // use if only want to take the first track
        track_iter = tracks.begin();
        double TKU_chiSquare = (*track_iter)->chi2();
        int TKU_ndof = (*track_iter)->ndf();
        Chi2file << TKU_chiSquare << "\n";
        doffile << TKU_ndof << "\n";
    }


//------------------------------------------------------------------------------------------------------
// Saving of events that pass all cuts
//------------------------------------------------------------------------------------------------------

//    std::cerr<<"**Event accepted"<<std::endl;
    _USMCset.append_instance(globalVars.UScentre_absorber);
    _DSMCset.append_instance(globalVars.DScentre_absorber);
    Collect_EndOfDS.append_instance(globalVars.USend_of_DStracker);
    Collect_StartOfDS.append_instance(globalVars.USstart_of_DStracker);
    Collect_USTOF0.append_instance(globalVars.USTOF0);
    Collect_USTOF1.append_instance(globalVars.USTOF1);
    Collect_StartofUST.append_instance(globalVars.UStrackerUS53);
    Collect_EndofUST.append_instance(globalVars.UStrackerUS11);
    Collect_DiffuserStart.append_instance(globalVars.DiffuserStart);
    Collect_DiffuserEnd.append_instance(globalVars.DiffuserEnd);

//   std::cerr<<"Cut added to MC data"<<std::endl;
    mccuts_accept->Fill("All Cuts",1);

    
  }
  _USMCset.save_collection(USRefname);
  _DSMCset.save_collection(DSRefname);
  Collect_EndOfDS.save_collection(EndofDStrackerrefname);
  Collect_StartOfDS.save_collection(StartofDStrackerrefname);
  Collect_USTOF0.save_collection(Collect_USTOF0refname);
  Collect_USTOF1.save_collection(Collect_USTOF1refname);
  Collect_StartofUST.save_collection(Collect_StartofUSTrefname);
  Collect_EndofUST.save_collection(Collect_EndofUSTrefname);
  Collect_DiffuserStart.save_collection(Collect_DiffuserStartrefname);
  Collect_DiffuserEnd.save_collection(Collect_DiffuserEndrefname);
  TOF0file.close();
  TOF1file.close();
  TOF2file.close();
  Chi2file.close();
  doffile.close();

  std::cerr<<"Reference: Size of USset: "<<_USMCset.N()<<std::endl;
  std::cerr<<"Reference: Size of DSset: "<<_DSMCset.N()<<std::endl;
  std::cerr<<"Reference: Size of Collect_EndOfDS: "<<Collect_EndOfDS.N()<<std::endl;
  std::cerr<<"Reference: Collect_StartOfDS: "<<Collect_StartOfDS.N()<<std::endl;
  std::cerr<<"Reference: Size of Collect_USTOF0: "<<Collect_USTOF0.N()<<std::endl;
  std::cerr<<"Reference: Size of collect_USTOF1: "<<Collect_USTOF1.N()<<std::endl;
  std::cerr<<"Reference: Size of Collect_StartofUST: "<<Collect_StartofUST.N()<<std::endl;
  std::cerr<<"Reference: Size of Collect_EndofUST: "<<Collect_EndofUST.N()<<std::endl;
  std::cerr<<"Reference: Size of Collect_DiffuserStart: "<<Collect_DiffuserStart.N()<<std::endl;
  std::cerr<<"Reference: Size of Collect_DiffuserEnd: "<<Collect_DiffuserEnd.N()<<std::endl;


  make_beam_histograms(USPostTKU, "Upstream, Data Reference", "ref_03US_PostTKU");
  make_beam_histograms(DSPostTKU, "Downstream, Data Reference", "ref_03DS_PostTKU"); 

}

void MCSAnalysis::generateMCSResponse(){

  chain->SetBranchAddress("TOFBranch", &tofevent);
  chain->SetBranchAddress("SciFiBranch", &scifievent);
  chain->SetBranchAddress("CkovBranch", &ckovevent);
  chain->SetBranchAddress("KLBranch", &klevent);
  chain->SetBranchAddress("EMRBranch", &emrevent);
  chain->SetBranchAddress("MCEvent", &mcevent);
  // Loop over the training sample.
  Collection USRecSet, DSRecSet, USVirtSet, DSVirtSet;
  Collection USAllTOF, DSAllTOF, USPreRadSel, DSPreRadSel;
  Collection USVAllTOF, DSVAllTOF, USVPreRadSel, DSVPreRadSel;
  // Collection USeRecSet, DSeRecSet, USmuRecSet, DSmuRecSet;
  int nMuAbsSel=0, nEAbsSel=0, nMuAbsAll=0, nEAbsAll=0;
  TH1D* mom_resUS = new TH1D("mom_resUS", ";p_{rec} - p_{MC} (MeV/c)",
			     2000, -100, 100);
  TH1D* mom_resDS = new TH1D("mom_resDS", ";p_{rec} - p_{MC} (MeV/c)",
			     2000, -100, 100);
  TH2D* mom_responseUS = new TH2D("mom_responseUS", ";p_{rec} (MeV/c); p_{MC} (MeV/c)",
				  300, 0, 300, 300, 0, 300);
  TH2D* mom_responseDS = new TH2D("mom_responseDS", ";p_{rec} (MeV/c); p_{MC} (MeV/c)",
				  300, 100, 300, 300, 0, 300);
  TH2D* mom_responseABS = new TH2D("mom_responseABS", ";p_{rec} (MeV/c); p_{MC} (MeV/c)",
				  300, 100, 300, 300, 0, 300);
  int ngood=0;
  int Nentries = chain->GetEntries();
  for (int j=0; j<Nentries; j++){
    chain->GetEntry(j);
    if (j%100000==0) std::cout<<"MC Event "<<j<<", selected "<<ngood<<" events.\n"; 
    // if(fabs(mcevent->GetPrimary()->GetParticleId()) != 13) continue;
    // Select events that produce virtual plane hits and pass the data selection cuts.
    Vars USAbsHit, DSAbsHit, USTrackerRefHit, DSTrackerRefHit;
    double pz = 0.;
    mccuts_accept->Fill("All Events",1);
    if ( !PIDSelection(false) ) continue; // event_ok=false;
    if( !findVirtualPlanes() ) continue;
    mccuts_accept->Fill("Found Virtual Planes",1);
    FillVarsVirtual(USAbsHit, USabsPlaneI);
    FillVarsVirtual(DSAbsHit, DSabsPlaneI);
    FillVarsVirtual(USTrackerRefHit, USrefplaneI);
    FillVarsVirtual(DSTrackerRefHit, DSrefplaneI);
    if (mcevent->GetVirtualHits()->at(USrefplaneI).GetParticleId()==-13 &&
	mcevent->GetVirtualHits()->at(DSrefplaneI).GetParticleId()==-11)
      nEAbsAll++;
    if (mcevent->GetVirtualHits()->at(USrefplaneI).GetParticleId()==-13)
      nMuAbsAll++;
    // Apply selection cuts as with the data
    bool event_ok=true;
    if ( !MatchUSDS() ) {
      if (jUS == -1 || kUS == -1){
	cuts_accept->Fill("US Tracker Found", 1);
	event_ok=false;
      } 
    }
    if (event_ok) mccuts_accept->Fill("US Track Found",1);

    if (event_ok) pz = MomentumFromTOF(false);
    if (event_ok) mccuts_accept->Fill("TOF Selection",1);
    // if (scifievent->scifitracks().size() != 2) continue;
    FillCollectionSciFi(USPreRadSel, jUS, kUS, pz, 0);
    FillCollectionSciFi(DSPreRadSel, jDS, kDS, pz, 1);
    USVPreRadSel.append_instance(USAbsHit);
    DSVPreRadSel.append_instance(DSAbsHit);
    if ( !RadialSelection(pz, 19948.8,meanp) ) event_ok=false;
    if (event_ok) mccuts_accept->Fill("Fiducial Selection",1);
    if ( !RadialSelection(pz,_sys["diffpos"],100) ) event_ok=false;
    if (event_ok) mctrue_mom->Fill(mcevent->GetVirtualHits()->at(USrefplaneI).GetMomentum().z());
    if (event_ok) mctrue_mom->Fill(mcevent->GetVirtualHits()->at(USrefplaneI).GetMomentum().z());
    Vars USSciFiRec, DSSciFiRec;
    if (event_ok){
      FillVarsSciFi(USSciFiRec, jUS, kUS, pz, 0);
      FillVarsSciFi(DSSciFiRec, jDS, kDS, pz, 1);
      USRecSet.append_instance(USSciFiRec);
      DSRecSet.append_instance(DSSciFiRec);
      if (mcevent->GetVirtualHits()->at(USrefplaneI).GetParticleId()==-13 &&
	  mcevent->GetVirtualHits()->at(USrefplaneI).GetParticleId()==-11)
	nEAbsSel++;
      if (mcevent->GetVirtualHits()->at(USrefplaneI).GetParticleId()==-13)
	nMuAbsSel++;
      
      mom_resUS->Fill(USSciFiRec.pz - USTrackerRefHit.pz);
      mom_resDS->Fill(DSSciFiRec.pz - DSTrackerRefHit.pz);
      mom_responseUS->Fill(USSciFiRec.pz, USTrackerRefHit.pz);
      mom_responseDS->Fill(DSSciFiRec.pz, DSTrackerRefHit.pz);
      mom_responseABS->Fill(DSSciFiRec.pz, USAbsHit.pz);
    }
    USVirtSet.append_instance(USAbsHit);
    DSVirtSet.append_instance(DSAbsHit);
    FillMCSResponse(event_ok, USSciFiRec, DSSciFiRec, USAbsHit, DSAbsHit);
    FillMuScattResponse(event_ok, USSciFiRec, DSSciFiRec, USAbsHit, DSAbsHit);
    ngood++;
  }
  std::cout<<"For all events simulated there are "<<nEAbsAll<<" decay electrons and "<<nMuAbsAll<<" muons.\n";
  std::cout<<"For the selected MC events there are "<<nEAbsSel<<" decay electrons and "<<nMuAbsSel<<" muons.\n";
  make_beam_histograms(USVirtSet, "Upstream, Data", "VirtMCUS");
  make_beam_histograms(DSVirtSet, "Downstream, Data", "VirtMCDS");
  //make_beam_histograms(USAllTOF, "Upstream, Reconstructed Simulation", "recMCUS_alltof");
  //make_beam_histograms(DSAllTOF, "Downstream, Reconstructed Simulation", "recMCDS_alltof");
  make_beam_histograms(USPreRadSel, "Upstream, Reconstructed Simulation", "recMCUS_prerad");
  make_beam_histograms(DSPreRadSel, "Downstream, Reconstructed Simulation", "recMCDS_prerad");
  make_scattering_acceptance_histograms(USVPreRadSel,
					DSVPreRadSel, DSPreRadSel,
					"Virtual Projection","VirtPreRad");

  make_beam_histograms(USRecSet, "Upstream, Reconstructed Simulation", "recMCUS");
  make_beam_histograms(DSRecSet, "Downstream, Reconstructed Simulation", "recMCDS");
  make_scattering_acceptance_histograms(USVirtSet,
					DSVirtSet,DSRecSet,
					"Virtual Projection","VirtProj");
  mom_resUS->Write();
  mom_resDS->Write();
  mom_responseUS->Write();
  mom_responseDS->Write();
  mom_responseABS->Write();
}

void MCSAnalysis::TruthData(int mode){

  mcchain->SetBranchAddress("MCEvent", &mcevent);
  mcemptychain->SetBranchAddress("MCEvent", &mcevent);
  //Collection USTruthSet, DSTruthSet;
  Collection USTruthzero, DSTruthzero;
//  Vars USTruthzeroHit, DSTruthzeroHit;
  int ngood=0;
//  bool event_ok=true;
  int Nentries;
  if (mode == 5) {
     Nentries = mcemptychain->GetEntries();
  }
  else {
     Nentries = mcchain->GetEntries();
  }
  for (int j=0; j<Nentries; j++){
//	  std::cout << j << std::endl;
//	  std::cout << ngood << std::endl;
    if (mode == 5) {
	    mcemptychain->GetEntry(j);
    }
    else {
	    mcchain->GetEntry(j);
    }

    if (j%100000==0) std::cout<<"MC Event "<<j<<", selected "<<ngood<<" events.\n"; 
    // if(fabs(mcevent->GetPrimary()->GetParticleId()) != 13) continue;
    // Select events that produce virtual plane hits and pass the data selection cuts.

    if( !findVirtualPlanes() ) continue;
    if ( !TruthMatchUSDS() ) {
      if (jUS == -1 || kUS == -1){
	continue;
      //  if (jDS == -1 || kDS == -1){
      //  	continue;
      //}
      }
    }

    double pz = 0.;
    if ( !TruthTime(false) ) continue; // event_ok=false;
    //std::cout << "USrefplaneI " << USrefplaneI << std::endl;
    //std::cout << "DSrefplaneI " << DSrefplaneI << std::endl;

    if ( !TruthRadialSelection(pz, 19948.8,meanp,USrefplaneI) ) continue;
    if ( !TruthRadialSelection(pz,_sys["diffpos"],100,USrefplaneI) ) continue;
    Vars USTruthAbsHit, DSTruthAbsHit;
    FillVarsVirtual(USTruthAbsHit, USabsPlaneI);
    FillVarsVirtual(DSTruthAbsHit, DSabsPlaneI);
    std::cout << "mcevent->GetVirtualHits()->at(USabsPlaneI).GetMomentum().x() " << mcevent->GetVirtualHits()->at(USabsPlaneI).GetMomentum().x() << std::endl;
    USTruthSet.append_instance(USTruthAbsHit);
    DSTruthSet.append_instance(DSTruthAbsHit);
    std::cout << USTruthSet.N() << std::endl;
	    std::cout << "USTruthAbsHit.px " << USTruthAbsHit.px << std::endl;
	    std::cout << "USTruthSet.E(USTruthSet.N()-1).px " << USTruthSet.E(USTruthSet.N()-1).px << std::endl;
    /*
    double thetaXMC = atan(mcevent->GetVirtualHits()->at(USabsPlaneI).GetMomentum().x()/mcevent->GetVirtualHits()->at(USabsPlaneI).GetMomentum().z()) - atan(mcevent->GetVirtualHits()->at(DSabsPlaneI).GetMomentum().x()/mcevent->GetVirtualHits()->at(DSabsPlaneI).GetMomentum().z());
    if (thetaXMC>=0 && thetaXMC<=0.004){
          std::cout << "atan(mcevent->GetVirtualHits()->at(USabsPlaneI).GetMomentum().x()/mcevent->GetVirtualHits()->at(USabsPlaneI).GetMomentum().z()) " << atan(mcevent->GetVirtualHits()->at(USabsPlaneI).GetMomentum().x()/mcevent->GetVirtualHits()->at(USabsPlaneI).GetMomentum().z()) << std::endl;
          std::cout << "atan(mcevent->GetVirtualHits()->at(DSabsPlaneI).GetMomentum().x()/mcevent->GetVirtualHits()->at(DSabsPlaneI).GetMomentum().z()) " << atan(mcevent->GetVirtualHits()->at(DSabsPlaneI).GetMomentum().x()/mcevent->GetVirtualHits()->at(DSabsPlaneI).GetMomentum().z()) << std::endl;
          std::cout << "mcevent->GetVirtualHits()->at(USabsPlaneI).GetMomentum().x()/mcevent->GetVirtualHits()->at(USabsPlaneI).GetMomentum().z() " << mcevent->GetVirtualHits()->at(USabsPlaneI).GetMomentum().x()/mcevent->GetVirtualHits()->at(USabsPlaneI).GetMomentum().z() << std::endl;
          std::cout << "mcevent->GetVirtualHits()->at(DSabsPlaneI).GetMomentum().x()/mcevent->GetVirtualHits()->at(DSabsPlaneI).GetMomentum().z() " << mcevent->GetVirtualHits()->at(DSabsPlaneI).GetMomentum().x()/mcevent->GetVirtualHits()->at(DSabsPlaneI).GetMomentum().z() << std::endl;
	  FillVarsVirtual(USTruthzeroHit, USabsPlaneI);
	  FillVarsVirtual(DSTruthzeroHit, DSabsPlaneI);
	  USTruthzero.append_instance(USTruthzeroHit);
	  DSTruthzero.append_instance(DSTruthzeroHit);
    }
    */
  }
  TruthGraph(USTruthSet, DSTruthSet);
  /*
    for (int i=0; i<USTruthSet.N(); i++){
	    std::cout << "USTruthSet.E(i).px " << USTruthSet.E(i).px << std::endl;
    }
    */
  make_beam_histograms(USTruthSet, "Upstream, Truth", "TruthUS");
  make_beam_histograms(DSTruthSet, "Downstream, Truth", "TruthDS");
  //make_beam_histograms(USTruthzero, "Upstream, Truth zero", "TruthUSzero");
  //make_beam_histograms(DSTruthzero, "Downstream, Truth zero", "TruthDSzero");
}

bool MCSAnalysis::MatchUSDS(){

  bool trackmatched = true;
  jUS=-1;
  jDS=-1;
  kUS=-1;
  kDS=-1;
 
  int notrackUS = 0;
  int notrackDS = 0;
  for ( size_t j=0; j<scifievent->scifitracks().size(); j++){
//      double maxUS=0.0, minDS=44000;
      int tracker = scifievent->scifitracks()[j]->tracker();
      if(tracker==0){
	  notrackUS += 1;
      }
      if(tracker==1){
	  notrackDS += 1;
      }
  }
  trackno->Fill(notrackUS,notrackDS);
  
  if( scifievent->scifitracks().size() == 1 || 
      scifievent->scifitracks().size() == 2){
    for ( size_t j=0; j<scifievent->scifitracks().size(); j++){
      int npoints = scifievent->scifitracks()[j]->scifitrackpoints().size();
      double maxUS=0.0, minDS=44000;
      int tracker = scifievent->scifitracks()[j]->tracker();
      for ( int k=0; k<npoints; k++){
	double zpos = 
	  scifievent->scifitracks()[j]->scifitrackpoints()[k]->pos().z();
	if(tracker==0 && zpos > maxUS){
	  maxUS = zpos;
	  kUS = k;
	  jUS = j;
	}
	if(tracker==1 && zpos < minDS){      
	  minDS = zpos;
	  kDS = k;
	  jDS = j;
	}
      }
      if (jUS != -1 && kUS != -1 && jDS != -1 && kDS != 1) {
	USrefplaneZ = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->pos().z();
	DSrefplaneZ = minDS;
	// std::cout<<"Identified US track by ["<<jUS<<", "<<kUS<<"] at z = "<<USrefplaneZ<<" mm and DS track by ["<<jDS<<", "<<kDS<<"] at z = "<<DSrefplaneZ<<" mm."<<std::endl; 
      }
    }
    if (jUS == -1 || kUS == -1) {
      //std::cout<<"Failed US track by ["<<jUS<<", "<<kUS<<"] and DS track by ["<<jDS<<", "<<kDS<<"]"<<std::endl; 
      trackmatched = false;
    }
  }
  
  else trackmatched = false;
  return trackmatched;  
}

bool MCSAnalysis::TruthMatchUSDS(){

  bool trackmatched = true;
  jUS=-1;
  jDS=-1;
  kUS=-1;
  kDS=-1;

  for (unsigned int jj=0; jj<mcevent->GetVirtualHitsSize();jj++) {
   int station_id = mcevent->GetVirtualHits()->at(jj).GetStationId();
   if (station_id == 47){
      kUS = 1;
      jUS = 1;
   }

   if (station_id == 49){
	  kDS = 1;
	  jDS = 1;
	}
      }
    if (jUS == -1 || kUS == -1) {
      trackmatched = false;
    }
  return trackmatched;  
}
bool MCSAnalysis::PIDSelection(bool isdata=true){

  // Consider the TOF hits first 
  double rawTOF0HitTime = -1., rawTOF1HitTime = -1., rawTOF2HitTime = -1.;
  /*
  if ( tofevent->GetTOFEventSlabHit().GetTOF0SlabHitArraySize() == 2 ) 
    rawTOF0HitTime  = ( tofevent->GetTOFEventSlabHit().GetTOF0SlabHitArray()[0].GetRawTime() +
			tofevent->GetTOFEventSlabHit().GetTOF0SlabHitArray()[1].GetRawTime() ) / 2.;
  else
    return false;
  if ( tofevent->GetTOFEventSlabHit().GetTOF1SlabHitArraySize() == 2 ) 
    rawTOF1HitTime  = ( tofevent->GetTOFEventSlabHit().GetTOF1SlabHitArray()[0].GetRawTime() +
			tofevent->GetTOFEventSlabHit().GetTOF1SlabHitArray()[1].GetRawTime() ) / 2.;
  else
    return false;
  if ( tofevent->GetTOFEventSlabHit().GetTOF2SlabHitArraySize() == 2 ) 
    rawTOF2HitTime  = ( tofevent->GetTOFEventSlabHit().GetTOF2SlabHitArray()[0].GetRawTime() +
			tofevent->GetTOFEventSlabHit().GetTOF2SlabHitArray()[1].GetRawTime() ) / 2.;
  else
    return false;
  */
  if ( tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArraySize() >= 1 )
    rawTOF0HitTime = tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[0].GetTime();
  else
    return false;
  if ( tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArraySize() >= 1 )
    rawTOF1HitTime = tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime();
  else
    return false;
  if (  tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArraySize() >= 1 )
    rawTOF2HitTime = tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray()[0].GetTime();
  else // Allow for TOF2 to not be hit
    rawTOF2HitTime = rawTOF1HitTime + 100.0 * 8.22475 / 0.299792458; // ns. 
  // return false;
  if(isdata){
    tof10->Fill(rawTOF1HitTime - rawTOF0HitTime);
    tof21->Fill(rawTOF2HitTime - rawTOF1HitTime);
    if ( rawTOF1HitTime - rawTOF0HitTime < TOF_lower_limit ||
	 rawTOF1HitTime - rawTOF0HitTime > TOF_upper_limit) return false;
    // if ( rawTOF2HitTime - rawTOF1HitTime < TOF_lower_limit ||
    //      rawTOF2HitTime - rawTOF1HitTime > TOF_upper_limit) return false;
    tof10_sel->Fill(rawTOF1HitTime - rawTOF0HitTime);
    tof21_sel->Fill(rawTOF2HitTime - rawTOF1HitTime);
  }
  else {
    mctof10->Fill(rawTOF1HitTime - rawTOF0HitTime);
    mctof21->Fill(rawTOF2HitTime - rawTOF1HitTime);

    if ( rawTOF1HitTime - rawTOF0HitTime < TOF_lower_limit ||
	 rawTOF1HitTime - rawTOF0HitTime > TOF_upper_limit) return false;
    // if ( rawTOF2HitTime - rawTOF1HitTime < TOF_lower_limit ||
    //      rawTOF2HitTime - rawTOF1HitTime > TOF_upper_limit) return false;
    mctof10_sel->Fill(rawTOF1HitTime - rawTOF0HitTime);
    mctof21_sel->Fill(rawTOF2HitTime - rawTOF1HitTime);
  }
  return true;
  
}

bool MCSAnalysis::TruthTime(bool isdata=true){

  MAUS::TOFHitArray* tofhit = mcevent->GetTOFHits();
  float tof0;
  float tof1;
  for (unsigned int i=0;i<tofhit->size();i++) {
     if (tofhit->at(i).GetPosition().Z() - 5287 < 20 and tofhit->at(i).GetPosition().Z() - 5287 > -20){  
        tof0 = tofhit->at(i).GetTime();   
     }
     if (tofhit->at(i).GetPosition().Z() - 12929 < 20 and tofhit->at(i).GetPosition().Z() - 12929 > -20){        tof1 = tofhit->at(i).GetTime();
     }
  }
  //std::cout << "tof1-tof0 " << tof1-tof0 << std::endl;
  if ((tof1-tof0) > TOF_upper_limit || (tof1-tof0) < TOF_lower_limit) return false; 
  return true;
}

bool MCSAnalysis::RadialSelection(double pz, double pos, double radius){
  bool selected = true;
  if (jUS == -1 || kUS == -1) 
    selected = false;
  else {
    Vars USplane;
    USplane.X = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->pos().x() + _sys["alXUS"];
    USplane.Y = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->pos().y() + _sys["alYUS"];
    USplane.Z = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->pos().z();
    // if ( sqrt(xpos*xpos + ypos*ypos) > meanp) selected = false;
    USplane.dXdz = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->gradient().x() +
      tan(_sys["thXUS"] * atan(1.)/45.0);
    USplane.dYdz = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->gradient().y() +
      tan(_sys["thYUS"] * atan(1.)/45.0);
    USplane.pz   = pz; // scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->mom().z();
    USplane.px   = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->mom().x() +
      tan(_sys["thXUS"] * atan(1.)/45.0) * USplane.pz;
    USplane.py   = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->mom().y() +
      tan(_sys["thXUS"] * atan(1.)/45.0) * USplane.pz;
    //      scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->mom().z();
    // if ( sqrt(dXdz*dXdz + dYdz*dYdz) > sigmap) selected = false;
//    double abspos = _sys["abspos"];
    double phi = atan2(USplane.dYdz, USplane.dXdz);
//    double zdiff = 2*fabs(USplane.Z - pos);
    // double xabs  = dXdz > 0 ? (dXdz + sigmap*cos(phi)) * zdiff + xpos : (dXdz - sigmap*cos(phi)) * zdiff + xpos;
    // double yabs  = dYdz > 0 ? (dYdz + sigmap*sin(phi)) * zdiff + ypos : (dYdz - sigmap*sin(phi)) * zdiff + ypos;
    Vars DSproj;
    if (pos == _sys["diffpos"]) {
	    DSproj = PropagateVarsMu(USplane, pos);
    }
    else {
    DSproj = PropagateVarsMu(USplane, _sys["abspos"]);
    USplane.dXdz += sigmap*cos(phi);
    USplane.dYdz += sigmap*sin(phi);
    USplane.px   += sigmap*cos(phi)*USplane.pz;
    USplane.py   += sigmap*sin(phi)*USplane.pz;
    // double xabs  = (dXdz + sigmap*cos(phi)) * zdiff + xpos;
    // double yabs  = (dYdz + sigmap*sin(phi)) * zdiff + ypos;
    DSproj = PropagateVarsMu(USplane, pos);
    DSproj.X += 5.402;
    DSproj.dXdz = DSproj.dXdz*cos(-0.002384)-sin(-0.002384);
    }
    if ( sqrt(DSproj.X*DSproj.X + DSproj.Y*DSproj.Y) > radius) selected = false;
    if (pos == 19948.8) projradius->Fill(sqrt(DSproj.X*DSproj.X + DSproj.Y*DSproj.Y));
    if (pos == _sys["diffpos"]) diffradius->Fill(sqrt(DSproj.X*DSproj.X + DSproj.Y*DSproj.Y));
  }
  return selected;
}

bool MCSAnalysis::TruthRadialSelection(double pz, double pos, double radius, int j){
  bool selected = true;
  if (jUS == -1 || kUS == -1) 
    selected = false;
  else {
    Vars USplane;
    USplane.X = mcevent->GetVirtualHits()->at(j).GetPosition().x();
    USplane.Y = mcevent->GetVirtualHits()->at(j).GetPosition().y();
    USplane.Z = mcevent->GetVirtualHits()->at(j).GetPosition().z();
    USplane.dXdz = mcevent->GetVirtualHits()->at(j).GetMomentum().x()/mcevent->GetVirtualHits()->at(j).GetMomentum().z();
    USplane.dYdz = mcevent->GetVirtualHits()->at(j).GetMomentum().y()/mcevent->GetVirtualHits()->at(j).GetMomentum().z();
    USplane.pz   = mcevent->GetVirtualHits()->at(j).GetMomentum().z();
    USplane.px   = mcevent->GetVirtualHits()->at(j).GetMomentum().x();
    USplane.py   = mcevent->GetVirtualHits()->at(j).GetMomentum().y();
    double phi = atan2(USplane.dYdz, USplane.dXdz);
//    double zdiff = 2*fabs(USplane.Z - pos);
    USplane.dXdz += sigmap*cos(phi);
    USplane.dYdz += sigmap*sin(phi);
    USplane.px   += sigmap*cos(phi)*USplane.pz;
    USplane.py   += sigmap*sin(phi)*USplane.pz;
    Vars DSproj = PropagateVarsMu(USplane, pos);
    if ( sqrt(DSproj.X*DSproj.X + DSproj.Y*DSproj.Y) > radius) selected = false;
  }
  return selected;
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

std::vector<double> MCSAnalysis::RotDefineProjectionAngles(Vars US, Vars DS,int i){

  std::vector<double> projTheta;
  double USnorm = 1./sqrt(1 + US.dXdz*US.dXdz + US.dYdz*US.dYdz);
  TVector3 u(US.dXdz*USnorm, US.dYdz*USnorm, USnorm);
  //double y = -cos(i*3.14/180);
//  double x = sin(i*3.14/180);
  TVector3 s(0, -1, 0);
  s.RotateZ(i*TMath::Pi()/180);
  double snorm = 1./sqrt(s[2]*s[2] + s[1]*s[1] + s[0]*s[0]);
  TVector3 s_(s[0]*snorm,s[1]*snorm,s[2]*snorm);
  TVector3 v = s.Cross(u);
  double vnorm = 1./sqrt(v[2]*v[2] + v[1]*v[1] + v[0]*v[0]); 
  TVector3 v_(v[0]*vnorm,v[1]*vnorm,v[2]*vnorm);
  TVector3 w = v_.Cross(u);
  double Wnorm  = 1./sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
  w[0] *= Wnorm;
  w[1] *= Wnorm;
  w[2] *= Wnorm;
  double DSnorm = 1./sqrt(1 + DS.dXdz*DS.dXdz + DS.dYdz*DS.dYdz);
  TVector3 d(DS.dXdz*DSnorm, DS.dYdz*DSnorm, DSnorm);
  //projTheta.push_back(atan(d.Dot(w)/d.Mag()*w.Mag()));
  //projTheta.push_back(atan(d.Dot(v_)/v_.Mag()*d.Mag()));
  projTheta.push_back(atan(d.Dot(w)/d.Dot(u)));
  projTheta.push_back(atan(d.Dot(v_)/d.Dot(u)));

  //if (i<60 && i>35){
	  //std::cout << "i " << i << std::endl;
	  //std::cout << "y " << y << std::endl;
	  //std::cout << "x " << x << std::endl;
  //}
	  //std::cout << "d.Dot(w) " << d.Dot(w) << std::endl;
	  //std::cout << "d.Dot(v_) " << d.Dot(v_) << std::endl;

  projTheta.push_back( acos( ( (1 + US.dXdz * DS.dXdz + US.dYdz * DS.dYdz )/
			       sqrt(1 + US.dXdz*US.dXdz + US.dYdz*US.dYdz)/
			       sqrt(1 + DS.dXdz*DS.dXdz + DS.dYdz*DS.dYdz))) );
  projTheta.push_back(d.Dot(w)); 
  projTheta.push_back(d.Dot(v_)); 
  TVector3 yvec(0,-1,0);
  //std::cout << "yvec.Angle(d) " << yvec.Angle(d) << std::endl;
  //std::cout << "w.Angle(d) " << w.Angle(d) << std::endl;
  projTheta.push_back(yvec.Angle(d)); 
  projTheta.push_back(w.Angle(d)); 
  projTheta.push_back(v_.Angle(d)); 
  projTheta.push_back(d[0]);
  projTheta.push_back(d[1]);
  //std::cout << "d[0] " << d[0] << std::endl;
  //std::cout << "d[1] " << d[1] << std::endl;
  return projTheta;
}

double MCSAnalysis::PathLengthInLH2(double pz){
    
    if (jUS != -1 && kUS != -1 && jDS != -1 && kDS != -1){
    float x_up;
    float y_up;
    float z_up;
    float x_down;
    float y_down;
    float z_down;
    float x_up_it;
    float y_up_it;
//    float z_up_it;
    float x_down_it;
    float y_down_it;
//    float z_down_it;
    float rad_up = 0;
    float curvey = 0;
    float curvey_down = 0;
    float rad_down = 0;
    int i = 0;
    
    double poly[5];
    poly[0] = 5.00929921969443e-6;
    poly[1] = -0.0006657044;
    poly[2] = 0.0338057657;
    poly[3] = -0.8148399317;
    poly[4] = 10.8671889437;

    std::vector<double> path_length;
    std::vector<double> vpath_length;

    // Path length for muon between TOF1 and absorber
    Vars USplane;
    USplane.X = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->pos().x() + _sys["alXUS"];
    USplane.Y = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->pos().y() + _sys["alYUS"];
    USplane.Z = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->pos().z();
    USplane.dXdz = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->gradient().x() +
      tan(_sys["thXUS"] * atan(1.)/45.0);
    USplane.dYdz = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->gradient().y() +
      tan(_sys["thYUS"] * atan(1.)/45.0);
    USplane.pz   = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->mom().z();
    USplane.px   = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->mom().x() +
      tan(_sys["thXUS"] * atan(1.)/45.0) * USplane.pz;
    USplane.py   = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->mom().y() +
      tan(_sys["thXUS"] * atan(1.)/45.0) * USplane.pz;
//    Vars USabsoproj = PropagateVarsMu(USplane, _sys["abspos"]);
    Vars USabsfront = PropagateVarsMu(USplane, 16729.03);
    
    //Path length for muon between absorber and TOF2
    //vpath_length.clear();
    Vars DSplane;
    DSplane.X = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->pos().x() + _sys["alXDS"];
    DSplane.Y = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->pos().y() + _sys["alYDS"];
    DSplane.Z = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->pos().z();
    DSplane.dXdz = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->gradient().x() +
      tan(_sys["thXDS"] * atan(1.)/45.0);
    DSplane.dYdz = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->gradient().y() +
      tan(_sys["thYDS"] * atan(1.)/45.0);
    DSplane.pz   = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->mom().z();
    DSplane.px   = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->mom().x() +
      tan(_sys["thXDS"] * atan(1.)/45.0) * DSplane.pz;
    DSplane.py   = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->mom().y() +
      tan(_sys["thXDS"] * atan(1.)/45.0) * DSplane.pz;
//    Vars DSabsoproj = PropagateVarsMu(DSplane, _sys["abspos"]);
    Vars DSabsback = PropagateVarsMu(DSplane, 17179.33);
   
                                x_up = USabsfront.X;
                                y_up = USabsfront.Y;
                                z_up = 16779.18;
    rad_up = sqrt(pow(USabsfront.X,2)+pow(USabsfront.Y,2));
    while (rad_up >= curvey && i < 50) {
                      curvey = (poly[0] * pow(i,5)) + (poly[1] * pow(i,4)) + (poly[2] * pow(i,3)) + (poly[3] * pow(i,2)) + (poly[4] * i);
                      x_up_it = USabsfront.X + (USabsfront.dXdz * (i));
                      y_up_it = USabsfront.Y + (USabsfront.dYdz * (i));
                      rad_up = sqrt(pow(x_up_it,2) + pow(y_up_it,2));
		      if (rad_up < curvey) {
                                x_up = x_up_it;
                                y_up = y_up_it;
                                z_up = i+16779.18;
                      }
                      i++;
              }
             
                                x_down = DSabsback.X;
                                y_down = DSabsback.Y;
                                z_down = 17129.18;
              i = 0;
	      rad_down = sqrt(pow(DSabsback.X,2)+pow(DSabsback.Y,2));
              while (rad_down >= curvey_down && i < 50) {
                      curvey_down = sqrt(pow(-(poly[0] * pow(i,5)) - (poly[1] * pow(i,4)) - (poly[2] * pow(i,3)) - (poly[3] * pow(i,2)) - (poly[4] * i),2));
                      x_down_it = DSabsback.X + (DSabsback.dXdz * (-i));
                      y_down_it = DSabsback.Y + (DSabsback.dYdz * (-i));
                      rad_down = sqrt(pow(x_down_it,2) + pow(y_down_it,2));
                      if (rad_down < curvey_down) {
                                x_down = x_down_it;
                                y_down = y_down_it;
                                z_down = 17129.18-i;
                      }
                      i++;
              }

                float dist_absorb = sqrt(pow(x_down-x_up,2)+pow(y_down-y_up,2)+pow(z_down-z_up,2)); 
		/*
    std::cout << "x_down " << x_down << std::endl;
    std::cout << "x_up " << x_up << std::endl;
    std::cout << "y_down " << y_down << std::endl;
    std::cout << "y_up " << y_up << std::endl;
    std::cout << "z_down " << z_down << std::endl;
    std::cout << "z_up " << z_up << std::endl;
                std::cout << "dist_absorb " << dist_absorb << std::endl;
		*/
    pathlengthabs->Fill(dist_absorb);
    }
}

std::vector<double> MCSAnalysis::CalculatePathLength(double pz){
        
    std::vector<double> path_length;
    std::vector<double> vpath_length;

    // Path length for muon between TOF1 and absorber
    Vars USplane;
    /*
    std::cout << "jUS " << jUS << std::endl;
    std::cout << "kUS " << kUS << std::endl;
    */
    USplane.X = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->pos().x() + _sys["alXUS"];
    USplane.Y = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->pos().y() + _sys["alYUS"];
    USplane.Z = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->pos().z();
    USplane.dXdz = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->gradient().x() +
      tan(_sys["thXUS"] * atan(1.)/45.0);
    USplane.dYdz = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->gradient().y() +
      tan(_sys["thYUS"] * atan(1.)/45.0);
    USplane.pz   = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->mom().z();
    USplane.px   = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->mom().x() +
      tan(_sys["thXUS"] * atan(1.)/45.0) * USplane.pz;
    USplane.py   = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->mom().y() +
      tan(_sys["thXUS"] * atan(1.)/45.0) * USplane.pz;
    Vars USabsoproj = PropagateVarsMu(USplane, _sys["abspos"]);
    Vars USTOF1proj = PropagateVarsMu(USplane, _sys["TOF1_z"]);
    /*
    std::cout << USabsoproj.X << std::endl;
    std::cout << USplane.X << std::endl;
    std::cout << USTOF1proj.X << std::endl;
    */
    vpath_length.push_back(USabsoproj.X - USTOF1proj.X);
    vpath_length.push_back(USabsoproj.Y - USTOF1proj.Y);
    vpath_length.push_back(USabsoproj.Z - USTOF1proj.Z);
    path_length.push_back(sqrt(vpath_length[0]*vpath_length[0] + vpath_length[1]*vpath_length[1] + vpath_length[2]*vpath_length[2]));

    //Path length for muon between absorber and TOF2
    vpath_length.clear();
    Vars DSplane;
    /*
    std::cout << "jDS " << jDS << std::endl;
    std::cout << "kDS " << kDS << std::endl;
    */
    DSplane.X = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->pos().x() + _sys["alXDS"];
    DSplane.Y = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->pos().y() + _sys["alYDS"];
    DSplane.Z = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->pos().z();
    DSplane.dXdz = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->gradient().x() +
      tan(_sys["thXDS"] * atan(1.)/45.0);
    DSplane.dYdz = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->gradient().y() +
      tan(_sys["thYDS"] * atan(1.)/45.0);
    DSplane.pz   = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->mom().z();
    DSplane.px   = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->mom().x() +
      tan(_sys["thXDS"] * atan(1.)/45.0) * DSplane.pz;
    DSplane.py   = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->mom().y() +
      tan(_sys["thXDS"] * atan(1.)/45.0) * DSplane.pz;
    Vars DSabsoproj = PropagateVarsMu(DSplane, _sys["abspos"]);
    Vars DSTOF2proj = PropagateVarsMu(DSplane, _sys["TOF2_z"]);
    vpath_length.push_back(DSTOF2proj.X - DSabsoproj.X);
    vpath_length.push_back(DSTOF2proj.Y - DSabsoproj.Y);
    vpath_length.push_back(DSTOF2proj.Z - DSabsoproj.Z);
    path_length.push_back(sqrt(vpath_length[0]*vpath_length[0] + vpath_length[1]*vpath_length[1] + vpath_length[2]*vpath_length[2]));

    return path_length;

}
std::vector<double> MCSAnalysis::rCalculatePathLength(double pz){
        
    std::vector<double> path_length;
    std::vector<double> vpath_length;

    // Path length for muon between TOF1 and absorber
    Vars TOF0;
    Vars TOF1;
    /*
    std::cout << "jUS " << jUS << std::endl;
    std::cout << "kUS " << kUS << std::endl;
    */    
    if( int(tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray().size()) > 0) {
      TOF0.X = tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[0].GetGlobalPosX();
      TOF0.Y = tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[0].GetGlobalPosY();
      TOF0.Z = tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[0].GetGlobalPosZ();
    }
    if( int(tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()) > 0) {
      TOF1.X = tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetGlobalPosX();
      TOF1.Y = tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetGlobalPosY();
      TOF1.Z = tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetGlobalPosZ();
    }
    /*
    std::cout << USabsoproj.X << std::endl;
    std::cout << USplane.X << std::endl;
    std::cout << USTOF1proj.X << std::endl;
    */
    vpath_length.push_back(TOF1.X - TOF0.X);
    vpath_length.push_back(TOF1.Y - TOF0.Y);
    vpath_length.push_back(TOF1.Z - TOF0.Z);
    path_length.push_back(sqrt(vpath_length[0]*vpath_length[0] + vpath_length[1]*vpath_length[1] + vpath_length[2]*vpath_length[2]));

    //Path length for muon between absorber and TOF2
    vpath_length.clear();

    return path_length;

}
double MCSAnalysis::CorMomFromTOF(double pz, double mat, double diff){

   double t_BB = TimeofFlight();
   double dt0 = (_sys["TOF1_z"] - _sys["TOF0_z"]) / 0.299792458 / 1000.;
   double pzBB = 105.65/sqrt(pow(t_BB,2)/pow(dt0,2)-1);
   double Eini = sqrt(pzBB*pzBB + 105*105);

   double Iair = 85.7e-6;
   double Iscin = 64.7e-6;
   double IAl = 166e-6;
   double IHe = 41.6e-6;
   double ILH2 = 21.8e-6;
   double ILiH = 36.5e-6;
   double ICu = 322e-6;
   double IW = 727e-6;

   double zair = 897.16;
   double zscin = 5.09;
   double zAl = 16e-4;
   double zHe = 11.3;
   double zLH2 = 17.5;
   double zLiH = 3.25;
   double zCu = 0.891;
   double zW = 0.84;
   
   double Zair = 0.49919;
   double Zscin = 0.54141;
   double ZAl = 13;
   double ZHe = 2;
   double ZLH2 = 1;
   double ZLiH = 0.50321;
   double ZCu = 29;
   double ZW = 74;
   
   double Aair = 1;
   double Ascin = 1;
   double AAl = 26.9815385;
   double AHe = 4.002602;
   double ALH2 = 1.008;
   double ALiH = 1;
   double ACu = 64;
   double AW = 184;

   double Rair = 1.205e-3;
   double Rscin = 1.032;
   double RAl = 2.699;
   double RHe = 1.663e-4;
   double RLH2 = 0.0708;
   double RLiH = 0.82;
   double RCu = 8.96;
   double RW = 19.3;
   
   double hwair = 0.71e-6;
   double hwscin = 21.54e-6;
   double hwAl = 32.86e-6;
   double hwHe = 0.26e-6;
   double hwLH2 = 7.64e-6;
   double hwLiH = 18.51e-6;
   double hwCu = 58.27e-6;
   double hwW = 80.32e-6;
   
   double mostprobBBair = MostProbBB(pzBB,Iair,Zair,Aair,hwair,Rair,zair);
   double mostprobBBscin = MostProbBB(pzBB,Iscin,Zscin,Ascin,hwscin,Rscin,zscin);
   double mostprobBBAl = MostProbBB(pzBB,IAl,ZAl,AAl,hwAl,RAl,zAl);
   double mostprobBBHe = MostProbBB(pzBB,IHe,ZHe,AHe,hwHe,RHe,zHe);
   double mostprobBBLH2 = MostProbBB(pzBB,ILH2,ZLH2,ALH2,hwLH2,RLH2,zLH2);
   double mostprobBBLiH = MostProbBB(pzBB,ILiH,ZLiH,ALiH,hwLiH,RLiH,zLiH);
   double mostprobBBCu = MostProbBB(pzBB,ICu,ZCu,ACu,hwCu,RCu,zCu);
   double mostprobBBW = MostProbBB(pzBB,IW,ZW,AW,hwW,RW,zW);
   
   /*
   std::cout << "mostprobBBair " << mostprobBBair << std::endl;
   std::cout << "mostprobBBscin " << mostprobBBscin << std::endl;
   std::cout << "mostprobBBAl " << mostprobBBAl << std::endl;
   std::cout << "mostprobBBHe " << mostprobBBHe << std::endl;
   std::cout << "mostprobBBLH2 " << mostprobBBLH2 << std::endl;
   std::cout << "mostprobBBLiH " << mostprobBBLiH << std::endl;
   std::cout << "mostprobBBCu " << mostprobBBCu << std::endl;
   std::cout << "mostprobBBW " << mostprobBBW << std::endl;
   */

//  double BBair = BetheBloch(pzBB,Iair,Zair,Aair,hwair);
//   double BBscin = BetheBloch(pzBB,Iscin,Zscin,Ascin,hwscin);
//   double BBAl = BetheBloch(pzBB,IAl,ZAl,AAl,hwAl);
//   double BBHe = BetheBloch(pzBB,IHe,ZHe,AHe,hwHe);
//   double BBLH2 = BetheBloch(pzBB,ILH2,ZLH2,ALH2,hwLH2);
   /*
   std::cout << "BBair " << BBair << std::endl;
   std::cout << "BBscin " << BBscin << std::endl;
   std::cout << "BBAl " << BBAl << std::endl;
   std::cout << "BBHe " << BBHe << std::endl;
   std::cout << "BBLH2 " << BBLH2 << std::endl;
   */

//   double ELair = BBair*zair*Rair;
//   double ELscin = BBscin*zscin*Rscin;
//   double ELAl = BBAl*zAl*RAl;
//   double ELHe = BBHe*zHe*RHe;
//   double ELLH2 = BBLH2*zLH2*RLH2;
   /*
   std::cout << "ELair " << ELair << std::endl;
   std::cout << "ELscin " << ELscin << std::endl;
   std::cout << "ELAl " << ELAl << std::endl;
   std::cout << "ELHe " << ELHe << std::endl;
   std::cout << "ELLH2 " << ELLH2 << std::endl;
   */
  
   //std::cout << diff << std::endl;
   double pcor;
   if (material=="LH2" && mat==0) {
	   double totEL = mostprobBBair+mostprobBBscin+mostprobBBAl+mostprobBBHe+mostprobBBLH2;
	   if (diff==1) totEL += mostprobBBCu + mostprobBBW;
	   double Efin = Eini - totEL;
	   pcor = sqrt(pow(Efin,2)-pow(105,2));
           cor_mom->Fill(pcor);
   }
   if (material=="LiH" && mat==0) {
	   double totEL = mostprobBBair+mostprobBBscin+mostprobBBAl+mostprobBBHe+mostprobBBLiH;
	   //if (diff==1) totEL += mostprobBBCu + mostprobBBW;
	   double Efin = Eini - totEL;
	   pcor = sqrt(pow(Efin,2)-pow(105,2));
           cor_mom->Fill(pcor);
   }
   if (mat==1) {
	   double totEL = mostprobBBair+mostprobBBscin+mostprobBBAl+mostprobBBHe;
	   //if (diff==1) totEL += mostprobBBCu + mostprobBBW;
	   double Efin = Eini - totEL;
	   pcor = sqrt(pow(Efin,2)-pow(105,2));
   }

   //std::cout << "pzBB " << pzBB << std::endl;
   //std::cout << "pcor " << pcor << std::endl;


   
   // Initialise and collect initial TOF and pz
   /*
   double t_initial = TimeofFlight();
   double pz_cor = 0;
 
   
   //if (t_initial != 100) {
   t_initial = t_initial*0.299792458;
   //
   //std::vector<double> path_length = CalculatePathLength(pz);
   //double s1 = path_length.at(0)/1000;
   //if (path_length.size() == 2) {
   //double s2 = path_length.at(1)/1000;
   //}
   //

   std::vector<double> rpath_length = rCalculatePathLength(pz);
   double deltas = rpath_length.at(0)/1000;
   double s_1 = 12.929-16.952;
   double s_2 = 21.139-16.952;
   double c = 299792458;
   double pz01 = 105.65*(7.64)/sqrt(pow(t_initial,2)-pow((7.64),2));
   t_initial = t_initial/0.299792458;
   t_initial = t_initial/1000000000;

   // Paul 6th short order TOF01 
   
   double dEdx = 2.5;
   s_1 = 5.287-16.95;
   s_2 = 12.929-16.95;
   double term1 = pow(dEdx,2)*pow((s_1+s_2),2)/(4*pow(105.65,2));
   double term2 = -(dEdx*c*t_initial*(s_1+s_2))/((12.929-5.287)*105);
   double term3 = -1;
   double term4 = pow(c*t_initial/(12.929-5.287),2)-1;
   std::cout << "term1 " << term1 << std::endl;
   std::cout << "term2 " << term2 << std::endl;
   std::cout << "term3 " << term3 << std::endl;
   std::cout << "term4 " << term4 << std::endl;
   double x = 105/pz01;
   My2Function1D myf4;
   myf4.a = term1;
   myf4.b = term2;
   myf4.cp = term3;
   myf4.d = term4;
   ROOT::Math::GradFunctor1D  f4(myf4); 
   ROOT::Math::RootFinder rfn4(ROOT::Math::RootFinder::kGSL_NEWTON);
   rfn4.SetFunction(f4, x);
   if (isfinite(myf4.operator()(x)) && isfinite(myf4.Derivative(x))) {
   rfn4.Solve();
   }
   double p01shortPaul6thfor = 105.5/rfn4.Root();
   
   double pz_an = p01shortPaul6thfor;
   */
/*
   // Paul 6th order TOF12 
  
   double pz12 = 0;
   double p12Paul6thfor = 0;
   double t_ds =TimeofFlight12();
   if (t_ds != 100) {
   t_ds = t_ds*0.299792458;
   pz12 = 105.65*(8.209)/sqrt(pow(t_ds,2)-pow((8.209),2));
   double t_12 = TimeofFlight12()*0.299792458;
   t_12 = t_12/0.299792458;
   t_12 = t_12/1000000000;
   dEdx = 0.5;
   term1 = pow(dEdx,2)*pow((s_1+s_2),2)/(4*pow(105.65,2));
   term2 = -(dEdx*c*t_12*(s_1+s_2))/(8.21*105);
   term3 = -1;
   term4 = pow(c*t_12/(8.21),2)-1;

   std::cout << "term1 " << term1 << std::endl;
   std::cout << "term2 " << term2 << std::endl;
   std::cout << "term3 " << term3 << std::endl;
   std::cout << "term4 " << term4 << std::endl;

   x = 105/pz12;
   My2Function1D myf2;
   myf2.a = term1;
   myf2.b = term2;
   myf2.cp = term3;
   myf2.d = term4;
   //std::cout << "myf2.operator()(pz) " << myf2.operator()(x) << std::endl;
   //std::cout << "myf2.Derivative(pz) " << myf2.Derivative(x) << std::endl;
   ROOT::Math::GradFunctor1D  f2(myf2); 
   ROOT::Math::RootFinder rfn2(ROOT::Math::RootFinder::kGSL_NEWTON);
   rfn2.SetFunction(f2, x);
   if (isfinite(myf2.operator()(x)) && isfinite(myf2.Derivative(x))) {
   rfn2.Solve();
   }
   p12Paul6thfor = 105.5/rfn2.Root();
  
   pz_cor = p12Paul6thfor;
   }
   */
   double MCTruth_pz_mid = 0;
   double MCTruth_pz_up = 0;
   double MCTruth_pz_down = 0;
   double TOF0E;
   double TOF1E;
   double TOF2E;
   for ( size_t j=0; j < mcevent->GetVirtualHits()->size(); j++){
	   energyloss->Fill(mcevent->GetVirtualHits()->at(j).GetPosition().z(),sqrt(pow(mcevent->GetVirtualHits()->at(j).GetMomentum().z(),2)+pow(105,2)));
	   if (mcevent->GetVirtualHits()->at(j).GetPosition().z()-17000<50 && mcevent->GetVirtualHits()->at(j).GetPosition().z()-17000>-50) {
		   MCTruth_pz_mid = mcevent->GetVirtualHits()->at(j).GetMomentum().z();
	   }
	   if (mcevent->GetVirtualHits()->at(j).GetPosition().z()-16803.7<50 && mcevent->GetVirtualHits()->at(j).GetPosition().z()-16803.7>-50) {
		   MCTruth_pz_up = mcevent->GetVirtualHits()->at(j).GetMomentum().z();
	   }
	   if (mcevent->GetVirtualHits()->at(j).GetPosition().z()-17101.3<50 && mcevent->GetVirtualHits()->at(j).GetPosition().z()-17101.3>-50) {
		   MCTruth_pz_down = mcevent->GetVirtualHits()->at(j).GetMomentum().z();
	   }
	   if (mcevent->GetVirtualHits()->at(j).GetPosition().z()-5000<50 && mcevent->GetVirtualHits()->at(j).GetPosition().z()-5000>-50) {
		   TOF0E = sqrt(pow(mcevent->GetVirtualHits()->at(j).GetMomentum().z(),2)+pow(105,2));
	   }
	   if (mcevent->GetVirtualHits()->at(j).GetPosition().z()-12500<50 && mcevent->GetVirtualHits()->at(j).GetPosition().z()-12500>-50) {
		   TOF1E = sqrt(pow(mcevent->GetVirtualHits()->at(j).GetMomentum().z(),2)+pow(105,2));
	   }
	   if (mcevent->GetVirtualHits()->at(j).GetPosition().z()-21039<50 && mcevent->GetVirtualHits()->at(j).GetPosition().z()-21039>-50) {
		   TOF2E = sqrt(pow(mcevent->GetVirtualHits()->at(j).GetMomentum().z(),2)+pow(105,2));
	   }

   }
   TOF0Energy->Fill(TOF0E);
   TOF1Energy->Fill(TOF1E);
   TOF2Energy->Fill(TOF2E);
   if (MCTruth_pz_up != 0 && MCTruth_pz_down != 0){
//   double true_delta = MCTruth_pz_up - MCTruth_pz_down;	   
//   double MCTruth_pz = (MCTruth_pz_up+MCTruth_pz_down)/2;
   double res01 = pzBB - MCTruth_pz_mid;
   //double res12 = pz12 - MCTruth_pz_mid;
   //double res12p6 = p12Paul6thfor - MCTruth_pz_mid;
   double res01short = pcor - MCTruth_pz_mid;
   residual01->Fill(res01);
   //residual12->Fill(res12);
   //residual12p6->Fill(res12p6);
   residual01short->Fill(res01short);

   //TOFcom->Fill(pz_an,pcor);
   TOF01vsMCTruth->Fill(pzBB,MCTruth_pz_mid);
   //TOF12vsMCTruth->Fill(pz12,MCTruth_pz_mid);
   //TOF12Paul6thforvsMCTruth->Fill(p12Paul6thfor,MCTruth_pz_mid);
   
   TOF01shortPaul6thforvsMCTruth->Fill(pcor,MCTruth_pz_mid);
   //std::cout << "MCTruth_pz_mid " << MCTruth_pz_mid << std::endl;
   	}

   //}
   double t_ds =TimeofFlight12();
   if (t_ds != 100) {
   t_ds = t_ds*0.299792458;
   pcor = 105.65*(8.209)/sqrt(pow(t_ds,2)-pow((8.209),2));
   }

   //std::cout << "pcor 12 " << pcor << std::endl;
   PathLengthInLH2(pcor);
   return pcor;
}

double MCSAnalysis::TimeofFlight(){
  double rawTOF1HitTime = -1., rawTOF0HitTime = -1.;
  if( int(tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()) > 0)
    rawTOF1HitTime  = 
      tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime();
  if( int(tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray().size()) > 0)
    rawTOF0HitTime = 
      tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[0].GetTime();
  double dt  = 100.0; // something rediculously large as a default number
  if ( rawTOF1HitTime != -1 && rawTOF0HitTime != -1 ){
    dt  = rawTOF1HitTime - rawTOF0HitTime; 
  }
  return dt;
}

double MCSAnalysis::TimeofFlight12(){
  double rawTOF1HitTime = -1., rawTOF2HitTime = -1.;
  if( int(tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()) > 0)
    rawTOF1HitTime  = 
      tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime();
  if( int(tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray().size()) > 0)
    rawTOF2HitTime = 
      tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray()[0].GetTime();
  double dt  = 100.0; // something rediculously large as a default number
  if ( rawTOF1HitTime != -1. && rawTOF2HitTime != -1. ){
    dt  = rawTOF2HitTime - rawTOF1HitTime; 
  }
  return dt;
}

double MCSAnalysis::BetheBloch(double pz, double Imat, double Z, double A, double hw){

   double beta = pow(pz,2)/(pow(105.65,2)+pow(pz,2));
   double gamma = 1/sqrt(1-pow(beta,2));
   double W = 2*0.511*pow(beta,2)*pow(gamma,2)/(1+2*gamma*0.511/105.65+pow(0.511/105.65,2));
   double I = Imat;
   double density = log(hw/I)+log(beta*gamma)-1/2;
   double dEdxpre = 0.307075*Z/(pow(beta,2)*A);   
   double dEdxterm1 = 0.5*log(2*0.511*pow(beta,2)*pow(gamma,2)*W/pow(I,2)) - pow(beta,2) - density;
   double dEdx = dEdxpre * dEdxterm1;
   /*
   std::cout << "beta " << beta << std::endl;
   std::cout << "gamma " << gamma << std::endl;
   std::cout << "W " << W << std::endl;
   std::cout << "density " << density << std::endl;
   std::cout << "dEdxpre " << dEdxpre << std::endl;
   std::cout << "dEdxterm1 " << dEdxterm1 << std::endl;
   std::cout << "dEdx " << dEdx << std::endl;
   */
   return dEdx;
}

double MCSAnalysis::MostProbBB(double pz, double Imat, double Z, double A, double hw, double R, double z){

   double beta = pow(pz,2)/(pow(105.65,2)+pow(pz,2));
   double gamma = 1/sqrt(1-pow(beta,2));
   double I = Imat;
   double density = log(hw/I)+log(beta*gamma)-1/2;
   double E = 0.307075*Z*R*z/(2*pow(beta,2)*A);   
   double mostprobBB = E*(log(2*0.511*pow(beta,2)*pow(gamma,2)/I)+log(E/I)+0.2-pow(beta,2)-density);

   return mostprobBB;
}

TH1D* MCSAnalysis::trkreffix(TH1D* h1){

	TFile *f=new TFile(trkreffiname.c_str());
	TGraphAsymmErrors* efficiency_scat_x = (TGraphAsymmErrors*)f->Get("Effx_graph");
	double x;
	double y;
//	TCanvas* c1 = new TCanvas();
	for (int i=1; i<47; i++) {
		efficiency_scat_x->GetPoint(i,x,y);
		double ey = efficiency_scat_x->GetErrorYlow(i);
		if  (y!=0){
			h1->Sumw2();
			h1->SetBinContent(i,h1->GetBinContent(i)/y);
			h1->SetBinError(i,sqrt(h1->GetBinError(i)*h1->GetBinError(i)/h1->GetBinContent(i)*h1->GetBinContent(i)+ey*ey/y*y));
		}
		else {
			h1->Sumw2();
			h1->SetBinContent(i,h1->GetBinContent(i));
		}
	}
	f->Close();
        return h1;
}

TH1D* MCSAnalysis::trkreffiy(TH1D* h1){

	TFile *f=new TFile(trkreffiname.c_str());
	TGraphAsymmErrors* efficiency_scat_y = (TGraphAsymmErrors*)f->Get("Effy_graph");
	double x;
	double y;
//	TCanvas* c1 = new TCanvas();
	for (int i=1; i<47; i++) {
		efficiency_scat_y->GetPoint(i,x,y);
		double ey = efficiency_scat_y->GetErrorYlow(i);
		if  (y!=0){
			h1->Sumw2();
			h1->SetBinContent(i,h1->GetBinContent(i)/y);
			h1->SetBinError(i,sqrt(h1->GetBinError(i)*h1->GetBinError(i)/h1->GetBinContent(i)*h1->GetBinContent(i)+ey*ey/y*y));
		}
		else {
			h1->Sumw2();
			h1->SetBinContent(i,h1->GetBinContent(i));
		}
	}
	f->Close();
	
        return h1;
}

TH1D* MCSAnalysis::trkreffiscatt(TH1D* h1){

	TFile *f=new TFile(trkreffiname.c_str());
	TGraphAsymmErrors* efficiency_scat_scatt = (TGraphAsymmErrors*)f->Get("Effscatt_graph");
	double x;
	double y;
//	TCanvas* c1 = new TCanvas();
	for (int i=1; i<47; i++) {
		efficiency_scat_scatt->GetPoint(i,x,y);
		double ey = efficiency_scat_scatt->GetErrorYlow(i);
		if (y!=0){
			h1->Sumw2();
			h1->SetBinContent(i,h1->GetBinContent(i)/y);
			h1->SetBinError(i,sqrt(h1->GetBinError(i)*h1->GetBinError(i)/h1->GetBinContent(i)*h1->GetBinContent(i)+ey*ey/y*y));
		}
		else {
			h1->Sumw2();
			h1->SetBinContent(i,h1->GetBinContent(i));
		}
	}
	f->Close();
	
        return h1;
}

TH1D* MCSAnalysis::trkreffi2scatt(TH1D* h1){

	TFile *f=new TFile(trkreffiname.c_str());
	TGraphAsymmErrors* efficiency_scat_2scatt = (TGraphAsymmErrors*)f->Get("Eff2scatt_graph");
	double x;
	double y;
//	TCanvas* c1 = new TCanvas();
	for (int i=1; i<47; i++) {
		efficiency_scat_2scatt->GetPoint(i,x,y);
		double ey = efficiency_scat_2scatt->GetErrorYlow(i);
		if (y!=0){
			h1->Sumw2();
			h1->SetBinContent(i,h1->GetBinContent(i)/y);
			h1->SetBinError(i,sqrt(h1->GetBinError(i)*h1->GetBinError(i)/h1->GetBinContent(i)*h1->GetBinContent(i)+ey*ey/y*y));
		}
		else {
			h1->Sumw2();
			h1->SetBinContent(i,h1->GetBinContent(i));
		}
	}
	f->Close();
	
        return h1;
}

double MCSAnalysis::MomentumFromTOF(bool isdata=true){
    // Cuts remove events where the following statements do not make sense so we proceed without cuts.
  double rawTOF0HitTime = -1., rawTOF1HitTime = -1., rawTOF2HitTime = -1.;
  if( int(tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray().size()) > 0)
    rawTOF0HitTime  = 
      tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[0].GetTime();
  if( int(tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()) > 0)
    rawTOF1HitTime  = 
      tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime();
  if( int(tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray().size()) > 0)
    rawTOF2HitTime = 
      tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray()[0].GetTime();
  // tofevent->GetTOFEventSlabHit().GetTOF2SlabHitArray()[1].GetRawTime() ) / 2.;
  // need to hard code a few things here (unfortunately) pertaining to the geometry and the time of flight.
  double dt0 = (_sys["TOF1_z"] - _sys["TOF0_z"]) / 0.299792458 / 1000.; // ns. 
  double dt  = 10.0 * dt0; // something rediculously large as a default number
  double pz  = 105.65 / sqrt(dt*dt/dt0/dt0 - 1.0);
  if ( rawTOF1HitTime != -1 && rawTOF0HitTime != -1 ){
    dt  = rawTOF1HitTime - rawTOF0HitTime; 
    pz  = 105.65 / sqrt(dt*dt/dt0/dt0 - 1.0) - 36.1;
  }
  if ( rawTOF1HitTime != -1 && rawTOF2HitTime != -1 ){
    // Better estimate of the longitudinal momentum
    dt0 = (_sys["TOF2_z"] - _sys["TOF1_z"]) / 0.299792458 / 1000.; // ns. 
    dt  = rawTOF2HitTime - rawTOF1HitTime; 
//    double pz1 = pz;
    pz  = 105.65 / sqrt(dt*dt/dt0/dt0 - 1.0);
    //TOFcom->Fill(pz1,pz);
  }
  if(isdata)
    calc_mom->Fill(pz);
  else
    mccalc_mom->Fill(pz);
  return pz;
}
void MCSAnalysis::ConvolveWithInputDistribution(std::string distname){
  int failcount = 0;
  int passcount = 0;
//  int isfirst = 0;
  bool isGEANT;
  bool isCobb;
  bool isMoliere;
  if (distname.find(modelname1.c_str()) != std::string::npos)
    isGEANT = true;
  if (distname.find(modelname2.c_str()) != std::string::npos)
    isCobb = true;
  if (distname.find(modelname3.c_str()) != std::string::npos)
    isMoliere = true;

  std::cerr<<"Model file: "<<modelfile<<std::endl;
  TFile* infile = new TFile(modelfile.c_str());
  
  TH1D* hiswX = new TH1D("hiswX","", 1000, -5, 5);
  TH1D* hiswY = new TH1D("hiswY","", 1000, -5, 5);

  // Efficiency plots
//  TH1D* scatx = 
//    new TH1D("scatx","Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",
//	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

//  TH1D* scaty = 
//    new TH1D("scaty","Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",
//	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

//  TH1D* scatscat = 
//    new TH1D("scatscat","Scattering Angle between Momentum Vectors;#theta_{Scatt}; Events per mrad",
//	     _histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);

//  TH1D* scat2scatt = 
//    new TH1D("scat2scatt","Scattering Angle between Momentum Vectors;#theta^{2}_{Scatt}; Events per mrad",
// _histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);
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
      std::vector<double> nominalTheta = RotDefineProjectionAngles(_USMCset.E(i), _DSMCset.E(i), angdef);
//      double nomthetaX = nominalTheta[0];
//      double nomthetaY = nominalTheta[1];
//      double nomthetascat = nominalTheta[2];
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
      std::vector<double> projDTheta = RotDefineProjectionAngles(tmpvar, _DSMCset.E(i), angdef);
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
      
      if( sqrt(projvar.X*projvar.X + projvar.Y*projvar.Y) > fiducial_cut_radius ){
        failcount++;
//        std::cerr<<"ConvolveWithInputDistribution - event failed projvar check. Passed: "<< passcount<<" Failed: "<<failcount<<std::endl;
	tmpvar.X = 0.0;
	tmpvar.Y = 0.0;
	tmpvar.Z = 0.0;
	tmpvar.dXdz = 1./2.;
	tmpvar.dYdz = 1./2.;
      } else {passcount++;}
      std::vector<double> projTheta = RotDefineProjectionAngles(_USMCset.E(i), tmpvar, angdef);
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
    
      // Apply efficiency correction via weighting 
      /*
      double binx = scatx->Fill(thetaX);
      //std::cout << "binx " << binx << std::endl;
      TFile *f=new TFile(trkreffiname.c_str());
      TGraphAsymmErrors* efficiency_scat_x = (TGraphAsymmErrors*)f->Get("Effx_graph");
      double x;
      double y;
      efficiency_scat_x->GetPoint(binx,x,y);
      double effipred = y; 
      binx = scatx->Fill(nomthetaX);
      TFile *fempty=new TFile(trkreffiemptyname.c_str());
      TGraphAsymmErrors* efficiency_empty_scat_x = (TGraphAsymmErrors*)fempty->Get("Effx_graph");
      double a;
      double b;
      efficiency_empty_scat_x->GetPoint(binx,a,b);
      double effiempdata = b;
      if (effiempdata==0) effiempdata = 0.1;
      //std::cout << "effipred " << effipred << std::endl;
      //std::cout << "binx " << binx << std::endl;
      //std::cout << "effiempdata " << effiempdata << std::endl;
      double weightX = effipred/effiempdata; 
      //std::cout << "weightX " << weightX << std::endl;
      double biny = scaty->Fill(thetaY);
      //std::cout << "biny " << biny << std::endl;
      TGraphAsymmErrors* efficiency_scat_y = (TGraphAsymmErrors*)f->Get("Effy_graph");
      double xy;
      double yy;
      efficiency_scat_y->GetPoint(biny,xy,yy);
      double effipredy = yy; 
      biny = scaty->Fill(nomthetaY);
      TGraphAsymmErrors* efficiency_empty_scat_y = (TGraphAsymmErrors*)fempty->Get("Effy_graph");
      double ay;
      double by;
      efficiency_empty_scat_y->GetPoint(biny,ay,by);
      double effiempdatay = by;
      if (effiempdatay==0) effiempdatay = 0.1;
      double weightY = effipredy/effiempdatay; 
      //std::cout << "effipredy " << effipredy << std::endl;
      //std::cout << "biny " << biny << std::endl;
      //std::cout << "effiempdatay " << effiempdatay << std::endl;
      //std::cout << "weightY " << weightY << std::endl;
      double binscatt = scatscat->Fill(thetaScatt);
      TGraphAsymmErrors* efficiency_scat_scatt = (TGraphAsymmErrors*)f->Get("Effscatt_graph");
      double xscatt;
      double yscatt;
      efficiency_scat_scatt->GetPoint(binscatt,xscatt,yscatt);
      double effipredscatt = yscatt; 
      binscatt = scatscat->Fill(nomthetascat);
      TGraphAsymmErrors* efficiency_empty_scat_scatt = (TGraphAsymmErrors*)fempty->Get("Effscatt_graph");
      double ascatt;
      double bscatt;
      efficiency_empty_scat_scatt->GetPoint(binscatt,ascatt,bscatt);
      double effiempdatascatt = bscatt;
      if (effiempdatascatt==0) effiempdatascatt = 0.1;
      double weightscatt = effipredscatt/effiempdatascatt; 
      //std::cout << "binscatt " << binscatt << std::endl;
      //std::cout << "effipredscatt " << effipredscatt << std::endl;
      //std::cout << "effiempdatascatt " << effiempdatascatt << std::endl;
      //std::cout << "weightscatt " << weightscatt << std::endl;
      double bin2scatt = scat2scatt->Fill(thetaScatt*thetaScatt);
      TGraphAsymmErrors* efficiency_scat_2scatt = (TGraphAsymmErrors*)f->Get("Eff2scatt_graph");
      double x2scatt;
      double y2scatt;
      efficiency_scat_2scatt->GetPoint(bin2scatt,x2scatt,y2scatt);
      double effipred2scatt = y2scatt; 
      bin2scatt = scat2scatt->Fill(nomthetascat*nomthetascat);
      TGraphAsymmErrors* efficiency_empty_scat_2scatt = (TGraphAsymmErrors*)fempty->Get("Eff2scatt_graph");
      double a2scatt;
      double b2scatt;
      efficiency_empty_scat_2scatt->GetPoint(bin2scatt,a2scatt,b2scatt);
      double effiempdata2scatt = b2scatt;
      if (effiempdata2scatt==0) effiempdata2scatt = 0.1;
      if (effipred2scatt==0) effipred2scatt = 1;
      double weight2scatt = effipred2scatt/effiempdata2scatt; 
      f->Close();
      fempty->Close();
      //std::cout << "bin2scatt " << bin2scatt << std::endl;
      //std::cout << "effipred2scatt " << effipred2scatt << std::endl;
      //std::cout << "effiempdata2scatt " << effiempdata2scatt << std::endl;
      //std::cout << "weight2scatt " << weight2scatt << std::endl;
      //std::cout << "weightX " << weightX << std::endl;
      //std::cout << "weightY " << weightY << std::endl;
      //std::cout << "weightscatt " << weightscatt << std::endl;
      //std::cout << "weight2scatt " << weight2scatt << std::endl;
     
      hiswX->Fill(weightX);
      hiswY->Fill(weightY);
      if (isGEANT) {
    	  resp_thetaX.Fill(thetaX, d_thetaX, weightX);
    	  resp_thetaY.Fill(thetaY, d_thetaY, weightY);
    	  resp_thetaScatt.Fill(thetaScatt, dthetaScatt, weightscatt);
    	  resp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt, weight2scatt);
      }
      if (isCobb) {
    	  tresp_thetaX.Fill(thetaX, d_thetaX, weightX);
    	  tresp_thetaY.Fill(thetaY, d_thetaY, weightY);
    	  tresp_thetaScatt.Fill(thetaScatt, dthetaScatt, weightscatt);
    	  tresp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt, weight2scatt);
      }
      if (isMoliere) {
    	  mresp_thetaX.Fill(thetaX, d_thetaX, weightX);
    	  mresp_thetaY.Fill(thetaY, d_thetaY, weightY);
    	  mresp_thetaScatt.Fill(thetaScatt, dthetaScatt, weightscatt);
    	  mresp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt, weight2scatt);
      }
      */
      if (isGEANT) {
    	  resp_thetaX.Fill(thetaX, d_thetaX);
    	  resp_thetaY.Fill(thetaY, d_thetaY);
    	  resp_thetaScatt.Fill(thetaScatt, dthetaScatt);
    	  resp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt);
      }
      if (isCobb) {
    	  tresp_thetaX.Fill(thetaX, d_thetaX);
    	  tresp_thetaY.Fill(thetaY, d_thetaY);
    	  tresp_thetaScatt.Fill(thetaScatt, dthetaScatt);
    	  tresp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt);
      }
      if (isMoliere) {
    	  mresp_thetaX.Fill(thetaX, d_thetaX);
    	  mresp_thetaY.Fill(thetaY, d_thetaY);
    	  mresp_thetaScatt.Fill(thetaScatt, dthetaScatt);
    	  mresp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt);
      }
       
      /*
      isfirst == 1 ? tresp_thetaX.Fill(thetaX, d_thetaX) : resp_thetaX.Fill(thetaX, d_thetaX);
      std::cout << "thetaX " << thetaX << std::endl;
      isfirst == 1 ? tresp_thetaY.Fill(thetaY, d_thetaY) : resp_thetaY.Fill(thetaY, d_thetaY);
      std::cout << "thetaY " << thetaY << std::endl;
      isfirst == 1 ? tresp_thetaScatt.Fill(thetaScatt, dthetaScatt) : resp_thetaScatt.Fill(thetaScatt, dthetaScatt); 
      std::cout << "thetaScatt " << thetaScatt << std::endl;
      isfirst == 1 ? tresp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt) : resp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt); 
      */
    }
  }

  TCanvas *c1 = new TCanvas();
  hiswX->Draw();
  c1->SaveAs("hiswX.pdf");
  c1->Clear();
  hiswY->Draw();
  c1->SaveAs("hiswY.pdf");
  delete c1;
  outfile->cd();
  thetaXUS_thetaXDS->Write();
  thetaYUS_thetaYDS->Write();
  thetaX_refconv->Write();
  thetaY_refconv->Write();
  thetaScatt_refconv->Write();
  theta2Scatt_refconv->Write();
  thetaScatt_refconv_vp->Write();

}

void MCSAnalysis::ConvolveWithVirtualInputDistribution(std::string distname){
//  int isfirst = 0;
  bool isGEANT;
  bool isCobb;
  bool isMoliere;
  if (distname.find(modelname1.c_str()) != std::string::npos)
    isGEANT = true;
  if (distname.find(modelname2.c_str()) != std::string::npos)
    isCobb = true;
  if (distname.find(modelname3.c_str()) != std::string::npos)
    isMoliere = true;

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

  std::cout<<"Convolution with "<<tmpname<<" from "<<modelfile<<std::endl;
  TH1D* hx = (TH1D*)infile->Get(tmpname.c_str());
  TH1D* hy = (TH1D*)infile->Get(tmpname.c_str());

  // Collection DSConvSet;
  // for (int l=0; l<10; l++){
  for (int i=0; i<USTruthSet.N(); i++){
      std::vector<double> nominalTheta = RotDefineProjectionAngles(USTruthSet.E(i), DSTruthSet.E(i), angdef);
//      double nomthetaX = nominalTheta[0];
//      double nomthetaY = nominalTheta[1];
    for (int j=0; j<20; j++){
      double dthetaX = hx->GetRandom() * _sys["resX"];
      double dthetaY = hy->GetRandom() * _sys["resY"];
      // First project the upstream track to the absorber 
      double zabspos = _sys["abspos"] + 0.0;
      Vars projvarAbs = PropagateVarsMu(USTruthSet.E(i), zabspos);
      double xabs = projvarAbs.X;  /// _USMCset.E(i).X + _USMCset.E(i).dXdz * dzabsUS;
      double yabs = projvarAbs.Y;  /// _USMCset.E(i).Y + _USMCset.E(i).dYdz * dzabsUS;
      // Now add the angle from the model to the downstream measurement.
      double dXdz_abs = DSTruthSet.E(i).dXdz + tan(dthetaY);
      double dYdz_abs = DSTruthSet.E(i).dYdz + tan(dthetaX);
      // double d_thetaY    = atan(dXdz_abs) - atan(_DSMCset.E(i).dXdz);
      // double d_thetaX    = atan(dYdz_abs) - atan(_DSMCset.E(i).dYdz);
      
      // Project the track into the downstream reference plane
      // double xref = xabs + dXdz_abs * dzabsDS;
      // double yref = yabs + dYdz_abs * dzabsDS;
      Vars tmpvar = DSTruthSet.E(i);
      tmpvar.X = xabs;
      tmpvar.Y = yabs;
      tmpvar.Z = _sys["abspos"];
      tmpvar.dXdz = dXdz_abs;
      tmpvar.dYdz = dYdz_abs;
      tmpvar.px   = dXdz_abs * USTruthSet.E(i).pz;
      tmpvar.py   = dYdz_abs * USTruthSet.E(i).pz;
      tmpvar.pz   = USTruthSet.E(i).pz;
      tmpvar.TOF12= USTruthSet.E(i).TOF12;
      tmpvar.TOF01= USTruthSet.E(i).TOF01;
      std::vector<double> projDTheta = RotDefineProjectionAngles(tmpvar, DSTruthSet.E(i), angdef);
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
	tmpvar.X = 0.0;
	tmpvar.Y = 0.0;
	tmpvar.Z = 0.0;
	tmpvar.dXdz = 1./2.;
	tmpvar.dYdz = 1./2.;
      }
      std::vector<double> projTheta = RotDefineProjectionAngles(USTruthSet.E(i), tmpvar, angdef);
      double thetaY = projTheta[1];   /// atan(tmpvar.dXdz) - atan(_USMCset.E(i).dXdz);
      double thetaX = projTheta[0];   /// atan(tmpvar.dYdz) - atan(_USMCset.E(i).dYdz);
      // double cosScatt = ( (1 + _USMCset.E(i).dXdz * tmpvar.dXdz +
      //		   _USMCset.E(i).dYdz * tmpvar.dYdz )/
      //		  sqrt(1 + _USMCset.E(i).dXdz*_USMCset.E(i).dXdz +
      //		       _USMCset.E(i).dYdz*_USMCset.E(i).dYdz)/
      //		  sqrt(1 + tmpvar.dXdz*tmpvar.dXdz +
      //		       tmpvar.dYdz*tmpvar.dYdz));
      double thetaScatt = projTheta[2];  /// acos(cosScatt);
      thetaXUS_thetaXDS->Fill(atan(USTruthSet.E(i).dXdz), atan(tmpvar.dXdz));
      thetaYUS_thetaYDS->Fill(atan(USTruthSet.E(i).dYdz), atan(tmpvar.dYdz));
      thetaX_refconv->Fill(thetaX);
      thetaY_refconv->Fill(thetaY);
      thetaScatt_refconv->Fill(thetaScatt);
      theta2Scatt_refconv->Fill(thetaScatt*thetaScatt);
      thetaScatt_refconv_vp->Fill(USTruthSet.E(i).pz, thetaScatt);
    
      if (isGEANT) {
    	  resp_thetaX.Fill(thetaX, d_thetaX);
    	  resp_thetaY.Fill(thetaY, d_thetaY);
    	  resp_thetaScatt.Fill(thetaScatt, dthetaScatt);
    	  resp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt);
      }
      if (isCobb) {
    	  tresp_thetaX.Fill(thetaX, d_thetaX);
    	  tresp_thetaY.Fill(thetaY, d_thetaY);
    	  tresp_thetaScatt.Fill(thetaScatt, dthetaScatt);
    	  tresp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt);
      }
      if (isMoliere) {
    	  mresp_thetaX.Fill(thetaX, d_thetaX);
    	  mresp_thetaY.Fill(thetaY, d_thetaY);
    	  mresp_thetaScatt.Fill(thetaScatt, dthetaScatt);
    	  mresp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt);
      }
      /*
      isfirst == 1 ? tresp_thetaX.Fill(thetaX, d_thetaX) : resp_thetaX.Fill(thetaX, d_thetaX);
      std::cout << "thetaX " << thetaX << std::endl;
      isfirst == 1 ? tresp_thetaY.Fill(thetaY, d_thetaY) : resp_thetaY.Fill(thetaY, d_thetaY);
      std::cout << "thetaY " << thetaY << std::endl;
      isfirst == 1 ? tresp_thetaScatt.Fill(thetaScatt, dthetaScatt) : resp_thetaScatt.Fill(thetaScatt, dthetaScatt); 
      std::cout << "thetaScatt " << thetaScatt << std::endl;
      isfirst == 1 ? tresp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt) : resp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt); 
      */
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
    std::vector<double> projTheta = RotDefineProjectionAngles(_USMCset.E(i), _DSMCset.E(i), angdef);
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

  /*
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
  */
}


void MCSAnalysis::DoDeconvolution(std::string model, int n_sel=1){
  // The basic methods associated with the RooUnfolding package
  // First generate a histogram of the measured data.
  
//  int isfirst = 0;
//  if (model.find(modelname1.c_str()) != std::string::npos)
//    isfirst = 1;

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
  
  //1D Histogram Settings
  TH1D *defineHist2(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup);
  TH1D *Histo[180];
  TH1D *Histoy[180];
  TH1D *Histodotw[180];
  TH1D *Histodotv[180];
  TH1D *Histoyd[180];
  TH1D *Histowd[180];
  TH1D *Histowx[180];
  TH1D *Histowy[180];

  for (int l=0;l<180;l++){
      std::string  tmpname = "thetaX_measdata_";
      tmpname += l;
      Histo[l] = defineHist2(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
      tmpname = "thetaY_measdata_";
      tmpname += l;
      Histoy[l] = defineHist2(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
      tmpname = "dotw_";
      tmpname += l;
      Histodotw[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),200, -1, 1);
      tmpname = "dotv_";
      tmpname += l;
      Histodotv[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),200, -1, 1);
      tmpname = "yd_";
      tmpname += l;
      Histoyd[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),100, 1, 3);
      tmpname = "wd_";
      tmpname += l;
      Histowd[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),100, 1, 3);
      tmpname = "wx_";
      tmpname += l;
      Histowx[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),100, -2, 2);
      tmpname = "wy_";
      tmpname += l;
      Histowy[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),100, -2, 2);
 
  }
  
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
    std::vector<double> projTheta = RotDefineProjectionAngles(_USset.E(i), _DSset.E(i), angdef);
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
    for (int l=0;l<180;l++){
     if (_DSset.E(i).dXdz<0.2 && _DSset.E(i).dYdz<0.2) {
	    std::vector<double> RotprojTheta = RotDefineProjectionAngles(_USset.E(i), _DSset.E(i),l);
            double RotthetaY = RotprojTheta[1];
            double RotthetaX = RotprojTheta[0];
	    //std::cout << "RotprojTheta[1] " << RotprojTheta[1] << std::endl;
	    Histo[l]->Fill(RotthetaX);
            Histoy[l]->Fill(RotthetaY);
            //if (RotprojTheta[3]>0.2) Histodotw[l]->Fill(RotprojTheta[3]); 
            //if (RotprojTheta[3]<-0.2) Histodotw[l]->Fill(RotprojTheta[3]); 
            //if (RotprojTheta[4]>0.2) Histodotv[l]->Fill(RotprojTheta[4]); 
            //if (RotprojTheta[4]<-0.2) Histodotv[l]->Fill(RotprojTheta[4]); 
            Histodotw[l]->Fill(RotprojTheta[3]); 
            Histodotw[l]->Fill(RotprojTheta[3]); 
            Histodotv[l]->Fill(RotprojTheta[4]); 
            Histodotv[l]->Fill(RotprojTheta[4]); 
            Histoyd[l]->Fill(RotprojTheta[5]); 
            Histowd[l]->Fill(RotprojTheta[6]); 
            Histowx[l]->Fill(RotprojTheta[8]); 
            Histowy[l]->Fill(RotprojTheta[9]); 
            //Histov_x[l]->Fill(RotprojTheta[8[0]]); 
            //Histov_y[l]->Fill(RotprojTheta[8[0]]); 
  }
    }
  }
    
  if (model.find(modelname1.c_str()) != std::string::npos) {
	  fresp_thetaX = resp_thetaX;
	  fresp_thetaY = resp_thetaY;
	  fresp_thetaScatt = resp_thetaScatt;
	  fresp_theta2Scatt = resp_theta2Scatt;
  }
  if (model.find(modelname2.c_str()) != std::string::npos) {
	  fresp_thetaX = tresp_thetaX;
	  fresp_thetaY = tresp_thetaY;
	  fresp_thetaScatt = tresp_thetaScatt;
	  fresp_theta2Scatt = tresp_theta2Scatt;
  }
  if (model.find(modelname3.c_str()) != std::string::npos) {
	  fresp_thetaX = mresp_thetaX;
	  fresp_thetaY = mresp_thetaY;
	  fresp_thetaScatt = mresp_thetaScatt;
	  fresp_theta2Scatt = mresp_theta2Scatt;
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
  TCanvas* c2 = new TCanvas("c2","A Simple Graph Example",200,10,1400,1000);
  TGraph* convergence = new TGraph(100,xaxis,chi2);
  convergence->GetXaxis()->SetTitle("iterations");
  convergence->GetYaxis()->SetTitle("#chi^{2} MC reco & Truth");
  convergence->Draw();
  c2->SaveAs("convergence.pdf");
  c2->SaveAs("convergence.root");
  c2->Clear();
  TGraph* convergenceeffi = new TGraph(100,xaxis,chi2effi);
  convergenceeffi->GetXaxis()->SetTitle("iterations");
  convergenceeffi->GetYaxis()->SetTitle("#chi^{2} MC reco & Truth");
  convergenceeffi->Draw();
  c2->SaveAs("convergenceeffi.pdf");
  c2->SaveAs("convergenceeffi.root");
  delete c2;
  delete convergence;
  }
  */
  /*
  thetaX_measdata = trkreffix(thetaX_measdata);
  thetaY_measdata = trkreffiy(thetaY_measdata);
  thetaScatt_measdata = trkreffiscatt(thetaScatt_measdata);
  theta2Scatt_measdata = trkreffi2scatt(theta2Scatt_measdata);
  */

  RooUnfoldBayes unfold_thetaX(&fresp_thetaX, thetaX_measdata, int(_sys["niter"]));
  unfold_thetaX.Print();
  RooUnfoldBayes unfold_thetaY(&fresp_thetaY, thetaY_measdata, int(_sys["niter"]));
  unfold_thetaY.Print();
  RooUnfoldBayes 
    unfold_thetaScatt(&fresp_thetaScatt, thetaScatt_measdata, int(_sys["niter"]));
  unfold_thetaScatt.Print();
  RooUnfoldBayes 
    unfold_theta2Scatt(&fresp_theta2Scatt, theta2Scatt_measdata, int(_sys["niter"]));
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
  //thetaX_reco = trkreffix(thetaX_reco);
  TH1D* thetaX_reco_noeffi = (TH1D*)unfold_thetaX.Hreco();
  tmpname = "thetaX_reco_noeffi";
  tmpname += model;
  thetaX_reco_noeffi->SetName(tmpname.c_str());
  TH1D* thetaY_reco = (TH1D*)unfold_thetaY.Hreco();
  tmpname = "thetaY_reco";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaY_reco->SetName(tmpname.c_str());
  thetaY_reco->SetTitle(";#Delta #theta_{Y}");
  tmpname = "thetaY_reco_noeffi";
  TH1D* thetaY_reco_noeffi = (TH1D*)unfold_thetaY.Hreco();
  tmpname += model;
  thetaY_reco_noeffi->SetName(tmpname.c_str());
  //thetaY_reco = trkreffiy(thetaY_reco);
  TH1D* thetaScatt_reco = (TH1D*)unfold_thetaScatt.Hreco();
  tmpname = "thetaScatt_reco";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaScatt_reco->SetName(tmpname.c_str());
  thetaScatt_reco->SetTitle(";#theta_{Scatt}");
  //thetaScatt_reco = trkreffiscatt(thetaScatt_reco);
  TH1D* theta2Scatt_reco = (TH1D*)unfold_theta2Scatt.Hreco();
  tmpname = "theta2Scatt_reco";
  tmpname += model;
  // if(j>0) tmpname += j;
  theta2Scatt_reco->SetName(tmpname.c_str());
  theta2Scatt_reco->SetTitle(";#theta^{2}_{Scatt}");
  //theta2Scatt_reco = trkreffi2scatt(theta2Scatt_reco);
  
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
 

  TH2D* dotw = new TH2D("dotw","dotw", 180, 0, 180, 200, -1, 1);
  TH2D* dotv = new TH2D("dotv","dotv", 180, 0, 180, 200, -1, 1);
  TH2D* histoyd = new TH2D("histoyd","RMS with plane definition", 180, 0, 180, 100, 1, 3);
  TH2D* histowd = new TH2D("histowd","RMS with plane definition", 180, 0, 180, 100, 1, 3);
  TH2D* histowx = new TH2D("histowx","RMS with plane definition", 180, 0, 180, 100, -2, 2);
  TH2D* histowy = new TH2D("histowy","RMS with plane definition", 180, 0, 180, 100, -2, 2);
  outfile->cd();
  double x[180];
  double xer[180];
  double y[180];
  double yer[180];
  double xaxis[180];
  double xaxiser[180];
  for (int l=0;l<180;l++){
     //std::cout << "Histo[l]->GetRMS() " << Histo[l]->GetRMS() << std::endl;
     //std::cout << "Histoy[l]->GetRMS() " << Histoy[l]->GetRMS() << std::endl;
     xaxis[l]=l;
     xaxiser[l]=0;
     x[l]=Histo[l]->GetRMS();
     xer[l]=0.00054;
     y[l]=Histoy[l]->GetRMS();
     yer[l]=0.0006;
     dotw->Fill(l,Histodotw[l]->GetMean());
     dotv->Fill(l,Histodotv[l]->GetMean());
     histoyd->Fill(l,Histoyd[l]->GetMean());
     histowd->Fill(l,Histowd[l]->GetMean());
     histowx->Fill(l,Histowx[l]->GetMean());
     histowy->Fill(l,Histowy[l]->GetMean());
     //Histo[l]->Write();
     //Histoy[l]->Write();
     //Histowx[l]->Write();
     //Histowy[l]->Write();
     //Histoyd[l]->Write(); 
     //Histowd[l]->Write(); 
     //Histodotw[l]->Write(); 
     //Histodotv[l]->Write(); 
  }
  //histowx->Write();
  //histowy->Write();
  //histoyd->Write();
  //histowd->Write();
  //dotw->Write();
  //dotv->Write();
  thetaXUS_thetaXDS->Write();
  thetaYUS_thetaYDS->Write();
  //thetaX_measdata = trkreffix(thetaX_measdata);
  outfile->cd();
  thetaX_measdata->Write();
  //thetaY_measdata = trkreffiy(thetaY_measdata);
  outfile->cd();
  thetaY_measdata->Write();
  //thetaScatt_measdata = trkreffiscatt(thetaScatt_measdata);
  outfile->cd();
  thetaScatt_measdata->Write();
  //theta2Scatt_measdata = trkreffi2scatt(theta2Scatt_measdata);
  outfile->cd();
  theta2Scatt_measdata->Write();
  thetaScatt_measdata_vp->Write();
  thetaX_reco->Write();
  thetaY_reco->Write();
  thetaX_reco_noeffi->Write();
  thetaY_reco_noeffi->Write();
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
  //convergence->Write();
  
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
  thetaX_reco_noeffi->Draw();
  c1->Print("thetaX_reco_noeffi.pdf");
  c1->Clear();
  thetaY_reco_noeffi->Draw();
  c1->Print("thetaY_reco_noeffi.pdf");
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
  TGraphErrors* RotDefHis = new TGraphErrors(180, xaxis,x,xaxiser,xer);
  TGraphErrors* RotDefHisy = new TGraphErrors(180, xaxis,y,xaxiser,yer);
  RotDefHis->SetFillColor(6);
  RotDefHis->SetFillStyle(3005);
  RotDefHis->SetMarkerStyle(34);
  RotDefHis->SetTitle("");
  RotDefHis->Draw("Pa3");
  RotDefHis->GetXaxis()->SetTitle("Angle around Z axis (#circ)");
  RotDefHis->GetYaxis()->SetTitle("RMS of scattering distribution");
  RotDefHis->GetYaxis()->SetTitleOffset(1.5);
  RotDefHis->Write();
  c1->SaveAs("RotDefHis.pdf");
  c1->Clear();
  RotDefHisy->SetFillColor(6);
  RotDefHisy->SetFillStyle(3005);
  RotDefHisy->SetMarkerStyle(34);
  RotDefHis->SetTitle("");
  RotDefHisy->Draw("Pa3");
  RotDefHisy->Write();
  RotDefHisy->GetXaxis()->SetTitle("Angle around Z axis (#circ)");
  RotDefHisy->GetYaxis()->SetTitle("RMS of scattering distribution");
  RotDefHisy->GetYaxis()->SetTitleOffset(1.5);
  c1->SaveAs("RotDefHisy.pdf");
  c1->Clear();
  c1->SetLogy();
  scattering_proj_x->Draw();
  c1->Print("scattering_proj_x.pdf");
  c1->Clear();
  scattering_proj_y->Draw();
  c1->Print("scattering_proj_y.pdf");
  /*
  c1->Clear();
  c1->SetLogy();
  gStyle->SetOptStat(0);
  convergence->Draw();
  //convergence->GetXaxis()->SetTitle("No. of iterations");
  //convergence->GetYaxis()->SetTitle("#chi^{2}_{n}-#chi^{2}_{n-1}");
  c1->SaveAs("convergence.pdf");
  c1->SaveAs("convergence.root");
  */
  
  delete c1;
  delete thetaX_measdata;
  delete thetaY_measdata;
  delete thetaScatt_measdata;
  delete thetaScatt_measdata_vp;
  delete thetaX_reco;
  delete thetaY_reco;
  //delete thetaX_reco_noeffi;
  //delete thetaY_reco_noeffi;
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
  //delete convergence;

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
  TH1D* scattering_proj_x = new TH1D("scattering_proj_x_DC","Change in Projected Angle (X);#Delta#theta_{X}; Events per radian", 
				     NBINS, scat_bin_array);
  TH1D* scattering_proj_y = new TH1D("scattering_proj_y_DC","Change in Projected Angle (Y);#Delta#theta_{Y}; Events per radian", 
				     NBINS, scat_bin_array);
  
  TH2D* projposUSDSdiff = 
    new TH2D("projposUSDSdiff_ref",
	     "Difference of US and DS projections at absorber;#Delta x_{DS-US} (mm); #Delta y_{DS-US}",
	     180, -400, 400, 180, -400, 400);

  for(int i=0; i<_DSset.N(); i++){
    std::vector<double> projTheta = RotDefineProjectionAngles(_USset.E(i), _DSset.E(i), angdef);
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
    scattering_proj_x->Fill(thetaX);
    scattering_proj_y->Fill(thetaY);
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
    if (USrefplaneI < mcevent->GetVirtualHits()->size() &&
	DSrefplaneI < mcevent->GetVirtualHits()->size() &&
	DSabsPlaneI < mcevent->GetVirtualHits()->size() &&
	USabsPlaneI < mcevent->GetVirtualHits()->size()){
      planesfound = true;
      return planesfound;
    }
  }
  if (1) { 
    double USmindiff = 9999, DSmindiff=9999;
    int USminI=-1, DSminI=-1;
    double USrefpos = _sys["abspos"] - 500.;
    double DSrefpos = _sys["abspos"] + 500.;
    for ( size_t j=0; j < mcevent->GetVirtualHits()->size(); j++){
      double ztest = mcevent->GetVirtualHits()->at(j).GetPosition().z();
      double USzdiff = fabs(ztest - USrefpos);
      double DSzdiff = fabs(ztest - DSrefpos);
      //std::cout<<j<<"\t"<<ztest<<"\t"<<USzdiff<<"\t"<<DSzdiff<<std::endl;
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
    USrefpos = _sys["USref"];
    DSrefpos = _sys["DSref"];
    
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
    //std::cout << "planesfound " << planesfound << std::endl;
    return planesfound;
  }
}
void MCSAnalysis::FillMuScattResponse(bool event_ok, Vars& US, Vars& DS, Vars& USMC, Vars& DSMC){
  double thetaYMC = atan(DSMC.dXdz) - atan(USMC.dXdz);
  double thetaXMC = atan(DSMC.dYdz) - atan(USMC.dYdz);
  double XMC = DSMC.X;
  double YMC = DSMC.Y;
  //theta_true_x_graph->Fill(thetaXMC);
  //theta_true_y_graph->Fill(thetaYMC);
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
    
    //theta_true_x_graph->Fill(thetaXMC);
    //theta_true_y_graph->Fill(thetaYMC);
  }
  else{
    resp_thetaX.Miss(thetaXMC);
    resp_thetaY.Miss(thetaYMC);
    resp_thetaScatt.Miss(thetaScattMC);
    resp_theta2Scatt.Miss(thetaScattMC * thetaScattMC);

    //theta_true_x_graph->Fill(thetaXMC);
    //theta_true_y_graph->Fill(thetaYMC);
  }
  
}
 
void MCSAnalysis::TruthGraph(Collection& USMCTruthset, Collection& DSMCTruthset){
 
  int curr_sel = int(DSMCTruthset.N());
  int curr_k = 0;
  for(int i=curr_k; i<curr_sel; i++){
  std::vector<double> projDTheta =RotDefineProjectionAngles(USMCTruthset.E(i), DSMCTruthset.E(i), angdef); 
  double thetaYMC = projDTheta[1];
  double thetaXMC = projDTheta[0];
  double thetascattMC = projDTheta[2];

	  //std::cout << "i " << i << std::endl;
	  //std::cout << "projDTheta[0] " << projDTheta[0] << std::endl;
  theta_true_x_graph->Fill(thetaXMC);
  theta_true_y_graph->Fill(thetaYMC);
  theta_true_scat_graph->Fill(thetascattMC);
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
  tmpvar.px = mcevent->GetVirtualHits()->at(j).GetMomentum().x();
  tmpvar.py = mcevent->GetVirtualHits()->at(j).GetMomentum().y();
//  double px = mcevent->GetVirtualHits()->at(j).GetMomentum().x();
//  double py = mcevent->GetVirtualHits()->at(j).GetMomentum().y();
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

void MCSAnalysis::FillCollectionSciFi(Collection& Set, int j, int k, double pz, int isDS, bool project){      

  if(j < int(scifievent->scifitracks().size()) && j != -1){
    if(k < int(scifievent->scifitracks()[j]->scifitrackpoints().size()) && k != -1){
      Vars tmpvar;
      FillVarsSciFi(tmpvar, j, k, pz, isDS);
      if (project){
	Vars newvar = PropagateVarsMu(tmpvar, _sys["abspos"]);
	tmpvar = newvar;
      }
      Set.append_instance(tmpvar);
    }
  }
  if( j == -1 || k == -1){
    Vars tmpvar;
    tmpvar.X = 0.0;
    tmpvar.Y = 0.0;
    tmpvar.Z = 0.0;
    double v1 = rand() % 100;
    if (v1/100 >0.5) {
    tmpvar.dXdz = -1./2.;
    tmpvar.dYdz = -1./2.;
    tmpvar.px   = -pz/2.;
    tmpvar.py   = -pz/2.;
    tmpvar.pz   = pz;
    }
    else {
    tmpvar.dXdz = 1./2.;
    tmpvar.dYdz = 1./2.;
    tmpvar.px   = pz/2.;
    tmpvar.py   = pz/2.;
    tmpvar.pz   = pz;
    }
    tmpvar.pid  = -13;
    tmpvar.isgood = false;

    if( int(tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()) > 0 &&
	int(tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray().size()) > 0)
      tmpvar.TOF12 = tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray()[0].GetTime() 
	- tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime();
    else
      tmpvar.TOF12 = 100.0 * 8.22475 / 0.299792458;
    if( int(tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray().size()) > 0 &&
	int(tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()) > 0)
      tmpvar.TOF01 = tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime() 
	- tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[0].GetTime();
    else
      tmpvar.TOF01 = 100.0 * 7.64186 / 0.299792458;
    
    Set.append_instance(tmpvar);
  }
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
    double v1 = rand() % 100;
    if (v1/100 >0.5) {
    tmpvar.dXdz = -1./2.;
    tmpvar.dYdz = -1./2.;
    tmpvar.px   = -pz/2.;
    tmpvar.py   = -pz/2.;
    tmpvar.pz   = pz;
    }
    else {
    tmpvar.dXdz = 1./2.;
    tmpvar.dYdz = 1./2.;
    tmpvar.px   = pz/2.;
    tmpvar.py   = pz/2.;
    tmpvar.pz   = pz;
    }
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
  TH2D* dXdzdYdz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -0.125, 0.125, 90, -0.125, 0.125);
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
  tmptitle = desc + ";TOF01 (ns)";
  tmpname  = suffix + "_TOF01";
  TH1D* TOF01 = new TH1D(tmpname.c_str(),tmptitle.c_str(), 100, 25.0, 35.0);
  
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
    TOF01->Fill(Set.E(i).TOF01);
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
  TOF01->GetXaxis()->SetLabelSize(0.05);
  TOF01->GetXaxis()->SetTitleSize(0.05);
  TOF01->GetYaxis()->SetLabelSize(0.05);
  TOF01->GetYaxis()->SetTitleSize(0.05);

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
  TOF01->Draw("colz");
  tmpfile = suffix + "_TOF01.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();

  outfile->cd();
  XY->Write();
  XdXdz->Write();
  XdYdz->Write();
  XTOF12->Write();
  YdXdz->Write();
  YdYdz->Write();
  YTOF12->Write();
  dXdzdYdz->Write();
  dXdzTOF12->Write();
  dYdzTOF12->Write();
  TOF01->Write();

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
  delete TOF01;
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

  tmptitle = desc + ";|#hat{r}^{proj to LiH}_{US} - #hat{r}^{proj to LiH}_{DS} | (mm); No of events";
  tmpname  = suffix + "_appRef";
  TH1D* appRef = new TH1D(tmpname.c_str(), tmptitle.c_str(), 200, 0, 200);
  tmptitle = desc + ";X_{US}^{proj to LiH} - X_{DS}^{proj to LiH} (mm);Y_{US}^{proj to LiH} - Y_{DS}^{proj to LiH} (mm)";
  tmpname  = suffix + "_projRefDiff";
  TH2D* projRefDiff = new TH2D(tmpname.c_str(), tmptitle.c_str(), 200, -100, 100, 200, -100, 100);
  tmptitle = desc + ";dXdz_{US}^{proj to LiH} - dXdz_{DS}^{proj to LiH} (mm);dYdz_{US}^{proj to LiH} - dYdz_{DS}^{proj to LiH} (mm)";
  tmpname  = suffix + "_projRefDiffddz";
  TH2D* projRefDiffddz = new TH2D(tmpname.c_str(), tmptitle.c_str(), 90, -0.125, 0.125, 90, -0.125, 0.125);

  double latchZ = 0.0;
  Vars USproj;
  Vars DSproj;
  for (int i=0; i<USset.N(); i++){
	  /*
    if (DSset.E(i).Z != 0 && latchZ == 0.0){
      latchZ = DSset.E(i).Z;
    } 
    if (latchZ != 0.0 && fabs(DSset.E(i).Z - latchZ) < 5){
//      std::cerr<<"make_acceptance_histograms - Propagating to DSset.E(i).Z"<<std::endl;
      DSproj = PropagateVarsMu(USset.E(i), DSset.E(i).Z);
      */
    //} else if(latchZ != 0.0){
      USproj = PropagateVarsMu(USset.E(i), latchZ);
      DSproj = PropagateVarsMu(DSset.E(i), latchZ);
    //} else {
    //  continue;
    //}
    
    Vars predDiff = USproj - DSproj;
    double approach = sqrt(predDiff.X * predDiff.X + 
			   predDiff.Y * predDiff.Y +
			   predDiff.Z * predDiff.Z);
    appRef->Fill(approach); // can be very large if the ds track is lost
    projRefDiff->Fill(predDiff.X, predDiff.Y);
    posXY->Fill(DSproj.X, DSproj.Y);
    posDXDY->Fill(DSproj.dXdz, DSproj.dYdz);
    projRefDiffddz->Fill(predDiff.dXdz, predDiff.dYdz);
    
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

  TCanvas *c1 = new TCanvas();
  appRef->Draw("P");
  c1->SaveAs("appRef.pdf");
  c1->Clear();
  projRefDiff->Draw("P");
  c1->SaveAs("projRefDiff.pdf");
  c1->Clear();
  projRefDiffddz->Draw("P");
  c1->SaveAs("projRefDiffddz.pdf");
  delete c1;
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
  projRefDiffddz->Write();

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
    std::vector<double> projTheta = RotDefineProjectionAngles(USset.E(i), DSset.E(i),angdef);
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

void MCSAnalysis::ReturnWords(std::string inputstring,std::vector<std::string> &v){
  std::string word="";
  v.clear();

//  std::cerr<<"Input String: "<<inputstring<<std::endl;
  for (auto x : inputstring){
    if (x == ' '){
      v.push_back(word);
//      std::cerr<<"Word being added: "<<word<<std::endl;
      word="";
    } else {
      word=word+x;
//      std::cerr<<"Still constructing word: "<<word<<std::endl;
    }
  }
  v.push_back(word);
}

double MCSAnalysis::MaxRadius(double bz,double startx, double starty, double startz, double startpx, double startpy, double endx, double endy, double endz, double endpx, double endpy, double endpz){
//     std::cerr<<"Axial magnetic field  : "<<bz<<std::endl;

//------------------------------------------------------------------------------------------------------
//    Set the maximum r2 on each track point between tp and the next tp that
//    the track excurts from the axis. Estimate the region between adjacent track
//    points by extrapolating assuming uniform field and calculate the excursion
//    analytically. Sets the max_r2 parameter on the track.
//    Get azimuthal angles
//    - tp_0: track point with lower z
//    - tp_1: track point with higher z
//    - bz: magnetic field
//    Returns azimuthal angles, sorted for direction of rotation. If more than one
//    rotation, add the delta phi on (so that we get extra turns).

     double q = 1;
     double phi_start, phi_end;
     if (q*bz > 0){
        phi_start = atan2(startpx, startpy);
        phi_end = atan2(endpx, endpy);
     }
     else if (q*bz < 0){
        phi_start = atan2(endpx, endpy);
        phi_end = atan2(startpx, startpy);
     }
     double k0 = q*bz*c_light/endpz;
     while (phi_start < 0) { //phi_start is always in domain 0 < 2 pi
        phi_start += 2*pi;
     }
     while (phi_end < phi_start) { // phi_end is always > phi_start
        phi_end += 2*pi;
     }
    // for larger steps, we can make more than one rotation
     double delta_phi = k0*(startz - endz);
     double n_rot = abs(delta_phi)/pi/2.;
//    std::cerr<<"delta_phi: "<<delta_phi<<" n_rot: "<<n_rot<<std::endl;
     phi_end += floor(n_rot)*2.*pi;
/*
    Get the azimuthal angle at which the trajectory has maximum excursion from
    the axis
    - tp_0: track point with lower z
    - tp_1: track point with higher z
    - bz: magnetic field
    - phi_start: phi at which we start
    Returns azimuthal angles, sorted for direction of rotation. We add on 2pi so
    that it is always greater than phi_start.
*/
     double tp_start_x, tp_start_y, tp_start_pt;//tp_start_z,
     if (q*bz > 0){
        tp_start_x = startx;
        tp_start_y = starty;
//        tp_start_z = startz;
        tp_start_pt = sqrt(startpx*startpx+startpy*startpy);
     }
     else if (q*bz < 0){
        tp_start_x = endx;
        tp_start_y = endy;
//        tp_start_z = endz;
        tp_start_pt = sqrt(endpx*endpx+endpy*endpy);
     }
    // coordinates of helix centre
     double pt = tp_start_pt;
     double x0 = tp_start_x + pt*cos(phi_start)/q/bz/c_light;
     double y0 = tp_start_y - pt*sin(phi_start)/q/bz/c_light;
     double r0 = abs(pt/q/bz/c_light); // wavenumber is qBm/pz dz
//    std::cerr<<"pt: "<<pt<<" x0: "<<x0<<" y0: "<<y0<<" r0: "<<r0<<std::endl;
    //print 'x0:', x0, 'y0:', y0, 'r0:', r0, 'k0:', k0
    //azimuthal angle of maximum extent (when r is in line with x0, y0)
     double direction = abs(q*bz)/q/bz;
     double phi_max = atan2(direction*y0, -direction*x0);
     while (phi_max < phi_start){// phi_max is always > phi_start
        phi_max += 2*pi;
     }
     double max_r2 = (sqrt(x0*x0 + y0*y0)+r0); // vector to centre + radius of track
     double r2_start = startx*startx+starty*starty;
     double r2_end = endx*endx+endy*endy;
//    std::cerr<<"max_r2: "<<max_r2<<" r2_start: "<<r2_start<<" r2_end: "<<r2_end<<std::endl;
     double fiducial_radius;
     if ((phi_max > phi_start) and (phi_end > phi_max)){ // phi_max is between phi_end and phi_start
        fiducial_radius = max_r2;
     } else {
        fiducial_radius = max(r2_start, r2_end);
     }
     return fiducial_radius;
}
        //print "   tp0 p:", tp_0["hit"]["px"], tp_0["hit"]["py"], "tp1 p:", tp_1["hit"]["px"], tp_1["hit"]["py"]
        //print "   ", str(round(phi_start, 1)).rjust(8), str(round(phi_end, 1)).rjust(8), str(round(phi_max, 1)).rjust(8), str(round(r2_start, 1)).rjust(8), str(round(max_r2, 1)).rjust(8), str(round(r2_end, 1)).rjust(8)



multiVars MCSAnalysis::read_globals(double DS11, double DS53, double US53, double US11, double diffuserlow, double diffuserhigh, double TOF0pos, double TOF1pos){
// This is the main function to read the data from the root files
  double dz = 0.2;
  multiVars tempMV;
  tempMV.UScentre_absorber=reset_Vars();
  tempMV.USend_of_DStracker=reset_Vars();
  tempMV.USstart_of_DStracker=reset_Vars();
  tempMV.USTOF0=reset_Vars();
  tempMV.USTOF1=reset_Vars();
  tempMV.DScentre_absorber=reset_Vars();
  tempMV.UStrackerUS11=reset_Vars();
  tempMV.UStrackerUS53=reset_Vars();
  tempMV.DiffuserStart=reset_Vars();
  tempMV.DiffuserEnd=reset_Vars();
  tempMV.UScentre_absorber.Z=_sys["abspos"];
  tempMV.USend_of_DStracker.Z= DS53;
  tempMV.USstart_of_DStracker.Z= DS11;
  tempMV.USTOF0.Z= TOF0pos;
  tempMV.USTOF1.Z= TOF1pos;
  tempMV.DScentre_absorber.Z=_sys["abspos"];
  tempMV.UStrackerUS11.Z= US11;
  tempMV.UStrackerUS53.Z= US53;
  tempMV.DiffuserStart.Z= diffuserlow;
  tempMV.DiffuserEnd.Z= diffuserhigh;    
    
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
//        std::cerr<<"Chain type: "<<chain->get_chain_type()<<std::endl;

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
                            if(a_track_pos.Z() >= tempMV.UScentre_absorber.Z-dz && a_track_pos.Z() <= tempMV.UScentre_absorber.Z+dz){
                                tempMV.UScentre_absorber=FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.USend_of_DStracker.Z-dz && a_track_pos.Z() <= tempMV.USend_of_DStracker.Z+dz){
                                tempMV.USend_of_DStracker = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.USstart_of_DStracker.Z-dz && a_track_pos.Z() <= tempMV.USstart_of_DStracker.Z+dz){
                                tempMV.USstart_of_DStracker = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.USTOF0.Z-dz && a_track_pos.Z() <= tempMV.USTOF0.Z+dz){
                                tempMV.USTOF0 = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.USTOF1.Z-dz && a_track_pos.Z() <= tempMV.USTOF1.Z+dz){
                                tempMV.USTOF1 = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.UStrackerUS11.Z-dz && a_track_pos.Z() <= tempMV.UStrackerUS11.Z+dz){
                                tempMV.UStrackerUS11 = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.UStrackerUS53.Z-dz && a_track_pos.Z() <= tempMV.UStrackerUS53.Z+dz){
                                tempMV.UStrackerUS53 = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.DiffuserStart.Z-dz && a_track_pos.Z() <= tempMV.DiffuserStart.Z+dz){
                                tempMV.DiffuserStart = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.DiffuserEnd.Z-dz && a_track_pos.Z() <= tempMV.DiffuserEnd.Z+dz){
                                tempMV.DiffuserEnd = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            }
                        }
                    }
                  }
            }

//Chain5 - Get Downstream data from DSDaughter

            MAUS::DataStructure::Global::PrimaryChain* DSchain=chain->GetDSDaughter();
            some_tracks = DSchain->GetMatchedTracks();
            for(some_track_iterator = some_tracks.begin(); some_track_iterator < some_tracks.end(); ++some_track_iterator){
                MAUS::DataStructure::Global::Track* a_track = (*some_track_iterator);
                if(a_track != NULL){
                    if(a_track->get_pid() == -13 || a_track->get_pid() == 13){
                        std::vector< const MAUS::DataStructure::Global::TrackPoint * > track_points = a_track->GetTrackPoints(virtual_plane);
                        for(unsigned int p = 0; p < track_points.size(); p++){
                            TLorentzVector a_track_mom = track_points.at(p)->get_momentum();
                            TLorentzVector a_track_pos = track_points.at(p)->get_position();
                            if((a_track_pos.Z() >= DS11-dz && a_track_pos.Z() <= DS11+dz)&&tempMV.DScentre_absorber.isgood==false){
                                tempMV.DScentre_absorber=FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                                tempMV.DScentre_absorber=PropagateVarsMu(tempMV.DScentre_absorber,_sys["abspos"]);
                            }else if((a_track_pos.Z() >= DS11+0.65-dz && a_track_pos.Z() <= DS11+0.65+dz)&&tempMV.DScentre_absorber.isgood==false){
                                tempMV.DScentre_absorber=FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                                tempMV.DScentre_absorber=PropagateVarsMu(tempMV.DScentre_absorber,_sys["abspos"]);
                            }else if((a_track_pos.Z() >= DS11+1.3-dz && a_track_pos.Z() <= DS11+1.3+dz)&&tempMV.DScentre_absorber.isgood==false){
                                tempMV.DScentre_absorber=FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                                tempMV.DScentre_absorber=PropagateVarsMu(tempMV.DScentre_absorber,_sys["abspos"]);
                            }
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
                            if(a_track_pos.Z() >= tempMV.UScentre_absorber.Z-dz && a_track_pos.Z() <= tempMV.UScentre_absorber.Z+dz){
                                tempMV.UScentre_absorber=FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.USend_of_DStracker.Z-dz && a_track_pos.Z() <= tempMV.USend_of_DStracker.Z+dz){
                                tempMV.USend_of_DStracker = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.USstart_of_DStracker.Z-dz && a_track_pos.Z() <= tempMV.USstart_of_DStracker.Z+dz){
                                tempMV.USstart_of_DStracker = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.USTOF0.Z-dz && a_track_pos.Z() <= tempMV.USTOF0.Z+dz){
                                tempMV.USTOF0 = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.USTOF1.Z-dz && a_track_pos.Z() <= tempMV.USTOF1.Z+dz){
                                tempMV.USTOF1 = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.UStrackerUS11.Z-dz && a_track_pos.Z() <= tempMV.UStrackerUS11.Z+dz){
                                tempMV.UStrackerUS11 = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.UStrackerUS53.Z-dz && a_track_pos.Z() <= tempMV.UStrackerUS53.Z+dz){
                                tempMV.UStrackerUS53 = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.DiffuserStart.Z-dz && a_track_pos.Z() <= tempMV.DiffuserStart.Z+dz){
                                tempMV.DiffuserStart = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            } else if(a_track_pos.Z() >= tempMV.DiffuserEnd.Z-dz && a_track_pos.Z() <= tempMV.DiffuserEnd.Z+dz){
                                tempMV.DiffuserEnd = FillVars(a_track_mom,a_track_pos,a_track->get_pid());
                            }
                            if (tempMV.UScentre_absorber.isgood==true && tempMV.USend_of_DStracker.isgood==true   && tempMV.USstart_of_DStracker.isgood==true
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
     if(tempMV.USstart_of_DStracker.isgood==false){
        tempMV.USstart_of_DStracker=PropagateVarsMu(tempMV.UStrackerUS11, DS11);
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
     if(tempMV.DiffuserStart.isgood==false){
        tempMV.DiffuserStart=PropagateVarsMu(tempMV.UStrackerUS11, diffuserlow);
     }              
     if(tempMV.DiffuserEnd.isgood==false){
        tempMV.DiffuserEnd=PropagateVarsMu(tempMV.UStrackerUS11, diffuserhigh);
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

TH1D *defineHist2(const char* fname, const char* ftitle, Int_t fnbinsx, Double_t fxlow, Double_t fxup)
{
  TH1D *fhis1D = new TH1D(fname, ftitle, fnbinsx, fxlow, fxup);
  fhis1D->SetMinimum(0.001);
  fhis1D->GetXaxis()->SetTitle(ftitle);
  Double_t binning = ( (fhis1D->GetXaxis()->GetXmax()) - (fhis1D->GetXaxis()->GetXmin()) ) / (fhis1D->GetNbinsX());
  fhis1D->GetYaxis()->SetTitle( Form("Events / (%.2f [Mev/c^{2}])", binning) );
  fhis1D->GetYaxis()->SetTitleOffset(1.6);
  fhis1D->GetYaxis()->SetLabelSize(0.035);
  return fhis1D;
}




