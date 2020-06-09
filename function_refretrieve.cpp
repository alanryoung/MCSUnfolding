/////////////////////////////////////////////////////////////////////
// Function Refretrieve
/////////////////////////////////////////////////////////////////////

void MCSAnalysis::refretrieve(){
  Collection USretrieved, DSretrieved, EndofDS, StartofDS, USPostMom, DSPostMom, DSsetUncorrected;
  Collection Collect_USTOF0;
  USretrieved.load_collection(USRefname);
  DSretrieved.load_collection(DSRefname);
  EndofDS.load_collection(EndofDStrackerrefname);
  StartofDS.load_collection(StartofDStrackerrefname);
  Collect_USTOF0.load_collection(Collect_USTOF0refname);

  std::cerr<<"Entries in ref USretrieved: "<< USretrieved.N()<<std::endl;
  std::cerr<<"Entries in ref DSretrieved: "<< DSretrieved.N()<<std::endl;
  std::cerr<<"Entries in EndofDSref: "<< EndofDS.N()<<std::endl;
  std::cerr<<"Entries in StartofDSref: "<< StartofDS.N()<<std::endl;
  std::cerr<<"Entries in Collect_USTOF0: "<< Collect_USTOF0.N()<<std::endl;

  outfile->cd();
  BTField* field = dynamic_cast<BTField*>(MAUS::Globals::GetMCFieldConstructor());
  double pi = 3.14159265358979323846;
  double c_light=299.792458;
  double x = 0; double y = 0; double t = 0;
//Plot axial magnetic field
  for(int i=0;i<1000;i++){
     double z=(_sys["dataDS11"]+((i+0.5)*(_sys["dataDS53"]-_sys["dataDS11"])/1000));
     std::vector<double> magfield=field->GetFieldValue2(x,y,z,t);
     Bfield_DStrackerX->Fill(z,magfield[0]);
     Bfield_DStrackerY->Fill(z,magfield[1]);
     Bfield_DStrackerZ->Fill(z,magfield[2]);
  }


  for(int i=0;i<USretrieved.N();i++){
//  for(int i=0;i<100000;i++){
     double endx=EndofDS.E(i).X;
     double endy=EndofDS.E(i).Y;
     double endz=EndofDS.E(i).Z;
     double endpx=EndofDS.E(i).px;
     double endpy=EndofDS.E(i).py;
     double endpz=EndofDS.E(i).pz;
     double startx=StartofDS.E(i).X;
     double starty=StartofDS.E(i).Y;
     double startz=StartofDS.E(i).Z;
     double startpx=StartofDS.E(i).px;
     double startpy=StartofDS.E(i).py;
     double startpz=StartofDS.E(i).pz;
     std::vector<double> magfield=field->GetFieldValue2(x,y,startz,t);
     double bz=magfield[2];


//------------------------------------------------------------------------------------------------------
// Fiducial Cut
//------------------------------------------------------------------------------------------------------
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
//     std::cerr<<"q: "<<q<<" bz: "<<bz<<" c_light: "<<c_light<<" endpz: "<<endpz<<" k0: "<<k0<<std::endl;
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
    double tp_start_x, tp_start_y, tp_start_z, tp_start_pt;
    if (q*bz > 0){
        tp_start_x = startx;
        tp_start_y = starty;
        tp_start_z = startz;
        tp_start_pt = sqrt(startpx*startpx+startpy*startpy);
    }
    else if (q*bz < 0){
        tp_start_x = endx;
        tp_start_y = endy;
        tp_start_z = endz;
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
        //print "   tp0 p:", tp_0["hit"]["px"], tp_0["hit"]["py"], "tp1 p:", tp_1["hit"]["px"], tp_1["hit"]["py"]
        //print "   ", str(round(phi_start, 1)).rjust(8), str(round(phi_end, 1)).rjust(8), str(round(phi_max, 1)).rjust(8), str(round(r2_start, 1)).rjust(8), str(round(max_r2, 1)).rjust(8), str(round(r2_end, 1)).rjust(8)  
//------------------------------------------------------------------------------------------------------
    
    
//     double fiducial_radius = sqrt(EndofDS.E(i).X*EndofDS.E(i).X + EndofDS.E(i).Y*EndofDS.E(i).Y);
     mcpass_Fid->Fill(fiducial_radius);
//     mcpass_Fid_wide->Fill(fiducial_radius);
     if ( fiducial_radius > fiducial_cut_radius){
        if (cut_fiducial) continue;
     } 

//---------------------------------------------------------------------------------------
// Momentum at absorber cut. 
//---------------------------------------------------------------------------------------
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

