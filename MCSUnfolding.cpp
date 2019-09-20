#include "MCSAnalysis.h"

// xml parser
#include "libxml/tree.h"
#include "libxml/parser.h"
#include "libxml/xpath.h"
#include "libxml/xpathInternals.h"

#if defined(LIBXML_XPATH_ENABLED) && defined(LIBXML_SAX1_ENABLED)

struct specvals {
  std::string outfile;
  std::string outfilename;
  std::string dataname;
  std::string emptydataname;
  std::string trainname;
  std::string modelname;
  std::string modelname2;
  std::string geometryname;
  std::string USdatafile;
  std::string DSdatafile;
  std::string USreffile;
  std::string DSreffile;
  std::string EndofDStrackerrefname;
  std::string EndofDStrackerdataname;
  std::string StartofDStrackerrefname;
  std::string StartofDStrackerdataname;

  std::string model1;			       
  std::string model2;
  std::string model3;
  std::string material;
  std::string trkreffiname;
  std::string trkreffiemptyname;
  std::map<std::string, double> sys;
  std::map<std::string, double> histlimits;
  double angdef;
  double TOF_ll;
  double TOF_ul;
  double mom_ll;
  double mom_ul;
  double fid_grad;
  double fid_rad;
  int    mode;
  bool cut_TOF1_1sp;
  bool cut_TOF0_1sp;
  bool cut_TKU1track;
  bool cut_CHI2DOF;
  bool cut_TKURadial;
  bool cut_DiffuserRadial;
  bool cut_Fiducial;
  bool cut_TOF01;
  bool cut_TOF01ext;
  bool cut_AbsorberMomentum;
  bool cut_ExtrudeTKUTOF0;
  
}; 

static void print_element_names(xmlNode * a_node, specvals& spec)
{
  xmlNode *cur_node = NULL;
  for (cur_node= a_node; cur_node; cur_node = cur_node->next) { 
    if (cur_node->type == XML_ELEMENT_NODE) {
      const xmlChar* id = xmlCharStrdup("id");
      const xmlChar* nm = xmlCharStrdup("name");
      const xmlChar* vl = xmlCharStrdup("value");
      char* val1 = (char*)xmlGetProp(cur_node, id);
      char* val2 = (char*)xmlGetProp(cur_node, nm);
      char* val3 = (char*)xmlGetProp(cur_node, vl);
      printf("node type: Element, name: %s, with prop id: %s, prop name: %s, prop value: %s\n", 
	     cur_node->name, val1, val2, val3);
      if ( xmlStrEqual(cur_node->name, xmlCharStrdup("file")) ) {
	if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("outfile")) ){
	  spec.outfile = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("datafile")) ){
	  spec.dataname = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("outfilename")) ){
	  spec.outfilename = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("trainfile")) ){
	  spec.trainname = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("modelfile")) ){
	  spec.modelname = (char*)xmlGetProp(cur_node, nm);
        } else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("modelfile2")) ){
	  spec.modelname = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("USdatafile")) ){
	  spec.USdatafile = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("DSdatafile")) ){
	  spec.DSdatafile = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("USreffile")) ){
	  spec.USreffile = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("DSreffile")) ){
	  spec.DSreffile = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("EndofDStrackerrefname")) ){
	  spec.EndofDStrackerrefname = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("EndofDStrackerdataname")) ){
	  spec.EndofDStrackerdataname = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("StartofDStrackerrefname")) ){
	  spec.StartofDStrackerrefname = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("StartofDStrackerdataname")) ){
	  spec.StartofDStrackerdataname = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("model1")) ) {
	  spec.model1 = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("model2")) ) {
	  spec.model2 = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("model3")) ) {
	  spec.model3 = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("material")) ) {
	  spec.material = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("geofile")) ) {
	  spec.geometryname = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("trkreffiname")) ) {
	  spec.trkreffiname = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("trkreffiemptyname")) ) {
	  spec.trkreffiemptyname = (char*)xmlGetProp(cur_node, nm);
	}else {
	  std::cout<<"File designation not recognized. Will ignore "
		   << xmlGetProp(cur_node, id) 
		   << std::endl;
	}
      }

      if (xmlStrEqual(cur_node->name, xmlCharStrdup("cuts")) ) {
	if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("TOF_ll")) ){
	  spec.TOF_ll = std::atof((char*)xmlGetProp(cur_node, vl));
	} else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("TOF_ul")) ){
	  spec.TOF_ul = std::atof((char*)xmlGetProp(cur_node, vl));
	} else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("mom_ll")) ){
	  spec.mom_ll = std::atof((char*)xmlGetProp(cur_node, vl));
	} else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("mom_ul")) ){
	  spec.mom_ul = std::atof((char*)xmlGetProp(cur_node, vl));
	} else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("fid_grad")) ){
	  spec.fid_grad = std::atof((char*)xmlGetProp(cur_node, vl));
	} else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("fid_rad")) ){
	  spec.fid_rad = std::atof((char*)xmlGetProp(cur_node, vl));
	} else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("mode")) ){
	  spec.mode = std::atof((char*)xmlGetProp(cur_node, vl));
        } else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("cut_TOF1_1sp"))){
          spec.cut_TOF1_1sp = std::atoi((char*)xmlGetProp(cur_node,vl));
        } else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("cut_TOF0_1sp"))){
          spec.cut_TOF0_1sp = std::atoi((char*)xmlGetProp(cur_node,vl));
        } else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("cut_TKU1track"))){
          spec.cut_TKU1track = std::atoi((char*)xmlGetProp(cur_node,vl));
        } else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("cut_CHI2DOF"))){
          spec.cut_CHI2DOF = std::atoi((char*)xmlGetProp(cur_node,vl));
        } else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("cut_TKURadial"))){
          spec.cut_TKURadial = std::atoi((char*)xmlGetProp(cur_node,vl));
        } else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("cut_DiffuserRadial"))){
          spec.cut_DiffuserRadial = std::atoi((char*)xmlGetProp(cur_node,vl));
        } else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("cut_Fiducial"))){
          spec.cut_Fiducial = std::atoi((char*)xmlGetProp(cur_node,vl));
        } else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("cut_TOF01"))){
          spec.cut_TOF01 = std::atoi((char*)xmlGetProp(cur_node,vl));
        } else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("cut_TOF01ext"))){
          spec.cut_TOF01ext = std::atoi((char*)xmlGetProp(cur_node,vl));
        } else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("cut_AbsorberMomentum"))){
          spec.cut_AbsorberMomentum = std::atoi((char*)xmlGetProp(cur_node,vl));
        } else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("cut_ExtrudeTKUTOF0"))){
          spec.cut_ExtrudeTKUTOF0 = std::atoi((char*)xmlGetProp(cur_node,vl));
	} else {
	  std::cout<<"Cut designation not recognized. Will ignore "
		   << xmlGetProp(cur_node, nm) 
		   << std::endl;
	}
      }
      if (xmlStrEqual(cur_node->name, xmlCharStrdup("sys")) ){
	if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("angdef")) ){
	  spec.angdef = std::atof((char*)xmlGetProp(cur_node, vl));
	} else spec.sys[(char*)xmlGetProp(cur_node, nm)] = 
	  std::atof((char*)xmlGetProp(cur_node, vl));
      }
      if (xmlStrEqual(cur_node->name, xmlCharStrdup("histlimits")) ){
	spec.histlimits[(char*)xmlGetProp(cur_node, nm)] = 
	  std::atof((char*)xmlGetProp(cur_node, vl));
      }
    }

    print_element_names(cur_node->children, spec);
  }
}

int main(int argc, char* argv[]) { 
  // char* mrd = getenv("MAUS_ROOT_DIR");
  
  std::cout<<"Executing MCSUnfolding\n";
  std::cout<<"Algorithm to analyse reduced MAUS root\n";
  std::cout<<"trees for multiple scattering measurements\n";
  std::cout<<"and to unfold the distribution\n";
  std::cout<<"\nReceived "<<argc<<" arguements.\n";
  if( argc != 9 && argc != 2){
    std::cout<<"2 arguements are required:\n\n";
    std::cout<<"${exec} ${spec file}\n";
    return -1;
  }

  std::string mode_tree = "";

  specvals spec;

  spec.outfile = "";
  spec.outfilename = "";
  spec.dataname  = "";
  spec.trainname = "";
  spec.modelname = "";
  spec.modelname2 = "";
  spec.geometryname = "";
  spec.USdatafile = "";
  spec.DSdatafile = "";
  spec.USreffile = "";
  spec.DSreffile = "";
  spec.EndofDStrackerrefname = "";
  spec.EndofDStrackerdataname = "";
  spec.StartofDStrackerrefname = "";
  spec.StartofDStrackerdataname = "";
  spec.TOF_ll = 27.0;
  spec.TOF_ul = 42.0;
  spec.mom_ll = 140;
  spec.mom_ul = 300;
  spec.fid_grad = 0.01;
  spec.fid_rad  = 150.0;
  spec.mode     = 0;

  if ( argc == 2 ){
    xmlInitParser();
    xmlDocPtr doc=NULL;


    assert(argv[1]);

    doc = xmlReadFile(argv[1], NULL, 0);
    if (doc == NULL) { 
      printf("error: could not parse file %s\n", argv[1]);
    }
    xmlNode * root_element = xmlDocGetRootElement(doc);

    print_element_names(root_element, spec);

    xmlFreeDoc(doc);

  }
  else if ( argc == 12 ){
    spec.outfile = argv[1];
    spec.dataname  = argv[2];
    spec.trainname = argv[3];
    spec.TOF_ll    = std::atof(argv[4]);
    spec.TOF_ul    = std::atof(argv[5]);
    spec.mom_ll	   = std::atof(argv[6]);
    spec.mom_ul	   = std::atof(argv[7]);
    spec.fid_rad   = std::atof(argv[8]);
    spec.modelname = argv[9];
    spec.modelname2 = argv[9];
    spec.mode      = std::atoi(argv[10]);
    spec.angdef    = std::atof(argv[11]);
    spec.outfilename = argv[12];
  }
    
  if (spec.mode >= 0){
    mode_tree = "recon_reduced_tree";
  } else {
    mode_tree = "recon_reduced_MCtree";
  }
  
  std::cout<<"Writing outfile to "<<spec.outfile<<std::endl; 
//  std::cout<<"Reading outfilename at "<<spec.outfilename<<std::endl; 
  std::cout<<"Reading data from "<<spec.dataname<<std::endl; 
  std::cout<<"Reading reference data from "<<spec.trainname<<std::endl; 
  std::cout<<"Reading models from "<<spec.modelname<<std::endl;
  std::cout<<"Reading models from "<<spec.modelname2<<std::endl; 
  std::cout<<"Reading Upstream data "<<spec.USdatafile<<std::endl;
  std::cout<<"Reading Downstream data "<<spec.DSdatafile<<std::endl;
  std::cout<<"Reading Upstream Ref "<<spec.USreffile<<std::endl;
  std::cout<<"Reading Downstream Ref "<<spec.DSreffile<<std::endl;
  std::cout<<"Reading End of Downstream tracker Ref "<<spec.EndofDStrackerrefname<<std::endl;
  std::cout<<"Reading End of Downstream tracker data "<<spec.EndofDStrackerdataname<<std::endl;
  std::cout<<"Reading Start of Downstream tracker Ref "<<spec.StartofDStrackerrefname<<std::endl;
  std::cout<<"Reading Start of Downstream tracker data "<<spec.StartofDStrackerdataname<<std::endl;


  std::cout<<"Use "<<spec.geometryname<<" for propagation\n";
  std::cout<<"\n";
  std::cout<<"angdef " <<spec.angdef<<std::endl;
  std::cout<<"TOF selection between "
	   <<spec.TOF_ll<<" and "<<spec.TOF_ul<<std::endl;
  std::cout<<"Momentum selection between " 
           <<spec.mom_ll<<" and "<<spec.mom_ul<<std::endl;
  std::cout<<"Fiducial selection contained by radius "
	   <<spec.fid_rad<<" with scattering "<<spec.fid_grad<<std::endl;
  std::cout<<"\n";
  std::cout<<"Cuts Selected\n";
  std::cout<<"Require exactly 1 TOF1 Spacepoint: "<<spec.cut_TOF1_1sp<<std::endl;
  std::cout<<"Require exactly 1 TOF0 Spacepoint: "<<spec.cut_TOF0_1sp<<std::endl;
  std::cout<<"Require exactly 1 track in Upstream tracker: "<<spec.cut_TKU1track<<std::endl;
  std::cout<<"Chi-squared per degree of freedom: "<<spec.cut_CHI2DOF<<std::endl;
  std::cout<<"Maximum radius of muon in upstream tracker: "<<spec.cut_TKURadial<<std::endl;
  std::cout<<"Maximum radius of muon in diffuser: "<<spec.cut_DiffuserRadial<<std::endl;
  std::cout<<"Fiducial: "<<spec.cut_Fiducial<<std::endl;
  std::cout<<"TOF01 time consistent with Muon: "<<spec.cut_TOF01<<std::endl;
  std::cout<<"Extruded TOF01 time consistent with Muon: "<<spec.cut_TOF01ext<<std::endl;
  std::cout<<"Momentum at Absorber: "<<spec.cut_AbsorberMomentum<<std::endl;
  std::cout<<"Successfully extrude from Upstream tracker to TOF0: "<<spec.cut_ExtrudeTKUTOF0<<std::endl;
  std::cout<<"\n";
  std::cout<<"Attempting to create analysis object"<<std::endl;
//  MCSAnalysis anal("recon_reduced_tree", "Truth_reduced_tree", mode_tree, spec.outfile, spec.histlimits);
  MCSAnalysis anal("reduced_tree", "Truth_reduced_tree", "reduced_tree", spec.outfile, spec.histlimits);
  std::cout<<"Analysis object created"<<std::endl;
  
  TFile* data = new TFile(spec.dataname.c_str());
  if(data->IsZombie()){
    std::cout<<"Data file is a zombie. Aborting."<<std::endl;
    return -1;
  }
  TFile* training = new TFile(spec.trainname.c_str());
  if(training->IsZombie()){
    std::cout<<"Training file is a zombie. Aborting."<<std::endl;
    return -1;
  }
  
  anal.Setangdef(spec.angdef);
  anal.SetTOFLowerLimit(spec.TOF_ll);
  anal.SetTOFUpperLimit(spec.TOF_ul);
  anal.SetMomentumLowerLimit(spec.mom_ll);
  anal.SetMomentumUpperLimit(spec.mom_ul);
  anal.SetRadialLimit(spec.fid_rad);
  anal.SetGradientLimit(spec.fid_grad);
  std::cout<<"Check for segmentation fault 1"<<std::endl;
  anal.SetModelFileName(spec.modelname.c_str());
  anal.SetModelName1(spec.model1);
  anal.SetModelName2(spec.model2);
  anal.SetModelName3(spec.model3);
  anal.SetMaterial(spec.material);
  std::cout<<"Check for segmentation fault 2"<<std::endl;
  std::ostringstream strs;
  anal.SetUSDataname(spec.USdatafile);
  anal.SetDSDataname(spec.DSdatafile);
  anal.SetUSRefname(spec.USreffile);
  anal.SetDSRefname(spec.DSreffile);
  anal.SetEndofDStrackerrefname(spec.EndofDStrackerrefname);
  anal.SetEndofDStrackerdataname(spec.EndofDStrackerdataname);
  anal.SetStartofDStrackerrefname(spec.StartofDStrackerrefname);
  anal.SetStartofDStrackerdataname(spec.StartofDStrackerdataname);
  anal.SetTrkrEffiName(spec.trkreffiname.c_str());
  anal.SetTrkrEffiEmptyName(spec.trkreffiemptyname.c_str());
  anal.SetParentGeometryFile(spec.geometryname.c_str());
  anal.SetFileName(spec.outfilename.c_str());
  std::cout<<"Check for segmentation fault 3"<<std::endl;

  anal.SetCutTOF1SP1(spec.cut_TOF1_1sp);
  anal.SetCutTOF0SP1(spec.cut_TOF0_1sp);
  anal.SetCutTKU1track(spec.cut_TKU1track);
  anal.SetCutCHI2DOF(spec.cut_CHI2DOF);
  anal.SetCutTKURadial(spec.cut_TKURadial);
  anal.SetCutDiffuserRadial(spec.cut_DiffuserRadial);
  anal.SetCutFiducial(spec.cut_Fiducial);
  anal.SetCutTOF01(spec.cut_TOF01);
  anal.SetCutTOF01ext(spec.cut_TOF01ext);
  anal.SetCutAbsorberMomentum(spec.cut_AbsorberMomentum);
  anal.SetCutExtrudeTKUTOF0(spec.cut_ExtrudeTKUTOF0);
  std::cout<<"Check for segmentation fault 4"<<std::endl;
/*  
  void* dir = gSystem->OpenDirectory(spec.trainname.c_str());
  std::cout<<"Main Training File Directory: "<< dir <<std::endl;
  const char *ent;
  TString fn = spec.trainname.c_str();
  std::cout<<"Check for segmentation fault 5"<<std::endl;
  std::cout<<"Main Training File: "<< fn <<std::endl;
  while (ent = gSystem->GetDirEntry(dir)) {
      fn = fn.Append(ent);
      std::cout<<"Loop Main Training File: "<< fn <<std::endl;
          if (fn.Contains("ZeroAbs")) {
          if (fn.Contains(".root")) {
	          anal.GetRefTree()->Add(fn);
	  }
      }
      fn = spec.trainname.c_str();
  }
  gSystem->FreeDirectory(dir);

  dir = gSystem->OpenDirectory(spec.dataname.c_str());
  std::cout<<"Main Data File Directory: "<< dir <<std::endl;
  fn = spec.dataname.c_str();
  std::cout<<"Check for segmentation fault 6"<<std::endl;
  std::cout<<"Main Data File: "<< fn <<std::endl;
  while (ent = gSystem->GetDirEntry(dir)) {
      fn = fn.Append(ent);
      std::cout<<"Loop Main Data File: "<< fn <<std::endl;
          if (fn.Contains("LiH")) {
          if (fn.Contains(".root")) {
	          anal.GetTree()->Add(fn);
	  }
      }
      fn = spec.dataname.c_str();
  }
  std::cout<<"Check for segmentation fault 7"<<std::endl;

  gSystem->FreeDirectory(dir);

  dir = gSystem->OpenDirectory(spec.dataname.c_str());
  fn = spec.dataname.c_str();
  std::cout<<"MC Data File: "<< fn <<std::endl;
  while (ent = gSystem->GetDirEntry(dir)) {
      fn = fn.Append(ent);
      std::cout<<"Loop MC Data File: "<< fn <<std::endl;
          if (fn.Contains("LiH")) {
          if (fn.Contains(".root")) {

	          anal.GetMCTree()->Add(fn);
	  }
      }
      fn = spec.dataname.c_str();
  }
  std::cout<<"Check for segmentation fault 8"<<std::endl;

  gSystem->FreeDirectory(dir);
  
  dir = gSystem->OpenDirectory(spec.trainname.c_str());
  fn = spec.trainname.c_str();
  std::cout<<"MC Training File: "<< fn <<std::endl;
  while (ent = gSystem->GetDirEntry(dir)) {
      fn = fn.Append(ent);
      std::cout<<"Loop MC Training File: "<< fn <<std::endl;
          if (fn.Contains("ZeroAbs")) {
          if (fn.Contains(".root")) {
	          anal.GetMCEmptyTree()->Add(fn);
	  }
      }
      fn = spec.trainname.c_str();
  }
  std::cout<<"Check for segmentation fault 10"<<std::endl;

  gSystem->FreeDirectory(dir);

*/
  anal.GetRefTree()->Add(spec.trainname.c_str());

  anal.GetTree()->Add(spec.dataname.c_str());
  
  for(std::map<std::string, double>::iterator it=spec.sys.begin();
      it != spec.sys.end(); ++it){
    anal.SetSysOffset(it->first, it->second);
    std::cout<<"Adding systematic effect "<<it->first
	     <<" with value "<<it->second<<std::endl;
  }
  std::cout<<"Executing MCSAnalysis\n";
  anal.Execute(spec.mode);

  anal.Write();

  return 0;

}

#else
int main(void) {
  fprintf(stderr, "XPath support not compiled in\n");
  exit(1);
}
#endif
