#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

// #include "TMatrixD.h"

class Vars {
 public:
  double X;
  double Y;
  double Z;
  double dXdz;
  double dYdz;
  double px;
  double py;
  double pz;
  double TOF12;
  double TOF01;
  bool   isgood;
  int    pid;
  double TOFX;

  Vars operator+(const Vars& right);
  Vars operator-(const Vars& right);
  Vars operator*(const Vars& right);
  Vars operator/(const Vars& right);
  Vars operator*(const double right);
  Vars operator/(const double right);
  void Zero();
};

class Collection {
 public:
  Collection() {};
  ~Collection() {};

 private:
  std::vector<Vars> _Set;
  Vars _mean;
  // TMatrixD _cov;
  std::vector<std::vector<double> > _cov;
  double _emittance;
  double _mass;

 public:
  std::vector<Vars> Set() { return _Set; }
  int N() { return _Set.size();}
  Vars E(int i) {return _Set.at(i); }
  void append_instance(Vars instance) {
    _Set.push_back(instance);
  }
  void init_mean();
  void calc_mean();
  void init_cov();
  void calc_cov();
  double Determinant(std::vector<std::vector<double> > a,int n);
  void calc_emittance(int n, double mass);
  Vars mean() { 
    _Set.size()>0 ? calc_mean() : init_mean();
    return _mean;
  }
  std::vector<std::vector<double> > cov() {
    _Set.size()>1 ? calc_cov() : init_cov();
    return _cov;
  }
  double emittance(int n, double mass) {
    if ( _Set.size() > 1 ){
      calc_cov();
      calc_emittance(n, mass);
      return _emittance;
    }
    else {
      init_cov();
      return 0;
    }
  }
};

struct multiVars {
    Vars USTOF0;
    Vars USTOF1;
    Vars UScentre_absorber;
    Vars UStrackerUS11;
    Vars UStrackerUS53;
    Vars USend_of_DStracker;
    Vars DScentre_absorber;
};
