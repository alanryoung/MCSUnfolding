
#include "Collection.h"

Vars Vars::operator+(const Vars& right){
  Vars result;
  result.X = X + right.X;
  result.Y = Y + right.Y;
  result.Z = Z + right.Z;
  result.dXdz = dXdz + right.dXdz;
  result.dYdz = dYdz + right.dYdz;
  result.px   = px + right.px;
  result.py   = py + right.py;
  result.pz   = pz + right.pz;
  result.TOF12 = TOF12 + right.TOF12;
  result.TOF01 = TOF01 + right.TOF01;
  result.isgood = isgood || right.isgood;
  result.pid    = abs(pid) == 13 && abs(right.pid) == 13 ? pid : 0;
  return result;
}

Vars Vars::operator-(const Vars& right){
  Vars result;
  result.X = X - right.X;
  result.Y = Y - right.Y;
  result.Z = Z - right.Z;
  result.dXdz = dXdz - right.dXdz;
  result.dYdz = dYdz - right.dYdz;
  result.px   = px - right.px;
  result.py   = py - right.py;
  result.pz   = pz - right.pz;
  result.TOF12 = TOF12 - right.TOF12;
  result.TOF01 = TOF01 - right.TOF01;
  result.isgood = isgood || right.isgood;
  result.pid    = abs(pid) == 13 && abs(right.pid) == 13 ? pid : 0;
  return result;
}

Vars Vars::operator*(const Vars& right){
  Vars result;
  result.X = X * right.X;
  result.Y = Y * right.Y;
  result.Z = Z * right.Z;
  result.dXdz = dXdz * right.dXdz;
  result.dYdz = dYdz * right.dYdz;
  result.px   = px * right.px;
  result.py   = py * right.py;
  result.pz   = pz * right.pz;
  result.TOF12 = TOF12 * right.TOF12;
  result.TOF01 = TOF01 * right.TOF01;
  result.isgood = isgood && right.isgood;
  result.pid    = abs(pid) == 13 && abs(right.pid) == 13 ? pid : 0;
  return result;
}


Vars Vars::operator*(const double right){
  Vars result;
  result.X = X * right;
  result.Y = Y * right;
  result.Z = Z * right;
  result.dXdz = dXdz * right;
  result.dYdz = dYdz * right;
  result.px   = px * right;
  result.py   = py * right;
  result.pz   = pz * right;
  result.TOF12 = TOF12 * right;
  result.TOF01 = TOF01 * right;
  result.isgood = isgood;
  result.pid    = pid;
  return result;
}

Vars Vars::operator/(const Vars& right){
  Vars result;
  result.X = X / right.X;
  result.Y = Y / right.Y;
  result.Z = Z / right.Z;
  result.dXdz = dXdz / right.dXdz;
  result.dYdz = dYdz / right.dYdz;
  result.px   = px / right.px;
  result.py   = py / right.py;
  result.pz   = pz   / right.pz;
  result.TOF12 = TOF12 / right.TOF12;
  result.TOF01 = TOF01 / right.TOF01; 
  result.isgood = isgood && right.isgood;
  result.pid    = abs(pid) == 13 && abs(right.pid) == 13 ? pid : 0;
  return result;
}


Vars Vars::operator/(const double right){
  Vars result;
  result.X = X / right;
  result.Y = Y / right;
  result.Z = Z / right;
  result.dXdz = dXdz / right;
  result.dYdz = dYdz / right;
  result.px   = px / right;
  result.py   = py / right;
  result.pz   = pz   / right; 
  result.TOF12 = TOF12 / right;
  result.TOF01 = TOF01 / right;
  result.isgood = isgood;
  result.pid    = pid;
  return result;
}

void Vars::Zero(){
  X = 0.0;
  Y = 0.0;
  Z = 0.0;
  dXdz = 0.0;
  dYdz = 0.0;
  px   = 0.0;
  py   = 0.0;
  pz   = 0.0;
  TOF12 = 0.0;
  TOF01 = 0.0;
  isgood = false;
  pid   = 0;
}

void Collection::init_mean() {
  _mean.Zero();
}

void Collection::calc_mean() { 
  init_mean();
  std::vector<Vars>::iterator xit = _Set.begin();
  double N = double(_Set.size());
  do {
    if (!(*xit).isgood) continue;
    _mean = _mean + (*xit);
    xit++;
  } while(xit != _Set.end());
  _mean = _mean / N;
}

void Collection::init_cov() {
  // _cov.ResizeTo(a);
  for ( int i=0; i<4; i++){
    std::vector<double> tmp;
    for ( int j=0; j<4; j++)
      tmp.push_back(0.0);
    _cov.push_back(tmp);
  }
}

void Collection::calc_cov() {
  init_cov();
  Vars diag;
  diag.Zero();
  
  std::vector<Vars>::iterator xit = _Set.begin();
  double N = double(_Set.size());
  do {
    if (!(*xit).isgood) continue;
    diag = diag + ((*xit) - _mean)*((*xit) - _mean);
    _cov[0][1] += ((*xit).X - _mean.X)*((*xit).dXdz - _mean.dXdz);
    _cov[0][2] += ((*xit).X - _mean.X)*((*xit).Y - _mean.Y);
    _cov[0][3] += ((*xit).X - _mean.X)*((*xit).dYdz - _mean.dYdz);
    _cov[1][2] += ((*xit).dXdz - _mean.dXdz)*((*xit).Y - _mean.Y);
    _cov[1][3] += ((*xit).dXdz - _mean.dXdz)*((*xit).dYdz - _mean.dYdz);
    _cov[2][3] += ((*xit).Y - _mean.Y)*((*xit).dYdz - _mean.dYdz);
    xit++;
  } while ( xit != _Set.end() );
  diag = diag / double(N - 1);
  _cov[0][0] = diag.X;
  _cov[1][1] = diag.dXdz;
  _cov[2][2] = diag.Y;
  _cov[3][3] = diag.dYdz;
  _cov[0][1] /= double(N - 1);
  _cov[0][2] /= double(N - 1);
  _cov[0][3] /= double(N - 1);
  _cov[1][0] = _cov[0][1];
  _cov[1][2] /= double(N - 1);
  _cov[1][3] /= double(N - 1);
  _cov[2][0] = _cov[0][2];
  _cov[2][1] = _cov[1][2];
  _cov[2][3] /= double(N - 1);
  _cov[3][0] = _cov[0][3];
  _cov[3][1] = _cov[1][3];
  _cov[3][2] = _cov[2][3];
  
}

double Collection::Determinant(std::vector<std::vector<double> > a,int n)
{
  int i,j,j1,j2;
  double det = 0;
  std::vector<std::vector<double> > m;
  
  for (int i=0; i<n; i++){
    std::vector<double> temp;
    for (int j=0; j<4; j++)
      temp.push_back(0.0);
    m.push_back(temp);
  }

  if (n < 1) { /* Error */

  } else if (n == 1) { /* Shouldn't get used */
    det = a[0][0];
  } else if (n == 2) {
    det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
  } else {
    det = 0;
    for (j1=0;j1<n;j1++) {
      for (i=1;i<n;i++) {
	j2 = 0;
	for (j=0;j<n;j++) {
	  if (j == j1)
	    continue;
	  m[i-1][j2] = a[i][j];
	  j2++;
	}
      }
      det += pow(-1.0,1.0+j1+1.0) * a[0][j1] * Determinant(m,n-1);
    }
  }
  return(det);
}


void Collection::calc_emittance(int n, double mass){
  _emittance=0;
  _mass=mass;
  double det = Determinant(_cov, n);
  _emittance = pow(det, 1./double(n));
  _emittance /= mass;
  // return _emittance;
}

void Collection::save_collection(std::string filename){
   std::ofstream savefile;
   savefile.open(filename);
   for(int i=0; i<_Set.size(); i++){
     Vars xit =_Set.at(i);
     savefile << xit.X << " ";
     savefile << xit.Y << " ";;
     savefile << xit.Z << " ";;
     savefile << xit.dXdz << " ";;
     savefile << xit.dYdz << " ";;
     savefile << xit.px << " ";;
     savefile << xit.py << " ";;
     savefile << xit.pz << " ";;
     savefile << xit.TOF12 << " ";;
     savefile << xit.TOF01 << " ";;
     savefile << xit.isgood << " ";;
     savefile << xit.pid << " ";;
     savefile << xit.TOFX << "\n";;
   }
   savefile.close();
} 

 void Collection::ReturnWords(std::string inputstring,std::vector<std::string> &v){
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

void Collection::load_collection(std::string filename){
   std::ifstream savefile;
   std::vector<std::string> results;
   std::cerr<<"Running Load Collection"<<std::endl;
   Vars tempvars;
   savefile.open(filename);
   std::string line;
   if(savefile.is_open()){
     std::cerr<<"File is open for reading"<<std::endl;
     while(getline(savefile,line)){
//       std::cerr<<"Line retrieved for parsing"<<std::endl;
       if (line.find(' ') != std::string::npos){
//         std::cerr<<"Line contains spaces"<<std::endl;
//         std::cerr<<"Line: "<<line<<std::endl;
         ReturnWords(line, results);

//         std::istringstream(line);
/*         std::cerr<<"Attempting to create  vector containing data"<<std::endl;
         std::vector<std::string> results((std::istream_iterator<std::string>(line)),
                                          std::istream_iterator<std::string>());
*/
//         std::cerr<<"Vector created"<<std::endl;
/*         std::cerr <<"Size of returned vector: "<<results.size()<<std::endl;
         std::cerr<<"X: "<<results[0];
         std::cerr<<" Y: "<<results[1];
         std::cerr<<" Z: "<<results[2];
         std::cerr<<" dXdz: "<<results[3];
         std::cerr<<" dYdz: "<<results[4];
         std::cerr<<" px: "<<results[5];
         std::cerr<<" py: "<<results[6];
         std::cerr<<" pz: "<<results[7];
         std::cerr<<" TOF12: "<<results[8];
         std::cerr<<" TOF01: "<<results[9];
         std::cerr<<" isgood: "<<results[10];
         std::cerr<<" pid: "<<results[11];
         std::cerr<<" TOFX: "<<results[12]<<std::endl;
*/
         tempvars.X = stod(results[0]);
//         std::cerr<<"Successfully added X"<<std::endl;
         tempvars.Y = stod(results[1]);
         tempvars.Z = stod(results[2]);
         tempvars.dXdz = stod(results[3]);
         tempvars.dYdz = stod(results[4]);
         tempvars.px = stod(results[5]);
         tempvars.py = stod(results[6]);
         tempvars.pz = stod(results[7]);
         tempvars.TOF12 = stod(results[8]);
         tempvars.TOF01 = stod(results[9]);
         tempvars.isgood = stoi(results[10]);
         tempvars.pid = stoi(results[11]);
         tempvars.TOFX = stod(results[12]);
//         std::cerr<<"Line parsed and ready for input into Collection"<<std::endl;
         _Set.push_back(tempvars);
//         std::cerr<<"Data added to Collection"<<std::endl;
       }
     }
   }
   savefile.close();

}
