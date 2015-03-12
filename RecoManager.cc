///
/// \file DCRecoManager.cc
/// Pulse Reconstruction Manager 
///
/// \author Pau Novella - CIEMAT 2008
///


#include<RecoManager.h>

using namespace std;

//************************************************************
RecoManager::RecoManager(){
//************************************************************
  
  
  _algo = "slidingWindow";
  
  _algos["simpleWindow"] = new simpleWindow(true);

  _algos["slidingWindow"] = new slidingWindow(true);

  _algos["maxWindow"] = new maxWindow(true);
  
  _algos["peakWindow"] = new peakWindow(true);

  _algos["biPoPeakWindow"] = new biPolarPeakWindow(true);
  
  _selectedalgo = _algos[_algo];
}



//************************************************************
RecoManager::~RecoManager(){
//************************************************************
  
  map<std::string,Algo*>::iterator iter;   

  for( iter = _algos.begin(); iter != _algos.end(); iter++ ) {
    delete iter->second;
  }
 
}

//*****************************************************************
void RecoManager::setChargeAna(string name){
//*****************************************************************
  
  this->getAlgo()->setChargeAna(name);

}


//*****************************************************************
void RecoManager::config(std::string algo, bool pedEvt){
//*****************************************************************
  
 
  this->checkAlgo(algo);
  _algo=algo;
  _selectedalgo = _algos[_algo];
  
  this->getAlgo()->computeProfPedestal(pedEvt);

}


//*****************************************************************
void RecoManager::setPEDwindow(size_t tmin, size_t tmax){
//*****************************************************************
  
  this->getAlgo()->setPEDwindow(tmin,tmax);

}

//*****************************************************************
void RecoManager::setPedAlgo(string algo){
//*****************************************************************
  
  this->getAlgo()->setPedAlgo(algo);

}

//*****************************************************************
void RecoManager::setQmin(double q){
//*****************************************************************
   
  //_Qmin=q;
  map<std::string,Algo*>::iterator iter;   
  for( iter = _algos.begin(); iter != _algos.end(); iter++ ) {
  iter->second->setQmin(q);
  }
  
}

//*****************************************************************
void RecoManager::setImin(double I){
//*****************************************************************
  
  //_Imin = I;
  map<std::string,Algo*>::iterator iter;   
  for( iter = _algos.begin(); iter != _algos.end(); iter++ ) {
  iter->second->setImin(I);
  }
  
}

//*****************************************************************
void RecoManager::setWmin(size_t w){
//*****************************************************************
  
  //_Imin = I;
  map<std::string,Algo*>::iterator iter;   
  for( iter = _algos.begin(); iter != _algos.end(); iter++ ) {
  iter->second->setWmin(w);}
  
}

//*****************************************************************
void RecoManager::setAWRmax(double max){
//*****************************************************************
  
  map<std::string,Algo*>::iterator iter;   
  for( iter = _algos.begin(); iter != _algos.end(); iter++ ) {
  iter->second->setAWRmax(max);}
  
}

//*****************************************************************
bool RecoManager::isGoodPulse(){
//*****************************************************************
  
  return (this->getNpeaks());
  
}

//*****************************************************************
void RecoManager::setFilter(size_t nfilter){
//*****************************************************************
   
  map<std::string,Algo*>::iterator iter;   
  for( iter = _algos.begin(); iter != _algos.end(); iter++ ) {
    iter->second->setFilter(nfilter);
  }
  
}

//*****************************************************************
void RecoManager::setExtPedestal(double ped, double rms){
//*****************************************************************
  
  map<std::string,Algo*>::iterator iter;   
  for( iter = _algos.begin(); iter != _algos.end(); iter++ ) {
    iter->second->setExtPedestal(ped,rms);
  }
  
  //this->getAlgo()->setPedestal(ped,rms);

}

//*****************************************************************
bool RecoManager::sameChExtPedRMS(){
//*****************************************************************
  
  return this->getAlgo()->sameChExtPedRMS();

}

//*****************************************************************
bool RecoManager::useExtPedMean(){
//*****************************************************************
    
  return this->getAlgo()->useExtPedMean();

}

//*****************************************************************
bool RecoManager::useExtPedRMS(){
//*****************************************************************
  
  return this->getAlgo()->useExtPedRMS();
  
}

//*****************************************************************
void RecoManager::useExtPedMean(bool ok){
//*****************************************************************
  
  map<std::string,Algo*>::iterator iter;   
  for( iter = _algos.begin(); iter != _algos.end(); iter++ ) {
    iter->second->useExtPedMean(ok);
  }
  
}

//*****************************************************************
void RecoManager::useExtPedRMS(bool ok){
//*****************************************************************
  
  map<std::string,Algo*>::iterator iter;   
  for( iter = _algos.begin(); iter != _algos.end(); iter++ ) {
    iter->second->useExtPedRMS(ok);
  }
  
}

//*****************************************************************
double RecoManager::computePED(const vector<unsigned short>& p){
//*****************************************************************
  
  RPPulse rp = RPPulse(p); 

  return this->getAlgo()->computePED(rp);

}

//*****************************************************************
int RecoManager::getPedStatus(){
//*****************************************************************
  
  return this->getAlgo()->getPedStatus();

}

//*****************************************************************
bool RecoManager::isGoodPedMean(){
//*****************************************************************
  
  return this->getAlgo()->isGoodPedMean();

}



//*****************************************************************
bool RecoManager::isGoodPedRMS(){
//*****************************************************************
  
  return this->getAlgo()->isGoodPedRMS();

}

//*****************************************************************
bool RecoManager::isGoodExtPedMean(){
//*****************************************************************
  
  return this->getAlgo()->isGoodExtPedMean();

}


//*****************************************************************
bool RecoManager::isGoodExtPedRMS(){
//*****************************************************************
  
  return this->getAlgo()->isGoodExtPedRMS();

}

//*****************************************************************
bool RecoManager::isBiPolar(){
//*****************************************************************
  
  return this->getAlgo()->isBiPolar();

}

//*****************************************************************
bool RecoManager::isGoodPedestal(){
//*****************************************************************
  
  return this->getAlgo()->isGoodPedestal();

}

//*****************************************************************
bool RecoManager::isGoodExtPedestal(){
//*****************************************************************
  
  return this->getAlgo()->isGoodExtPedestal();

}


//*****************************************************************
void RecoManager::reco(const vector<unsigned short>& p){
//*****************************************************************
  
  RPPulse rp = RPPulse(p); 

  this->getAlgo()->reco(rp);

}

//******************************************************************
double RecoManager::computeQ(const vector<unsigned short>& p){ 
//******************************************************************
  
  RPPulse rp = RPPulse(p); 

  return this->getAlgo()->computeQ(rp);
 
}

//******************************************************************
int RecoManager::computeT0(const vector<unsigned short>& p){ 
//******************************************************************
  
  RPPulse rp = RPPulse(p); 

  return this->getAlgo()->computeT0(rp);
  
}

//******************************************************************
int RecoManager::computeT1(const vector<unsigned short>& p){ 
//******************************************************************

  RPPulse rp = RPPulse(p); 

  return this->getAlgo()->computeT1(rp);
  
}

//******************************************************************
int RecoManager::computeTstart(const vector<unsigned short>& p){ 
//******************************************************************

  RPPulse rp = RPPulse(p); 

  return this->getAlgo()->computeTstart(rp);
  
}

//******************************************************************
int RecoManager::computeTend(const vector<unsigned short>& p){ 
//******************************************************************

  RPPulse rp = RPPulse(p); 

  return this->getAlgo()->computeTend(rp);
  
}


//******************************************************************
int RecoManager::computeTmax(const vector<unsigned short>& p){ 
//******************************************************************
  
  RPPulse rp = RPPulse(p); 

  return this->getAlgo()->computeTmax(rp);
  
}

//******************************************************************
int RecoManager::computeTrise(const vector<unsigned short>& p){ 
//******************************************************************

  RPPulse rp = RPPulse(p); 

  return this->getAlgo()->computeTrise(rp);
  
}

//******************************************************************
int RecoManager::computeTfall(const vector<unsigned short>& p){ 
//******************************************************************
  
  RPPulse rp = RPPulse(p);   
  
  return this->getAlgo()->computeTfall(rp);
  
}


//******************************************************************
//void RecoManager::setCalibOn(bool ok){ 
//******************************************************************
  
//  return this->getAlgo()->setCalibOn(ok);
  
//}

//******************************************************************
//void RecoManager::setTcalFile(string file){ 
//******************************************************************
  
//  return this->getAlgo()->setTcalFile(file);
  
//}

//******************************************************************
//void RecoManager::setQcalFile(string file){ 
//******************************************************************
  
//  return this->getAlgo()->setQcalFile(file);
  
//}


//******************************************************************
void RecoManager::setNsigOverPed(double nsigma){ 
//******************************************************************
  
  return this->getAlgo()->setNsigOverPed(nsigma);
  
}

//******************************************************************
void RecoManager::setMaxPEDMean(double max){
//******************************************************************
  
  this->getAlgo()->setMaxPEDMean(max);

}

//******************************************************************
void RecoManager::setMinPEDMean(double min){
//******************************************************************
  
  this->getAlgo()->setMinPEDMean(min);

}

//******************************************************************
void RecoManager::setMaxPEDRMS(double max){
//******************************************************************
  
  this->getAlgo()->setMaxPEDRMS(max);

}

//******************************************************************
void RecoManager::setMinPedRMS(double min){
//******************************************************************
  
  this->getAlgo()->setMinPedRMS(min);

}


//******************************************************************
void RecoManager::computeProfPedestal(bool ok){ 
//******************************************************************
  
  /*
    Algorithm will compute pedestal for current waveform
    instead of using predefined values
   */

  map<std::string,Algo*>::iterator iter;   
  
  for( iter = _algos.begin(); iter != _algos.end(); iter++ ) {
    iter->second->computeProfPedestal(ok);}

}


//******************************************************************
double RecoManager::getPedMean(){ 
//******************************************************************
  
  return this->getAlgo()->getPedMean();
  
}

//******************************************************************
double RecoManager::getPedRMS(){ 
//******************************************************************
  
  return this->getAlgo()->getPedRMS();
  
}

//******************************************************************
double RecoManager::getChPedMean(){ 
//******************************************************************
  
  return this->getAlgo()->getChPedMean();
  
}

//******************************************************************
double RecoManager::getChPedRMS(){ 
//******************************************************************
  
  return this->getAlgo()->getChPedRMS();
  
}

//******************************************************************
double RecoManager::getExtPedMean(){ 
//******************************************************************
  
  return this->getAlgo()->getExtPedMean();
  
}

//******************************************************************
double RecoManager::getExtPedRMS(){ 
//******************************************************************
  
  return this->getAlgo()->getExtPedRMS();
  
}

//******************************************************************
int RecoManager::getT0(){ 
//******************************************************************
  
  return this->getAlgo()->getT0();
  
}

//******************************************************************
int RecoManager::getT1(){ 
//******************************************************************
  
  return this->getAlgo()->getT1();
  
}

//******************************************************************
int RecoManager::getTstart(){ 
//******************************************************************
  
  return this->getAlgo()->getTstart();
  
}

//******************************************************************
int RecoManager::getTstartError(){ 
//******************************************************************
  
  return this->getAlgo()->getTstartError();
  
}

//******************************************************************
int RecoManager::getTend(){ 
//******************************************************************
  
  return this->getAlgo()->getTend();
  
}

//******************************************************************
int RecoManager::getTrise(){ 
//******************************************************************
  
  return this->getAlgo()->getTrise();
  
}

//******************************************************************
int RecoManager::getTfall(){ 
//******************************************************************
  
  return this->getAlgo()->getTfall();
  
}

//******************************************************************
int RecoManager::getTmax(){ 
//******************************************************************
  
  return this->getAlgo()->getTmax();
  
}

//******************************************************************
double RecoManager::getQ(){ 
//******************************************************************
  
  return this->getAlgo()->getQ();
  
}

//******************************************************************
double RecoManager::getQerror(){ 
//******************************************************************
  
  return this->getAlgo()->getQerror();
  
}

//******************************************************************
double RecoManager::getImax(){ 
//******************************************************************
  
  return this->getAlgo()->getImax();
  
}

//******************************************************************
size_t RecoManager::getNpeaks(){ 
//******************************************************************
  
  return this->getAlgo()->getPeakQs().size();
  
}

//******************************************************************
vector<int> RecoManager::getPeakIntTs(){ 
//******************************************************************
  
  return this->getAlgo()->getPeakT0s();
  
}

//******************************************************************
vector<int> RecoManager::getPeakIntTends(){ 
//******************************************************************
  
  return this->getAlgo()->getPeakT1s();
  
}

//******************************************************************
vector<double> RecoManager::getPeakQs(){ 
//******************************************************************
  
  return this->getAlgo()->getPeakQs();
  
}

//******************************************************************
vector<int> RecoManager::getPeakTs(){ 
//******************************************************************
  
  return this->getAlgo()->getPeakTs();
  
}

//******************************************************************
vector<int> RecoManager::getPeakTends(){ 
//******************************************************************
  
  return this->getAlgo()->getPeakTends();
  
}

//******************************************************************
vector<int> RecoManager::getPeakTmaxs(){ 
//******************************************************************
  
  return this->getAlgo()->getPeakTmaxs();
  
}

//******************************************************************
vector<double> RecoManager::getPeakImaxs(){ 
//******************************************************************
  
  return this->getAlgo()->getPeakImaxs();
  
}

//******************************************************************
void RecoManager::reset(){ 
//******************************************************************
  
  return this->getAlgo()->reset();
  
}

//******************************************************************
size_t RecoManager::getPedWindowSize(){ 
//******************************************************************
  
  return this->getAlgo()->getPedWindowSize();
  
}

//******************************************************************
size_t RecoManager::getPedWindowMin(){ 
//******************************************************************
  
  return this->getAlgo()->getPedWindowMin();
  
}

//******************************************************************
size_t RecoManager::getPedWindowMax(){ 
//******************************************************************
  
  return this->getAlgo()->getPedWindowMax();
  
}


//*****************************************************************
bool RecoManager::find(std::string algo) const{
//*****************************************************************
  
   //! find an element in the store by name
         
  bool gotcha = false;
  
  map<std::string, Algo*>::const_iterator pi;
      
  for ( pi = _algos.begin(); 
	pi!= _algos.end(); ++pi)
    {
      if ( (pi->first) == algo) {
	gotcha = true;
	break;
      }
    } 

  return gotcha;
  
}

//*****************************************************************
void RecoManager::checkAlgo(std::string algo) const{
//*****************************************************************

  bool ok=this->find(algo);
  
  if (!ok){ 
    //Message::MSG("RecoPulse",kMFATAL,"Algorithm not found: "+algo);
  }

}
