///
/// \file DCRecoPulsePed.cc
/// Pedestal analysis algorithm 
///
/// \author Pau Novella - CIEMAT 2008
///

#include<math.h>

#include<RecoPulsePed.h>
#include <algorithm>

using namespace std;

//*************************************************************
RecoPulsePed::RecoPulsePed(){
//*************************************************************
  
  //Message::SetLevelMSG("RecoPulsePed",Message::GetLevelMSG());

  //Message::MSG("RecoPulsePed",kMINFO,"<<Construction>>");

  this->config(5000,1,3,0,50);

  _pedAlgoisfull=false;

  _pedTol = 0.5; 

}

//*************************************************************
void RecoPulsePed::config(double maxM,double minM,double mRMS,
		     size_t tmin,size_t tmax){
//*************************************************************
  
  //this->setVlevel(vlevel);

  _mean=NOPED; _RMS=NOPED;

  _extMean=NOPED; _extRMS=NOPED;

  this->setMaxMean(maxM); 
  
  //Message::MSG("RecoPulsePed",kMINFO,
  //       Form("Maximum pedestal mean set to: %f",maxM));

  this->setMinMean(minM); 
  
  //Message::MSG("RecoPulsePed",kMINFO,
  //	       Form("Minimum pedestal mean set to: %f",minM));
  
  this->setMaxRMS(mRMS);
  
  //Message::MSG("RecoPulsePed",kMINFO,
  //	       Form("Maximum pedestal RMS set to: %f",mRMS));

  this->setWindow(tmin,tmax);
  
  this->setNsigOverPed(5);
    
}

//*************************************************************
bool RecoPulsePed::isGoodPedMean() const{
//*************************************************************
  
  return (_mean<_maxMean && _mean>_minMean);

}

//*************************************************************
bool RecoPulsePed::isGoodPedRMS() const{
//*************************************************************
  
  return (_RMS<_maxRMS);

}

//*************************************************************
bool RecoPulsePed::isGoodExtPedMean() const{
//*************************************************************
  
  return (_extMean<_maxMean && _extMean>_minMean);

}

//*************************************************************
bool RecoPulsePed::isGoodExtPedRMS() const{
//*************************************************************
  
  return (_extRMS<_maxRMS);

}

// //*************************************************************
// bool RecoPulsePed::sameChExtPeds() const{
// //*************************************************************
  
//   bool ok = true;

//   if ( _extMean < _mean - _mean * _pedTol ) ok = false;

//   else if ( _extMean >  _mean + _mean * _pedTol) ok = false;
  
//   return ok;

// }

//*************************************************************
bool RecoPulsePed::sameChExtPedMean() const{
//*************************************************************
  
  return ( fabs(_mean - _extMean)< _pedTol );

}

//*************************************************************
bool RecoPulsePed::sameChExtPedRMS() const{
//*************************************************************
  
  return ( fabs(_RMS - _extRMS)< _pedTol );

}

//*************************************************************
bool RecoPulsePed::isGoodPedestal() const{
//*************************************************************

  bool ok1,ok2,ok3,ok4;
  
  ok1 = (_mean!=NOPED && _RMS!=NOPED );
  
  //if (!ok1) Message::MSG("RecoPulsePed",kMINFO,"PEDESTAL not defined");

  ok2 = this->isGoodPedMean();

  //if (!ok2) Message::MSG("RecoPulsePed",kMINFO,
  //			 Form("PEDESTAL mean out of bounds: %f",_mean));

  ok3 = this->isGoodPedRMS();
  
//if (!ok3) Message::MSG("RecoPulsePed",kMINFO,
//			 Form("PEDESTAL RMS too high: %f",_RMS));
  
  ok4 = true;

  //if (this->isGoodExtPedestal()){ok4 = this->sameChExtPeds();}

  return (ok1 && ok2 && ok3 && ok4);

}

//*************************************************************
bool RecoPulsePed::isGoodExtPedestal() const{
//*************************************************************

  bool ok1,ok2,ok3;

  ok1 = (_extMean!=NOPED && _extRMS!=NOPED );
  
  //if (!ok1) Message::MSG("RecoPulsePed",kMINFO,"EXT PEDESTAL not defined");

  ok2 = this->isGoodExtPedMean();

  //if (!ok2) Message::MSG("RecoPulsePed",kMINFO,
  //			 Form("EXT PEDESTAL mean out of bounds: %f",_mean));

  ok3 = this->isGoodExtPedRMS();
  
  //if (!ok3) Message::MSG("RecoPulsePed",kMINFO,
  //			 Form("EXT PEDESTAL RMS too high: %f",_RMS));
  
  

  return (ok1 && ok2 && ok3);

}

//*************************************************************
int RecoPulsePed::getPedStatus() const{
//*************************************************************

  int pedstatus;

  if (_mean==NOPED || _RMS==NOPED ) pedstatus = 0;

  else if (this->isGoodPedestal()) pedstatus = 1;

  else pedstatus = -1;

  return pedstatus;

}

//*************************************************************
int RecoPulsePed::getExtPedStatus() const{
//*************************************************************

  int pedstatus;

  if (_extMean==NOPED || _extRMS==NOPED ) pedstatus = 0;

  else if (this->isGoodExtPedestal()) pedstatus = 1;

  else pedstatus = -1;

  return pedstatus;

}

//**************************************************************
bool RecoPulsePed::computePedestal(RPPulse& ipulse,
				   size_t tmin, size_t tmax){
//*************************************************************

  bool ok=false;

  if(!_pedAlgoisfull)
    
    ok = this->computeOffsetPedestal(ipulse,tmin,tmax);
  
  else 

    ok = this->computeFullPedestal(ipulse);
  
  //Message::MSG("RecoPulsePed",kMDEBUG,Form("Pedestal mean: %f",_mean));

  //Message::MSG("RecoPulsePed",kMDEBUG,Form("Pedestal RMS: %f",_RMS));
  
  return ok;

}

//**************************************************************
bool RecoPulsePed::computeOffsetPedestal(RPPulse& ipulse,
					 size_t tmin, size_t tmax){
//*************************************************************

  const vector<unsigned short>& pulse = ipulse.getProfile();

  if (tmin || tmax) this->setWindow(tmin,tmax);

  if (!_tmax) _tmax = pulse.size();  
  
  double pedInt=0;
  
  for (size_t i=_tmin; i<_tmax; i++){ pedInt += pulse[i]; }
    
  _mean = pedInt / (_tmax-_tmin) ;
  
  //_mean = round(_mean);

  double qdev = 0;

  for (size_t i=_tmin; i<_tmax; i++){ qdev+=pow(pulse[i]-_mean,2); }
  
  _RMS = sqrt(qdev / (_tmax-_tmin) );
  
  _wsize = _tmax - _tmin;
  
  //return this->isGoodPedestal();
  return true;

}

//**************************************************************
bool RecoPulsePed::computeFullPedestal(RPPulse& ipulse){
//*************************************************************
  
  vector<unsigned short> pulse = ipulse.getProfile();
  
  std::sort(pulse.begin(), pulse.end());
  
  bool ok=false; double mean = 0;
  
  while(!ok){
    
    mean = 0; 

    for(size_t i=0; i<pulse.size();i++){mean+=pulse[i]*1./pulse.size();}

    double posMaxDev = fabs(mean-pulse[0]);

    double negMaxDev = fabs(mean-pulse[pulse.size()-1]);
    
    ok = ( fabs(posMaxDev - negMaxDev) < 1);
    
    if (negMaxDev < posMaxDev) pulse.erase(pulse.begin());

    else pulse.erase(pulse.end()-1);

  } //end of while
  
  //_mean = round(mean);
  
  _mean = mean;
  
  double qdev = 0;

  for (size_t i=0; i<pulse.size(); i++){ qdev+=pow(pulse[i]-_mean,2); }

  _RMS = sqrt(qdev *1. / (pulse.size()) );
  
  _wsize = pulse.size();
  
  double fullPedRMS = _RMS;

  if (fullPedRMS>_maxRMS){ // compute ped from waveform offset

  this->computeOffsetPedestal(ipulse,_tmin,_tmax);}

  if (fullPedRMS < _RMS) _RMS = fullPedRMS;
  
  return true;
}


//**************************************************************
void RecoPulsePed::setWindow(size_t tmin,size_t tmax){
//*************************************************************

  this->checkWindow(tmin,tmax);

  _tmin=tmin; _tmax=tmax;
  
}

//**************************************************************
void RecoPulsePed::checkWindow(size_t tmin,size_t tmax) const{
//*************************************************************
  
  if (tmin>tmax && tmax!=0){
    
    //Message::MSG("RecoPulsePed",kMFATAL,"PEDESTAL window not well defined");
    
  }
  
}

//*****************************************************************
double RecoPulsePed::computeThreshold(double mean, double RMS){
//*****************************************************************
  
  //_threshold = mean - _nsigma * RMS; 

  double threshold = mean - _nsigma * RMS; 

  return threshold;
  
}

// //*****************************************************************
// bool RecoPulsePed::isOverThreshold(double amp){
// //*****************************************************************
  
//   this->checkStatus();

//   bool ok=false;
  
//   ok = (this->getThreshold() > amp);

//   return ok;
  
// }

//*****************************************************************
void RecoPulsePed::checkStatus() const{
//*****************************************************************
  
  if (!this->status()){
      
    //Message::MSG("RecoPulsePed",kMFATAL,"PEDESTAL not defined");
  }

}
