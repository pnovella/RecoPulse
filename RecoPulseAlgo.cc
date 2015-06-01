///
/// \file DCRecoPulseAlgo.cc
/// Pulse Reconstruction Algorithms 
///
/// \author Pau Novella - CIEMAT 2008
///


#include<RecoPulseAlgo.h>
#include<algorithm>

//============================ Algo ===============================//

using namespace std;

//******************************************************************
Algo::Algo(bool pedEvt){
//******************************************************************
  
  //Message::SetLevelMSG("RecoPulseAlgo",Message::GetLevelMSG());

  //Message::MSG("RecoPulseAlgo",gate::VERBOSE,"<<Construction>>");

  this->computeProfPedestal(pedEvt);
  
  ped = new RecoPulsePed();

  timer = new RecoPulseTime();

  charge = new RecoPulseQ();
  
  //calib = new RecoPulseCal(); 

  //qsat = new RecoPulseSat(vlevel);

  _Q=0; _QError = 0; 
  
  _Imax=0;

  _T0=NOTIME;_T1=NOTIME;

  _Tstart=NOTIME;_TstartError=0;_Tend=NOTIME; 

  _Tmax=NOTIME;_Trise=NOTIME;_Tfall=NOTIME;
  
  _Imin=5; //
  
  _Qmin=0; _Wmin=0; _AWRmax=1e9;

  _ped=NOPED;_pedRMS=NOPED;
  
  _pedStatus = 0;
  
  _pedRMS_min = 0.4; // default value !!!

  _nFilter=0; // no filter by default

  //_calibOn = false;
  
  _isBiPolar = false;

  _useExtPedRMS = false;
  
  _useExtPedMean = false;

}

//******************************************************************
void Algo::reset(){
//******************************************************************
  
  _Q=0; 

  _T0=NOTIME;_T1=NOTIME;
   
  _Tstart=NOTIME;_Tend=NOTIME;

  _Tmax=NOTIME;_Trise=NOTIME;_Tfall=NOTIME;
  
  _QError=0; _TstartError=0;

  _Qs.clear();_Ts.clear();_Tends.clear(); _Tmaxs.clear(); _Imaxs.clear(); 

  _T0s.clear(); _T1s.clear(); 

  _TsError.clear();_QsError.clear();
  
  _nFilter=0; // no filter by default
  
  _Imax=0;
  
  _isBiPolar = false;

  if (_pedEvt){ 
    
    this->getPedAna()->reset();

    _ped=NOPED; _pedRMS=NOPED;
    
    _pedStatus = 0;
  }

  this->getTimer()->reset();

  this->getChargeAna()->reset();

}

//******************************************************************
Algo::~Algo(){
//******************************************************************

  delete ped; delete timer; delete charge; //delete calib;

}

//******************************************************************
void Algo::setPEDwindow(size_t tmin, size_t tmax){
//******************************************************************

  this->getPedAna()->setWindow(tmin,tmax);

}

//******************************************************************
void Algo::setMaxPEDMean(double max){
//******************************************************************

  this->getPedAna()->setMaxMean(max);

}

//******************************************************************
void Algo::setMinPEDMean(double min){
//******************************************************************

  this->getPedAna()->setMinMean(min);

}

//******************************************************************
void Algo::setMaxPEDRMS(double max){
//******************************************************************

  this->getPedAna()->setMaxRMS(max);

}

//******************************************************************
//void Algo::setMinPedWidth(double min){
//******************************************************************

//  this->getPedAna()->setMinPedWidth(min);

//}


//******************************************************************
void Algo::setNsigOverPed(double thr){ 
//******************************************************************
  
  this->getPedAna()->setNsigOverPed(thr);
  
  //double threshold = this->getPedAna()->getThreshold();/??????????
  
  //this->getTimer()->getTimeDef()->setThreshold(threshold);/???????????

}


//******************************************************************
size_t Algo::getPedWindowSize(){
//******************************************************************

  return this->getPedAna()->getWindowSize();

}

//******************************************************************
void Algo::setPedAlgo(string algo){
//******************************************************************

  this->getPedAna()->setPedAlgo(algo);

}

//******************************************************************
void Algo::setTstartAlgo(string algo){
//******************************************************************

  this->getTimer()->getTimeDef()->setTstartAlgo(algo);

}

//******************************************************************
//void Algo::setQcalFile(string file){
//******************************************************************
  
//  this->getCalibrator()->readCalFile(file,"Q");

//}

//******************************************************************
//void Algo::setTcalFile(string file){
//******************************************************************
  
//  this->getCalibrator()->readCalFile(file,"T");

//}

//******************************************************************
//void Algo::calibrateQ(double rawq, double& q, double& sq){
//******************************************************************
  
//  this->getCalibrator()->calibrateQ(rawq,q,sq);
  

//}

//******************************************************************
//void Algo::calibrateT(double rawq, int& t, int& st){
//******************************************************************
  
//  this->getCalibrator()->calibrateT(rawq,t,st);
  
//}

//******************************************************************
double Algo::computeImax(RPPulse& ipulse){
//******************************************************************
  
  if (_Tmax==NOTIME) this->computeTmax(ipulse);
  
  if (_Tmax==NOTIME) return 0.0;

  const vector<unsigned short>& pulse = ipulse.getProfile();
  
  _Imax = _ped - pulse[(int) _Tmax];
  
  for (size_t i=0; i<_Tmaxs.size();i++){ 
    
    _Imaxs.push_back(_ped - pulse[(int) _Tmaxs[i]]);}

  return _Imax;
  
}

//******************************************************************
void Algo::extendWindow(RPPulse& ipulse){
//******************************************************************
  
  /* extend window on both right and left sides
 
  until amplitude goes below threshold */

  //size_t nsigma = abs( this->getPedAna()->getNsigOverPed() );
  
  const vector<unsigned short>& pulse = ipulse.getProfile();

  while (true){
    
    if (_T1>= (int) pulse.size()-1 || _T1==NOTIME) break;
    
    if (ipulse.isIntegrated(_T1)) break;

    double Qsample = (round(_ped) - pulse[_T1]); 
    
    if (Qsample<=0) break;

    //if (Qsample<1) break;
    
    double Q_ped = 1. * _pedRMS * 1.;  
    
    if (Qsample > Q_ped){ _T1++; _Q += (_ped - pulse[_T1]); }
      
    else break;
      
  }
  
  while (true){
    
    if (_T0>= (int)pulse.size()-1 || _T0<=0) break;
    
    if (ipulse.isIntegrated(--_T0)) break;

    double Qsample = (round(_ped) - pulse[_T0]); 
    
    if (Qsample<=0) break;

    //if (Qsample<1) break;
    
    double Q_ped = 1. * _pedRMS * 1.;  
    
    if (Qsample > Q_ped){ _Q += (_ped - pulse[_T0]); }
      
    else {_T0++; break;}

  }
  

}

//******************************************************************
void Algo::reco(RPPulse& ipulse){
//******************************************************************

  if (_pedEvt){ _ped = this->computePED(ipulse);}

  _Q = this->computeQ(ipulse);
  
  _T0 = this->computeT0(ipulse);

  _T1 = this->computeT1(ipulse);

  _Tstart = this->computeTstart(ipulse);

  _Tend = this->computeTend(ipulse);

  _Tmax = this->computeTmax(ipulse);
  
  _Trise = this->computeTrise(ipulse);

  _Tfall = this->computeTfall(ipulse);
  
  this->applyUserCuts(ipulse);

  //if (_calibOn && _Q) this->calibrate();
    
    

}

// //******************************************************************
// void Algo::calibrate(){
// //******************************************************************
  
//   double rawq = _Q;

//   this->calibrateQ(rawq,_Q,_QError); // TO REVIEW!!!!
  
//   this->calibrateT(rawq,_Tstart,_TstartError);
  
//   for (size_t i = 0; i<_Qs.size();i++){

//     rawq = _Qs[i]; 

//     _QsError.push_back(0);_TsError.push_back(0);

//     this->calibrateQ(rawq,_Qs[i],_QsError[i]);

//     this->calibrateT(rawq,_Ts[i],_TsError[i]);

//   }

// }

//******************************************************************
double Algo::getQmin(){
//******************************************************************
  
  /*
    Compute minimum charge that can be distinguished 
    from pedestal integration
   */
  
 
  bool status = (_T0 == NOTIME || _T1 == NOTIME || _T0 > _T1);

  gate::Assert(status == false,__FILE__,__LINE__,
	       gate::internal_logic("Q window not well defined!"));

  int wsize = _T1 - _T0; 

  double nsigma = abs( this->getPedAna()->getNsigOverPed() );
  
  double noise = max(_pedRMS,_pedRMS_min); 

  double Q_ped = nsigma * noise * sqrt((double)wsize);  
  
  double Qmin=Q_ped;

  return Qmin;
  

}

//******************************************************************
void  Algo::applyUserCuts(RPPulse& ipulse){
//******************************************************************
  
  bool ok = true;
  
  if (!_Qmin && !_Imin && !_Wmin) return; // no cuts
  
  //Message::MSG("RecoPulseAlgo",gate::DUMP,"Applying user cuts...");

  vector<double> Qs;
  
  vector<int> Ts,Tends,T0s,T1s;
  
  vector<int> Tmaxs; vector<double> Imaxs; 
  
  //double thr = this->getPedAna()->getThreshold();
  
  for (size_t i=0; i<_Qs.size();i++){
    
    bool ok1 = true;

    double Twidth;
    
    if (fabs(_Qs[i])<_Qmin) ok1 = false; 
    
    else if (fabs(_Imaxs[i])<_Imin) ok1 = false; 

    //else if (ipulse.getPeakWidth(_ped,thr,_T0s[i],_T1s[i])<_Wmin) ok1 = false;
    
    else if (ipulse.getContSamples(_ped,_T0s[i],_T1s[i])<_Wmin) ok1 = false; 
    
    else if ((Twidth = (_Tends[i]-_Ts[i])) && _Imaxs[i]/Twidth > _AWRmax) ok1 = false;

    if (!ok1){ 
      //Message::MSG("RecoPulseAlgo",gate::VERBOSE,
		   
      //	   Form("Bad quality peak: [%i,%i] ns",_Ts[i],_Tends[i]));
      ok = false;
      continue;
    }
        
    Qs.push_back(_Qs[i]); Ts.push_back(_Ts[i]); Tends.push_back(_Tends[i]); 
    
    Tmaxs.push_back(_Tmaxs[i]); Imaxs.push_back(_Imaxs[i]);
    
    T0s.push_back(_T0s[i]); T1s.push_back(_T1s[i]);;

  }
  
  if (ok) return;
  
  _Qs = Qs; _Ts = Ts; _Tends = Tends; _T0s = T0s; _T1s = T1s; 

  _Tmaxs = Tmaxs; _Imaxs = Imaxs; 
  
  // change global Q, Tstart, Tmax, Imax, Tend
  
  _Q = 0; _Tmax=0; _Imax=0; _Tstart=10000; _T0=10000; _T1=0;
  
  _isBiPolar = false;

  for (size_t i=0; i<_Qs.size();i++){
    
    if (_Qs[i]>0) _Q += _Qs[i];
    
    else _isBiPolar = true;

    if (_Ts[i]<_Tstart) _Tstart = _Ts[i];

    if (_Tmaxs[i]>_Tmax) _Tmax = _Tmaxs[i];

    if (_Imaxs[i]>_Imax) _Imax = _Imaxs[i];
     
    if (_T0s[i]<_T0) _T0 = _T0s[i];

    if (_T1s[i]>_T1) _T1 = _T1s[i];

  }

  // copute new Tend
  
  if(_Qs.size()){ 

    this->getTimer()->getTimeDef()->setT0(_T0s[_T0s.size()-1]); 

    this->getTimer()->getTimeDef()->setT1(_T1s[_T1s.size()-1]); 

    this->getTimer()->getTimeDef()->setTmax(_Tmaxs[_T1s.size()-1]); 
  
    _Tend = this->getTimer()->computeTend(ipulse); //last peak
  
  }
  else{ // all pulses were rejected 
        
    _Tend = NOTIME;

  }

}

//========================= simpleWindow ==========================//


//******************************************************************
simpleWindow::simpleWindow(bool pedEvt) : Algo(pedEvt){ 
//******************************************************************
  

  //Message::MSG("RecoPulseAlgo",gate::VERBOSE,"<<SimpleWindow Construction>>");

  this->getTimer()->config("windowTime");

  this->getChargeAna()->config("fixedTimeQWindow");

  this->setT0(0); this->setT1(0); // analyze full window 
  
  _T0s.push_back(0);  _T1s.push_back(0); 
}

//******************************************************************
double simpleWindow::computePED(RPPulse& ipulse){ 
//******************************************************************
  
  this->getPedAna()->computePedestal(ipulse);
  
  int status = this->getPedAna()->getPedStatus();

  int ext_status = this->getPedAna()->getExtPedStatus();
  
  //if (status<0) 
    
    //Message::MSG("RecoPulseAlgo",gate::VERBOSE,"Channel Pedestal not OK");
  
  double mean; double RMS;
  
  //------------------------------------------------------------

  double chMean =  this->getPedAna()->getMean();

  double chRMS =  this->getPedAna()->getRMS();

  double extMean =  this->getPedAna()->getExtMean();

  double extRMS =  this->getPedAna()->getExtRMS();


  if (ext_status==1){ // external ped estimation is available
    
    if (chRMS == 0 || this->getPedAna()->sameChExtPedRMS()){ // good ch RMS (= ext. RMS) 

      if (_useExtPedRMS){ RMS =  extRMS;}
    
      else RMS =  chRMS;
    
      if (_useExtPedMean){ mean =  extMean; }
      
      else{ mean = chMean; }
          
    }
    else{ // bad Ch RMS: use external trigger
      
      mean =  extMean; RMS =  extRMS;}
  }
  
  else{ // no ext. ped: use channel estimation
   
    mean =  this->getPedAna()->getMean();
    
    RMS =  this->getPedAna()->getRMS();}

  //----------------------------------------------

  this->setPedestal(mean,RMS);

  this->setPedStatus(status);

  if (_nFilter){ ipulse.filter(mean,mean,_nFilter); }
    
  return mean; 
  
}

//******************************************************************
void simpleWindow::setPedestal(double mean, double rms){ 
//******************************************************************

  this->setPedMean(mean); this->setPedRMS(rms);

  //this->getPedAna()->setMean(mean); //?????????????????????

  //this->getPedAna()->setRMS(rms); //???????????????????????
  
  rms =  max(rms,_pedRMS_min);

  double threshold = this->getPedAna()->computeThreshold(mean,rms);
  
  this->getTimer()->getTimeDef()->setThreshold(threshold);

  this->getTimer()->getTimeDef()->setPedestal(mean);

  this->getChargeAna()->setPedestal(mean);
  
}

//******************************************************************
void simpleWindow::setExtPedestal(double mean, double rms){ 
//******************************************************************
  
  this->getPedAna()->setExtMean(mean);

  this->getPedAna()->setExtRMS(rms);

}

//******************************************************************
double simpleWindow::computeQ(RPPulse& ipulse){ 
//******************************************************************
 

  _Qs.clear();

  const vector<unsigned short>& pulse = ipulse.getProfile();
  
  int t1 = this->getT1();

  if (!t1) t1 = pulse.size();
  
  //else if (t1 > (int)pulse.size()){}
      
  gate::Assert(t1 < (int)pulse.size(),__FILE__,__LINE__,
		 
	       gate::internal_logic(Form("T1 (%i)  larger than pulse size (%i)",
				   
				   t1,(int)pulse.size())));
  

  this->getChargeAna()->setWindow(this->getT0(),t1);

  this->getChargeAna()->setWindow(this->getT0(),this->getT1());
  
  _Q = this->getChargeAna()->computeQ(ipulse);

  //if ( _Q < this->getQmin() ) _Q = 0; // pedestal integration 
  
  //if (_Q)
  //Message::MSG("RecoPulseAlgo",gate::DUMP,Form("Pulse Charge: %f",_Q));
  
  _Qs.push_back(_Q);

  return _Q; 
  
 }

//******************************************************************
int simpleWindow::computeT0(RPPulse& ipulse){ 
//******************************************************************
  
  ipulse.size();//dummy stuff

  //Message::MSG("RecoPulseAlgo",gate::NORMAL,"T0 set by user!");

  
  return _T0;
  
}

//******************************************************************
int simpleWindow::computeT1(RPPulse& ipulse){ 
//******************************************************************
  
  ipulse.size();//dummy stuff
  
  //Message::MSG("RecoPulseAlgo",gate::NORMAL,"T0 set by user!");

  return _T1;
  
}

//******************************************************************
int simpleWindow::computeTstart(RPPulse& ipulse){ 
//******************************************************************

  _Ts.clear();
  
  int t1 = _T1;

  if (!t1) t1 = ipulse.size();

  this->getTimer()->getTimeDef()->setT0(_T0);

  this->getTimer()->getTimeDef()->setT1(t1);
  
  _Tstart = this->getTimer()->computeTstart(ipulse);
  
  _Ts.push_back(_Tstart);

  return _Tstart;
  
}

//******************************************************************
int simpleWindow::computeTend(RPPulse& ipulse){ 
//******************************************************************
  
  _Tends.clear();

  int t1 = _T1;

  if (!t1) t1 = ipulse.size();

  this->getTimer()->getTimeDef()->setT0(_T0);

  this->getTimer()->getTimeDef()->setT1(t1);

  _Tend = this->getTimer()->computeTend(ipulse);
  
  _Tends.push_back(_Tend);

  return _Tend;
  
}

//******************************************************************
int simpleWindow::computeTmax(RPPulse& ipulse){ 
//******************************************************************
  
  _Tmaxs.clear();
  
  int t1 = _T1;

  if (!t1) t1 = ipulse.size();

  this->getTimer()->getTimeDef()->setT0(_T0);

  this->getTimer()->getTimeDef()->setT1(t1);

  _Tmax = this->getTimer()->computeTmax(ipulse);
  
  _Tmaxs.push_back(_Tmax);

  return _Tmax; 
  
}

//******************************************************************
int simpleWindow::computeTrise(RPPulse& ipulse){ 
//******************************************************************

  _Trise = this->getTimer()->computeTrise(ipulse);

  return _Trise;
  
}

//******************************************************************
int simpleWindow::computeTfall(RPPulse& ipulse){ 
//******************************************************************
  
  _Tfall = this->getTimer()->computeTfall(ipulse);

  return _Tfall;
  
}

//******************************************************************
void simpleWindow::reco(RPPulse& ipulse){
//******************************************************************
  
  if (_pedEvt) _ped = this->computePED(ipulse);
  
  this->computeQ(ipulse);
  
  //_T0 = this->computeT0(ipulse);

  //_T1 = this->computeT1(ipulse);

  this->computeTstart(ipulse);

  this->computeTend(ipulse);

  this->computeTmax(ipulse);

  this->computeImax(ipulse);
  
  this->computeTrise(ipulse);

  this->computeTfall(ipulse);
  
  //this->applyUserCuts(ipulse);
  
  //if (_calibOn && _Q) this->calibrate();
    
   
}

//******************************************************************
void simpleWindow::setWindow(size_t tmin, size_t tmax){
//******************************************************************
  
  _T0s.clear();_T1s.clear();

  this->setT0((int)tmin);

  this->setT1((int)tmax);
  
  _T0s.push_back(_T0);

  _T1s.push_back(_T1);

}


//******************************************************************
void simpleWindow::reset(){
//******************************************************************
  
  _Q=0;

  //_T0=0;_T1=0;

  _Imax=0;
   
  _Tstart=NOTIME;_Tend=NOTIME;

  _Tmax=NOTIME;_Trise=NOTIME;_Tfall=NOTIME;
  
  _QError=0; _TstartError=0;

  if (_pedEvt){ 
    
    this->getPedAna()->reset();

    _ped=NOPED;_pedRMS=NOPED;

  }

  this->getTimer()->reset();
  
  //this->getChargeAna()->reset();

}

//========================= thresholdWindow =======================//


//******************************************************************
thresholdWindow::thresholdWindow(bool pedEvt) 
  : simpleWindow(pedEvt){ 
//******************************************************************

  //Message::MSG("RecoPulseAlgo",gate::VERBOSE,"<<thresholdWindow Construction>>");


  this->getTimer()->config("thresholdTime");

  this->getChargeAna()->config("fixedTimeQWindow");

}

//******************************************************************
int thresholdWindow::computeT0(RPPulse& ipulse){ 
//******************************************************************
  
  _T0 = this->getTimer()->computeT0(ipulse);
  
  return _T0;
  
}

//******************************************************************
int thresholdWindow::computeT1(RPPulse& ipulse){ 
//******************************************************************
  
  _T1 = this->getTimer()->computeT1(ipulse);
  
  return _T1;
  
}

//******************************************************************
double thresholdWindow::computeQ(RPPulse& ipulse){ 
//******************************************************************
    
  int tmin = _T0 = this->computeTstart(ipulse);
  
  int tmax = _T1 = this->computeTend(ipulse);
  
  if (_T0==NOTIME) {_Q=0; return _Q;}

  if (_T1==NOTIME) tmax = _T1 = ipulse.size();

  this->setT0(tmin); this->setT1(tmax); 
  
  this->getChargeAna()->setWindow(tmin,tmax);
  
  _Q = this->getChargeAna()->computeQ(ipulse);
  
  if ( _Q < this->getQmin() ) _Q = 0; // pedestal integration
  
  //if (_Q)
    //Message::MSG("RecoPulseAlgo",gate::DUMP,Form("Pulse Charge: %f",_Q));

  return _Q;
  
}

//******************************************************************
void thresholdWindow::reco(RPPulse& ipulse){
//******************************************************************

  if (_pedEvt) _ped = this->computePED(ipulse);

  this->computeQ(ipulse);
  
  this->computeTmax(ipulse);

  this->computeImax(ipulse);

  this->computeTrise(ipulse);

  this->computeTfall(ipulse);

  this->applyUserCuts(ipulse);

  //if (_calibOn && _Q) this->calibrate();

}

//******************************************************************
void thresholdWindow::reset(){
//******************************************************************
  
  _Q=0;

  _T0=NOTIME;_T1=NOTIME;
   
  _Tstart=NOTIME;_Tend=NOTIME;

  _Tmax=NOTIME;_Trise=NOTIME;_Tfall=NOTIME;
  
  _QError=0; _TstartError=0;

  _Imax=0;

  if (_pedEvt){ 
    
    this->getPedAna()->reset();

    _ped=NOPED;_pedRMS=NOPED;

  }

  this->getTimer()->reset();

  //this->getChargeAna()->reset();

}

//========================= slidingWindow =======================//

//******************************************************************
slidingWindow::slidingWindow(bool pedEvt) 
  : simpleWindow(pedEvt){ 
//******************************************************************
  
  
  //Message::MSG("RecoPulseAlgo",gate::VERBOSE,"<<SlidingWindow Construction>>");


  this->getTimer()->config("windowTime");
  
  this->getChargeAna()->config("slidingQWindow");
  
  // default window size
  
  this->setWindowSize(56);
 
}


//******************************************************************
void slidingWindow::setChargeAna(string name){ 
//******************************************************************
  
  
  bool status = (name!="slidingQWindow" && name!="slidingSPQWindow");

  gate::Assert(status == false,__FILE__,__LINE__,
	     gate::internal_logic("Not a valid charge analyzer: "+name));

  this->getChargeAna()->config(name);

}

//******************************************************************
void slidingWindow::setWindowSize(size_t size){ 
//******************************************************************
  
  this->getChargeAna()->setWindowSize(size);
  
}




//******************************************************************
double slidingWindow::computeQ(RPPulse& ipulse){ 
//******************************************************************
    
  //Message::MSG("RecoPulseAlgo",gate::DUMP,"computeQ function");

  _T1s.clear();_Ts.clear();_Tends.clear();_Qs.clear();

  _TsError.clear();_QsError.clear();
  
  vector<double> Qs;
    
  vector<int> T0s; vector<int> T1s; 
 
  double Qp=1; double Q = 0; double Qfactor=1e6;

  //RPPulse tpulse = RPPulse(ipulse);
  
  while (Qp) { //loop over peaks

    this->getChargeAna()->setWindow(0,ipulse.size());
    
    Qp = this->computePeakQ(ipulse); Q += Qp;  
    
    if (!Qp) break; // no signal found

    //tpulse.smoothPed(round(_ped),_T0,_T1);//   
    ipulse.integrated(_T0,_T1);
    
    Qs.push_back(Qp/Qfactor+_T0); // nasty trick to sort Qs according to time
    
    T0s.push_back(_T0); T1s.push_back(_T1);
  
  }

  //peak counting
  
  if (Qs.size()){
    
  std::sort(T0s.begin(),T0s.end());
  
  std::sort(T1s.begin(),T1s.end());

  std::sort(Qs.begin(),Qs.end());
  
  _T0s.push_back(T0s[0]); _T1s.push_back(T1s[0]); 

  _Qs.push_back((Qs[0]-T0s[0])*Qfactor); // first peak
  
  for (size_t i=1; i<Qs.size();i++){
  
    if (_T1s[_T1s.size()-1]>=T0s[i]){ //same peak
     
	_T1s[_T1s.size()-1]=T1s[i]; _Qs[_Qs.size()-1]+=(Qs[i]-T0s[i])*Qfactor;
    }
    else{ 

      _T0s.push_back(T0s[i]); _T1s.push_back(T1s[i]); // new peak
    
      _Qs.push_back((Qs[i]-T0s[i])*Qfactor);}
    
  }}
  
  // set global info
 
  _Q = Q; 
  
  if (_T0s.size()) {_T0 = _T0s[0]; _T1 = _T1s[_T1s.size()-1];}

  else{_T0 = 0; _T1 = ipulse.size();}

  this->getTimer()->getTimeDef()->setT0(_T0);

  this->getTimer()->getTimeDef()->setT1(_T1);
  
  //if (_Q)
    
    //Message::MSG("RecoPulseAlgo",gate::DUMP,Form("Pulse Charge: %f",_Q));
  
  return _Q;
  
}


//******************************************************************
double slidingWindow::computePeakQ(RPPulse& ipulse){ 
//******************************************************************
  
  _Q = this->getChargeAna()->computeQ(ipulse);
  
  // T0 and T1 computed along with Q

  _T0= this->getChargeAna()->getTmin();
  
  _T1= this->getChargeAna()->getTmax();
  
 
  //------- extend window. 
  
  this->extendWindow(ipulse);
  
  //-------- end of extended window approach

  //set T0 and T1 in window timmer
  
  this->getTimer()->getTimeDef()->setT0(_T0);

  this->getTimer()->getTimeDef()->setT1(_T1);
  
 
  if ( _Q < this->getQmin() ){ 
    
    //Message::MSG("RecoPulseAlgo",gate::DUMP,Form("Peak within [%i,%i]",_T0,_T1));

    //Message::MSG("RecoPulseAlgo",gate::DUMP,
    //		 Form("... not over threshold: %0.1f>%0.1f",this->getQmin(),_Q));
    
    _Q = 0; // pedestal integration
    
    return _Q;

  }
  
  //Message::MSG("RecoPulseAlgo",gate::DUMP,
	       
  //	       Form("Peak within [%i,%i]: %0.2f DUQ",_T0,_T1,_Q));

  return _Q; 
  
}

//******************************************************************
int slidingWindow::computeT0(RPPulse& ipulse){ 
//******************************************************************
  
  if (_T0 == NOTIME) this->computeQ(ipulse); 
  
  return  _T0;
  
}

//******************************************************************
int slidingWindow::computeT1(RPPulse& ipulse){ 
//******************************************************************
  
  if (_T1 == NOTIME) this->computeQ(ipulse); 

  return _T1;
  
}

//******************************************************************
int slidingWindow::computeTmax(RPPulse& ipulse){ 
//******************************************************************
  

  if (_T0 == NOTIME) this->computeQ(ipulse); 
 
  if (!_Qs.size()) return NOTIME;
  

  _Tmaxs.clear(); _Imaxs.clear();
  
  _Tmax=NOTIME;
   
  for (size_t i=0; i < _Qs.size(); i++){
    
    this->getTimer()->getTimeDef()->setT0(_T0s[i]);

    this->getTimer()->getTimeDef()->setT1(_T1s[i]);

    int Tmax = this->getTimer()->computeTmax(ipulse);
  
    _Tmaxs.push_back(Tmax);
    
    //Message::MSG("RecoPulseAlgo",gate::DUMP,
    //		 Form("Peak within [%i,%i]",_T0s[i],_T1s[i]));

    //Message::MSG("RecoPulseAlgo",gate::DUMP,Form("... Maximum: %i",Tmax));
    
    const vector<unsigned short>& pulse = ipulse.getProfile();
    
    if (Tmax!=NOTIME){ 

      if (_Tmax==NOTIME || pulse[Tmax] < pulse[_Tmax ]) _Tmax=Tmax;}    
  }

  this->getTimer()->getTimeDef()->setT0(_T0);

  this->getTimer()->getTimeDef()->setT1(_T1);


  //if (_Tmax)
    //Message::MSG("RecoPulseAlgo",gate::DUMP,Form("Pulse Maximum: %i",_Tmax));

  return _Tmax; 
  
}

//******************************************************************
int slidingWindow::computeTstart(RPPulse& ipulse){ 
//******************************************************************
  
  if (_T0 == NOTIME) this->computeQ(ipulse); 

  if (_Tmax == NOTIME) this->computeTmax(ipulse); 
  
  if (!_Qs.size()) return NOTIME;

  //compute tstart of each peak
  
  _Ts.clear();
  
  for (size_t i=0; i < _Qs.size(); i++){
    
    this->getTimer()->getTimeDef()->setT0(_T0s[i]);

    this->getTimer()->getTimeDef()->setT1(_T1s[i]);

    this->getTimer()->getTimeDef()->setTmax(_Tmaxs[i]);

    int Tpeak = this->getTimer()->computeTstart(ipulse);
  
    _Ts.push_back(Tpeak);
    
    //Message::MSG("RecoPulseAlgo",gate::DUMP,
    //		 Form("Peak within [%i,%i]",_T0s[i],_T1s[i]));

    //Message::MSG("RecoPulseAlgo",gate::DUMP,Form("... start time: %i",Tpeak));
   
  }
  
  //compute  global tstart

  this->getTimer()->getTimeDef()->setT0(_T0); //restore global T0

  this->getTimer()->getTimeDef()->setT1(_T1); //restore global T1

  this->getTimer()->getTimeDef()->setTmax(_Tmax); //restore global Tmax
  
  //_Tstart = this->getTimer()->computeTstart(ipulse);//global
  
  _Tstart = _Ts[0]; //first peak

  this->getTimer()->getTimeDef()->setTstart(_Tstart);
  
  //Message::MSG("RecoPulseAlgo",gate::DUMP,Form("Pulse start time: %i",_Tstart));

  return _Tstart;
  
}


//******************************************************************
int slidingWindow::computeTend(RPPulse& ipulse){ 
//******************************************************************
  
  if (_T1 == NOTIME) this->computeQ(ipulse); 

  if (_Tmax == NOTIME) this->computeTmax(ipulse); 

  if (!_Qs.size()) return NOTIME;
  
  
  //---------------------------

  //compute tend of each peak
  
  _Tends.clear();
  
  for (size_t i=0; i < _Qs.size(); i++){
    
    this->getTimer()->getTimeDef()->setT0(_T0s[i]);

    this->getTimer()->getTimeDef()->setT1(_T1s[i]);

    this->getTimer()->getTimeDef()->setTmax(_Tmaxs[i]);

    int Tpeak = this->getTimer()->computeTend(ipulse);
  
    _Tends.push_back(Tpeak);
    
    //Message::MSG("RecoPulseAlgo",gate::DUMP,
    //Form("Peak within [%i,%i]",_T0s[i],_T1s[i]));

  //Message::MSG("RecoPulseAlgo",gate::DUMP,Form("... end time: %i",Tpeak));
   
  }
  //_------------------------------------

  _Tend = _Tends[_Tends.size()-1];

  this->getTimer()->getTimeDef()->setTend(_Tend);

  this->getTimer()->getTimeDef()->setT0(_T0); //restore global T0

  this->getTimer()->getTimeDef()->setT1(_T1); //restore global T1

  this->getTimer()->getTimeDef()->setTmax(_Tmax); //restore global Tmax
  
//Message::MSG("RecoPulseAlgo",gate::DUMP,Form("Pulse end time: %i",_Tend));

  return _Tend;
  
}


//******************************************************************
void slidingWindow::reco(RPPulse& ipulse){
//******************************************************************

  if (_pedEvt){ 
    
    _ped = this->computePED(ipulse);
        
  }
  
  //Check saturation
  //ipulse = this->getSatAna()->smoothSat(ipulse);
  //

  _Q = this->computeQ(ipulse);
  
 
  _Tmax = this->computeTmax(ipulse);
  
  _Imax = this->computeImax(ipulse);

  _Tstart = this->computeTstart(ipulse);

  _Tend = this->computeTend(ipulse);

  _Trise = this->computeTrise(ipulse);

  _Tfall = this->computeTfall(ipulse);

  this->applyUserCuts(ipulse);

  //if (_calibOn && _Q) this->calibrate();
 
 

}


//******************************************************************
void slidingWindow::checkPed(){ 
//******************************************************************

  
  bool status = (this->getPedMean()==NOPED);

  gate::Assert(status == false,__FILE__,__LINE__,
	       gate::internal_logic("Pedestal not defined!"));

}

//******************************************************************
void slidingWindow::reset(){
//******************************************************************
  
  _Q=0;

  _T0=NOTIME;_T1=NOTIME;
   
  _Tstart=NOTIME;_Tend=NOTIME;

  _Tmax=NOTIME;_Trise=NOTIME;_Tfall=NOTIME;
  
  _QError=0; _TstartError=0;

  _Qs.clear(); _T0s.clear(); _T1s.clear(); 
  
  _Ts.clear();_Tends.clear(); _Tmaxs.clear();_Imaxs.clear();
  
  _TsError.clear();_QsError.clear();
  
  _Imax=0;

  if (_pedEvt){ 
    
    this->getPedAna()->reset();

    _ped=NOPED;_pedRMS=NOPED;

  }

  this->getTimer()->reset();

  this->getChargeAna()->reset();
  
}

//========================== maxWindow =========================//

//******************************************************************
maxWindow::maxWindow(bool pedEvt) 
  : simpleWindow(pedEvt){ 
//******************************************************************
  
  
  //Message::MSG("RecoPulseAlgo",gate::VERBOSE,"<<MaxWindow Construction>>");

  this->getTimer()->config("maxAmpTime");

  this->getChargeAna()->config("fixedTimeQWindow");

  this->setWindow(2,8);

  _T0 =_T1 = NOTIME;

}

//******************************************************************
void maxWindow::setChargeAna(string name){ 
//******************************************************************
  
 
  bool status = (name!="fixedTimeQWindow" && name!="splineQ");

  gate::Assert(status == false,__FILE__,__LINE__,
	       gate::internal_logic("Not a valid charge analyzer: "+name));

  this->getChargeAna()->config(name);

}

//******************************************************************
double maxWindow::computeQ(RPPulse& ipulse){ 
//******************************************************************

  _T1s.clear();_Ts.clear(); _Tends.clear();_Qs.clear();_Tmaxs.clear();_Imaxs.clear();

  _TsError.clear();_QsError.clear();

  vector<double> Qs;
    
  vector<int> T0s; vector<int> T1s; vector<int> Tmaxs; 
  
  double Qp=1; double Q = 0; double Qfactor=1e6;
  
  //RPPulse tpulse = RPPulse(ipulse);

  while (Qp) {
    
    this->getTimer()->reset();
    
    Qp = this->computePeakQ(ipulse); Q += Qp;  

    if (!Qp) break;

    //tpulse.smoothPed(_ped,_T0,_T1);
    ipulse.integrated(_T0,_T1);

    Qs.push_back(Qp/Qfactor+_T0); 
    
    T0s.push_back(_T0); T1s.push_back(_T1);

    Tmaxs.push_back(_Tmax);
    
  }
   

  //peak counting
  
  if (Qs.size()){
  
  std::sort(T0s.begin(),T0s.end());
  
  std::sort(T1s.begin(),T1s.end());

  std::sort(Tmaxs.begin(),Tmaxs.end());

  std::sort(Qs.begin(),Qs.end());
  
  _T0s.push_back(T0s[0]); _T1s.push_back(T1s[0]); 

  _Tmaxs.push_back(Tmaxs[0]); _Qs.push_back((Qs[0]-T0s[0])*1e6);
  
  for (size_t i=1; i<Qs.size();i++){
    
    if (_T1s[_T1s.size()-1]>=T0s[i]){//same peak
      
      _T1s[_T1s.size()-1]=T1s[i]; _Qs[0]+=(Qs[i]-T0s[i])*Qfactor;
      
      if (Tmaxs[_T1s.size()-1]<=Tmaxs[i]) Tmaxs[_T1s.size()-1]=Tmaxs[i];

      }

    else{ //new peak found
      
      _T0s.push_back(T0s[i]); _T1s.push_back(T1s[i]); 
    
      _Tmaxs.push_back(Tmaxs[i]);_Qs.push_back((Qs[i]-T0s[i])*Qfactor);
    
    }}}
  
  // set global info
 
  _Q = Q; 

  if (_T0s.size()) {_T0 = _T0s[0]; _T1 = _T1s[_T1s.size()-1];}

  else{_T0 = 0; _T1 = ipulse.size();}

  this->getTimer()->getTimeDef()->setT0(_T0);

  this->getTimer()->getTimeDef()->setT1(_T1);  

  //if (_Q)

    //Message::MSG("RecoPulseAlgo",gate::DUMP,Form("Pulse Charge: %f",_Q));

  return _Q;
}
  
//******************************************************************
double maxWindow::computePeakQ(RPPulse& ipulse){ 
//******************************************************************
  
  int tmin = this->getTimer()->computeT0(ipulse);

  int tmax = this->getTimer()->computeT1(ipulse);
  
  if (tmin==NOTIME){ 
    
    //Message::MSG("RecoPulseAlgo",gate::NORMAL,"Could not compute T0!");

    //Message::MSG("RecoPulseAlgo",gate::NORMAL,"Integrating from first sample");
    
    tmin = 0; }
  
  if (tmax==NOTIME){ 

    //Message::MSG("RecoPulseAlgo",gate::NORMAL,"Could not compute T1!");
    
    //Message::MSG("RecoPulseAlgo",gate::NORMAL,"Integrating from last sample");

    tmax = ipulse.size(); }
  
  this->getChargeAna()->setWindow(tmin,tmax);
  
  //testing

  _T0 = tmin; _T1 = tmax;

  _Tmax = this->getTimer()->computeTmax(ipulse);
  
  _Q = this->getChargeAna()->computeQ(ipulse);
  
  if (_Q<0){  _Q=0; return _Q;}

  //------- extend window. 
  
  this->extendWindow(ipulse);
  
  //-------- end of extended window approach


  if ( _Q < this->getQmin() ) _Q = 0; // pedestal integration
  
  return _Q; 
  
}

//******************************************************************
int maxWindow::computeT0(RPPulse& ipulse){ 
//******************************************************************
  
  if (!_Q) this->computeQ(ipulse);

  //_T0 = this->getTimer()->computeT0(ipulse);

  return  _T0;
  
}

//******************************************************************
int maxWindow::computeT1(RPPulse& ipulse){ 
//******************************************************************
  
  if (!_Q) this->computeQ(ipulse);

  //_T1 = this->getTimer()->computeT1(ipulse); 
  
  return  _T1;

}



//******************************************************************
int maxWindow::computeTstart(RPPulse& ipulse){ 
//******************************************************************
  
  if (_T0 == NOTIME) this->computeQ(ipulse); 
  
  if (!_Qs.size()) return NOTIME;

  //compute tstart of each peak
  
  _Ts.clear();
  
  for (size_t i=0; i < _Qs.size(); i++){
        
    this->getTimer()->getTimeDef()->setTmax(_Tmaxs[i]);
    
    this->getTimer()->computeT0(ipulse);

    this->getTimer()->computeT1(ipulse);

    int Tpeak = this->getTimer()->computeTstart(ipulse);
    
    _Ts.push_back(Tpeak);
    
    //Message::MSG("RecoPulseAlgo",gate::DUMP,
    //		 Form("Peak within [%i,%i]",_T0s[i],_T1s[i]));

    //Message::MSG("RecoPulseAlgo",gate::DUMP,Form("... start time: %i",Tpeak));

  }
  
 
  _Tstart = _Ts[0]; //first peak
  
  this->getTimer()->getTimeDef()->setTstart(_Tstart);
  
  //Message::MSG("RecoPulseAlgo",gate::DUMP,Form("Pulse arrival time: %i",_Tstart));

  this->getTimer()->getTimeDef()->setT0(_T0); //restore global T0

  this->getTimer()->getTimeDef()->setT1(_T1); //restore global T1
  
  this->getTimer()->getTimeDef()->setTmax(_Tmax); //restore global Tmax

  

  return _Tstart;
  
}



//******************************************************************
int maxWindow::computeTend(RPPulse& ipulse){ 
//******************************************************************
  
  if (_T1 == NOTIME) this->computeQ(ipulse); 
  
  if (!_Qs.size()) return NOTIME;

  
  //compute tend of each peak
  
  _Tends.clear();
  
  for (size_t i=0; i < _Qs.size(); i++){
        
    this->getTimer()->getTimeDef()->setTmax(_Tmaxs[i]);
    
    this->getTimer()->computeT0(ipulse);

    this->getTimer()->computeT1(ipulse);

    int Tpeak = this->getTimer()->computeTend(ipulse);
    
    _Tends.push_back(Tpeak);
    
    //Message::MSG("RecoPulseAlgo",gate::DUMP,
    //Form("Peak within [%i,%i]",_T0s[i],_T1s[i]));

  //Message::MSG("RecoPulseAlgo",gate::DUMP,Form("... end time: %i",Tpeak));

  }
  
 
  _Tend = _Tends[_Tends.size()-1]; //first peak


  this->getTimer()->getTimeDef()->setTend(_Tend);

  this->getTimer()->getTimeDef()->setT0(_T0); //restore global T0

  this->getTimer()->getTimeDef()->setT1(_T1); //restore global T1
  
  this->getTimer()->getTimeDef()->setTmax(_Tmax); //restore global T1
      
//Message::MSG("RecoPulseAlgo",gate::DUMP,Form("Pulse end time: %i",_Tend));

  return _Tend;
  
}




//******************************************************************
void maxWindow::setWindow(size_t left, size_t right){ 
//******************************************************************
  
  this->getTimer()->getTimeDef()->setWindow(left,right);
  
}

//******************************************************************
void maxWindow::reco(RPPulse& ipulse){
//******************************************************************
  


  if (_pedEvt) _ped = this->computePED(ipulse);
  
  this->computeQ(ipulse);
  
  this->computeTmax(ipulse);

  this->computeImax(ipulse);

  this->computeTstart(ipulse);
  
  this->computeTend(ipulse);
  
  this->computeTrise(ipulse);

  this->computeTfall(ipulse);

  this->computeT0(ipulse);

  this->computeT1(ipulse);
  
  this->applyUserCuts(ipulse);

  //if (_calibOn && _Q) this->calibrate();

 
 
}

//******************************************************************
void maxWindow::reset(){
//******************************************************************
  
  _Q=0;

  _T0=NOTIME;_T1=NOTIME;
   
  _Tstart=NOTIME;_Tend=NOTIME;

  _Tmax=NOTIME;_Trise=NOTIME;_Tfall=NOTIME;
  
  _QError=0; _TstartError=0;

  _Qs.clear(); _T0s.clear(); _T1s.clear(); 
  
  _Ts.clear();_Tends.clear(); _Imaxs.clear();

  _TsError.clear();_QsError.clear();
  
  _Imax=0;

  if (_pedEvt){ 
    
    this->getPedAna()->reset();

    _ped=NOPED;_pedRMS=NOPED;

  }

  this->getTimer()->reset();

  //this->getChargeAna()->reset();

}

//========================== peakWindow =========================//

//******************************************************************
peakWindow::peakWindow(bool pedEvt) : simpleWindow(pedEvt){ 
//******************************************************************
  
    
  this->getTimer()->config("windowTime");

  this->getChargeAna()->config("fixedTimeQWindow");
  
  _Wmin=0;

  _T0 =_T1 = NOTIME;

}

//******************************************************************
int peakWindow::computeT0(RPPulse& ipulse){ 
//******************************************************************
    
  if (_T0==NOTIME) this->computeQ(ipulse);

  return _T0;

}

//******************************************************************
int peakWindow::computeT1(RPPulse& ipulse){ 
//******************************************************************
  
  if (_T1==NOTIME) this->computeQ(ipulse);
  
  return _T1;

}

//******************************************************************
int peakWindow::computeTmax(RPPulse& ipulse){ 
//******************************************************************

   if (_T0 == NOTIME) this->computeQ(ipulse); 
 
  //compute tstart of each peak

   if (!_Qs.size()) return NOTIME;

  _Tmaxs.clear();_Imaxs.clear();

  for (size_t i=0; i < _Qs.size(); i++){
    
    this->getTimer()->getTimeDef()->setT0(_T0s[i]);

    this->getTimer()->getTimeDef()->setT1(_T1s[i]);

    int Tmax = this->getTimer()->computeTmax(ipulse);
  
    _Tmaxs.push_back(Tmax);
    
    //Message::MSG("RecoPulseAlgo",gate::DUMP,
    //		 Form("Peak within [%i,%i]",_T0s[i],_T1s[i]));

    //Message::MSG("RecoPulseAlgo",gate::DUMP,Form("... maximum at: %i",Tmax));

  }

  this->getTimer()->getTimeDef()->setT0(_T0);

  this->getTimer()->getTimeDef()->setT1(_T1);

  _Tmax = this->getTimer()->computeTmax(ipulse);
  
  //if (_Tmax)
    //Message::MSG("RecoPulseAlgo",gate::DUMP,Form("Pulse maximum at: %i",_Tmax));

  return _Tmax; 
  
}

//******************************************************************
int peakWindow::computeTstart(RPPulse& ipulse){ 
//******************************************************************
  
  if (_T0 == NOTIME) this->computeQ(ipulse); 
  
  if (!_Qs.size()) return NOTIME;

  //compute tstart of each peak
  
  _Ts.clear();
  
  for (size_t i=0; i < _Qs.size(); i++){
    
    this->getTimer()->getTimeDef()->setT0(_T0s[i]);

    this->getTimer()->getTimeDef()->setT1(_T1s[i]);
    
    this->getTimer()->getTimeDef()->setTmax(_Tmaxs[i]);

    int Tpeak = this->getTimer()->computeTstart(ipulse);
  
    _Ts.push_back(Tpeak);

    //Message::MSG("RecoPulseAlgo",gate::DUMP,
    //		 Form("Peak within [%i,%i]",_T0s[i],_T1s[i]));

    //Message::MSG("RecoPulseAlgo",gate::DUMP,Form("... start time: %i",Tpeak));
  }
  
  //compute  tstar of first peak
  
  //_Tstart = this->getTimer()->computeTstart(ipulse);

  _Tstart = _Ts[0];

  this->getTimer()->getTimeDef()->setT0(_T0); //restore global T0

  this->getTimer()->getTimeDef()->setT1(_T1); //restore global T1

  this->getTimer()->getTimeDef()->setTmax(_Tmax); //restore global Tmax

  this->getTimer()->getTimeDef()->setTstart(_Tstart); //restore global Tstart
 
  //Message::MSG("RecoPulseAlgo",gate::DUMP,Form("Pulse start time: %i",_Tstart));
 
  return _Tstart;
  
}

//******************************************************************
int peakWindow::computeTend(RPPulse& ipulse){ 
//******************************************************************
  
  if (_T0 == NOTIME) this->computeQ(ipulse); 
  
  if (!_Qs.size()) return NOTIME;
  
   //compute tstart of each peak
  
  _Tends.clear();
  
  for (size_t i=0; i < _Qs.size(); i++){
    
    this->getTimer()->getTimeDef()->setT0(_T0s[i]);

    this->getTimer()->getTimeDef()->setT1(_T1s[i]);
    
    this->getTimer()->getTimeDef()->setTmax(_Tmaxs[i]);

    int Tpeak = this->getTimer()->computeTend(ipulse);
  
    _Tends.push_back(Tpeak);

    //Message::MSG("RecoPulseAlgo",gate::DUMP,
    //		 Form("Peak within [%i,%i]",_T0s[i],_T1s[i]));

    //Message::MSG("RecoPulseAlgo",gate::DUMP,Form("... end time: %i",Tpeak));
  }
  
  

  _Tend = _Tends[_Tends.size()-1];

  
  this->getTimer()->getTimeDef()->setTend(_Tend); //set global end

  this->getTimer()->getTimeDef()->setT0(_T0); //restore global T0

  this->getTimer()->getTimeDef()->setT1(_T1); //restore global T1
  
  this->getTimer()->getTimeDef()->setTmax(_Tmax); //restore global Tmax

  
  //Message::MSG("RecoPulseAlgo",gate::DUMP,Form("Pulse end time: %i",_Tend));

  return _Tend;
  
}


//******************************************************************
double peakWindow::computeQ(RPPulse& ipulse){ 
//******************************************************************
  
  _T1s.clear();_T1s.clear();_Qs.clear();

  _TsError.clear();_QsError.clear();

  this->checkPed(); 
  
  double thr = _ped - this->getPedAna()->getNsigOverPed()*_pedRMS;
  
  if (!ipulse.getPeakGatherer().status()) ipulse.getPeaks(thr,30);

  //if (!ipulse.getPeakGatherer().status()) ipulse.getPeaks(_ped,30);
  
  vector<int> T0s = ipulse.getPeakGatherer().getT0s();
  
  vector<int> T1s = ipulse.getPeakGatherer().getT1s();

  for (size_t i = 0; i<T0s.size();i++){
    
    _T0 = T0s[i]; _T1 = T1s[i];
    
    this->getChargeAna()->setWindow(_T0,_T1);
    
    double Q = this->getChargeAna()->computeQ(ipulse);
    
    if ( Q < this->getQmin() ) {
      
      //Message::MSG("RecoPulseAlgo",
      //	   gate::DUMP,Form("Peak within [%i,%i]",_T0,_T1));

      //Message::MSG("RecoPulseAlgo",gate::DUMP,
      //	   Form("... not over threshold: %0.1f>%0.1f",
      //		this->getQmin(),Q));
      Q = 0; }
    
  
    else{ _Q+=Q; _Qs.push_back(Q);
    
      _T0s.push_back(_T0);_T1s.push_back(_T1); }
    
  }
  
  if (!_T0s.size()) {_T0s.push_back(0); _T1s.push_back(0);}
  
  _T0 = _T0s[0]; _T1 = _T1s[_T1s.size()-1];
  

  this->getTimer()->getTimeDef()->setT0(_T0);

  this->getTimer()->getTimeDef()->setT1(_T1);
  
  //if (_Q)
    //Message::MSG("RecoPulseAlgo",gate::DUMP,Form("Pulse Charge: %f",_Q));

  return _Q;

}

//******************************************************************
void peakWindow::reco(RPPulse& ipulse){ 
//******************************************************************
  
  if (_pedEvt) _ped = this->computePED(ipulse);
  
  this->computeQ(ipulse);
  
  this->computeTmax(ipulse);

  this->computeImax(ipulse);

  this->computeTstart(ipulse);

  this->computeTend(ipulse);

  this->computeTrise(ipulse);
  
  this->computeTfall(ipulse);
  
  this->applyUserCuts(ipulse);
  
  //if (_calibOn && _Q) this->calibrate();

  

}

//******************************************************************
void peakWindow::checkPed(){ 
//******************************************************************

  
  bool status = (this->getPedMean()==NOPED);
  
  gate::Assert(status == false,__FILE__,__LINE__,
	       gate::internal_logic("Pedestal not defined!"));
 

}

//******************************************************************
void peakWindow::reset(){
//******************************************************************
  
  _Q=0;

  _T0=NOTIME;_T1=NOTIME;
   
  _Tstart=NOTIME;_Tend=NOTIME;

  _Tmax=NOTIME;_Trise=NOTIME;_Tfall=NOTIME;
  
  _QError=0; _TstartError=0;

  _Qs.clear(); _T0s.clear(); _T1s.clear();

  _Tmaxs.clear();_Imaxs.clear();

  _Ts.clear();_Tends.clear();
  
  _TsError.clear();_QsError.clear();
  
  _Imax=0;

  if (_pedEvt){ 
    
    this->getPedAna()->reset();

    _ped=NOPED;_pedRMS=NOPED;

  }

  this->getTimer()->reset();

  this->getChargeAna()->reset();

  _isBiPolar = false;
  
}

//========================== biPolarPeakWindow =========================//

//******************************************************************
biPolarPeakWindow::biPolarPeakWindow(bool pedEvt) : peakWindow(pedEvt){ 
//******************************************************************
  

  //Message::MSG("RecoPulseAlgo",gate::VERBOSE,
  //	       "<<BiPolarPeakWindow Construction>>");

  this->getTimer()->config("windowTime");

  this->getTimer()->getTimeDef()->
    
    setTstartAlgo("max");// both pos and neg pulses

  this->getChargeAna()->config("fixedTimeQWindow");
  
  _T0 =_T1 = NOTIME;

}

//******************************************************************
double biPolarPeakWindow::computeQ(RPPulse& ipulse){ 
//******************************************************************
  
  _T1s.clear();_T1s.clear();_Qs.clear();

  _TsError.clear();_QsError.clear();

  this->checkPed(); 

  if (!ipulse.getPeakGatherer().status()) ipulse.getBPPeaks(_ped);
  
  vector<int> T0s = ipulse.getPeakGatherer().getT0s();
  
  vector<int> T1s = ipulse.getPeakGatherer().getT1s();

  for (size_t i = 0; i<T0s.size();i++){
    
    _T0 = T0s[i]; _T1 = T1s[i];
    
    this->getChargeAna()->setWindow(_T0,_T1);
    
    double Q = this->getChargeAna()->computeQ(ipulse);
    
    if ( fabs(Q) < this->getQmin() ) {
      
      //Message::MSG("RecoPulseBPAlgo",
      //	   gate::DUMP,Form("Peak within [%i,%i]",_T0,_T1));

      //Message::MSG("RecoPulseBPAlgo",gate::DUMP,
      //	   Form("... not over threshold: %0.1f>%0.1f",
      //		this->getQmin(),Q));
      Q = 0; }
    
    
    else{ 
      
      if(Q>0){ _Q += Q; }// total charge only computed with neg pulses 
      
      else{ _isBiPolar = true;}

      _Qs.push_back(Q);
      
      _T0s.push_back(_T0); _T1s.push_back(_T1); 
      
    }
    
  }
  
  if (!_T0s.size()) {_T0s.push_back(0); _T1s.push_back(0);}
  
  _T0 = _T0s[0]; _T1 = _T1s[_T1s.size()-1];
  

  this->getTimer()->getTimeDef()->setT0(_T0);

  this->getTimer()->getTimeDef()->setT1(_T1);
  
  //if (_Q)
    //Message::MSG("RecoPulseAlgo",gate::DUMP,Form("Pulse Charge: %f",_Q));
  
  return _Q;

}
