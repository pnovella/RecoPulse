///
/// \file DCRecoPulseTime.cc
/// Time reconstruction algorithms 
///
/// \author Pau Novella - CIEMAT 2008
///


#include<RecoPulseTime.h>
#include<math.h>
#include <stdlib.h>

using namespace std;

//*****************************************************************
timeDef::timeDef(){
//*****************************************************************

  this->reset();

}



//*****************************************************************
int timeDef::computeTrise(RPPulse& ipulse){
//*****************************************************************
  
  if (this->getTstart()==NOTIME) this->computeTstart(ipulse);

  if (this->getTmax()==NOTIME) this->computeTmax(ipulse);

  int trise =  this->getTmax() - this->getTstart();

  this->setTrise(trise);
  
  return trise;

}

//*****************************************************************
int timeDef::computeTfall(RPPulse& ipulse){
//*****************************************************************
  

  if (this->getTend()==NOTIME) this->computeTend(ipulse);

  if (this->getTmax()==NOTIME) this->computeTmax(ipulse);

  int tfall =  this->getTend() - this->getTmax();

  this->setTfall(tfall);
  
  return tfall;

}



//======================= thresholdTime ==========================//


//*****************************************************************
thresholdTime::thresholdTime(double thr) : timeDef(){
//*****************************************************************
    
  _threshold=thr;
  
}

//*****************************************************************
int thresholdTime::computeT0(RPPulse& ipulse){
//*****************************************************************

  /*
    Returns sample index of pulse defining t0
  */
  
  const vector<unsigned short>& pulse = ipulse.getProfile();

  int isample=0; bool found = false;
  
  vector<unsigned short>::const_iterator iamp;
	  
  for( iamp = pulse.begin(); iamp != pulse.end(); iamp++ ){ 

    const  double Amp = (*iamp);
    
    if (this->isOverThreshold(Amp)) { found=true; break; }
    
    isample++;

  }
  
  if (!found) isample = NOTIME;
  
  this->setT0(isample);

  return isample;

}

//*****************************************************************
int thresholdTime::computeT1(RPPulse& ipulse){
//*****************************************************************

  /*
    Returns sample index of pulse defining t1=t_end
  */
  
  
  const vector<unsigned short>& pulse = ipulse.getProfile();

  if (this->getT0()==NOTIME) this->computeT0(ipulse);
  
  //if (this->getT0()==NOTIME) this->setT0(ipulse.size());//no pulse found

  int isample=this->getT0(); bool found = false;
  
  vector<unsigned short>::const_iterator iamp;
  
  for( iamp = pulse.begin()+this->getT0(); iamp != pulse.end(); iamp++ ){ 
    
    const  double Amp = (*iamp);
    
    if (!this->isOverThreshold(Amp)) { found = true; break; }
    
    isample++;

  } if (!found) isample = NOTIME;
  
  this->setT1(isample);

  return isample;

}

//*****************************************************************
int thresholdTime::computeTstart(RPPulse& ipulse){
//*****************************************************************
  
  int tstart = this->computeT0(ipulse);

  this->setTstart(tstart);

  return tstart;

}

//*****************************************************************
int thresholdTime::computeTend(RPPulse& ipulse){
//*****************************************************************
  
  int tend = this->computeT1(ipulse);

  this->setTend(tend);

  return tend;
}

//*****************************************************************
int thresholdTime::computeTmax(RPPulse& ipulse){
//*****************************************************************

  /*
    Returns sample index of pulse with max amp in [Tstart,Tend]
  */

  const vector<unsigned short>& pulse = ipulse.getProfile();

  if (this->getTstart()==NOTIME) this->computeTstart(ipulse);

  if (this->getTend()==NOTIME) this->computeTend(ipulse);

  int fsample=this->getTstart(); int lsample=this->getTend();
  
  if (fsample==NOTIME) fsample = 0;
  
  if (lsample==NOTIME) lsample = pulse.size();

  //int fsample=this->getT0(); int lsample=this->getT1();
  
  int msample = fsample; int isample = fsample; 

  //double ped = pulse[fsample-1];
  double ped = pulse[0];

  double max = fabs(pulse[fsample]-ped);
  //double max = 0;

  vector<unsigned short>::const_iterator iamp;
  
  for( iamp = pulse.begin() + fsample; 
       iamp != pulse.begin() + lsample; 
       iamp++ ){ 
    
    const  double Amp = (*iamp );
    
    if ( fabs( Amp - ped ) > max ) { 
      
      max = fabs(Amp-ped); 

      msample = isample; } 

    isample++;

  }

  this->setTmax(msample);

  return msample;

}

//*****************************************************************
bool thresholdTime::isOverThreshold(double amp){
//*****************************************************************
  
  this->checkStatus();

  bool ok=false;
  
  //if (_threshold<0) ok = (fabs(_threshold) > amp);
  
  //else ok = (_threshold < amp);
  
  ok = (_threshold > amp);

  return ok;
  
}

//*****************************************************************
void thresholdTime::checkStatus() const{
//*****************************************************************
  
  if (!this->status()){
    
    //Message::MSG("RecoPulseAlgo",kMFATAL,"Threshold not defined");
  }

}

//========================= windowTime =========================//

//*****************************************************************
windowTime::windowTime() : timeDef(){
//*****************************************************************
  
 
  
  _pedestal = NOPED;

  _threshold = 5;

  _startPoint = 0.1; // % MAX 

  _endPoint = 0.1; // % MAX 
  
  _timeAlgo = "threshold";

}

//*****************************************************************
int windowTime::computeTstart(RPPulse& ipulse){
//*****************************************************************

  int T = NOTIME;

  if (_timeAlgo=="max") T = computeMaxTstart(ipulse);
  
  else if (_timeAlgo=="weight") T = computeWTstart(ipulse);

  else if (_timeAlgo=="line") T = computeSLTstart(ipulse);

  else if (_timeAlgo=="threshold") T = computeThTstart(ipulse);

  else{
    
    gate::Assert(false,__FILE__,__LINE__,
		 
		 gate::internal_logic("Not a valid Tstart algo: "+_timeAlgo));
  }
    
  return T;
  
}

//*****************************************************************
int windowTime::computeTend(RPPulse& ipulse){
//*****************************************************************

  int T = NOTIME;

  if (_timeAlgo=="max") T = computeMaxTend(ipulse);
  
  else if (_timeAlgo=="threshold") T = computeThTend(ipulse);

  else{
    
    gate::Assert(false,__FILE__,__LINE__,
		 
		 gate::internal_logic("Not a valid Tend algo: "+_timeAlgo));
  }
    
  return T;
  
}

//*****************************************************************
int windowTime::computeWTstart(RPPulse& ipulse){
//*****************************************************************
  
  /*
    Compute weighted time within [T0,T1]

          Sum ( amplitude x time )
    t = ---------------------------
             Sum ( amplitude ) 
  */

  this->checkWindow();

  this->checkPedestal();

  const vector<unsigned short>& pulse = ipulse.getProfile();

  int tstart = 0;
  
  double num=0; double den=0;

  for (size_t i=(size_t)this->getT0(); i<(size_t)this->getT1(); i++){
    
    num += ( (pulse[i]-this->getPedestal()) * i ) ; 
    
    den +=  ( pulse[i] - this->getPedestal() );

  }

  tstart = (int)(num / den) ; 

  this->setTstart(tstart);

  return tstart;

}

//*****************************************************************
int windowTime::computeSLTstart(RPPulse& ipulse){
//*****************************************************************
   
  if (this->getTmax()==NOTIME) this->computeTmax(ipulse);
  
  const vector<unsigned short>& pulse = ipulse.getProfile();

  int x0 = this->getTmax();
  
  this->setTstart(x0);

  double ped = this->getPedestal();
 
  if (x0){
    
    int x1 = x0-1;
    
    double y0 = fabs( pulse[x0]-ped );
    
    double y1 = fabs( pulse[x1]-ped );
    
    double tstart = x0 - y0*(x1-x0)/(y1-y0);
    
    this->setTstart((int)tstart);

  }
  
  return this->getTstart();

}

//*****************************************************************
int windowTime::computeMaxTstart(RPPulse& ipulse){
//*****************************************************************
  
  if (this->getTmax()==NOTIME) this->computeTmax(ipulse);
  
  const vector<unsigned short>& pulse = ipulse.getProfile();

  size_t fsample = this->getTmax();
  
  this->setTstart(fsample);

  double ped = this->getPedestal();

  double max = fabs( pulse[fsample]-ped );
 
  for (size_t i=fsample; i>0; i--){
        
    //if ( fabs(pulse[i]) > (1+_startPoint) * max){
    
    if ( fabs(pulse[i]-ped) < max * _startPoint){
      
      this->setTstart(i); break; }  
  }
  
  return this->getTstart();

}

//*****************************************************************
int windowTime::computeThTstart(RPPulse& ipulse){
//*****************************************************************
   
  this->checkWindow();
  
  vector<unsigned short> pulse = ipulse.getProfile();

  int isample=this->getT0(); bool found = false;
  
  vector<unsigned short>::iterator iamp;
  
  for( iamp = pulse.begin()+this->getT0(); 
       iamp != pulse.begin()+this->getT1(); iamp++ ){ 
    
    const  double Amp = (*iamp); bool ok;
   
    ok = (this->getThreshold() > Amp);
  
    if (ok) { found = true; break; }
    
    isample++;

  } 
  if (!found) isample = this->getT0();
  
  this->setTstart(isample);

  return isample;

  
}



//*****************************************************************
int windowTime::computeMaxTend(RPPulse& ipulse){
//*****************************************************************
  
  if (this->getTmax()==NOTIME) this->computeTmax(ipulse);
  
  const vector<unsigned short>& pulse = ipulse.getProfile();

  int fsample = this->getTmax();
  
  if (fsample<1) return this->getTend();
  
  double ped = this->getPedestal();
  
  double max = fabs( pulse[fsample]-ped );

  const int t1 = this->getT1();

  this->setTend(t1);
  
  if (t1==0 || t1==NOTIME) return this->getTend();
  
  for (int i=t1-1; i> fsample ; i--){
    
    if ( fabs(pulse[i]-ped) > max * _endPoint){

      this->setTend(i+1); break; } 
  }
  

  return this->getTend();
  
}

//*****************************************************************
int windowTime::computeThTend(RPPulse& ipulse){
//*****************************************************************
   
  this->checkWindow();
  
  vector<unsigned short> pulse = ipulse.getProfile();

  int isample=this->getT1(); bool found = false;
  
  vector<unsigned short>::iterator iamp;
  
  for( iamp = pulse.begin()+this->getT1(); 
       iamp != pulse.begin()+this->getT0(); iamp-- ){ 
    
    const  double Amp = (*iamp); bool ok;
   

    ok = (this->getThreshold() > Amp);
  
    if (ok) { found = true; break; }
    
    isample--;

  } 
  if (!found) isample = this->getT0();
  
  this->setTstart(isample);

  return isample;

  
}




//*****************************************************************
int windowTime::computeTmax(RPPulse& ipulse){
//*****************************************************************

  /*
    Returns sample index of pulse with max amp
  */
  

 
  this->checkWindow();

  this->checkPedestal();

  int fsample=this->getT0(); int lsample=this->getT1();
  
  int msample = fsample; int isample = fsample; 
 
  const vector<unsigned short>& pulse = ipulse.getProfile();
 
  double ped = this->getPedestal();

  double max = fabs(pulse[fsample]-ped);
 
  vector<unsigned short>::const_iterator iamp;
 

  for( iamp = pulse.begin() + fsample; 
       iamp != pulse.begin() + lsample; 
       iamp++ ){ 
    
    const  double Amp = (*iamp );
   
    if ( fabs( Amp - ped ) > max ) { 
      
      max = fabs(Amp-ped); 

      msample = isample; } 

    isample++;

  }

  this->setTmax(msample);
  

  return msample;

}


//*****************************************************************
void windowTime::checkWindow(){
//*****************************************************************
    
  bool status = (this->getT0()==NOTIME);
  
  gate::Assert(status==false,__FILE__,__LINE__,
		 
	       gate::internal_logic("Window (T0) not defined"));
  
  status = (this->getT1()==NOTIME);
  
  gate::Assert(status==false,__FILE__,__LINE__,
		 
	       gate::internal_logic("Window (T1) not defined"));
  
}

//*****************************************************************
void windowTime::checkPedestal(){
//*****************************************************************
  
  bool status = (this->getPedestal()==NOPED);
  
  gate::Assert(status==false,__FILE__,__LINE__,
		 
	       gate::internal_logic("Pedestal not defined"));
}



//========================= maxAmpTime =========================//

//*****************************************************************
maxAmpTime::maxAmpTime() : windowTime(){
//*****************************************************************
    
    
  _lWindowSize = 0; _rWindowSize = 0;

}



//*****************************************************************
int maxAmpTime::computeTmax(RPPulse& ipulse){
//*****************************************************************
  
  const vector<unsigned short>& pulse = ipulse.getProfile();
    
  int fsample = 0; int lsample = pulse.size();

  int msample = fsample; int isample = fsample; 


  double ped = this->getPedestal();

  double max = fabs(pulse[fsample]-ped);

  vector<unsigned short>::const_iterator iamp;
  
  for( iamp = pulse.begin() + fsample; 
       iamp != pulse.begin() + lsample; 
       iamp++ ){ 
    
    const  double Amp = (*iamp );
    
    if ( fabs( Amp - ped ) > max ) { 
      
      max = fabs(Amp-ped); 

      msample = isample; } 

    isample++;

  }

  this->setTmax(msample);
  
  return msample;

}



//*****************************************************************
int maxAmpTime::computeT1(RPPulse& ipulse){
//*****************************************************************
  
 
  if (this->getTmax()==NOTIME) this->computeTmax(ipulse);

  //if (this->getTend()==NOTIME) this->computeTend(ipulse);
    
  int t1 = this->getTmax()+_rWindowSize;
  //int t1 = max((int)(this->getTmax()+_rWindowSize),this->getTend());
  
  if ( t1 > (int)ipulse.size() ) t1 = ipulse.size();

  this->setT1(t1);
  
  return this->getT1();

}

//*****************************************************************
int maxAmpTime::computeT0(RPPulse& ipulse){
//*****************************************************************
  
  
  if (this->getTmax()==NOTIME) this->computeTmax(ipulse);
  
  int t0 = this->getTmax()-_lWindowSize;
  
  if (t0<0) t0 = 0;

  this->setT0(t0);
  
  return this->getT0();

  
}

//*****************************************************************
void maxAmpTime::setWindow(int left,int right){
//*****************************************************************
  
  /*
    Set window size: t = [max-left,max+right]
   */
  
  _lWindowSize = left; _rWindowSize = right;
  
}


//======================= RecoPulseTime ========================//


//*****************************************************************
RecoPulseTime::RecoPulseTime(){
//*****************************************************************

  _timeDef="windowTime";

  _tools["thresholdTime"]= new thresholdTime();

  _tools["windowTime"]= new windowTime();

  _tools["maxAmpTime"]= new maxAmpTime();

  //_tools["splineTime"]= new splineTime();
  
  _ptrtimeDef = dynamic_cast<timeDef*>(_tools[_timeDef]);
}

//*****************************************************************
RecoPulseTime::~RecoPulseTime(){
//*****************************************************************
  
  map<std::string,timeDef*>::iterator iter;   

  for( iter = _tools.begin(); iter != _tools.end(); iter++ ) {
    delete iter->second;
  }
}


//*****************************************************************
void RecoPulseTime::config(std::string timedef){
//*****************************************************************
  
    
  this->checkTimeDef(timedef);
  
  _timeDef=timedef;

  _ptrtimeDef = dynamic_cast<timeDef*>(_tools[_timeDef]);
}

//*****************************************************************
timeDef* RecoPulseTime::getTimeDef(){
//*****************************************************************
  
  return _ptrtimeDef;

}

//*****************************************************************
int RecoPulseTime::computeT0(RPPulse& ipulse){
//*****************************************************************
  
  return this->getTimeDef()->computeT0(ipulse);
  
}

//*****************************************************************
int RecoPulseTime::computeT1(RPPulse& ipulse){
//*****************************************************************
  
  return  this->getTimeDef()->computeT1(ipulse);
   
}

//*****************************************************************
int RecoPulseTime::computeTstart(RPPulse& ipulse){
//*****************************************************************
  
  return this->getTimeDef()->computeTstart(ipulse);
  
}

//*****************************************************************
int RecoPulseTime::computeTend(RPPulse& ipulse){
//*****************************************************************
  
  return  this->getTimeDef()->computeTend(ipulse);
   
}

//*****************************************************************
int RecoPulseTime::computeTmax(RPPulse& ipulse){
//*****************************************************************
  
  return  this->getTimeDef()->computeTmax(ipulse);
   
}

//*****************************************************************
int RecoPulseTime::computeTrise(RPPulse& ipulse){
//*****************************************************************
  
  return  this->getTimeDef()->computeTrise(ipulse);
   
}

//*****************************************************************
int RecoPulseTime::computeTfall(RPPulse& ipulse){
//*****************************************************************
  
  return  this->getTimeDef()->computeTfall(ipulse);
   
}

//*****************************************************************
void RecoPulseTime::reset(){
//*****************************************************************
  
  this->getTimeDef()->reset();
   
}


//*****************************************************************
bool RecoPulseTime::find(std::string timedef) const{
//*****************************************************************
  
   //! find an element in the store by name
         
  bool gotcha = false;
  
  map<std::string, timeDef*>::const_iterator pi;
      
  for ( pi = _tools.begin(); 
	pi!= _tools.end(); ++pi)
    {
      if ( (pi->first) == timedef) {
	gotcha = true;
	break;
      }
    } 

  return gotcha;
  
}

//*****************************************************************
void RecoPulseTime::checkTimeDef(std::string timedef) const{
//*****************************************************************

  bool ok=this->find(timedef);
   
  gate::Assert(ok==true,__FILE__,__LINE__,
		 
	       gate::internal_logic("TimeDef not found: "+timedef));
}


