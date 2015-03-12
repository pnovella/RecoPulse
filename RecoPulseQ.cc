///
/// \file DCRecoPulseQ.cc
/// Charge reconstruction algorithms 
///
/// \author Pau Novella - CIEMAT 2008
///

#include<RecoPulseQ.h>
#include<math.h>
#include <stdlib.h>

using namespace std;

//==========================Qint=========================//

//**************************************************************
Qint::Qint(){
//**************************************************************
  

  this->reset();

}


//**************************************************************
void Qint::setWindow(size_t tmin,size_t tmax) {
//*************************************************************
  
  this->checkWindow(tmin,tmax);

  _tmin=tmin;_tmax=tmax;

}

//**************************************************************
void Qint::checkWindow(size_t tmin,size_t tmax) const{
//*************************************************************
  
  if (tmin>tmax && tmax!=0){
   
    //Message::MSG("RecoPulseAlgo",kMFATAL,"Q sample window not well defined");
       
  }
  
  if ((int)tmin==NOTIME || (int)tmax==NOTIME){
    
    //Message::MSG("RecoPulseAlgo",kMFATAL,"Q sample window not defined");

  }

}

//**************************************************************
inline void Qint::checkPedestal() const{
//*************************************************************
  
  if (this->getPedestal()==NOPED){
    
    //Message::MSG("RecoPulseAlgo",kMFATAL,"Pedestal not defined");
  }
  

}


//========================= fixedQTimeWindow =======================//

//*****************************************************************
fixedTimeQWindow::fixedTimeQWindow():Qint(){
//*****************************************************************
  
 
}


//*****************************************************************
double fixedTimeQWindow::computeQ(RPPulse& ipulse){
//*****************************************************************
  
  //vector<double> pulse = ipulse.getProfile();

  return this->computeQ(ipulse,this->tmin(),this->tmax());

}

//*****************************************************************
double fixedTimeQWindow::computeQ(RPPulse& ipulse,size_t tmin,size_t tmax){
//*****************************************************************
  
  this->checkPedestal();

  if (tmin || tmax) this->setWindow(tmin,tmax);
  
  size_t max = this->tmax();

  size_t min = this->tmin();
  
  const vector<unsigned short>& pulse = ipulse.getProfile();

  if (!max) max = pulse.size();  
  
  double Int=0;
  
  for (size_t i=min; i<max; i++){ 

    if (ipulse.isIntegrated(i)) continue;

    Int += ( this->getPedestal() - pulse[i] ); 
    
  }  
  
  double Q = Int; //if (Q<0) Q = 0; !!!
  
  this->setQ(Q);

  return Q;

}

//========================= slidingQWindow =======================//


//*****************************************************************
slidingQWindow::slidingQWindow(): Qint(){
//*****************************************************************
  
  _wsize = NOTIME; _step = 1; 

 
}

//*****************************************************************
double slidingQWindow::computeQ(RPPulse& ipulse){
//*****************************************************************

  this->checkStatus();
  
  const int maxWindowMin = (int)ipulse.size()-_wsize;
  int tmin = this->tmin(), tminMax = this->tmin();

  int tmax = tmin + _wsize;

  // Find the Q of the first window position. We allow rawQ to be
  // negative and carry through negative values so that the integrated
  // charge of the shifted windows is not a function of the path we used
  // to get to them. However, we never report a negative Q.
  double rawQ = this->computeWindowQ(ipulse, tmin, tmax);
  double Q = rawQ > 0? rawQ: 0;
  double Qmax = Q;

  tmin += _step; 
  tmax += _step; 
  
  // We've done one full computation of rawQ for the first window.
  // Now find the Qs of the shifted windows by a cheaper method:
  // subtract from one end and add to the other.
  const vector<unsigned short>& pulse = ipulse.getProfile();

  const double ped = this->getPedestal();

  while(tmin <= maxWindowMin){
 
    // Keep track of the difference in the number of samples masked
    // out because they are already used by another pulse.  We don't
    // add or subtract the value of the waveform for these, and we 
    // have to correct the amount of pedestal subtracted.
    int peddiff = 0;

    for(int i = 0; i < _step; i++){
      if(ipulse.isIntegrated(tmin - _step + i)) peddiff++;
      else  rawQ  +=  pulse [tmin - _step + i];

      if(ipulse.isIntegrated(tmax - _step + i)) peddiff--;
      else  rawQ  -=  pulse [tmax - _step + i];
    }
    rawQ += peddiff * ped; 

    Q = rawQ > 0? rawQ: 0;

    if(Q > Qmax){
      Qmax = Q;
      tminMax = tmin;
    } 
    tmin += _step; 

    tmax += _step; 
  }
  
  this->setQ(Qmax);

  this->setWindow(tminMax,tminMax+_wsize);

  return this->getQ();
  
}

//*****************************************************************
double slidingQWindow::computeWindowQ(RPPulse& ipulse,
					size_t tmin,size_t tmax){
//*****************************************************************
  
  this->checkPedestal();
   
  const vector<unsigned short>& pulse = ipulse.getProfile();

  double pulsesum = 0;
  size_t integratedcount = 0;

  for (size_t i=tmin; i<tmax; i++){ 
    // Don't count samples already used by another reconstructed pulse
    if(ipulse.isIntegrated(i)){
      integratedcount++;
      continue;
    }
    pulsesum -= pulse[i]; 
  }  


  return pulsesum + (tmax-tmin - integratedcount)*this->getPedestal();
}

//*****************************************************************
void slidingQWindow::checkStatus(){
//*****************************************************************
  
  if (_wsize==NOTIME){

    //Message::MSG("RecoPulseAlgo",kMFATAL,"Window size not defined!");

  }

}

//**************************************************************
void slidingQWindow::setWindow(size_t tmin,size_t tmax) {
//*************************************************************
  
  this->checkWindow(tmin,tmax);

  this->setTmin(tmin); this->setTmax(tmax);

}

//========================= SlidingSPQWindow =======================//


//*****************************************************************
slidingSPQWindow::slidingSPQWindow(): slidingQWindow(){
//*****************************************************************
  

}

//*****************************************************************
double slidingSPQWindow::computeWindowQ(RPPulse& ipulse,
			       size_t tmin,size_t tmax){
//*****************************************************************
  
 
  this->checkPedestal();
  
  if (!ipulse.getSpline().status()) ipulse.genSpline();
  
  if (!tmin) tmin = 1;

  double qint = ipulse.getSpline().integrate(tmin-1,tmax-1);
  
  double qped = this->getPedestal() * (tmax-tmin);
  
  double Q = qped - qint; if (Q<0) Q = 0;
  
  this->setQ(Q);
 
  return Q;
  
}



//=========================== splineQ =========================//


//*****************************************************************
splineQ::splineQ(): Qint(){
//*****************************************************************
   
}

//*****************************************************************
double splineQ::computeQ(RPPulse& ipulse){
//*****************************************************************
  
  if (!ipulse.getSpline().status()) ipulse.genSpline();
  
  //if (!this->tmax()) this->setTmax(ipulse.size());//??
  
  //if (this->tmax()==ipulse.size()) this->setTmax(ipulse.size()-1);//??
  
  int tmin = this->tmin(); if (!tmin) this->setTmin(1);

  double qint = ipulse.getSpline().integrate(this->tmin()-1,this->tmax()-1);
  
  double qped = this->getPedestal() * (this->getWindowSize());
  
  double Q = qped - qint; if (Q<0) Q = 0;

  this->setQ(Q);
 
  return Q;
  
}


//=========================== RecoPulseQ =========================//


//*****************************************************************
RecoPulseQ::RecoPulseQ(){
//*****************************************************************
  

  _Qint = "fixedTimeQWindow";

  _tools["fixedTimeQWindow"]= new fixedTimeQWindow();
  
  _tools["slidingQWindow"]= new slidingQWindow();

  //_tools["slidingSPQWindow"]= new slidingSPQWindow();

  //_tools["splineQ"]= new splineQ();
  
  _ptrQint = dynamic_cast<Qint*>(_tools[_Qint]);
}

//*****************************************************************
RecoPulseQ::~RecoPulseQ(){
//*****************************************************************
  
  map<std::string,Qint*>::iterator iter;   

  for( iter = _tools.begin(); iter != _tools.end(); iter++ ) {
    delete iter->second;
  }
}



//*****************************************************************
void RecoPulseQ::reset(){
//*****************************************************************
  
  return this->getQint()->reset();

}

//*****************************************************************
void RecoPulseQ::config(std::string qint){
//*****************************************************************
  
 
  this->checkQint(qint);
  
  _Qint=qint;

  _ptrQint = dynamic_cast<Qint*>(_tools[_Qint]);
}



//*****************************************************************
Qint* RecoPulseQ::getQint(){
//*****************************************************************

  return _ptrQint; 

}

//*****************************************************************
void RecoPulseQ::setPedestal(double p){
//*****************************************************************
  
  map<std::string,Qint*>::iterator iter;   

  for( iter = _tools.begin(); iter != _tools.end(); iter++ ) {
    iter->second->setPedestal(p);
  }
  //this->getQint()->setPedestal(p);
  
}

//*****************************************************************
void RecoPulseQ::setWindow(size_t tmin,size_t tmax){
//*****************************************************************
  
  this->getQint()->setWindow(tmin,tmax);
  
}

//*****************************************************************
void RecoPulseQ::setWindowSize(size_t size){
//*****************************************************************
  
  map<std::string,Qint*>::iterator iter;   

  for( iter = _tools.begin(); iter != _tools.end(); iter++ ) {
    iter->second->setWindowSize(size);
  }
  //this->getQint()->setWindowSize(size);
  
}

//*****************************************************************
double RecoPulseQ::computeQ(RPPulse& ipulse){
//*****************************************************************
  
  return this->getQint()->computeQ(ipulse);

}

//*****************************************************************
size_t RecoPulseQ::getWindowSize(){
//*****************************************************************
  
  return this->getQint()->getWindowSize();

}

//*****************************************************************
size_t RecoPulseQ::getTmin() {
//*****************************************************************
  
  return this->getQint()->tmin();

}

//*****************************************************************
size_t RecoPulseQ::getTmax() {
//*****************************************************************
  
  return this->getQint()->tmax();

}


//*****************************************************************
bool RecoPulseQ::find(std::string timedef) const{
//*****************************************************************
  
   //! find an element in the store by name
         
  bool gotcha = false;
  
  map<std::string, Qint*>::const_iterator pi;
      
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
void RecoPulseQ::checkQint(std::string qint) const{
//*****************************************************************

  bool ok=this->find(qint);
  
  if (!ok){ 
    
    //Message::MSG("RecoPulseAlgo",kMFATAL,"Qint not found: "+qint);
    
  }
}
 
