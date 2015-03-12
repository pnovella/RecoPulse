

#include<RecoPulseUtils.h>

using namespace std;

//****************************************************
spline3::spline3(){
//****************************************************

  _nodes.resize(0);
  
  _values.resize(0);

  _status = false;

}

//****************************************************
spline3::spline3(const vector<double>& y){
//****************************************************
  
  vector<double> x;
  
  for (size_t i = 0; i< y.size(); i++) x.push_back(i);

  _nodes = x;

  _values = y;
  
  _status = false;

  //this->gen();

}

//****************************************************
spline3::spline3(const vector<unsigned short>& y){
//****************************************************
  
  vector<double> new_y; 

  new_y = vector<double>(y.begin(),y.end());;
  
  vector<double> x;
  
  for (size_t i = 0; i< new_y.size(); i++) x.push_back(i);

  _nodes = x;

  _values = new_y;
  
  _status = false;

}



//****************************************************
spline3::spline3(const vector<double>& x, const vector<double>& y){
//****************************************************
  
  _nodes = x;
 
  _values = y;
  
  _status = false;
  
  //this->gen();

}



//*****************************************************
void spline3::gen(){
//*****************************************************
  
  /*
    Natural splines implementation
   */
  
  //if (!_status) return; 
  
 
  _acoefs = _values; //natural

  size_t n = _nodes.size();                               
  
  size_t nint = n - 1 ;
  
  size_t neq = nint - 1; 
  
  //resize coef vectors

  _bcoefs.resize(nint); _ccoefs.resize(nint); _dcoefs.resize(nint);
  
  //upper diagonal vector
  
  vector<double> h(nint);
  
  for(size_t i=0; i<nint; i++) h[i]=_nodes[i+1]-_nodes[i];

  //diagonal vector
  
  vector<double> u(nint), v(nint);  //u(0), v(0) not defined
  for( size_t i=1; i<nint; i++) u[i]=2*(h[i-1]+h[i]);
  for( size_t i=0; i<neq; i++)
    v[i+1]=3*( (_values[i+2]-_values[i+1])/h[i+1] - 
	       (_values[i+1]-_values[i])/h[i] );

  //Gauss
  for( size_t i=1; i<neq; i++)
    {
      u[i+1]=u[i+1]-h[i]*h[i]/u[i];
      v[i+1]=v[i+1]-h[i]*v[i]/u[i];
    }

  //Substitution
  _ccoefs[nint-1]=v[nint-1]/u[nint-1];
  _dcoefs[nint-1]= -_ccoefs[nint-1]/(3*h[nint-1]);
  for( size_t k=nint-2; k>0; k--)
    {
      _ccoefs[k]=(v[k]-h[k]*_ccoefs[k+1])/u[k];
      _dcoefs[k]=(_ccoefs[k+1]-_ccoefs[k])/(3*h[k]);
    }

  //values at _nodes[0]
  //_ccoefs[0]= T(0);
  _ccoefs[0]= 0;
  _dcoefs[0]=_ccoefs[1]/(3*h[0]);
  _bcoefs[0]=(_values[1]-_values[0])/h[0]-_ccoefs[1]*h[0]/3;

  //compute b[k]
  for( size_t k=0; k<nint-1; k++)  
    _bcoefs[k+1]=_bcoefs[k]+(_ccoefs[k]+_ccoefs[k+1])*h[k];
  
  
  
  _status = true; 
}


//*****************************************************
int spline3::getIndex(double x){
//*****************************************************

  //search index k
  
  int n = _nodes.size();
  int k = 0;
  int kinf = 0;
  int ksup = n - 1;

  while(ksup-kinf>1)
    {
      k=(ksup+kinf)/2;
      if (_nodes[k]>x) ksup=k;
      else kinf=k;
    }
  k=kinf;

  return k;

}


//*****************************************************
double spline3::eval(double x){
//*****************************************************

  if (!_status) return 0; 

  int k = this->getIndex(x);
  
  x=x-_nodes[k];
  
  return _acoefs[k]+ x*( _bcoefs[k]+ x*( _ccoefs[k] +_dcoefs[k]*x ) );

}

//*****************************************************
double spline3::integral(double x, size_t index){
//*****************************************************
  
  x=x-_nodes[index];

  return _acoefs[index]* x + 
    _bcoefs[index] * pow(x,2)/2 + 
    _ccoefs[index] * pow(x,3)/3 + 
    _dcoefs[index] * pow(x,4)/4;
  
}

//*****************************************************
double spline3::integrate(size_t min,size_t max){
//*****************************************************
  
  if (!_status) return -1; 
  
  if (!max) max = _nodes.size()-1; 
  
  else{ if (max>_nodes.size()-1) max = _nodes.size()-1;}
  
  double Int = 0;
  
  for (size_t i = min; i < max; i += 1 ){
    
    double intX0 = this->integral(_nodes[i],i);

    double intX1 = this->integral(_nodes[i+1],i);
    
    double qi = (intX1- intX0);

    Int +=qi; 
    
  }
 
  return Int;

}



//*****************************************************
void peakGatherer::getPeaks(double thr,size_t nsamp){
//*****************************************************
  
  _T0s.clear(); _T1s.clear();

  size_t nsamppeak = 0; 
  
  //cout<<thr<<" "<<nsamp<<endl;

  for (size_t i = 0; i < _prof.size(); i++){
    
    if ( (thr - _prof[i]) > 0 ){ // new peak ?
      
      nsamppeak++; 
            
    }
    
    else{ // no peak or end of peak
      
      if (nsamppeak >= nsamp){ // end of peak
	
	_T0s.push_back(i-nsamppeak); // set T0
	
	_T1s.push_back(i); // set T1

      }

      nsamppeak = 0;
    }
    
  }
  
  if (nsamppeak >= nsamp){ // unfinished peak
  
    	_T0s.push_back(_prof.size()-nsamppeak); // set T0
	
	_T1s.push_back(_prof.size()); // set T1
  }
}

//*****************************************************
void peakGatherer::getBPPeaks(double thr,size_t nsamp){
//*****************************************************
  
  _T0s.clear(); _T1s.clear();

  size_t nsamppeak_pos = 0; 

  size_t nsamppeak_neg = 0; 
  
  for (size_t i = 0; i < _prof.size(); i++){
    
    if ( (round(thr) - _prof[i]) > 0 ){ // new peak positive peak?
      
      nsamppeak_pos++; 
      
      // comes from a negative pulse?

      if (nsamppeak_neg >= nsamp){ // end of peak
	
	_T0s.push_back(i-nsamppeak_neg); // set T0
	
	_T1s.push_back(i);
      }
      nsamppeak_neg=0; 
    }
    
    else if ( (_prof[i] - round(thr)) > 0 ){ // new negative peak ?
      
      nsamppeak_neg++; 
      
      // comes from positive pulse?

      if (nsamppeak_pos >= nsamp){ // end of peak
	
	_T0s.push_back(i-nsamppeak_pos); // set T0
	
	_T1s.push_back(i);
      }
    
      nsamppeak_pos=0; 
    }
    
    else { // on the threshold: look if any peak ended here
      
      if ( (nsamppeak_pos >= nsamp) ){
	
	_T0s.push_back(i-nsamppeak_pos); // set T0
	
	_T1s.push_back(i); // set T1

      }
      else if (nsamppeak_neg >= nsamp){
	
	_T0s.push_back(i-nsamppeak_neg); // set T0
	
	_T1s.push_back(i); // set T1
      }
      
      nsamppeak_neg = 0; nsamppeak_pos = 0;
    }
    
  }// end of waveform loop
  
  if (nsamppeak_pos >= nsamp){ // unfinished peak
  
    	_T0s.push_back(_prof.size()-nsamppeak_pos); // set T0
	
	_T1s.push_back(_prof.size()); // set T1
  }

  else if (nsamppeak_neg >= nsamp){ // unfinished peak
  
    	_T0s.push_back(_prof.size()-nsamppeak_neg); // set T0
	
	_T1s.push_back(_prof.size()); // set T1
  }

}

//*****************************************************
const vector<unsigned short>& RPPulse::getSubProfile(double ped){
//*****************************************************
  
  _subProf.resize(0);

  for (size_t i = 0; i < this->size(); i++){
    
    _subProf.push_back((int)round(ped) - _prof[i]);

  }
  
  return _subProf;

}

//*****************************************************
void RPPulse::smoothPed(double ped,size_t min, size_t max){
//*****************************************************
  
  /*
    Set to pedestal level a given region of profile.
   */


  vector<unsigned short> newProf;
  
  for (size_t i=0; i< this->size(); i++){ 
    
    if (i>=min && i<max) newProf.push_back((int)round(ped)); 
  
    else newProf.push_back(_prof[i]); 

  }
  
  _prof = newProf;

  //_spline = spline3(_prof);

  _peaks = peakGatherer(_prof);

}

//*****************************************************
void RPPulse::integrated(size_t min, size_t max){
//*****************************************************
  
  for (size_t i=min; i<max;i++)

    {_integrated[i] = true;}
  
  return;

}


//*****************************************************
void RPPulse::filter(double ped,double thr, size_t nsamples){
//*****************************************************
  
  /*
    Remove from waveform spureous spikes
   */

  
  _peaks.getPeaks(thr,nsamples);
  
  vector<int> T0s = _peaks.getT0s();

  vector<int> T1s = _peaks.getT1s();
  
  // no peaks
  if (!T0s.size()) return this->smoothPed(ped,0,_prof.size());
  
  this->smoothPed(ped,0,T0s[0]);// before first peak
  
  this->smoothPed(ped,T1s[T1s.size()-1],_prof.size()); //after last peak

  for(size_t i=0;  i<T0s.size()-1; i++){ // between peaks
    
    this->smoothPed(ped,T1s[i],T0s[i+1]);
  }

}

//*****************************************************
bool RPPulse::hasPeak(double ped, double thr, size_t nsamples,
			   int t0, int t1){
//*****************************************************
  
  
  thr = ped - thr; 

  
  size_t sc = 0;

  for (int i=t0; i<t1; i++){
    
    double sub = (ped - _prof[i]);
    
    if (sub>thr){ 
    
      sc++;
      
      if (sc>=nsamples) return true;
    }

    else sc=0;

  }
    
  return false;

}

//*****************************************************
size_t RPPulse::getContSamples(double ped,size_t t0, size_t t1){
//*****************************************************
  
  /* compute maximum number of contiguous samples over pedestal */

  if (!t1) t1 = _prof.size();

  size_t nCS = 0; size_t nCSmax = 0;

  const double rped = round(ped);
  for (size_t i=t0; i<t1; i++){
    
    if (rped - _prof[i]>0){ nCS++; }

    else{ if (nCS>nCSmax) nCSmax=nCS ; nCS=0;}    
  }
  if (nCS) {if (nCS>nCSmax) nCSmax=nCS;}

  return nCSmax;
  
}

//*****************************************************
size_t RPPulse::getPeakWidth(double ped, double thr,int t0, int t1){
//*****************************************************
  
  /* compute maximum number of contiguous samples over threshold */
  
  /* TO REVIEW !!!!!!!!!!!!!!!!!!!!!*/

  if (!t0) t1 = _prof.size();
  
  thr = ped - thr;

  size_t nCS = 0; size_t nCSmax = 0;

  for (int i=t0; i<t1; i++){
    
    double sub = (ped - _prof[i]);
    
    if (sub>thr){ nCS++; }

    else{ if (nCS>nCSmax) nCSmax=nCS ; nCS=0;}
  }
  
  return nCSmax;
  
}
