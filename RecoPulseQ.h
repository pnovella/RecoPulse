/**
 * @file DCRecoPulseQ.hh
 * Charge reconstruction algorithms
 * 
 * @class RecoPulseQ
 *
 * @ingroup RecoPulse
 *
 * @brief Charge reconstruction algorithms
 *
 * Main Algorithm classes for pulse charge reconstrution.  
 * RecoPulseQ class is an storage of an arbitrary number of
 * charge measuremt techniques (also classes). 
 *
 * @author Pau Novella, Carmen Palomares $Author: pnovella $
 *
 * @version $Revision: 1.5 $
 *
 * @date $Date: 2009-10-28 08:01:17 $
 *
 * Contact: pau.novella@ciemat.es
 *
 */

#ifndef DCRecoPulseQ_hh
#define DCRecoPulseQ_hh

#include<map>
#include<string>
#include<vector>

#include<Messenger.h>
#include<RecoPulseUtils.h>


class Qint {

private:
  
  size_t _tmin, _tmax;
  
  double _pedestal;
  
  double _Q;
  

protected:


public:

  Qint();

  virtual ~Qint(){};
  
  virtual double computeQ(RPPulse&) = 0;

  virtual void setWindow(size_t tmin,size_t tmax);
  
  virtual void reset(){_Q=0;_tmin=0;_tmax=0;}
 

  void setPedestal(double ped){_pedestal=ped;}
  
  void setTmin(size_t tmin) {_tmin=(int)tmin;}

  void setTmax(size_t tmax) {_tmax=(int)tmax;}
  
  void setQ(double Q){ _Q = Q;}

  virtual void setWindowSize(size_t){};

  double getQ() const{return _Q;}
  
  double getPedestal()const{return _pedestal;}
  
  virtual bool status() const {return (_pedestal!=NOPED);}
  
  size_t tmin() const {return _tmin;}

  size_t tmax() const {return _tmax;}
  
  size_t getWindowSize() const {return _tmax - _tmin;}
  

protected:
  
  void checkWindow(size_t,size_t) const;
  
  void checkPedestal() const;

 
};


class fixedTimeQWindow : public Qint {
    
public:
  
  fixedTimeQWindow();
  
  double computeQ(RPPulse&);

  double computeQ(RPPulse&,size_t tmin,size_t tmax);
  
};


class slidingQWindow : public Qint {
  

private:

  int _wsize;

  int _step;
  
  int _fstep;
  
  

public:
  
  slidingQWindow();
  
  double computeQ(RPPulse&);
  
  virtual double computeWindowQ(RPPulse&,size_t tmin,size_t tmax);
  
  void setWindowSize(size_t wsize) {_wsize=(int)wsize;}

  void setWindowStep(size_t step) {_step=(int)step;}
  
  

private:
  
  void checkStatus();

  void setWindow(size_t tmin,size_t tmax);
 
};

class slidingSPQWindow : public slidingQWindow {
  
public:
  
  slidingSPQWindow();

  double computeWindowQ(RPPulse&,size_t tmin,size_t tmax);

};

class splineQ : public Qint {

public:
  
  splineQ();

  double computeQ(RPPulse&);

}; 


class RecoPulseQ {

private:

  std::map<std::string, Qint*> _tools; 
  
  std::string _Qint;

  Qint * _ptrQint;
  

public:
  
  RecoPulseQ();

  virtual ~RecoPulseQ();

  void config(std::string qint);

  double computeQ(RPPulse&);

  void setPedestal(double p);

  void setWindow(size_t,size_t);

  void setWindowSize(size_t);

  void reset();

  Qint* getQint();
  
  bool find(std::string) const;

  size_t getWindowSize(); 
  
  size_t getTmin();

  size_t getTmax();

private:

  void checkQint(std::string) const;

};

#endif
