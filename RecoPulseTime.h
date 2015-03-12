/**
 * @file DCRecoPulseTime.hh
 * Time reconstruction algorithms
 *
 * @class RecoPulseTime
 *
 * @ingroup RecoPulse
 *
 * @brief Time reconstruction algorithms
 *
 * Main Algorithm classes for pulse time reconstrution.  
 * RecoPulseTime class is an storage of an arbitrary number of
 * time measuremt techniques (also classes). 
 *
 * @author Pau Novella, Carmen Palomares $Author: pnovella $
 *
 * @version $Revision: 1.8 $
 *
 * @date $Date: 2009-10-28 08:01:17 $
 *
 * Contact: pau.novella@ciemat.es
 *
 */

#ifndef DCRecoPulseTime_hh
#define DCRecoPulseTime_hh

#include<map>
#include<string>
#include<vector>
#include<Messenger.h>
#include<RecoPulseUtils.h>
#include<Error.h>

class timeDef {
  
private:

  int _T0,_T1; // window bounds

  int _Tstart,_Tend,_Tmax; 
  
  int _Trise, _Tfall; 
  
 
protected:
  
  std::string _timeAlgo;

public:
  
  
  timeDef();
  
  void reset(){_T0=NOTIME;_T1=NOTIME;_Tmax=NOTIME;
    _Tstart=NOTIME;_Tend=NOTIME;_Trise=NOTIME;_Tfall=NOTIME;}

  virtual ~timeDef(){};

  //mandatory
  virtual int computeTstart(RPPulse&) = 0;
  
  virtual int computeTend(RPPulse&) = 0 ;
  //

  virtual int computeT0(RPPulse&) {return NOTIME;}

  virtual int computeT1(RPPulse&) {return NOTIME;}

  virtual int computeTmax(RPPulse&){return NOTIME;}

  virtual int computeTrise(RPPulse&);

  virtual int computeTfall(RPPulse&);
  
  //
  virtual void setThreshold(double){};
  
  virtual void setPedestal(double){};

  

  void setT0(int T0){_T0=T0;}

  void setT1(int T1){_T1=T1;}
  
  virtual void setWindow(int tmin,int tmax){
    this->setT0(tmin); this->setT1(tmax);}

  void setTstart(int Tstart){_Tstart=Tstart;}

  void setTend(int Tend){_Tend=Tend;}

  void setTmax(int T){_Tmax=T;}

  void setTrise(int T){_Trise=T;}

  void setTfall(int T){_Tfall=T;}

  void setTstartAlgo(std::string algo){_timeAlgo=algo;}
  
  int getT0() const {return _T0;}
  
  int getT1() const {return _T1;}

  int getTstart() const {return _Tstart;}
  
  int getTend() const {return _Tend;}

  int getTmax() const {return _Tmax;}

  int getTise() const {return _Trise;}

  int getTfall() const {return _Tfall;}
  
  std::string getTstartAlgo()const {return _timeAlgo;}

};


class thresholdTime: public timeDef{

  /*
    T0 = first sample over threshold

    T1 = first sample >T0 under threshold
    
    Tstart = T0

    Tend = T1

   */

private:
  
  double _threshold;
 
public: 
  
  thresholdTime(double thr=5);
    
  virtual ~thresholdTime(){};

  int computeT0(RPPulse&);

  int computeT1(RPPulse&);

  virtual int computeTstart(RPPulse&);

  virtual int computeTend(RPPulse&);

  int computeTmax(RPPulse&);

  bool isOverThreshold(double amp);

  double getThreshold()const{return _threshold;}

  bool status() const {return (_threshold!=NOPED);}
  
  void setThreshold(double thr){_threshold=thr;}


private:
  
  void checkStatus() const;

};


class windowTime : public timeDef{

  /*
    Time definition according to a given window [T0,T1] 
    where pulse has been found
   */
  
private:

  double _pedestal;
  
  double _threshold;

  double _endPoint;

  double _startPoint;

public:

  windowTime();

  int computeTstart(RPPulse&);

  int computeWTstart(RPPulse&);

  int computeSLTstart(RPPulse&);

  int computeMaxTstart(RPPulse&);
  
  int computeThTstart(RPPulse&);
  
  int computeTend(RPPulse&);
  
  int computeMaxTend(RPPulse&);

  int computeThTend(RPPulse&);

  int computeTmax(RPPulse&);

  void setPedestal(double mean){_pedestal=mean;}

  double getPedestal() const{return _pedestal;}

  double getThreshold() const{return _threshold;}

  void setThreshold(double thr){_threshold=thr;}

private:

  void checkWindow();

  void checkPedestal();

  //void checkThreshold();
 
};

class maxAmpTime : public windowTime {
  
  /*
    Time estimation starting from pulse maximum.
   */

private:

  double _startPoint;

  double _endPoint;
  
  size_t _lWindowSize, _rWindowSize;

public:
  
  maxAmpTime();
  
  int computeTmax(RPPulse&); 
  
  int computeT1(RPPulse&); 

  int computeT0(RPPulse&); 
  
  void setWindow(int left,int right);

  void setStartPoint(double percent) {_startPoint=percent;}

  void setEndPoint(double percent) {_endPoint=percent;}

  double getStartPoint() {return _startPoint;}

  double getEndPoint() {return _endPoint;}

};

// class splineTime : public timeDef {
  
// public:
 
//   splineTime(size_t vlevel=1);

//   int computeT0(RPPulse&); //min node

//   int computeT1(RPPulse&); //max node

//   int computeTstart(RPPulse&); //max
  
//   int computeTend(RPPulse&); // find zero
  
//   int computeTmax(RPPulse&); 

// };

class RecoPulseTime{

private:
  
  std::map<std::string, timeDef*> _tools; 
  
  std::string _timeDef;
 
  timeDef * _ptrtimeDef; 
  
  

public:
  
  RecoPulseTime();

  virtual ~RecoPulseTime();
  
  void reset();
  
  void config(std::string);
  
  
  timeDef* getTimeDef();

  int computeT0(RPPulse&);
  
  int computeT1(RPPulse&);

  int computeTstart(RPPulse&);
  
  int computeTend(RPPulse&);

  int computeTmax(RPPulse&);

  int computeTrise(RPPulse&);

  int computeTfall(RPPulse&);
  
  //bool computeTimes(pulse);

  bool find(std::string) const;
  
private:

  void checkTimeDef(std::string) const;

  
  
};



#endif
