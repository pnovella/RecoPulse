/**
 * @file DCRecoPulseAlgo.hh
 * @brief Pulse Reconstruction Algorithms
 * 
 * @class RecoPulseAlgo
 *
 * @ingroup RecoPulse
 *
 * Main Algorithm classes for pulse reconstrution. 
 * Current version provides three general reconstruction algorithms:
 *
 * 1) peakWindow: searchs for peaks in the full read out
 * 2) maxAlgo: serachs for amplitude relative maximums
 * 3) slidingWindow: searchs for time windows with maximum integral 
 *
 * Each one contains tools for charge, time and pedestal measurement.
 *
 * Class Algo defines the basic structure for a reconstruction algorithm 
 *
 * @author Pau Novella <pau.novella@ific.uv.es>
 *
 * @version $Revision: 1.9 $
 *
 * @date $Date: 2009-10-28 08:01:17 $
 *
 */

#ifndef RecoPulseAlgo_hh
#define RecoPulseAlgo_hh

#include<RecoPulseQ.h>
#include<RecoPulseTime.h>
#include<RecoPulsePed.h>
#include<RecoPulseUtils.h>

#include<TSystem.h> // Form function. To be fixed!!!

#include<GATE/Error.h>


class Algo {

protected:
  
  /// pulse start time
  int _Tstart;
  
  /// pulse start time error
  int _TstartError;

  /// pulse end time
  int _Tend;
  
  /// maximum amplitude time
  int _Tmax;

  /// pulse fall time
  int _Tfall;
  
  /// pulse rise time
  int _Trise;
  
  /// lower bound of time integration window
  int _T0;
  
  /// upper bound of time integration window
  int _T1;
  
  /// reconstructed charge in full readout
  double _Q;

  /// error on reconstructed charge
  double _QError;
  
  /// maximum pulse amplitude
  double _Imax;

  /// Q threshold for charge reconstruction
  double _Qmin; 

  /// Amplitude threshold for charge reconstruction
  double _Imin;
  
  /// pulse width threshold for charge reconstruction
  size_t _Wmin;

  /// Amax over Width ratio threshold for charge reconstruction
  double _AWRmax;

  /// pedestal mean
  double _ped;
  
  /// pedestal RMS
  double _pedRMS;
  
  /// minimum noise to compute Q_threshold
  double _pedRMS_min;
  
  /// pedestal status
  int _pedStatus;
  
  /// use external ped mean
  bool _useExtPedMean;

  /// use external ped RMS
  bool _useExtPedRMS;
 
  /// true if bipolar pulse
  bool _isBiPolar;

  /// reconstructed charge of each peak found
  std::vector<double> _Qs;
  
  /// error on reconstructed charge of each peak found
  std::vector<double> _QsError;
  
  /// lower bound of time window for each peak found
  std::vector<int> _T0s; 
  
  /// upper bound of time window for each peak found
  std::vector<int> _T1s; 
  
  /// start time of each peak found
  std::vector<int> _Ts; 

  /// end time of each peak found
  std::vector<int> _Tends; 
  
  /// Error on start time of each peak found
  std::vector<int> _TsError; 
  
  /// maximum amplitude time of each peak found
  std::vector<int>_Tmaxs;

  /// maximum amplitude of each peak found
  std::vector<double>_Imaxs;
  
  /// filter parameter
  size_t _nFilter;
  
private:
  

  /// time reconstruction algorithm
  RecoPulseTime* timer;
  
  /// Pedestal analysis algorithm
  RecoPulsePed* ped;
  
  /// charge reconstruction algorithm
  RecoPulseQ* charge;
  
public:
  
  Algo(bool pedEvt=true);

  virtual ~Algo();

  virtual double computeQ(RPPulse&) = 0;

  virtual int computeT0(RPPulse&) = 0;

  virtual int computeT1(RPPulse&) = 0;

  virtual int computeTstart(RPPulse&) = 0;

  virtual int computeTend(RPPulse&) = 0;

  virtual int computeTmax(RPPulse&) = 0;

  virtual int computeTrise(RPPulse&) = 0;

  virtual int computeTfall(RPPulse&) = 0;
  
  virtual double computePED(RPPulse&) = 0;
  
  /// compute maximum amplitude
  double computeImax(RPPulse&);

  /// set pedestal mean and RMS
  virtual void setPedestal(double,double)=0;

  /// set external pedestal mean and RMS
  virtual void setExtPedestal(double,double)=0;
  
  /// set filter (waveform smoothing)
  virtual void setFilter(size_t nfilter){_nFilter=nfilter;}
  
  /// compute calibrated Q and error
  //virtual void calibrateQ(double,double&,double&);

  /// compute calibrated T and error
  //virtual void calibrateT(double,int&,int&); 
 
  
  // extend window until pedestal level is reached
  void extendWindow(RPPulse& ipulse);

  virtual void reset();
  
  /// set time integration window
  virtual void setWindow(size_t t0,size_t t1){this->setT0(t0);this->setT1(t1);}
    
  /// set time window for pedestal analysis
  void setPEDwindow(size_t,size_t);
  
  /// set pedestal mean
  void setPedMean(double ped){_ped=ped;} 
  
  /// set pedestal RMS
  void setPedRMS(double rms) {_pedRMS=rms;}
  
  /// set minimum noise for Q_threshold
  void setMinPedRMS(double rms) {_pedRMS_min=rms;}
  
  /// set pedestal status
  void setPedStatus(int status) {_pedStatus=status;}
  
  /// turn on/off event pedestal analysis
  void computeProfPedestal(bool ok){_pedEvt=ok;}
  
  /// use ext ped Mean
  void useExtPedMean(bool ok){_useExtPedMean = ok;}

  /// use ext ped RMS
  void useExtPedRMS(bool ok){_useExtPedRMS = ok;}

  /// use ext ped Mean
  bool useExtPedMean() const {return _useExtPedMean;}

  /// use ext ped RMS
  bool useExtPedRMS() const {return _useExtPedRMS;}


  /// set lower bound of pulse window
  void setT0(size_t tmin){_T0=(int)tmin;}
  
  /// set upper bound of pulse window
  void setT1(size_t tmax){_T1=(int)tmax;}
  
  /// set threshold as a minimum charge 
  void setQmin(double q){_Qmin=q;}

  /// set threshold as a minimum amplitude 
  void setImin(double I){_Imin=I;}

  /// set minimum number of samples over ped 
  virtual void setWmin(size_t w){_Wmin=w;}
  
  /// set maximum Amax over Width ratio 
  void setAWRmax(double max){_AWRmax=max;}

  /// set integration window size, if needed depending on algorithm
  virtual void setWindowSize(size_t){
    //Message::MSG("RecoPulseAlgo",DC::kMWARNING,"NOT IMPLEMENTED!");
  }
    
  
  /// choose charge integration technique
  virtual void setChargeAna(std::string){
    //Message::MSG("RecoPulseAlgo",DC::kMWARNING,"CAN'T SET WINDOW!");
  }

  
  /// set threshold in number of pedestal sigmas
  void setNsigOverPed(double thr = 5); 
  
  /// set maximum allowed pedestal mean
  void setMaxPEDMean(double);

  /// set minimum allowed pedestal mean
  void setMinPEDMean(double);
  
  /// set maximum allowed pedestal RMS
  void setMaxPEDRMS(double);
  
  /// set minimum pedestal width
  //void setMinPedWidth(double);

  /// choose algorithm for pedestal analysis
  void setPedAlgo(std::string);

  /// choose algorithm for start time reconstruction
  void setTstartAlgo(std::string);

  /// turn on/off calibration
  //void setCalibOn(bool ok) {_calibOn = ok;}
  
  /// set T calibration file
  //void setTcalFile(std::string);

  /// set Q calibration file
  //void setQcalFile(std::string);
    
  double getQ() const {return _Q;}

  double getQerror() const {return _QError;}

  int getT0() const {return _T0;}

  int getT1() const {return _T1;}

  int getTstart() const {return _Tstart;}

  int getTstartError() const {return _TstartError;}

  int getTend() const {return _Tend;}

  int getTmax() const {return _Tmax;}
  
  int getTfall() const {return _Tfall;}

  int getTrise() const {return _Trise;}

  double getImax() const {return _Imax;}

  std::vector<int> getPeakT0s() const {return _T0s;}
  
  std::vector<int> getPeakT1s() const {return _T1s;}

  std::vector<double> getPeakQs() const {return _Qs;}
  
  std::vector<int> getPeakTs() const {return _Ts;}

  std::vector<int> getPeakTends() const {return _Tends;}

  std::vector<double> getPeakQsError() const {return _QsError;}
  
  std::vector<int> getPeakTsError() const {return _TsError;}

  std::vector<int> getPeakTmaxs() const {return _Tmaxs;}

  std::vector<double> getPeakImaxs() const {return _Imaxs;}

  double getPedMean() const {return _ped;}

  double getPedRMS() const {return _pedRMS;}

  double getExtPedMean() const {return ped->getExtMean();}

  double getExtPedRMS() const {return ped->getExtRMS();}

  double getChPedMean() const {return ped->getMean();}

  double getChPedRMS() const {return ped->getRMS();}

  int getPedStatus() const {return _pedStatus;}
  
  size_t getPedWindowMin() const {return ped->getWindowMin();}

  size_t getPedWindowMax() const {return ped->getWindowMax();}

  bool isGoodPedMean() const {return ped->isGoodPedMean();}

  bool isGoodPedRMS() const {return ped->isGoodPedRMS();}
  
  bool isGoodExtPedMean() const {return ped->isGoodExtPedMean();}

  bool isGoodExtPedRMS() const {return ped->isGoodExtPedRMS();}

  //bool sameChExtPeds() const {return ped->sameChExtPeds();}

  bool isBiPolar() const {return _isBiPolar;}

  bool isGoodPedestal() const {return ped->isGoodPedestal();}

  bool isGoodExtPedestal() const {return ped->isGoodExtPedestal();}
  
  size_t getPedWindowSize();
  
  bool sameChExtPedRMS(){return ped->sameChExtPedRMS();}

  /// compute and retrieve charge threshold
  double getQmin();
  
  /// apply user defined cuts
  void applyUserCuts(RPPulse& ipulse);

  /// apply full reconstruction
  virtual void reco(RPPulse&);
  
protected:
  
  bool _pedEvt;
  
  //bool _calibOn;

  //messenger _m;

protected:
  
  /// retrieve time algorithm
  RecoPulseTime* getTimer(){return timer;}
  
  /// retrieve pedestal algorithm
  RecoPulsePed* getPedAna(){return ped;}
  
  /// retrieve charge algorithm
  RecoPulseQ* getChargeAna(){return charge;}

  // retrieve calibration engine
  //RecoPulseCal* getCalibrator(){return calib;}
  
  /// calibrate and add error to Q and Tstart
  //void calibrate();

  /// retrieve saturation algorithm
  //RecoPulseSat* getSatAna(){return qsat;}
  
};

class simpleWindow : public Algo{

public: 
  
  simpleWindow(bool pedEvt=true);
  
  double computeQ(RPPulse&);

  int computeT0(RPPulse&);

  int computeT1(RPPulse&);

  int computeTstart(RPPulse&);

  int computeTend(RPPulse&);

  int computeTmax(RPPulse&);

  int computeTrise(RPPulse&);
  
  int computeTfall(RPPulse&);

  double computePED(RPPulse&);
  
  void setPedestal(double, double);

  void setExtPedestal(double, double);
  
  virtual void setWindow(size_t,size_t);

  virtual void reco(RPPulse&);
  
  virtual void reset();

};

class thresholdWindow : public simpleWindow{

public: 
  
  thresholdWindow(bool pedEvt=true);
  
  int computeT0(RPPulse&);
  
  int computeT1(RPPulse&);

  double computeQ(RPPulse&);
  
  void reco(RPPulse&);

  void reset();
 
private:

  void setWindow(size_t,size_t){
    //Message::MSG("RecoPulseAlgo",DC::kMWARNING,"CAN'T SET WINDOW!");
  }

};

class slidingWindow : public simpleWindow {
  
private: 


public:

  slidingWindow(bool pedEvt=true); 
  
  double computeQ(RPPulse&);

  int computeT0(RPPulse&);

  int computeT1(RPPulse&);

  int computeTmax(RPPulse&);

  int computeTstart(RPPulse&);

  int computeTend(RPPulse&);
  
  void setWindowSize(size_t);
  
  void checkPed(); 

  void reset();

  void reco(RPPulse&);

  void setWindow(size_t,size_t){
    //Message::MSG("RecoPulseAlgo",DC::kMWARNING,"CAN'T SET WINDOW!");
  }
  
  void setChargeAna(std::string);
  
  
private:

  double computePeakQ(RPPulse&);

};

class maxWindow : public simpleWindow{


public:
  
  maxWindow(bool pedEvt=true);
  
  double computeQ(RPPulse&);
  
  int computeT0(RPPulse&);

  int computeT1(RPPulse&);
  
  int computeTstart(RPPulse&);

  int computeTend(RPPulse&);

  void setWindow(size_t,size_t);

  void setChargeAna(std::string);

  void reco(RPPulse&);
  
  void reset();

private:
  
  double computePeakQ(RPPulse&);

};


class peakWindow : public simpleWindow{



protected:

  void checkPed();


public:
  
  peakWindow(bool pedEvt=true);
  
  double computeQ(RPPulse&);
    
  int computeT0(RPPulse&);

  int computeT1(RPPulse&);
  
  int computeTstart(RPPulse&);

  int computeTend(RPPulse&);

  int computeTmax(RPPulse&);
  
  void reco(RPPulse&);
  
  void reset();
  
  void setWmin(size_t w){_Wmin=w=0;}

};

class biPolarPeakWindow : public peakWindow{

public:
  
  biPolarPeakWindow(bool pedEvt=true);
  
  double computeQ(RPPulse&);

};

#endif


