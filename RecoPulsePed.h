/**
 * @file DCRecoPulsePed.hh
 * Pedestal analysis algorithm
 * 
 * @class RecoPulsePed
 *
 * @ingroup RecoPulse
 *
 * @brief Pedestal analysis algorithm
 *
 * Algorithm for channel pedestal measurement (mean and RMS), 
 * quality evaluation, and definition of signal threshold.  
 * Current version provides two different ways of pedestal estimations: 
 *
 * 1) Simple mean and RMS within a given window (no signal expected)
 * 2) full readout analysis (signal pulses might be there)
 *
 * @author Pau Novella, Carmen Palomares $Author: pnovella $
 *
 * @version $Revision: 1.7 $
 *
 * @date $Date: 2009-10-28 08:01:17 $
 *
 * Contact: pau.novella@ciemat.es
 *
 */

#ifndef DCRecoPulsePed_hh
#define DCRecoPulsePed_hh

#include<vector>
#include<Messenger.h>

#include<RecoPulseUtils.h>

class RecoPulsePed {

  
private:
  
  /// pedestal mean
  double _mean;
  
  /// pedestal RMS
  double _RMS;
  
  /// external (computed elsewhere) pedestal mean
  double _extMean;
  
  /// external (computed elsewhere) pedestal RMS
  double _extRMS;

  /// maximum allowed mean
  double _maxMean; 

  /// minimum allowed mean
  double _minMean; 
  
  /// maximum allowed RMS
  double _maxRMS;

  /// lower bound of window time
  size_t _tmin;
  
  /// upper bound of window time
  size_t _tmax;
  
  /// time window size
  size_t _wsize;
 
  /// sigmas over pedestal to define threshold
  double _nsigma; 
   
  ///! current threshold
  //double _threshold;

  /// name of pedestal algorithm applied
  bool _pedAlgoisfull;
  
  /// pedestal tolerance (%) to compare with ext trigger
  double _pedTol;

public:

  RecoPulsePed();
  
  virtual ~RecoPulsePed(){};

  /// configure pedestal analysis
  void config(double,double,double,size_t,size_t);

  void reset(){_mean=NOPED;_RMS=NOPED;_extMean=NOPED;_extRMS=NOPED;}

  /// compute pedestal according to specified technique
  bool computePedestal(RPPulse&,size_t tmin=0,size_t tmax=0);
  
  ///compute pedestal within an offset of the readout
  bool computeOffsetPedestal(RPPulse&,size_t tmin=0,size_t tmax=0);
  
  ///compute pedestal using full readout
  bool computeFullPedestal(RPPulse&);

  /// set maximum mean for a pedestal with status OK
  void setMaxMean(double max){_maxMean=max;}

  /// set minimum mean for a pedestal with status OK
  void setMinMean(double min){_minMean=min;}
  
  /// set maximum RMS for a pedestal with status OK
  void setMaxRMS(double max){_maxRMS=max;}

  /// set pedestal mean
  void setMean(double mean){_mean=mean;}
  
  /// set pedestal RMS
  void setRMS(double rms){_RMS=rms;}

  /// set external pedestal mean
  void setExtMean(double mean){_extMean=mean;}
  
  /// set external pedestal RMS
  void setExtRMS(double rms){_extRMS=rms;}
  
  /// set minim pedestal width. Will replace measured RMS if >RMS
  //void setMinPedWidth(double rms){_RMS_min = rms;}

  /// set window for pedetal analysis
  void setWindow(size_t tmin, size_t tmax);
  
  /// choose recosntruction technique
  void setPedAlgo(std::string algo){
    if(algo != "full" && algo != "offset") 
      //Message::MSG("RecoPulsePed",DC::kMFATAL,"Algo not defined: " + algo);
    _pedAlgoisfull = (algo == "full");
  }
  
  /// set number of sigmas over pedestal, thus defining threshold level
  void setNsigOverPed(double nsigma){_nsigma=nsigma;}
  
  /// retrieve threshold level
  //double getThreshold(){return _threshold;}
  
  /// compute threshod level
  double computeThreshold(double,double);
  
  /// true if amplitude is over threshold
  //bool isOverThreshold(double amp);
  
  /// true if pedestal is defined
  bool status() const {return (_mean!=NOPED);}
  
  /// retrieve pedestal mean
  double getMean() const {return _mean;}

  /// retrieve pedestal RMS
  double getRMS() const {return _RMS;}
  
  /// retrieve ext pedestal mean
  double getExtMean() const {return _extMean;}

  /// retrieve ext pedestal RMS
  double getExtRMS() const {return _extRMS;}
  

  /// retrieve number of sigmas defining the threshold
  double getNsigOverPed() const {return _nsigma;}
  
  /// get pedestal analysis window size
  size_t getWindowSize() const {return _wsize;}

  /// get pedestal analysis window: lower edge
  size_t getWindowMin() const {return _tmin;}

  /// get pedestal analysis window: upper edge
  size_t getWindowMax() const {return _tmax;}
  
  //! true if ped mean is withing bounds
  bool isGoodPedMean() const;

  //! true if ped RMS is withing bounds
  bool isGoodPedRMS() const;
    
   /// true if pedestal fulfills requirements
  bool isGoodPedestal() const;


  /// true if ext pedestal fulfills requirements
  bool isGoodExtPedestal() const;

  //! true if ext ped mean is withing bounds
  bool isGoodExtPedMean() const;

  //! true if ext ped RMS is withing bounds
  bool isGoodExtPedRMS() const;
  
  //! true if same channel and ext. trig. pedestals
  //bool sameChExtPeds() const;

  //! true if same channel and ext. trig. ped means
  bool sameChExtPedMean() const;

  //! true if same channel and ext. trig. ped RMS
  bool sameChExtPedRMS() const;

  /// retrieve pedestal status
  int getPedStatus() const;

  /// retrieve pedestal status
  int getExtPedStatus() const;
  
private:

  /// ensure analysis window is well defined
  void checkWindow(size_t,size_t) const;

  /// chech if pedestal is defined
  void checkStatus() const;

};
#endif
