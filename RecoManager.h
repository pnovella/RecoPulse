/**
 * @file DCRecoManager.hh
 * @brief Pulse Reconstruction Manager
 *
 * @class RecoManager
 *
 * @ingroup RecoPulse
 *
 * User-friendly user interface to set pulse reconstruction parameters 
 *
 * @author Pau Novella $Author: pnovella $
 *
 * @version $Revision: 1.7 $
 *
 * @date $Date: 2009-04-17 07:14:17 $
 *
 * Contact: pau.novella@ciemat.es
 *
 */

#ifndef DCRecoManager_hh
#define DCRecoManager_hh

#include <stdlib.h>
//#include <map>

//#include <DCRecoPulse-TypeDef.hh>

#include<RecoPulseAlgo.h>

class RecoManager {

 
private:

  std::map<std::string,Algo*> _algos; ///< algorithm storage

  Algo * _selectedalgo;
  
  std::string _algo; ///< selected reconstruction algorithm
  
private:

  /* cuts */

  //double _Qmin;
  //double _Imin;

public:
  
  /// constructor 
  RecoManager();
  
  /// default destructor
  virtual ~RecoManager();
  
  /// set algorithm and turn on/off channel pedestal analysis
  void config(std::string algo,bool pedEvt=true);
  
  /// retrieve selected algorithm
  Algo* getAlgo(){ return _selectedalgo;}
  
  /// find out if algorithm exists
  bool find(std::string) const;
  
  /// get algorithm name
  std::string getAlgoName() const {return _algo;}

  /// clear previous readout analysis
  void reset();
  

public:
  
  //-------- Algorithm alias functions ---------//
  

  /// perform full reconstruction
  void reco( const std::vector<unsigned short>&);
  
  //---- Charge ----//
  
  /// set charge measurament approach
  void setChargeAna(std::string);
  
  /// set cut on minimum Q
  void setQmin(double q);

  /// set cut on minimum maximum amplitude
  void setImin(double I);

  /// set cut on minimum pulse width 
  void setWmin(size_t w);
  
  /// set cut on maximum Amax/Width ratio 
  void setAWRmax(double r);

  /// compute pulse charge within window
  double computeQ( const std::vector<unsigned short>&);
  
  //---- time ----//

  /// compute start time for integration window
  int computeT0( const std::vector<unsigned short>&);
  
  /// compute end time for integration window
  int computeT1( const std::vector<unsigned short>&);
  
  /// compute arrival time
  int computeTstart( const std::vector<unsigned short>&);
  
  /// compute final time
  int computeTend( const std::vector<unsigned short>&);
  
  /// compute time for maximum amplitude
  int computeTmax( const std::vector<unsigned short>&);
  
  /// compute rise time 
  int computeTrise( const std::vector<unsigned short>&);

  /// compute fall time 
  int computeTfall( const std::vector<unsigned short>&);
  
  //---- pedestal ----//
  
   /// use ext ped Mean
  void useExtPedMean(bool ok);

  /// use ext ped RMS
  void useExtPedRMS(bool ok);

   /// use ext ped Mean
  bool useExtPedMean();

  /// use ext ped RMS
  bool useExtPedRMS();

  /// compute pedestal 
  double computePED( const std::vector<unsigned short>&);
  
  /// set window for pedestal estimation
  void setPEDwindow(size_t,size_t);
  
  /// set pedestal mean and RMS
  void setExtPedestal(double ped,double rms);  

  /// set minimum pedestal width
  void setMinPedRMS(double rms);
  
  /// turn on/off pedestal analysis
  void computeProfPedestal(bool ok);
  
  /// set threshold level 
  void setNsigOverPed(double nsigma);
  
  /// set maximum for good pedestal mean
  void setMaxPEDMean(double max);
  
  /// set minumum for good pedestal mean
  void setMinPEDMean(double min);
  
  /// set maximum for good pedestal RMS
  void setMaxPEDRMS(double max);
  
  /// set algorithm for pedestal estimation
  void setPedAlgo(std::string);
  
  /// enable calibration (bias correction)
  void setCalibOn(bool);

  /// set T calibration file
  void setTcalFile(std::string);

  /// set Q calibration file
  void setQcalFile(std::string);
  
  /// set filter (waveform smoothing)
  void setFilter(size_t);
  
  // reco varaibles
  
  /// get start time of integration window
  int getT0(); 
  
  /// get end time of integration window
  int getT1();
  
  /// get arrival time
  int getTstart(); 

  /// get arrival time error
  int getTstartError(); 

  /// get end time
  int getTend(); 

  /// get time with maximum amplitude
  int getTmax();
  
  /// get rise time
  int getTrise(); 
  
  /// get fall time
  int getTfall();
  
  /// get charge
  double getQ();

   /// get maximum amplitude
  double getImax();

  /// get charge error
  double getQerror();
  
  /// get number of peaks
  size_t getNpeaks();
  
  /// get start time of integration windows
  std::vector<int> getPeakIntTs();
  
  /// get ends time of integration windows
  std::vector<int> getPeakIntTends();
  
  /// get charge of each peak
  std::vector<double> getPeakQs();

  /// get start time of each peak
  std::vector<int> getPeakTs();

  /// get end time of each peak
  std::vector<int> getPeakTends();

  /// get max amplitude time of each peak
  std::vector<int> getPeakTmaxs();

  /// get max amplitude of each peak
  std::vector<double> getPeakImaxs();
  
  /// get used pedestal mean (ext or channel)
  double getPedMean();
  
  /// get used pedestal RMS (ext or channel)
  double getPedRMS();

  /// get channel pedestal mean
  double getChPedMean();
  
  /// get channel pedestal RMS
  double getChPedRMS();

  /// get ext pedestal mean
  double getExtPedMean();
  
  /// get ext pedestal RMS
  double getExtPedRMS();
  
  /// get pedestal integration window size
  size_t getPedWindowSize();

  /// get pedestal integration window lower edge
  size_t getPedWindowMin();

  /// get pedestal integration window upper edge
  size_t getPedWindowMax();
  
  /// get pedestal status
  int getPedStatus();

  /// true is pedestal mean is whithin bounds
  bool isGoodPedMean();

  /// true is pedestal mean is whithin bounds
  bool isGoodPedRMS();

   /// true is ext pedestal mean is whithin bounds
  bool isGoodExtPedMean();

  /// true is ext pedestal mean is whithin bounds
  bool isGoodExtPedRMS();
  
  /// true is ext. ped. RMS = ch. ped RMS
  bool sameChExtPedRMS();

  /// true if negative pulse
  bool isBiPolar();
  
  /// true is good pedestal
  bool isGoodPedestal();

  /// true is good ext. trig. pedestal
  bool isGoodExtPedestal();

  /// true if pulse fulfills requirements
  bool isGoodPulse();
  
  

  //-------------------------------------//

 
private:
  
  /// check if algorithm is in storage
  void checkAlgo(std::string) const;

};

#endif
