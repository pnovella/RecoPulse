#ifndef _RecoPulse__
#define _RecoPulse__

#include <GATE/Centella.h>
#include <RecoManager.h>

class RecoPulse : public gate::IAlgo {

 public:
  
  //! default contructor
  RecoPulse(gate::VLEVEL=gate::NORMAL,
	       std::string label="RecoPulseInstance");
  
  //! constructor with store with input parameters 
  RecoPulse(const gate::ParamStore& gs,
	       gate::VLEVEL=gate::NORMAL,
	       std::string label="RecoPulseInstance");
  
  //! destructor
  virtual ~RecoPulse(){};
  
  //! initialize algorithm
  bool initialize();        
  
  //! execute algorithm: process current event
  bool execute(gate::Event& evt);  
  
  //! finalize algorithm
  bool finalize();          
  
  //! set pulse label
  inline void SetPulseLabel(std::string lb) {_pulseLabel=lb;}

   //! set pulse label
  inline void ClearWF(bool ok) {_clearWF=ok;}

 private:
  
  //! set up RecoManager parameters
  void confRecoMan();
  
  //! initialize parameters
  void initParams();

 private:
  
  //! PMT pulse recosntruction manager
  RecoManager* _pmtRecoMan;

  //! SiPM pulse recosntruction manager
  RecoManager* _siRecoMan;
  
  //! PMT sampling width
  double _pmtSampWidth;
  
  //! SiPM sampling width
  double _sipmSampWidth;
    
 private:
  
  //! pulse label
  std::string _pulseLabel;
  
  //! clear waveform?
  bool _clearWF;
  
  //! minimum pulse amplitude
  double _Imin;

  //! number of sigmas above noise
  double _nSig;
  
  //! reconstruct SiPM waveforms?
  bool _noSiPM;
  
  //! reocntruct only this PMT channel
  int _PMTCh;

  ClassDef(RecoPulse,0)
    
};

#endif
