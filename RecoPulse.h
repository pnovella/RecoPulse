#ifndef _RecoPulse__
#define _RecoPulse__

#include <TSystem.h>

#include <Centella.h>
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
  
 private:
  
  //! pulse recosntruction manager
  RecoManager* _recoMan;

  ClassDef(RecoPulse,0)
    
};

#endif
