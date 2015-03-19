
#include<RecoPulse.h>

ClassImp(RecoPulse)

//==========================================================================
RecoPulse::RecoPulse(gate::VLEVEL vl, std::string label) : 
IAlgo(vl,"RecoPulse",0,label){
//==========================================================================
    
    _clearWF = false;

    _pulseLabel = "RecoPulse_v1";
}

//==========================================================================
RecoPulse::RecoPulse(const gate::ParamStore& gs, 
		     gate::VLEVEL vl, std::string label) :
  IAlgo(gs,vl,"RecoPulse",0,label){
//==========================================================================

    try{  _clearWF = gs.fetch_istore("CLEAR_WF");  }
    
    catch(exception& e) { _clearWF = false; }

    try{  _pulseLabel = gs.fetch_istore("PULSE_LABEL");  }

    catch(exception& e) { _pulseLabel = "RecoPulse_v1";}
    
}

//==========================================================================
bool RecoPulse::initialize(){
//==========================================================================

  _m.message("Intializing algorithm",this->getAlgoLabel(),gate::NORMAL);
  
  gate::Centella::instance()
    
    ->hman()->h1("Qch","Channel Q; Charge (ADC); Entries ",10,0,1000);
  
   gate::Centella::instance()
    
    ->hman()->h1("Qpeak","Pulse Q; Charge (ADC); Entries ",10,0,1000);
  
  _recoMan =  new RecoManager();
  
  _recoMan->config("peakWindow");

  _pmtSampWidth = 25;

  _sipmSampWidth = 1000;
  
  return true;

}

//==========================================================================
bool RecoPulse::execute(gate::Event& evt){
//==========================================================================

  _m.message("Executing algorithm",this->getAlgoLabel(),gate::VERBOSE);
  
  _m.message("Event number:",evt.GetEventID(),gate::VERBOSE);
  
  std::vector<gate::Hit*> pmts = evt.GetHits(gate::PMT);
  
  std::vector<gate::Hit*>::iterator ith;
  
  gate::Run* runInfo = &gate::Centella::instance()->getRun();

  for (ith = pmts.begin(); ith != pmts.end(); ++ith){
    
    _m.message("PMT",(*ith)->GetSensorID(),gate::VERBOSE);
    
    (*ith)->SetState(gate::RECOED);
    
    double gain = fabs(runInfo->GetSensor((*ith)->GetSensorID())->GetGain());

    gate::Waveform& wf = (*ith)->GetWaveform();
    
    const std::vector<std::pair<unsigned short,unsigned short> >& 
        
        pprof = wf.GetData(); 
    
    std::vector<short unsigned int> prof;
    
    std::vector<std::pair<unsigned short, unsigned short> >::const_iterator it;
    
    for (it = pprof.begin();it != pprof.end();++it){prof.push_back(it->second);}
    
    _recoMan->reset();
    
    _recoMan->reco(prof);
    
    for (size_t i=0; i< _recoMan->getNpeaks(); i++){
      
      gate::Pulse* pul = new gate::Pulse(); //check who is deleting!!!!
      
      pul->SetAmplitude(_recoMan->getPeakQs()[i]/gain);//TO CHECK: multiply per sample width?
      
      pul->SetStartTime(_recoMan->getPeakIntTs()[i]*_pmtSampWidth);
            
      pul->SetEndTime(_recoMan->getPeakIntTends()[i]*_pmtSampWidth );

      pul->SetMaxADC((int)_recoMan->getPeakImaxs()[i]);

      gate::Centella::instance()->hman()->fill("Qpeak",_recoMan->getPeakQs()[i]);
      
      pul->store("RecoPulse",_pulseLabel);

      wf.AddPulse(pul);

    }
    
    if (_clearWF) wf.ClearData();
    
    gate::Centella::instance()->hman()->fill("Qch",_recoMan->getQ());

  } //end of PMT loop
  
  return true;

}

//==========================================================================
bool RecoPulse::finalize(){
//==========================================================================

  _m.message("Finalising algorithm",this->getAlgoLabel(),gate::NORMAL);
  
  return true;

}
