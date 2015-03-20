
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
   
   this->confRecoMan();

  _pmtSampWidth = 25;

  _sipmSampWidth = 1000;
  
  return true;

}
//==========================================================================
void RecoPulse::confRecoMan(){
//==========================================================================
  
  _pmtRecoMan =  new RecoManager();
  
  _pmtRecoMan->config("peakWindow");
  
  _pmtRecoMan->setPEDwindow(0,100);
  
  _pmtRecoMan->setNsigOverPed(5);
  
  _pmtRecoMan->setImin(5);

  //---
  
  _siRecoMan =  new RecoManager();

  _siRecoMan->config("peakWindow");
  
  _siRecoMan->setPEDwindow(0,50);

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
    
    // TO FIX: this won't work with zero-suppression electronics!!!!!!

    for (it = pprof.begin();it != pprof.end();++it){prof.push_back(it->second);}
    
    _pmtRecoMan->reset();
    
    _pmtRecoMan->reco(prof);
    
    wf.SetBaseline(_pmtRecoMan->getChPedMean());

    wf.SetBaselineRMS(_pmtRecoMan->getChPedRMS());
    
    for (size_t i=0; i< _pmtRecoMan->getNpeaks(); i++){
      
      gate::Pulse* pul = new gate::Pulse(); //check who is deleting!!!!
      
      pul->SetSensorID((*ith)->GetSensorID());

      pul->SetAmplitude(_pmtRecoMan->getPeakQs()[i]/gain);
      
      pul->SetStartTime(_pmtRecoMan->getPeakIntTs()[i]*_pmtSampWidth);
            
      pul->SetEndTime(_pmtRecoMan->getPeakIntTends()[i]*_pmtSampWidth );

      pul->SetMaxADC((int)_pmtRecoMan->getPeakImaxs()[i]);

      gate::Centella::instance()->hman()
	
	->fill("Qpeak",_pmtRecoMan->getPeakQs()[i]);
      
      pul->store("RecoPulse",_pulseLabel);

      wf.AddPulse(pul);

    }
    
    if (_clearWF) wf.ClearData();
    
    gate::Centella::instance()->hman()->fill("Qch",_pmtRecoMan->getQ());

  } //end of PMT loop

  std::vector<gate::Hit*> sipms = evt.GetHits(gate::SIPM);

  std::vector<gate::Hit*>::iterator ith2;
  
  for (ith2 = sipms.begin(); ith2 != sipms.end(); ++ith2){
    
    /* just get baseline mean and RMS, do not reconstruct pulses */
    
    _m.message("SiPM",(*ith2)->GetSensorID(),gate::VERBOSE);
    
    (*ith2)->SetState(gate::RAW);
    
    gate::Waveform& wf = (*ith2)->GetWaveform();
    
    const std::vector<std::pair<unsigned short,unsigned short> >& 
        
        pprof = wf.GetData(); 
    
    std::vector<short unsigned int> prof;
    
    // TO FIX: this won't work with zero-suppression electronics!!!!!!

    std::vector<std::pair<unsigned short, unsigned short> >::const_iterator it;
    
    for (it = pprof.begin();it != pprof.end();++it){prof.push_back(it->second);}
    
    _siRecoMan->reset();
    
    _siRecoMan->computePED(prof);
    
    wf.SetBaseline(_siRecoMan->getChPedMean());

    wf.SetBaselineRMS(_siRecoMan->getChPedRMS());
    
  }

  return true;

}

//==========================================================================
bool RecoPulse::finalize(){
//==========================================================================

  _m.message("Finalising algorithm",this->getAlgoLabel(),gate::NORMAL);
  
  return true;

}
