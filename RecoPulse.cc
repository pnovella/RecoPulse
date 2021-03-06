
#include<RecoPulse.h>

ClassImp(RecoPulse)

//==========================================================================
RecoPulse::RecoPulse(gate::VLEVEL vl, std::string label) : 
IAlgo(vl,"RecoPulse",0,label){
//==========================================================================
        
    this->initParams();
}

//==========================================================================
RecoPulse::RecoPulse(const gate::ParamStore& gs, 
		     gate::VLEVEL vl, std::string label) :
  IAlgo(gs,vl,"RecoPulse",0,label){
//==========================================================================
  
    this->initParams();

    try{  _clearWF = gs.fetch_istore("CLEAR_WF");  }
    
    catch(exception& e) { }

    try{  _pulseLabel = gs.fetch_sstore("PULSE_LABEL");  }

    catch(exception& e) { }
    
    try{  _Imin = gs.fetch_dstore("MIN_AMPLITUDE");  }

    catch(exception& e) { }

    try{  _nSig = gs.fetch_dstore("N_SIGMA");  }

    catch(exception& e) { }
    
    try{  _PMTCh = gs.fetch_istore("PMT_ID");  }

    catch(exception& e) { }
    
    try{  _noSiPM = bool(gs.fetch_istore("NO_SIPM"));  }

    catch(exception& e) { }

}

//==========================================================================
void RecoPulse::initParams(){
//==========================================================================
  
  _clearWF = false;
   
  _pulseLabel = "RecoPulse_v1";
  
  _Imin = 5;

  _nSig = 5;
  
  _noSiPM = true;
  
  _PMTCh = -1;

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
  
  _pmtRecoMan->setNsigOverPed(_nSig);
  
  _pmtRecoMan->setImin(_Imin);

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
      
    int chID =  (*ith)->GetSensorID();
      
    _m.message("PMT",chID,gate::VERBOSE);
    
    if (_PMTCh >= 0 && _PMTCh!=chID){
        
        continue;  }

    (*ith)->SetState(gate::RECOED);
    
    //if ((*ith)->GetSensorID()<20) continue;
    
    double gain = fabs(runInfo->GetGeometry()->
                       GetSensor((*ith)->GetSensorID())->GetGain());
    
    //double gain = 1;

    gate::Waveform& wf = (*ith)->GetWaveform();
    
    const std::vector<std::pair<unsigned int,float> >& 
        
        pprof = wf.GetData(); 
    
    std::vector<unsigned int> prof;
    
    std::vector<std::pair<unsigned int, float> >::const_iterator it;
    
    // TO FIX: this won't work with zero-suppression electronics!!!!!!

    for (it = pprof.begin();it != pprof.end();++it){prof.push_back(it->second);}
    
    _pmtRecoMan->reset();
    
    _pmtRecoMan->reco(prof);
    
    wf.SetBaseline(_pmtRecoMan->getChPedMean());

    wf.SetBaselineSig(_pmtRecoMan->getChPedRMS());
    
    for (size_t i=0; i< _pmtRecoMan->getNpeaks(); i++){
      
      gate::Pulse* pul = new gate::Pulse(); //check who is deleting!!!!
      
      pul->SetSensorID((*ith)->GetSensorID());

      pul->SetCharge(_pmtRecoMan->getPeakQs()[i]/gain);
      
      pul->SetStartTime(_pmtRecoMan->getPeakIntTs()[i]*_pmtSampWidth);
            
      pul->SetEndTime(_pmtRecoMan->getPeakIntTends()[i]*_pmtSampWidth );

      pul->SetMaxTime(_pmtRecoMan->getPeakTmaxs()[i]*_pmtSampWidth);

      pul->SetMaxADC((int)_pmtRecoMan->getPeakImaxs()[i]);

      gate::Centella::instance()->hman()
	
	->fill("Qpeak",_pmtRecoMan->getPeakQs()[i]);
      
      pul->store("RecoPulse",_pulseLabel);

      wf.AddPulse(pul);

    }
    
    if (_clearWF) wf.ClearData();
    
    gate::Centella::instance()->hman()->fill("Qch",_pmtRecoMan->getQ());

  } //end of PMT loop
  
  if (_noSiPM) return true; // to be improved!!!

  std::vector<gate::Hit*> sipms = evt.GetHits(gate::SIPM);

  std::vector<gate::Hit*>::iterator ith2;
  
  for (ith2 = sipms.begin(); ith2 != sipms.end(); ++ith2){
    
    /* just get baseline mean and RMS, do not reconstruct pulses */
    
    _m.message("SiPM",(*ith2)->GetSensorID(),gate::VERBOSE);
    
    (*ith2)->SetState(gate::RAW);
    
    gate::Waveform& wf = (*ith2)->GetWaveform();
    
    const std::vector<std::pair<unsigned int,float> >& 
        
        pprof = wf.GetData(); 
    
    std::vector<unsigned int> prof;
    
    // TO FIX: this won't work with zero-suppression electronics!!!!!!

    std::vector<std::pair<unsigned int, float> >::const_iterator it;
    
    for (it = pprof.begin();it != pprof.end();++it){prof.push_back(it->second);}
    
    _siRecoMan->reset();
    
    _siRecoMan->computePED(prof);
    
    wf.SetBaseline(_siRecoMan->getChPedMean());

    wf.SetBaselineSig(_siRecoMan->getChPedRMS());

  }

  return true;

}

//==========================================================================
bool RecoPulse::finalize(){
//==========================================================================

  _m.message("Finalising algorithm",this->getAlgoLabel(),gate::NORMAL);
  
  return true;

}
