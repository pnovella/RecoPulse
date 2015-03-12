
#include<RecoPulse.h>

ClassImp(RecoPulse)

//==========================================================================
RecoPulse::RecoPulse(gate::VLEVEL vl, std::string label) : 
IAlgo(vl,"RecoPulse",0,label){
//==========================================================================


}

//==========================================================================
RecoPulse::RecoPulse(const gate::ParamStore& gs, 
		     gate::VLEVEL vl, std::string label) :
  IAlgo(gs,vl,"RecoPulse",0,label){
//==========================================================================


}

//==========================================================================
bool RecoPulse::initialize(){
//==========================================================================

  _m.message("Intializing algorithm",this->getAlgoLabel(),gate::NORMAL);
  
  gate::Centella::instance()
    
    ->hman()->h1("Qch","Channel Q; Charge (ADC); Entries ",10,0,100);
  
  _recoMan =  new RecoManager();
  
  _recoMan->config("peakWindow");

  return true;

}

//==========================================================================
bool RecoPulse::execute(gate::Event& evt){
//==========================================================================

  _m.message("Executing algorithm",this->getAlgoLabel(),gate::VERBOSE);
  
  _m.message("Event number:",evt.GetEventID(),gate::VERBOSE);

  std::vector<gate::Hit*> pmts = evt.GetHits(gate::PMT);
  
  std::vector<gate::Hit*>::iterator ith;
  
  for (ith = pmts.begin(); ith != pmts.end(); ++ith){
    
    _m.message("PMT",(*ith)->GetSensorID(),gate::NORMAL);
    
    gate::Waveform& wf = (*ith)->GetWaveform();
    
    const std::vector<std::pair<double,double> >& pprof = wf.GetData(); 
    
    std::vector<short unsigned int> prof;
    
    std::vector<double> times;

    std::vector<std::pair<double,double> >::const_iterator it;
    
    for (it = pprof.begin(); it != pprof.end(); ++it){
      
      times.push_back(it->first); prof.push_back(it->second);}
    
    _recoMan->reset();
    
    _recoMan->reco(prof);
    
    for (size_t i=0; i< _recoMan->getNpeaks(); i++){
      
      gate::Pulse* pul = new gate::Pulse(); //check who is deleting!!!!
      
      pul->SetAmplitude(_recoMan->getPeakQs()[i]);

      //pul->SetStartTime(times[_recoMan->getPeakTs()[i]]);
      pul->SetStartTime(times[_recoMan->getPeakIntTs()[i]]);
      
      //pul->SetEndTime(times[_recoMan->getPeakTends()[i]]);
      pul->SetEndTime(times[_recoMan->getPeakIntTends()[i]]);
      
      pul->store("RecoPulse",1);

      wf.AddPulse(pul);
      
    }
   
    gate::Centella::instance()->hman()->fill("Qch",_recoMan->getQ());
  
  }

  
  return true;

}

//==========================================================================
bool RecoPulse::finalize(){
//==========================================================================

  _m.message("Finalising algorithm",this->getAlgoLabel(),gate::NORMAL);
  
  return true;

}
