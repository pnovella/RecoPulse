{
  
  gSystem->Load("libCore.so");
  gSystem->Load("libRIO.so");
  gSystem->Load("libHist.so");

  gSystem->Load("$(GATE_DIR)/lib/libGATE.so");
  gSystem->Load("$(GATE_DIR)/lib/libGATEIO.so");
  gSystem->Load("$(GATE_DIR)/lib/libGATEUtils.so");

  gSystem->Load("$(GATE_DIR)/lib/libRecoPulse.so");

  RecoPulse* rp = new RecoPulse(gate::NORMAL,"RP");
  
  std::string ifname = "../../../DATA/DST_Na22_3200.root";
  std::string ofname = "../../../DATA/DST_Na22_3200_RP.root";

  gate::Centella::instance(gate::NORMAL);
  gate::Centella::instance()->addInputFile(ifname);
  gate::Centella::instance()->addOutputFile(ofname);
  gate::Centella::instance()->setNevents(1000000);
  gate::Centella::instance()->saveEvents(true);
  gate::Centella::instance()->saveHistos(true);
  gate::Centella::instance()->addAlgo("RecoPulse",rp);
    
  gate::Centella::instance()->run();
  
  gate::Centella::instance()->destroy();
  delete rp; 

}
