{
  
  gSystem->Load("libCore.so");
  gSystem->Load("libRIO.so");
  gSystem->Load("libHist.so");

  gSystem->Load("$(GATE_DIR)/lib/libGATE.so");
  gSystem->Load("$(GATE_DIR)/lib/libGATEIO.so");
  gSystem->Load("$(GATE_DIR)/lib/libGATEUtils.so");

  gSystem->Load("$(GATE_DIR)/lib/libRecoPulse.so");

  RecoPulse* rp = new RecoPulse(gate::NORMAL,"RP");

  //Message::SetLevelMSG("RecoPulseAlgo",gate::DUMP);

  gate::Centella::instance(gate::NORMAL);
  gate::Centella::instance()->addInputFile("../../../DATA/DST_Na22_3200.root");
  gate::Centella::instance()->addOutputFile("output_dst.root");
  gate::Centella::instance()->setNevents(1);
  gate::Centella::instance()->saveEvents(true);
  gate::Centella::instance()->saveHistos(true);
  gate::Centella::instance()->addAlgo("RecoPulse",rp);
    
  gate::Centella::instance()->run();
  
  gate::Centella::instance()->destroy();
  delete rp; 

}
