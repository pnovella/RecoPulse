{
  
  gSystem->Load("libCore.so");
  gSystem->Load("libRIO.so");
  gSystem->Load("libHist.so");

  gSystem->Load("$(GATE_DIR)/lib/libGATE.so");
  gSystem->Load("$(GATE_DIR)/lib/libGATEIO.so");
  gSystem->Load("$(GATE_DIR)/lib/libGATEUtils.so");

  gSystem->Load("$(GATE_DIR)/lib/libRecoPulse.so");

  RecoPulse* rp = new RecoPulse(gate::NORMAL,"RP");
  rp->SetPulseLabel("v0");
  rp->ClearWF(true);

  std::string ifname0 = "/data4/NEXT/users/pnovella/DATA/DST_3350_000.root";
  std::string ifname1 = "/data4/NEXT/users/pnovella/DATA/DST_3350_001.root";
  std::string ifname2 = "/data4/NEXT/users/pnovella/DATA/DST_3350_002.root";
  std::string ifname3 = "/data4/NEXT/users/pnovella/DATA/DST_3350_003.root";
  std::string ifname4 = "/data4/NEXT/users/pnovella/DATA/DST_3350_004.root";
  std::string ifname5 = "/data4/NEXT/users/pnovella/DATA/DST_3350_005.root";
  std::string ifname6 = "/data4/NEXT/users/pnovella/DATA/DST_3350_006.root";
  std::string ifname7 = "/data4/NEXT/users/pnovella/DATA/DST_3350_007.root";
  std::string ifname8 = "/data4/NEXT/users/pnovella/DATA/DST_3350_008.root";
  std::string ifname9 = "/data4/NEXT/users/pnovella/DATA/DST_3350_009.root";
  
  std::string ofname0 = "/data4/NEXT/users/pnovella/DATA/DST_3350_000_RP.root";
  std::string ofname1 = "/data4/NEXT/users/pnovella/DATA/DST_3350_001_RP.root";
  std::string ofname2 = "/data4/NEXT/users/pnovella/DATA/DST_3350_002_RP.root";
  std::string ofname3 = "/data4/NEXT/users/pnovella/DATA/DST_3350_003_RP.root";
  std::string ofname4 = "/data4/NEXT/users/pnovella/DATA/DST_3350_004_RP.root";
  std::string ofname5 = "/data4/NEXT/users/pnovella/DATA/DST_3350_005_RP.root";
  std::string ofname6 = "/data4/NEXT/users/pnovella/DATA/DST_3350_006_RP.root";
  std::string ofname7 = "/data4/NEXT/users/pnovella/DATA/DST_3350_007_RP.root";
  std::string ofname8 = "/data4/NEXT/users/pnovella/DATA/DST_3350_008_RP.root";
  std::string ofname9 = "/data4/NEXT/users/pnovella/DATA/DST_3350_009_RP.root";

  gate::Centella::instance(gate::NORMAL);

  gate::Centella::instance()->addInputFile(ifname0);
  gate::Centella::instance()->addInputFile(ifname1);
  gate::Centella::instance()->addInputFile(ifname2);
  gate::Centella::instance()->addInputFile(ifname3);
  gate::Centella::instance()->addInputFile(ifname4);
  gate::Centella::instance()->addInputFile(ifname5);
  gate::Centella::instance()->addInputFile(ifname6);
  gate::Centella::instance()->addInputFile(ifname7);
  gate::Centella::instance()->addInputFile(ifname8);
  gate::Centella::instance()->addInputFile(ifname9);

  gate::Centella::instance()->addOutputFile(ofname0);
  gate::Centella::instance()->addOutputFile(ofname1);
  gate::Centella::instance()->addOutputFile(ofname2);
  gate::Centella::instance()->addOutputFile(ofname3);
  gate::Centella::instance()->addOutputFile(ofname4);
  gate::Centella::instance()->addOutputFile(ofname5);
  gate::Centella::instance()->addOutputFile(ofname6);
  gate::Centella::instance()->addOutputFile(ofname7);
  gate::Centella::instance()->addOutputFile(ofname8);
  gate::Centella::instance()->addOutputFile(ofname9);

  gate::Centella::instance()->setNevents(1000000000);
  gate::Centella::instance()->saveEvents(true);
  gate::Centella::instance()->saveHistos(true);
  gate::Centella::instance()->addAlgo("RecoPulse",rp);
    
  gate::Centella::instance()->run();
  
  gate::Centella::instance()->destroy();
  delete rp; 

}
