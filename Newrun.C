#include "NewEventDisplay.C"

void Newrun(int num){
  gROOT->LoadMacro("NewEventDisplay.C");
  TChain *chain = new TChain("t_tap");
//  chain->Add("Jpsi_new3.root");
//  chain->Add("Zmumu_defaultMC16.root");
//  chain->Add("Zmumu_noRpcHit.root");
//  chain->Add("Zmumu_noRpcHitWide.root");
//  chain->Add("Zmumu_noRpcHitWideRoI.root");
//  chain->Add("Zmumu_noRpcHitWideRoI2.root");
//  chain->Add("Zmumu_noRpcHitWideRoI3.root");
//  chain->Add("Zmumu_noRpcHitWideRoI4.root");
//  chain->Add("Zmumu_noRpcHitWideRoISvc.root");
//  chain->Add("Zmumu_noRpcHitWideRoISvc12.root");
//  chain->Add("Zmumu_noRpcHitWideRoISvc12RoI.root");
//  chain->Add("Jpsi_default.root");
//  chain->Add("Jpsi_default2.root");
//  chain->Add("Jpsi_default_test.root");
//  chain->Add("Zmumu_MDTRegion.root");
//  chain->Add("Zmumu_outlier2.root");
//  chain->Add("Zmumu_new_Ib0386.root");
//  chain->Add("Zmumu_MyAOD2.root");
//  chain->Add("Jpsi_mc16_13TeV.root");
    chain->Add("/gpfs/fs6001/kayamash/dataset/Jpsi_mc16_13TeV1.root");
  //chain->Add("/gpfs/fs6001/kayamash/dataset/Zmumu300540_hadd.root");
  NewEventDisplay m(chain);
  m.Loop(num);
}

