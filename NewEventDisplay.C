#define NewEventDisplay_cxx
#include "NewEventDisplay.h"
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TBox.h>
#include <TArc.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMath.h>

#define NMEAMX 10
#define NCAND 6

void NewEventDisplay::Loop(int num)
{
//   In a ROOT session, you can do:
//      root> .L NewEventDisplay.C
//      root> NewEventDisplay t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   TGraph *tsegZR = new TGraph();
   tsegZR->SetMarkerStyle(22);
   tsegZR->SetMarkerColor(1);
   tsegZR->SetMarkerSize(0.7);
   tsegZR->SetName("tsegZR");

   TGraph *psegZR = new TGraph();
   psegZR->SetMarkerStyle(8);
   psegZR->SetMarkerColor(1);
   psegZR->SetMarkerSize(0.7);
   psegZR->SetName("psegZR");

   TGraph *mdtRegionLD_BI = new TGraph();
   mdtRegionLD_BI->SetMarkerStyle(2);
   mdtRegionLD_BI->SetMarkerColor(7);
   mdtRegionLD_BI->SetMarkerSize(0.5);
   mdtRegionLD_BI->SetName("mdtRegionLD_BI");
     
   TGraph *mdtRegionRD_BI = new TGraph();
   mdtRegionRD_BI->SetMarkerStyle(2);
   mdtRegionRD_BI->SetMarkerColor(7);
   mdtRegionRD_BI->SetMarkerSize(0.5);
   mdtRegionRD_BI->SetName("mdtRegionRD_BI");
     
   TGraph *mdtRegionLU_BI = new TGraph();
   mdtRegionLU_BI->SetMarkerStyle(2);
   mdtRegionLU_BI->SetMarkerColor(7);
   mdtRegionLU_BI->SetMarkerSize(0.5);
   mdtRegionLU_BI->SetName("mdtRegionLU_BI");
     
   TGraph *mdtRegionRU_BI = new TGraph();
   mdtRegionRU_BI->SetMarkerStyle(2);
   mdtRegionRU_BI->SetMarkerColor(7);
   mdtRegionRU_BI->SetMarkerSize(0.5);
   mdtRegionRU_BI->SetName("mdtRegionRU_BI");
     
   TGraph *mdtRegionLD_BM = new TGraph();
   mdtRegionLD_BM->SetMarkerStyle(2);
   mdtRegionLD_BM->SetMarkerColor(7);
   mdtRegionLD_BM->SetMarkerSize(0.5);
   mdtRegionLD_BM->SetName("mdtRegionLD_BM");
     
   TGraph *mdtRegionRD_BM = new TGraph();
   mdtRegionRD_BM->SetMarkerStyle(2);
   mdtRegionRD_BM->SetMarkerColor(7);
   mdtRegionRD_BM->SetMarkerSize(0.5);
   mdtRegionRD_BM->SetName("mdtRegionRD_BM");
     
   TGraph *mdtRegionLU_BM = new TGraph();
   mdtRegionLU_BM->SetMarkerStyle(2);
   mdtRegionLU_BM->SetMarkerColor(7);
   mdtRegionLU_BM->SetMarkerSize(0.5);
   mdtRegionLU_BM->SetName("mdtRegionLU_BM");
     
   TGraph *mdtRegionRU_BM = new TGraph();
   mdtRegionRU_BM->SetMarkerStyle(2);
   mdtRegionRU_BM->SetMarkerColor(7);
   mdtRegionRU_BM->SetMarkerSize(0.5);
   mdtRegionRU_BM->SetName("mdtRegionRU_BM");
     
   TGraph *mdtRegionLD_BO = new TGraph();
   mdtRegionLD_BO->SetMarkerStyle(2);
   mdtRegionLD_BO->SetMarkerColor(7);
   mdtRegionLD_BO->SetMarkerSize(0.5);
   mdtRegionLD_BO->SetName("mdtRegionLD_BO");
     
   TGraph *mdtRegionRD_BO = new TGraph();
   mdtRegionRD_BO->SetMarkerStyle(2);
   mdtRegionRD_BO->SetMarkerColor(7);
   mdtRegionRD_BO->SetMarkerSize(0.5);
   mdtRegionRD_BO->SetName("mdtRegionRD_BO");
     
   TGraph *mdtRegionLU_BO = new TGraph();
   mdtRegionLU_BO->SetMarkerStyle(2);
   mdtRegionLU_BO->SetMarkerColor(7);
   mdtRegionLU_BO->SetMarkerSize(0.5);
   mdtRegionLU_BO->SetName("mdtRegionLU_BO");
     
   TGraph *mdtRegionRU_BO = new TGraph();
   mdtRegionRU_BO->SetMarkerStyle(2);
   mdtRegionRU_BO->SetMarkerColor(7);
   mdtRegionRU_BO->SetMarkerSize(0.5);
   mdtRegionRU_BO->SetName("mdtRegionRU_BO");
     
   TGraph *rpcZR = new TGraph();
   rpcZR->SetMarkerStyle(24);
   rpcZR->SetMarkerColor(kBlue);
   rpcZR->SetMarkerSize(0.5);
   rpcZR->SetName("rpcZR");

   TGraph *mdtZR = new TGraph();
   mdtZR->SetMarkerStyle(24);
   mdtZR->SetMarkerColor(kRed);
   mdtZR->SetMarkerSize(0.5);
   mdtZR->SetName("mdtZR");

   TGraph *mdtZRIs1 = new TGraph();
   mdtZRIs1->SetMarkerStyle(30);
   mdtZRIs1->SetMarkerColor(kRed);
   mdtZRIs1->SetMarkerSize(1);
   mdtZRIs1->SetName("mdtZRIs1");

   TGraph *mdtZRIs2 = new TGraph();
   mdtZRIs2->SetMarkerStyle(32);
   mdtZRIs2->SetMarkerColor(kRed);
   mdtZRIs2->SetMarkerSize(1);
   mdtZRIs2->SetName("mdtZRIs2");

   TGraph *mdtZRIs3 = new TGraph();
   mdtZRIs3->SetMarkerStyle(5);
   mdtZRIs3->SetMarkerColor(kRed);
   mdtZRIs3->SetMarkerSize(1.0);
   mdtZRIs3->SetName("mdtZRIs3");

   TArc *mdtarc = new TArc();
   mdtarc->SetFillStyle(0);
   mdtarc->SetLineColor(kRed);
   mdtarc->SetLineWidth(1);

   TGraph *SPBIZR = new TGraph();
   SPBIZR->SetMarkerStyle(8);
   SPBIZR->SetMarkerColor(kPink);
   SPBIZR->SetMarkerSize(0.7);
   SPBIZR->SetName("SPBIZR");

   TGraph *SPBMZR = new TGraph();
   SPBMZR->SetMarkerStyle(8);
   SPBMZR->SetMarkerColor(kPink);
   SPBMZR->SetMarkerSize(0.7);
   SPBMZR->SetName("SPBMZR");

   TGraph *SPBOZR = new TGraph();
   SPBOZR->SetMarkerStyle(8);
   SPBOZR->SetMarkerColor(kPink);
   SPBOZR->SetMarkerSize(0.7);
   SPBOZR->SetName("SPBOZR");

   TGraph *SPEIZR = new TGraph();
   SPEIZR->SetMarkerStyle(8);
   SPEIZR->SetMarkerColor(kPink);
   SPEIZR->SetMarkerSize(0.7);
   SPEIZR->SetName("SPEIZR");

   TGraph *SPEMZR = new TGraph();
   SPEMZR->SetMarkerStyle(8);
   SPEMZR->SetMarkerColor(kPink);
   SPEMZR->SetMarkerSize(0.7);
   SPEMZR->SetName("SPEMZR");

   TGraph *SPEOZR = new TGraph();
   SPEOZR->SetMarkerStyle(8);
   SPEOZR->SetMarkerColor(kPink);
   SPEOZR->SetMarkerSize(0.7);
   SPEOZR->SetName("SPEOZR");

   TGraph *SPEEZR = new TGraph();
   SPEEZR->SetMarkerStyle(8);
   SPEEZR->SetMarkerColor(kPink);
   SPEEZR->SetMarkerSize(0.7);
   SPEEZR->SetName("SPEEZR");

   TGraph *SPCSCZR = new TGraph();
   SPCSCZR->SetMarkerStyle(8);
   SPCSCZR->SetMarkerColor(kPink);
   SPCSCZR->SetMarkerSize(0.7);
   SPCSCZR->SetName("SPCSCZR");

   TGraph *SPBEEZR = new TGraph();
   SPBEEZR->SetMarkerStyle(8);
   SPBEEZR->SetMarkerColor(kPink);
   SPBEEZR->SetMarkerSize(0.7);
   SPBEEZR->SetName("SPBEEZR");

   TGraph *SPBMEZR = new TGraph();
   SPBMEZR->SetMarkerStyle(8);
   SPBMEZR->SetMarkerColor(kPink);
   SPBMEZR->SetMarkerSize(0.7);
   SPBMEZR->SetName("SPBMEZR");

   TArc *SParc1 = new TArc();
   SParc1->SetFillStyle(0);
   SParc1->SetLineColor(2);
   SParc1->SetLineWidth(1);

   TArc *SParc2 = new TArc();
   SParc2->SetFillStyle(0);
   SParc2->SetLineColor(4);
   SParc2->SetLineWidth(1);

   TArc *SParc3 = new TArc();
   SParc3->SetFillStyle(0);
   SParc3->SetLineColor(3);
   SParc3->SetLineWidth(1);

   TArc *SParc4 = new TArc();
   SParc4->SetFillStyle(0);
   SParc4->SetLineColor(7);
   SParc4->SetLineWidth(1);


   for (Long64_t jentry = num; jentry < num+1 ;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

// set valiable //

      int tsegn = 0;
      double tsegX = 0, tsegY = 0, tsegZ = 0, tsegR = 0;
      double tsegPx = 0, tsegPy = 0, tsegPz = 0, tsegPr = 0;
      double tsegSlope = 0, tsegIntercept = 0;

      int psegn = 0;
      double psegX = 0, psegY = 0, psegZ = 0, psegR = 0;
      double psegPx = 0, psegPy = 0, psegPz = 0, psegPr = 0;
      double psegSlope = 0, psegIntercept = 0;

      double L1eta = 0, L1slope = 0;
      double L1etamin = 0, L1slopemin = 0;
      double L1etamax = 0, L1slopemax = 0;
     
      double Roadwidth;

      double RoadZmin[10], RoadZmax[10];
      double RoadAw[10],RoadBw[10],RoadIntercept[10];
      double mdtRegionZmin[10]; 
      double mdtRegionZmax[10]; 
      double mdtRegionRmin[10]; 
      double mdtRegionRmax[10]; 
      double mdtRegionEtamin[10]; 
      double mdtRegionEtamax[10]; 
      double mdtRegionSlopemin[10]; 
      double mdtRegionSlopemax[10]; 
      double frameZmin[10],frameZmax[10]; 
      double frameRmin[10],frameRmax[10]; 

      for(int station; station < 10; station++){
         RoadZmin[station] = 0;
         RoadZmax[station] = 0;
         RoadAw[station] = 0;
         RoadBw[station] = 0;
         mdtRegionZmin[station] = 0;
         mdtRegionZmax[station] = 0;
         mdtRegionRmin[station] = 0;
         mdtRegionRmax[station] = 0;
         mdtRegionEtamin[station] = 0;
         mdtRegionEtamax[station] = 0;
         mdtRegionSlopemin[station] = 0;
         mdtRegionSlopemax[station] = 0;
         frameZmin[station] = 0;
         frameZmax[station] = 0;
         frameRmin[station] = 0;
         frameRmax[station] = 0;
      }

      double rpcZ[500], rpcR[500];
      for(int lr = 0; lr < 500; lr++){
         rpcZ[lr] = 0;
         rpcR[lr] = 0;
      }

      double mdtZ[500], mdtR[500], mdtSpace[500],mdtSigma[500];
      for(int lm = 0; lm < 500; lm++){
         mdtZ[lm] = 0;
         mdtR[lm] = 0;
         mdtSpace[lm] = 0;
         mdtSigma[lm] = 0;
      }

      float mdtZ_BI[10], mdtR_BI[10], mdtSpace_BI[10],mdtSigma_BI[10],mdtWeight_BI[10];
      for(int lm = 0; lm < 10; lm++){
         mdtZ_BI[lm] = 0;
         mdtR_BI[lm] = 0;
         mdtSpace_BI[lm] = 0;
         mdtSigma_BI[lm] = 0;
         mdtWeight_BI[lm] = 0;
      }

      float mdtZ_BM[10], mdtR_BM[10], mdtSpace_BM[10],mdtSigma_BM[10],mdtWeight_BM[10];
      for(int lm = 0; lm < 10; lm++){
         mdtZ_BM[lm] = 0;
         mdtR_BM[lm] = 0;
         mdtSpace_BM[lm] = 0;
         mdtSigma_BM[lm] = 0;
         mdtWeight_BM[lm] = 0;
      }

      float mdtZ_BO[10], mdtR_BO[10], mdtSpace_BO[10],mdtSigma_BO[10],mdtWeight_BO[10];
      for(int lm = 0; lm < 10; lm++){
         mdtZ_BO[lm] = 0;
         mdtR_BO[lm] = 0;
         mdtSpace_BO[lm] = 0;
         mdtSigma_BO[lm] = 0;
         mdtWeight_BO[lm] = 0;
      }

      int nMdt_BI = 0, nMdt_BM = 0,nMdt_BO = 0;

      int IGLin_BI[10];
      for(int im = 0; im < 10; im++){
         IGLin_BI[im] = 3;
      }

      float ALin_BI = 0, BLin_BI = 0, DAB_BI[2][2], Chi2_BI = 0, PChi2_BI = 0;
      float SlopeCand_BI[NCAND], InterceptCand_BI[NCAND], Chi2Cand_BI[NCAND];

      double SP_Z_BI = 0, SP_R_BI = 0, SP_S_BI = 0, SP_I_BI = 0, SP_C_BI = 0;
      double SP_Z_BM = 0, SP_R_BM = 0, SP_S_BM = 0, SP_I_BM = 0, SP_C_BM = 0;
      double SP_Z_BO = 0, SP_R_BO = 0, SP_S_BO = 0, SP_I_BO = 0, SP_C_BO = 0;

      double SP_Z_EI = 0, SP_R_EI = 0, SP_S_EI = 0, SP_I_EI = 0, SP_C_EI = 0;
      double SP_Z_EM = 0, SP_R_EM = 0, SP_S_EM = 0, SP_I_EM = 0, SP_C_EM = 0;
      double SP_Z_EO = 0, SP_R_EO = 0, SP_S_EO = 0, SP_I_EO = 0, SP_C_EO = 0;

      double SP_Z_EE = 0, SP_R_EE = 0, SP_S_EE = 0, SP_I_EE = 0, SP_C_EE = 0;
      double SP_Z_CSC = 0, SP_R_CSC = 0, SP_S_CSC = 0, SP_I_CSC = 0, SP_C_CSC = 0;
      double SP_Z_BEE = 0, SP_R_BEE = 0, SP_S_BEE = 0, SP_I_BEE = 0, SP_C_BEE = 0;
      double SP_Z_BME = 0, SP_R_BME = 0, SP_S_BME = 0, SP_I_BME = 0, SP_C_BME = 0;

      float Cz1[1],Cr1[1],Radius1[1];
      float Cz2[1],Cr2[1],Radius2[1];
      float Cz3[1],Cr3[1],Radius3[1];
      float Cz4[1],Cr4[1],Radius4[1];

      int isL1pass = 0;

// set valiable end//

// fanction //

      TF1 *ftsegslope[10];
      TF1 *fpsegslope[10];

      TF1 *froicentor;
      TF1 *froimin;
      TF1 *froimax;

      TF1 *froad[10];
      TF1 *froadZmin[10];
      TF1 *froadZmax[10];

      TF1 *fmdtRegionRmin[10];
      TF1 *fmdtRegionRmax[10];
      TF1 *fmdtRegionEtamin[10];
      TF1 *fmdtRegionEtamax[10];

      TBox *bMdtRegion[10];

      TF1 *fSPBIslope;
      TF1 *fSPBMslope;
      TF1 *fSPBOslope;
      TF1 *fSPEIslope;
      TF1 *fSPEMslope;
      TF1 *fSPEOslope;
      TF1 *fSPEEslope;
      TF1 *fSPCSCslope;
      TF1 *fSPBEEslope;
      TF1 *fSPBMEslope;

      TF1 *fSPBIinterslope;

      TF1 *fSPBIslope2;
      TF1 *fSPBMslope2;
      TF1 *fSPBOslope2;

// fanction end//
 
      if(fabs(probe_eta) > 1.05){
         continue;
      }
      
// tag Offline segment //

      tsegn = tag_segment_n;

      for(int k = 0; k < tsegn; k++){
         if(tag_segment_x[k] == -99999 || tag_segment_x[k] == -88888) continue;
         tsegX = tag_segment_x[k];
         tsegY = tag_segment_y[k];
         tsegZ = tag_segment_z[k];
         tsegR = sqrt(tsegX * tsegX + tsegY * tsegY);
         tsegPx = tag_segment_px[k];
         tsegPy = tag_segment_py[k];
         tsegPz = tag_segment_pz[k];
         tsegPr = sqrt(tsegPx * tsegPx + tsegPy * tsegPy);

         tsegSlope = tsegPr/tsegPz;   
         tsegIntercept = tsegR - (tsegSlope * tsegZ);   

         tsegZR->SetPoint(k, tsegZ, tsegR);

         double length = 200 / ( 1 + sqrt(tsegPz * tsegPz + tsegPr * tsegPr) * sqrt(tsegPz * tsegPz + tsegPr * tsegPr));
//         std::cout << "length = " << length << std::endl;

         ftsegslope[k] = new TF1("f", "[0]*x+[1]", tsegZ-length, tsegZ+length);
         ftsegslope[k]->SetParameter(0, tsegSlope);
         ftsegslope[k]->SetParameter(1, tsegIntercept);
         ftsegslope[k]->SetLineColor(1);
         ftsegslope[k]->SetLineWidth(1);
         
      }  

// tag Offline segment //

// probe Offline segment //

      psegn = probe_segment_n;

      for(int k = 0; k < psegn; k++){
         if(probe_segment_x[k] == -99999 || probe_segment_x[k] == -88888) continue;
         psegX = probe_segment_x[k];
         psegY = probe_segment_y[k];
         psegZ = probe_segment_z[k];
         psegR = sqrt(psegX * psegX + psegY * psegY);
         psegPx = probe_segment_px[k];
         psegPy = probe_segment_py[k];
         psegPz = probe_segment_pz[k];
         psegPr = sqrt(psegPx * psegPx + psegPy * psegPy);

         psegSlope = psegPr/psegPz;   
         psegIntercept = psegR - (psegSlope * psegZ);   

         psegZR->SetPoint(k, psegZ, psegR);

         double length = 200 / ( 1 + sqrt(psegPz * psegPz + psegPr * psegPr) * sqrt(psegPz * psegPz + psegPr * psegPr));
//         std::cout << "length = " << length << std::endl;

         fpsegslope[k] = new TF1("f", "[0]*x+[1]", psegZ-length, psegZ+length);
         fpsegslope[k]->SetParameter(0, psegSlope);
         fpsegslope[k]->SetParameter(1, psegIntercept);
         fpsegslope[k]->SetLineColor(1);
         fpsegslope[k]->SetLineWidth(1);
         
      }  

// probe Offline segment end//

// trigger flag //

      for ( int j = 0; j < 25; j++){
  
         if(mes_name->at(j) == "mu4"){

// cout infomation //
            std::cout << "--------------------------------------------" << std::endl;
            std::cout << "probe_mesL1_pass = " << probe_mesL1_pass->at(j) <<std::endl;
            std::cout << "probe_mesL1_pt = " << probe_mesL1_pt->at(j) <<std::endl;
            std::cout << "probe_mesSA_pass = " << probe_mesSA_pass->at(j) <<std::endl;
            std::cout << "probe_mesSA_pt = " << probe_mesSA_pt->at(j) <<std::endl;

            std::cout << "--------------------------------------------" << std::endl;
// cout infomation end//

            if(probe_mesL1_pt->at(j) != -99999){

// roi centor //

              L1eta = probe_mesL1_eta->at(j);
              L1slope = tan(2*atan(exp(-L1eta)));

              froicentor = new TF1("froicentor", "[0]*x+[1]", -100000,100000);
              froicentor->SetParameter(0,L1slope);
              froicentor->SetParameter(1,0);
              froicentor->SetLineColor(kBlue);
              froicentor->SetLineWidth(1);
              froicentor->SetLineStyle(4);


              L1etamin = probe_mesL1_eta->at(j) - 0.05;
              L1slopemin = tan(2*atan(exp(-L1etamin)));

              froimin = new TF1("froimin", "[0]*x+[1]", -100000,100000);
              froimin->SetParameter(0,L1slopemin);
              froimin->SetParameter(1,0);
              froimin->SetLineColor(kBlue);
              froimin->SetLineWidth(1);
              froimin->SetLineStyle(8);


              L1etamax = probe_mesL1_eta->at(j) + 0.05;
              L1slopemax = tan(2*atan(exp(-L1etamax)));

              froimax = new TF1("froimax", "[0]*x+[1]", -100000,100000);
              froimax->SetParameter(0,L1slopemax);
              froimax->SetParameter(1,0);
              froimax->SetLineColor(kBlue);
              froimax->SetLineWidth(1);
              froimax->SetLineStyle(8);

// roi centor end //

// road //

              for( int is = 0; is < 10; is ++){

                 if (is == 0 ){
                    Roadwidth = 400; 
                 }
                 else if (is == 1){
                    Roadwidth = 200;
                 }
                 else if (is == 2){
                    Roadwidth = 400;
                 }
                 else{
                    Roadwidth = 0;
                 }

/*
                 if (is == 0 ){
                    Roadwidth = 500; 
                 }
                 else if (is == 1){
                    Roadwidth = 650;
                 }
                 else if (is == 2){
                    Roadwidth = 800;
                 }
                 else{
                    Roadwidth = 0;
                 }
*/

                 RoadAw[is] = probe_mesSA_roadAw->at(j)[is];
                 RoadBw[is] = probe_mesSA_roadBw->at(j)[is];
                 RoadIntercept[is] = Roadwidth * (sqrt(RoadAw[is]*RoadAw[is] + 1));

                 if(is == 0){
                 RoadZmin[is] = ( 4000 - RoadBw[is] ) / RoadAw[is];
                 RoadZmax[is] = ( 6000 - RoadBw[is] ) / RoadAw[is];
                 }
                 else if(is == 1){
                 RoadZmin[is] = ( 6000 - RoadBw[is] ) / RoadAw[is];
                 RoadZmax[is] = ( 9000 - RoadBw[is] ) / RoadAw[is];
                 }
                 else if(is == 2){
                 RoadZmin[is] = ( 9000 - RoadBw[is] ) / RoadAw[is];
                 RoadZmax[is] = ( 15000 - RoadBw[is] ) / RoadAw[is];
                 }

                 std::cout << "roadAw = " << RoadAw[is] << std::endl;

                 froad[is] = new TF1("froad","[0]*x+[1]",RoadZmin[is],RoadZmax[is]);
                 froad[is]->SetParameter(0,RoadAw[is]);
                 froad[is]->SetParameter(1,RoadBw[is]);
                 froad[is]->SetLineColor(8);
                 froad[is]->SetLineWidth(1);
                 froad[is]->SetLineStyle(3);

                 froadZmin[is] = new TF1("froadZmin","[0]*x+[1]-[2]",-100000,100000);
                 froadZmin[is]->SetParameter(0,RoadAw[is]);
                 froadZmin[is]->SetParameter(1,RoadBw[is]);
                 froadZmin[is]->SetParameter(2,RoadIntercept[is]);
                 froadZmin[is]->SetLineColor(8);
                 froadZmin[is]->SetLineWidth(2);
//                 froadZmin[is]->SetLineStyle(3);

                 froadZmax[is] = new TF1("froadZmax","[0]*x+[1]+[2]",-100000,100000);
                 froadZmax[is]->SetParameter(0,RoadAw[is]);
                 froadZmax[is]->SetParameter(1,RoadBw[is]);
                 froadZmax[is]->SetParameter(2,RoadIntercept[is]);
                 froadZmax[is]->SetLineColor(8);
                 froadZmax[is]->SetLineWidth(2);
//                 froadZmax[is]->SetLineStyle(3);
// road end //

// mdt region //

                 mdtRegionZmin[is] = probe_mesSA_zMin->at(j)[is];
                 mdtRegionZmax[is] = probe_mesSA_zMax->at(j)[is];
                 mdtRegionRmin[is] = probe_mesSA_rMin->at(j)[is];
                 mdtRegionRmax[is] = probe_mesSA_rMax->at(j)[is];
                 mdtRegionEtamin[is] = probe_mesSA_etaMin->at(j)[is];
                 mdtRegionEtamax[is] = probe_mesSA_etaMax->at(j)[is];
                 mdtRegionSlopemin[is] = tan(2*atan(exp(-mdtRegionEtamin[is])));
                 mdtRegionSlopemax[is] = tan(2*atan(exp(-mdtRegionEtamax[is])));

                 if( mdtRegionEtamax[is] > 0){

                    double Z_frameLD =  mdtRegionRmin[is]/mdtRegionSlopemin[is];
                    double Z_frameRD =  mdtRegionRmin[is]/mdtRegionSlopemax[is];
                    double Z_frameLU =  mdtRegionRmax[is]/mdtRegionSlopemin[is];
                    double Z_frameRU =  mdtRegionRmax[is]/mdtRegionSlopemax[is];

                    fmdtRegionEtamin[is] = new TF1("fSPBIslope", "[0]*x+[1]",  Z_frameLD,Z_frameLU);
                    fmdtRegionEtamin[is]->SetParameter(0,mdtRegionSlopemin[is] );
                    fmdtRegionEtamin[is]->SetParameter(1,0);
                    fmdtRegionEtamin[is]->SetLineColor(1);
                    fmdtRegionEtamin[is]->SetLineWidth(1);

                    fmdtRegionEtamax[is] = new TF1("fSPBIslope", "[0]*x+[1]",  Z_frameRD,Z_frameRU);
                    fmdtRegionEtamax[is]->SetParameter(0,mdtRegionSlopemax[is] );
                    fmdtRegionEtamax[is]->SetParameter(1,0);
                    fmdtRegionEtamax[is]->SetLineColor(1);
                    fmdtRegionEtamax[is]->SetLineWidth(1);                    

                    fmdtRegionRmin[is] = new TF1("fmdtRegionRmin", "[0]", Z_frameLD, Z_frameRD);
                    fmdtRegionRmin[is]->SetParameter(0,mdtRegionRmin[is] );
                    fmdtRegionRmin[is]->SetLineColor(1);
                    fmdtRegionRmin[is]->SetLineWidth(1);

                    fmdtRegionRmax[is] = new TF1("fmdtRegionRmax", "[0]", Z_frameLU, Z_frameRU);
                    fmdtRegionRmax[is]->SetParameter(0,mdtRegionRmax[is] );
                    fmdtRegionRmax[is]->SetLineColor(1);
                    fmdtRegionRmax[is]->SetLineWidth(1);

                 }
                 else {

                    double Z_frameLD =  mdtRegionRmin[is]/mdtRegionSlopemin[is];
                    double Z_frameRD =  mdtRegionRmin[is]/mdtRegionSlopemax[is];
                    double Z_frameLU =  mdtRegionRmax[is]/mdtRegionSlopemin[is];
                    double Z_frameRU =  mdtRegionRmax[is]/mdtRegionSlopemax[is];

                    fmdtRegionEtamin[is] = new TF1("fSPBIslope", "[0]*x+[1]",  Z_frameRD, Z_frameRU);
                    fmdtRegionEtamin[is]->SetParameter(0,(mdtRegionSlopemax[is]) );
                    fmdtRegionEtamin[is]->SetParameter(1,0);
                    fmdtRegionEtamin[is]->SetLineColor(1);
                    fmdtRegionEtamin[is]->SetLineWidth(1);

                    fmdtRegionEtamax[is] = new TF1("fSPBIslope", "[0]*x+[1]",  Z_frameLD,Z_frameLU);
                    fmdtRegionEtamax[is]->SetParameter(0,(mdtRegionSlopemin[is]) );
                    fmdtRegionEtamax[is]->SetParameter(1,0);
                    fmdtRegionEtamax[is]->SetLineColor(1);
                    fmdtRegionEtamax[is]->SetLineWidth(1);                    

                    fmdtRegionRmin[is] = new TF1("fmdtRegionRmin", "[0]", Z_frameLD, Z_frameRD);
                    fmdtRegionRmin[is]->SetParameter(0,mdtRegionRmin[is] );
                    fmdtRegionRmin[is]->SetLineColor(1);
                    fmdtRegionRmin[is]->SetLineWidth(1);

                    fmdtRegionRmax[is] = new TF1("fmdtRegionRmax", "[0]", Z_frameLU, Z_frameRU);
                    fmdtRegionRmax[is]->SetParameter(0,mdtRegionRmax[is] );
                    fmdtRegionRmax[is]->SetLineColor(1);
                    fmdtRegionRmax[is]->SetLineWidth(1);

                 }

                 if(probe_eta > 0){
                    frameZmin[is] = mdtRegionRmin[is]/tan(2*atan(exp(-mdtRegionEtamin[is]))); 
                    frameZmax[is] = mdtRegionRmax[is]/tan(2*atan(exp(-mdtRegionEtamax[is])));                  
                    frameRmin[is] = mdtRegionRmin[is]; 
                    frameRmax[is] = mdtRegionRmax[is];                  
                 }
                 else{
                    frameZmin[is] = mdtRegionRmax[is]/tan(2*atan(exp(-mdtRegionEtamin[is]))); 
                    frameZmax[is] = mdtRegionRmin[is]/tan(2*atan(exp(-mdtRegionEtamax[is]))); 
                    frameRmin[is] = mdtRegionRmin[is]; 
                    frameRmax[is] = mdtRegionRmax[is]; 
                 }
/*
                 if(mdtRegionZmin[is]<mdtRegionZmax[is]){
                    frameZmin[is] = mdtRegionZmin[is]; 
                    frameZmax[is] = mdtRegionZmax[is]; 
                    frameRmin[is] = mdtRegionRmin[is]; 
                    frameRmax[is] = mdtRegionRmax[is]; 
                 }
                 else{
                    frameZmin[is] = mdtRegionZmax[is]; 
                    frameZmax[is] = mdtRegionZmin[is]; 
                    frameRmin[is] = mdtRegionRmin[is]; 
                    frameRmax[is] = mdtRegionRmax[is]; 
                 }
*/

                 bMdtRegion[is] = new TBox(mdtRegionZmin[is], mdtRegionRmin[is], mdtRegionZmax[is], mdtRegionRmax[is]);
              }

              mdtRegionLD_BI->SetPoint(0, mdtRegionZmin[0], mdtRegionRmin[0] );
              mdtRegionRD_BI->SetPoint(0, mdtRegionZmax[0], mdtRegionRmin[0] );
              mdtRegionLU_BI->SetPoint(0, mdtRegionZmin[0], mdtRegionRmax[0] );
              mdtRegionRU_BI->SetPoint(0, mdtRegionZmax[0], mdtRegionRmax[0] );
             
              mdtRegionLD_BM->SetPoint(1, mdtRegionZmin[1], mdtRegionRmin[1] );
              mdtRegionRD_BM->SetPoint(1, mdtRegionZmax[1], mdtRegionRmin[1] );
              mdtRegionLU_BM->SetPoint(1, mdtRegionZmin[1], mdtRegionRmax[1] );
              mdtRegionRU_BM->SetPoint(1, mdtRegionZmax[1], mdtRegionRmax[1] );
             
              mdtRegionLD_BO->SetPoint(2, mdtRegionZmin[2], mdtRegionRmin[2] );
              mdtRegionRD_BO->SetPoint(2, mdtRegionZmax[2], mdtRegionRmin[2] );
              mdtRegionLU_BO->SetPoint(2, mdtRegionZmin[2], mdtRegionRmax[2] );
              mdtRegionRU_BO->SetPoint(2, mdtRegionZmax[2], mdtRegionRmax[2] );            

// mdt region end //

// rpc Hit //

               std::cout << "rpc mapping start " << std::endl;
         
               std::cout << "number of rpc Hit = " << probe_mesSA_rpcHitZ->front().size() << std::endl;
               for ( int ir = 0; ir < probe_mesSA_rpcHitZ->front().size(); ir++){
                  rpcZ[ir] = probe_mesSA_rpcHitZ->at(j)[ir];
                  rpcR[ir] = probe_mesSA_rpcHitR->at(j)[ir];
                  if( ( rpcZ[ir] < -0.001 || rpcZ[ir] > 0.001 ) && ( rpcR[ir] < -0.001 || rpcR[ir] > 0.001 ) ){
//                     rpcZR->SetPoint(ir,fabs(rpcZ[ir]), rpcR[ir]);
                     rpcZR->SetPoint(ir,rpcZ[ir], rpcR[ir]);
                  }
               }
               std::cout << "rpc mapping end " << std::endl;

// rpc Hit end//

// mdt Hit //

               std::cout << "mdt mapping start " << std::endl;
   
               std::cout << "number of mdt Hit = " << probe_mesSA_mdtHitZ->front().size() << std::endl;
               for ( int im = 0; im < probe_mesSA_mdtHitZ->front().size(); im++){

                  mdtZ[im] = probe_mesSA_mdtHitZ->at(j)[im];
                  mdtR[im] = probe_mesSA_mdtHitR->at(j)[im];
                  //mdtSpace[im] = probe_mesSA_mdtHitSpace->at(j)[im];
                  mdtSigma[im] = 0;
//                  mdtSigma[im] = probe_mesSA_mdtHitSigma->at(j)[im];

                  if( ( mdtZ[im] < -0.001 || mdtZ[im] > 0.001) &&  ( mdtR[im] < -0.001 || mdtR[im] > 0.001)  ){
//                     mdtZR->SetPoint(im,fabs(mdtZ[im]), mdtR[im]);

                     if(probe_mesSA_mdtHitIsOutlier->at(j)[im] == 0){
                        mdtZR->SetPoint(im,mdtZ[im], mdtR[im]);
                     }

                     else if(probe_mesSA_mdtHitIsOutlier->at(j)[im] == 1){
                        mdtZRIs1->SetPoint(im,mdtZ[im], mdtR[im]);
                     }

                     else if(probe_mesSA_mdtHitIsOutlier->at(j)[im] == 2){
                        mdtZRIs2->SetPoint(im,mdtZ[im], mdtR[im]);
                     }

                     else if(probe_mesSA_mdtHitIsOutlier->at(j)[im] == 3){
                        mdtZRIs3->SetPoint(im,mdtZ[im], mdtR[im]);
                     }

                     if(probe_mesSA_mdtHitIsOutlier->at(j)[im] == 0 && probe_mesSA_mdtHitR->at(j)[im] > 4000 && probe_mesSA_mdtHitR->at(j)[im] < 6000){
                        mdtZ_BI[im] = mdtZ[im];
                        mdtR_BI[im] = mdtR[im];
                        mdtSpace_BI[im] = mdtSpace[im];
                        mdtSigma_BI[im] = mdtSigma[im];
                        mdtWeight_BI[im] = 1 / (mdtSigma[im] * mdtSigma[im]);
                        nMdt_BI ++;
                     }

                     else if(probe_mesSA_mdtHitIsOutlier->at(j)[im] == 0 && probe_mesSA_mdtHitR->at(j)[im] > 6500 && probe_mesSA_mdtHitR->at(j)[im] < 9000){
                        mdtZ_BM[im - nMdt_BI] = mdtZ[im];
                        mdtR_BM[im - nMdt_BI] = mdtR[im];
                        mdtSpace_BM[im - nMdt_BI] = mdtSpace[im];
                        mdtSigma_BM[im - nMdt_BI] = mdtSigma[im];
                        mdtWeight_BM[im - nMdt_BI] = 1 / (mdtSigma[im] * mdtSigma[im]);
                        nMdt_BM ++;
                     }
                     else if(probe_mesSA_mdtHitIsOutlier->at(j)[im] == 0 && probe_mesSA_mdtHitR->at(j)[im] > 9000 && probe_mesSA_mdtHitR->at(j)[im] < 11000){
                        mdtZ_BO[im - nMdt_BI - nMdt_BM] = mdtZ[im];
                        mdtR_BO[im - nMdt_BI - nMdt_BM] = mdtR[im];
                        mdtSpace_BO[im - nMdt_BI - nMdt_BM] = mdtSpace[im];
                        mdtSigma_BO[im - nMdt_BI - nMdt_BM] = mdtSigma[im];
                        mdtWeight_BO[im - nMdt_BI - nMdt_BM] = 1 / (mdtSigma[im] * mdtSigma[im]);
                        nMdt_BO ++;
                     }
                  }
               }

               std::cout << "---number of IsOutlier0 MDT hits for each station---" << std::endl;
               std::cout << "inner station = " << nMdt_BI << std::endl;
               std::cout << "middle station = " << nMdt_BM << std::endl;
               std::cout << "outer station = " << nMdt_BO << std::endl;

               std::cout << "mdt mapping end " << std::endl;

               Circles ( nMdt_BI, mdtR_BI, mdtZ_BI, mdtSpace_BI, mdtSigma_BI, IGLin_BI, &ALin_BI, &BLin_BI, DAB_BI, &Chi2_BI, &PChi2_BI, SlopeCand_BI, InterceptCand_BI, Chi2Cand_BI );   
               
               std::cout << "---Circle fit result---"<< std::endl;
               std::cout << "inner station : SlopeCand_BI, InterceptCand_BI, Chi2Cand_BI = "<< SlopeCand_BI[0]<<", "<< InterceptCand_BI[0]<<", "<< Chi2Cand_BI[0]<<std::endl;
               std::cout << "inner station : SlopeCand_BI, InterceptCand_BI, Chi2Cand_BI = "<< SlopeCand_BI[1]<<", "<< InterceptCand_BI[1]<<", "<< Chi2Cand_BI[1]<<std::endl;
               std::cout << "inner station : SlopeCand_BI, InterceptCand_BI, Chi2Cand_BI = "<< SlopeCand_BI[2]<<", "<< InterceptCand_BI[2]<<", "<< Chi2Cand_BI[2]<<std::endl;


// mdt Hit end //

// MuonSA SuperPoint //

               std::cout << "SP mapping start " << std::endl;

               SP_Z_BI = probe_mesSA_superPointZ_BI->at(j);
               SP_Z_BM = probe_mesSA_superPointZ_BM->at(j);
               SP_Z_BO = probe_mesSA_superPointZ_BO->at(j);
               SP_Z_EI = probe_mesSA_superPointZ_EI->at(j);
               SP_Z_EM = probe_mesSA_superPointZ_EM->at(j);
               SP_Z_EO = probe_mesSA_superPointZ_EO->at(j);
               SP_Z_EE = probe_mesSA_superPointZ_EE->at(j);
               SP_Z_CSC = probe_mesSA_superPointZ_CSC->at(j);
               SP_Z_BEE = probe_mesSA_superPointZ_BEE->at(j);
               SP_Z_BME = probe_mesSA_superPointZ_BME->at(j);
   
               SP_R_BI = probe_mesSA_superPointR_BI->at(j);
               SP_R_BM = probe_mesSA_superPointR_BM->at(j);
               SP_R_BO = probe_mesSA_superPointR_BO->at(j);
               SP_R_EI = probe_mesSA_superPointR_EI->at(j);
               SP_R_EM = probe_mesSA_superPointR_EM->at(j);
               SP_R_EO = probe_mesSA_superPointR_EO->at(j);
               SP_R_EE = probe_mesSA_superPointR_EE->at(j);
               SP_R_CSC = probe_mesSA_superPointR_CSC->at(j);
               SP_R_BEE = probe_mesSA_superPointR_BEE->at(j);
               SP_R_BME = probe_mesSA_superPointR_BME->at(j);
   
               SP_S_BI = probe_mesSA_superPointSlope_BI->at(j);
               SP_S_BM = probe_mesSA_superPointSlope_BM->at(j);
               SP_S_BO = probe_mesSA_superPointSlope_BO->at(j);
               SP_S_EI = probe_mesSA_superPointSlope_EI->at(j);
               SP_S_EM = probe_mesSA_superPointSlope_EM->at(j);
               SP_S_EO = probe_mesSA_superPointSlope_EO->at(j);
               SP_S_EE = probe_mesSA_superPointSlope_EE->at(j);
               SP_S_CSC = probe_mesSA_superPointSlope_CSC->at(j);
               SP_S_BEE = probe_mesSA_superPointSlope_BEE->at(j);
               SP_S_BME = probe_mesSA_superPointSlope_BME->at(j);
   
               SP_I_BI = probe_mesSA_superPointIntercept_BI->at(j);
               SP_I_BM = probe_mesSA_superPointIntercept_BM->at(j);
               SP_I_BO = probe_mesSA_superPointIntercept_BO->at(j);
               SP_I_EI = probe_mesSA_superPointIntercept_EI->at(j);
               SP_I_EM = probe_mesSA_superPointIntercept_EM->at(j);
               SP_I_EO = probe_mesSA_superPointIntercept_EO->at(j);
               SP_I_EE = probe_mesSA_superPointIntercept_EE->at(j);
               SP_I_CSC = probe_mesSA_superPointIntercept_CSC->at(j);
               SP_I_BEE = probe_mesSA_superPointIntercept_BEE->at(j);
               SP_I_BME = probe_mesSA_superPointIntercept_BME->at(j);

               SP_C_BI = probe_mesSA_superPointChi2_BI->at(j);
               SP_C_BM = probe_mesSA_superPointChi2_BM->at(j);
               SP_C_BO = probe_mesSA_superPointChi2_BO->at(j);
               SP_C_EI = probe_mesSA_superPointChi2_EI->at(j);
               SP_C_EM = probe_mesSA_superPointChi2_EM->at(j);
               SP_C_EO = probe_mesSA_superPointChi2_EO->at(j);
               SP_C_EE = probe_mesSA_superPointChi2_EE->at(j);
               SP_C_CSC = probe_mesSA_superPointChi2_CSC->at(j);
               SP_C_BEE = probe_mesSA_superPointChi2_BEE->at(j);
               SP_C_BME = probe_mesSA_superPointChi2_BME->at(j);
   
               SPBIZR->SetPoint(0,SP_Z_BI, SP_R_BI);
               SPBMZR->SetPoint(0,SP_Z_BM, SP_R_BM);
               SPBOZR->SetPoint(0,SP_Z_BO, SP_R_BO);
               SPEIZR->SetPoint(0,SP_Z_EI, SP_R_EI);
               SPEMZR->SetPoint(0,SP_Z_EM, SP_R_EM);
               SPEOZR->SetPoint(0,SP_Z_EO, SP_R_EO);
               SPEEZR->SetPoint(0,SP_Z_EE, SP_R_EE);
               SPCSCZR->SetPoint(0,SP_Z_CSC, SP_R_CSC);
               SPBEEZR->SetPoint(0,SP_Z_BEE, SP_R_BEE);
               SPBMEZR->SetPoint(0,SP_Z_BME, SP_R_BME);

               std::cout << "SP mapping end " << std::endl;
        
  
               std::cout << "SPslope mapping start " << std::endl;
   
               double SPBIZwidth = 200/sqrt( 1 + (1/( SP_S_BI )) * (1/( SP_S_BI )) );
               std::cout << "SPBIZwidth = " << SPBIZwidth << std::endl;
//               fSPBIslope = new TF1("fSPBIslope", "[0]*x+[1]",  SP_Z_BI - SPBIZwidth, SP_Z_BI + SPBIZwidth);
               fSPBIslope = new TF1("fSPBIslope", "[0]*x+[1]",  SP_Z_BI - 100000, SP_Z_BI + 100000);
               fSPBIslope->SetParameter(0,1/SP_S_BI );
               fSPBIslope->SetParameter(1,SP_R_BI - (SP_Z_BI / SP_S_BI) );
               fSPBIslope->SetLineColor(kPink);
               fSPBIslope->SetLineWidth(1);
               fSPBIslope->SetLineStyle(2);
   
               double SPBMZwidth = 200/sqrt( 1 + (1/( SP_S_BM )) * (1/( SP_S_BM )) );
               std::cout << "SPBMZwidth = " << SPBMZwidth << std::endl;
//               fSPBMslope = new TF1("fSPBMslope", "[0]*x+[1]",  SP_Z_BM  - SPBMZwidth, SP_Z_BM + SPBMZwidth);
               fSPBMslope = new TF1("fSPBMslope", "[0]*x+[1]",  SP_Z_BM  - 100000, SP_Z_BM + 100000);
               fSPBMslope->SetParameter(0,1/SP_S_BM );
               fSPBMslope->SetParameter(1,SP_R_BM - (SP_Z_BM / SP_S_BM) );
               fSPBMslope->SetLineColor(kPink);
               fSPBMslope->SetLineWidth(1);
               fSPBMslope->SetLineStyle(2);
   
               double SPBOZwidth = 200/sqrt( 1 + ( 1/( SP_S_BO )) * (1/( SP_S_BO )) );
               std::cout << "SPBOZwidth = " << SPBOZwidth << std::endl;
//               fSPBOslope = new TF1("fSPBOslope", "[0]*x+[1]",  SP_Z_BO  - SPBOZwidth, SP_Z_BO + SPBOZwidth);
               fSPBOslope = new TF1("fSPBOslope", "[0]*x+[1]",  SP_Z_BO  - 100000, SP_Z_BO + 100000);
               fSPBOslope->SetParameter(0,1/SP_S_BO );
               fSPBOslope->SetParameter(1,SP_R_BO - (SP_Z_BO / SP_S_BO) );
               fSPBOslope->SetLineColor(kPink);
               fSPBOslope->SetLineWidth(1);
               fSPBOslope->SetLineStyle(2);
   
               std::cout << "SPslope mapping end " << std::endl;


               fSPBIinterslope = new TF1("fSPBIOslope", "[0]*x+[1]",  SP_Z_BI - 100000, SP_Z_BI + 100000);
               fSPBIinterslope->SetParameter(0,SP_R_BI/SP_Z_BI );
               fSPBIinterslope->SetParameter(1,0);
               fSPBIinterslope->SetLineColor(kPink);
               fSPBIinterslope->SetLineWidth(1);
               fSPBIinterslope->SetLineStyle(2);


               std::cout << "SPslope2 mapping" << std::endl;

               std::cout << "SP intercept = " << SP_I_BI << ", " << SP_I_BM << ", " << SP_I_BO << std::endl; 

               fSPBIslope2 = new TF1("fSPBIslope", "[0]*x+[1]", -100000, 100000);
               fSPBIslope2->SetParameter(0,1/SP_S_BI );
               fSPBIslope2->SetParameter(1,mdtR_BI[0] - ( (mdtZ_BI[0] + SP_I_BI)/SP_S_BI ) );
               std::cout << "SP absolute slope BI = " << 1/SP_S_BI << std::endl; 
               std::cout << "SP absolute intercept BI = " << mdtR_BI[0] - ( (mdtZ_BI[0] + SP_I_BI)/SP_S_BI ) << std::endl; 
               fSPBIslope2->SetLineColor(kRed);
               fSPBIslope2->SetLineWidth(2);
               fSPBIslope2->SetLineStyle(1);

               fSPBMslope2 = new TF1("fSPBMslope", "[0]*x+[1]", -100000, 100000);
               fSPBMslope2->SetParameter(0,1/SP_S_BM );
               fSPBMslope2->SetParameter(1,mdtR_BM[0] - ( (mdtZ_BM[0] + SP_I_BM)/SP_S_BM ));
               std::cout << "mdtR_BM[0] = " << mdtR_BM[0] << std::endl; 
               std::cout << "SP absolute slope BM = " << 1/SP_S_BM << std::endl; 
               std::cout << "SP absolute intercept BM = " << mdtR_BM[0] - ( (mdtZ_BM[0] + SP_I_BM)/SP_S_BM ) << std::endl; 
               fSPBMslope2->SetLineColor(kRed);
               fSPBMslope2->SetLineWidth(2);
               fSPBMslope2->SetLineStyle(1);

               fSPBOslope2 = new TF1("fSPBMslope", "[0]*x+[1]", -100000, 100000);
               fSPBOslope2->SetParameter(0,1/SP_S_BO );
               fSPBOslope2->SetParameter(1,mdtR_BO[0] - ( (mdtZ_BO[0] + SP_I_BO)/SP_S_BO ) );
               std::cout << "SP absolute slope BO = " << 1/SP_S_BO << std::endl; 
               std::cout << "SP absolute intercept BO = " << mdtR_BO[0] - ( (mdtZ_BO[0] + SP_I_BO)/SP_S_BO ) << std::endl; 
               fSPBOslope2->SetLineColor(kRed);
               fSPBOslope2->SetLineWidth(2);
               fSPBOslope2->SetLineStyle(1);


               std::cout << "SPslope2 mapping end" << std::endl;


// MuonSA SuperPoint end //
   
// SuperPoint Fitting //  

               if(fabs(SP_Z_BI) > 0.001 && fabs(SP_Z_BM) > 0.001 && fabs(SP_Z_BO) > 0.001){
                  Circle3point(SP_Z_BI, SP_R_BI, SP_Z_BM, SP_R_BM, SP_Z_BO, SP_R_BO, Cz1, Cr1, Radius1);
                  std::cout << "Cz1, Cr1, Radius1 = " << Cz1[0] << ", "<< Cr1[0] << ", " << Radius1[0] << std::endl;     
               }

               if(fabs(SP_Z_BI) > 0.001 && fabs(SP_Z_BM) > 0.001){
                  Circle2point1slope(SP_Z_BI, SP_R_BI, SP_Z_BM, SP_R_BM, 1/SP_S_BM, Cz2, Cr2, Radius2);
                  std::cout << "Cz2, Cr2, Radius2 = " << Cz2[0] << ", "<< Cr2[0] << ", " << Radius2[0] << std::endl;     
               }

               if(fabs(SP_Z_BM) > 0.001 && fabs(SP_Z_BO) > 0.001){
                  Circle2point1slope(SP_Z_BM, SP_R_BM, SP_Z_BO, SP_R_BO, 1/SP_S_BO, Cz3, Cr3, Radius3);
                  std::cout << "Cz3, Cr3, Radius3 = " << Cz3[0] << ", "<< Cr3[0] << ", " << Radius3[0] << std::endl;     
               }

               double SP_Slope_BIO = ( SP_R_BI - 0 ) / ( SP_Z_BI - 0);

               if(fabs(SP_Z_BI) > 0.001 && fabs(SP_Z_BM) > 0.001){
                  Circle2point1slope(SP_Z_BM, SP_R_BM, SP_Z_BI, SP_R_BI, SP_Slope_BIO, Cz4, Cr4, Radius4);
                  std::cout << "Cz4, Cr4, Radius4 = " << Cz4[0] << ", "<< Cr4[0] << ", " << Radius4[0] << std::endl;     
               }
// SuperPoint Fitting end //      

            }
            else{
               std::cout << "L1 not pass! " << std::endl;
               isL1pass = -1;
            }
         }
      }

   if(isL1pass != -1){      

      TCanvas *c1 = new TCanvas();
      c1->Divide(2, 2);

// inner staton //

      c1->cd(1);
      TH1F *frame1 = c1->DrawFrame( frameZmin[0] - 100, frameRmin[0] - 100, frameZmax[0] + 100,  frameRmax[0] + 100);  
      frame1->GetXaxis()->SetTitle("Z[mm]");
      frame1->GetYaxis()->SetTitle("R[mm]");
      frame1->SetTitle("inner");

//      bMdtRegion[0]->Draw("f");

      tsegZR->Draw("P,same");      
      for(int k = 0; k < tsegn; k++){
         ftsegslope[k]->Draw("same");
      }

      psegZR->Draw("P,same");      
      for(int k = 0; k < psegn; k++){
         fpsegslope[k]->Draw("same");
      }

      froicentor->Draw("same");
      froimin->Draw("same");
      froimax->Draw("same");

 
//      mdtRegionLD_BI->Draw("P,same");
//      mdtRegionRD_BI->Draw("P,same");
//      mdtRegionLU_BI->Draw("P,same");
//      mdtRegionRU_BI->Draw("P,same");


      fmdtRegionRmin[0]->Draw("l.same");
      fmdtRegionRmax[0]->Draw("l,same");
      fmdtRegionEtamin[0]->Draw("l.same");
      fmdtRegionEtamax[0]->Draw("l,same");

//      mdtZR->Draw("P,same");

      for ( int j = 0; j < 25; j++){
         if(mes_name->at(j) == "mu4"){
            for(int im = 0; im < probe_mesSA_mdtHitZ->front().size(); im++){
               if(probe_mesSA_mdtHitIsOutlier->at(j)[im] == 3){
                  mdtZRIs3->Draw("P,same");
               }
               if(probe_mesSA_mdtHitIsOutlier->at(j)[im] == 2){
                  mdtZRIs2->Draw("P,same");
               }
               if(probe_mesSA_mdtHitIsOutlier->at(j)[im] == 1){
                  mdtZRIs1->Draw("P,same");
               }
               if(probe_mesSA_mdtHitIsOutlier->at(j)[im] == 0){
                  mdtarc->DrawArc(mdtZ[im],mdtR[im],mdtSpace[im],0,360,"same");
               }
            }
         }
      }

      froad[0]->Draw("same");
      froadZmin[0]->Draw("same");
      froadZmax[0]->Draw("same");

      SPBIZR->Draw("P,same");
      fSPBIslope->Draw("same");
      fSPBIslope2->Draw("same");

// inner staton end //

//middle station //

      c1->cd(2);
      TH1F *frame2 = c1->DrawFrame( frameZmin[1] - 100, frameRmin[1] - 100, frameZmax[1] + 100,  frameRmax[1] + 100);  
      frame2->GetXaxis()->SetTitle("Z[mm]");
      frame2->GetYaxis()->SetTitle("R[mm]");
      frame2->SetTitle("middle");

//      bMdtRegion[1]->Draw("f");

      tsegZR->Draw("P,same");      
      for(int k = 0; k < tsegn; k++){
         ftsegslope[k]->Draw("same");
      }

      psegZR->Draw("P,same");      
      for(int k = 0; k < psegn; k++){
         fpsegslope[k]->Draw("same");
      }
 
      froicentor->Draw("same");
      froimin->Draw("same");
      froimax->Draw("same");

      rpcZR->Draw("P,same");


//      mdtRegionLD_BM->Draw("P,same");
//      mdtRegionRD_BM->Draw("P,same");
//      mdtRegionLU_BM->Draw("P,same");
//      mdtRegionRU_BM->Draw("P,same");


      fmdtRegionRmin[1]->Draw("l.same");
      fmdtRegionRmax[1]->Draw("l,same");
      fmdtRegionEtamin[1]->Draw("l,same");
      fmdtRegionEtamax[1]->Draw("l,same");

//      mdtZR->Draw("P,same");


      for ( int j = 0; j < 25; j++){
         if(mes_name->at(j) == "mu4"){
            for(int im = 0; im < probe_mesSA_mdtHitZ->front().size(); im++){
               if(probe_mesSA_mdtHitIsOutlier->at(j)[im] == 3){
                  mdtZRIs3->Draw("P,same");
               }
               if(probe_mesSA_mdtHitIsOutlier->at(j)[im] == 2){
                  mdtZRIs2->Draw("P,same");
               }
               if(probe_mesSA_mdtHitIsOutlier->at(j)[im] == 1){
                  mdtZRIs1->Draw("P,same");
               }
               if(probe_mesSA_mdtHitIsOutlier->at(j)[im] == 0){
                  mdtarc->DrawArc(mdtZ[im],mdtR[im],mdtSpace[im],0,360,"same");
               }
            }
         }
      }


      froad[1]->Draw("same");
      froadZmin[1]->Draw("same");
      froadZmax[1]->Draw("same");

      SPBMZR->Draw("P,same");
      fSPBMslope->Draw("same");
      fSPBMslope2->Draw("same");

//middle station end //

//outer station //

      c1->cd(3);
      TH1F *frame3 = c1->DrawFrame( frameZmin[2] - 100, frameRmin[2] - 100, frameZmax[2] + 100,  frameRmax[2] + 100);
      frame3->GetXaxis()->SetTitle("Z[mm]");
      frame3->GetYaxis()->SetTitle("R[mm]");
      frame3->SetTitle("outer");

//      bMdtRegion[2]->Draw("f");

      tsegZR->Draw("P,same");      
      for(int k = 0; k < tsegn; k++){
         ftsegslope[k]->Draw("same");
      }

      psegZR->Draw("P,same");      
      for(int k = 0; k < psegn; k++){
         fpsegslope[k]->Draw("same");
      }
 
      froicentor->Draw("same");
      froimin->Draw("same");
      froimax->Draw("same");

      rpcZR->Draw("P,same");


//      mdtRegionLD_BO->Draw("P,same");
//      mdtRegionRD_BO->Draw("P,same");
//      mdtRegionLU_BO->Draw("P,same");
//      mdtRegionRU_BO->Draw("P,same");


      fmdtRegionRmin[2]->Draw("l.same");
      fmdtRegionRmax[2]->Draw("l,same");
      fmdtRegionEtamin[2]->Draw("l,same");
      fmdtRegionEtamax[2]->Draw("l,same");

      for ( int j = 0; j < 25; j++){
         if(mes_name->at(j) == "mu4"){
            for(int im = 0; im < probe_mesSA_mdtHitZ->front().size(); im++){
               if(probe_mesSA_mdtHitIsOutlier->at(j)[im] == 3){
                  mdtZRIs3->Draw("P,same");
               }
               if(probe_mesSA_mdtHitIsOutlier->at(j)[im] == 2){
                  mdtZRIs2->Draw("P,same");
               }
               if(probe_mesSA_mdtHitIsOutlier->at(j)[im] == 1){
                  mdtZRIs1->Draw("P,same");
               }
               if(probe_mesSA_mdtHitIsOutlier->at(j)[im] == 0){
                  mdtarc->DrawArc(mdtZ[im],mdtR[im],mdtSpace[im],0,360,"same");

               }
            }
         }
      }

      froad[2]->Draw("same");
      froadZmin[2]->Draw("same");
      froadZmax[2]->Draw("same");

      SPBOZR->Draw("P,same");
      fSPBOslope->Draw("same");
      fSPBOslope2->Draw("same");

//outer station end //

//all station //

      c1->cd(4);


      TH1F *frame4 = c1->DrawFrame( -15000, 0, 15000, 15000);  
      frame4->GetXaxis()->SetTitle("Z[mm]");
      frame4->GetYaxis()->SetTitle("R[mm]");
      frame4->SetTitle("all");

/*
      if(frameZmin[0] < frameZmax[2]  ){
         TH1F *frame4 = c1->DrawFrame( frameZmin[0] - 500, frameRmin[0] - 500, frameZmax[2] + 500,  frameRmax[2] + 500);  
         frame4->GetXaxis()->SetTitle("Z[mm]");
         frame4->GetYaxis()->SetTitle("R[mm]");
         frame4->SetTitle("all");
      }
      else{
         TH1F *frame4 = c1->DrawFrame( frameZmin[2] - 500, frameRmin[0] - 500, frameZmax[0] + 500,  frameRmax[2] + 500);  
         frame4->GetXaxis()->SetTitle("Z[mm]");
         frame4->GetYaxis()->SetTitle("R[mm]");
         frame4->SetTitle("all");
      }
*/

//      bMdtRegion[0]->Draw("f");
//      bMdtRegion[1]->Draw("f,same");
//      bMdtRegion[2]->Draw("f,same");

      tsegZR->Draw("P,same");      

/*
      for(int k = 0; k < tsegn; k++){
         ftsegslope[k]->Draw("same");
      }
*/

      psegZR->Draw("P,same");
/*
      for(int k = 0; k < psegn; k++){
         fpsegslope[k]->Draw("same");
      }
*/

      froicentor->Draw("same");
      froimin->Draw("same");
      froimax->Draw("same");

/*
      mdtRegionLD_BI->Draw("P,same");
      mdtRegionRD_BI->Draw("P,same");
      mdtRegionLU_BI->Draw("P,same");
      mdtRegionRU_BI->Draw("P,same");

      mdtRegionLD_BM->Draw("P,same");
      mdtRegionRD_BM->Draw("P,same");
      mdtRegionLU_BM->Draw("P,same");
      mdtRegionRU_BM->Draw("P,same");

      mdtRegionLD_BO->Draw("P,same");
      mdtRegionRD_BO->Draw("P,same");
      mdtRegionLU_BO->Draw("P,same");
      mdtRegionRU_BO->Draw("P,same");
*/

      for(int is = 0; is < 3; is++){
         fmdtRegionRmin[is]->Draw("l.same");
         fmdtRegionRmax[is]->Draw("l,same");
         fmdtRegionEtamin[is]->Draw("l.same");
         fmdtRegionEtamax[is]->Draw("l,same");
      }

      rpcZR->Draw("P,same");

      froad[0]->Draw("same");
      froad[1]->Draw("same");
      froad[2]->Draw("same");

//      mdtZR->Draw("P,same");
//      mdtZRIs1->Draw("P,same");

      for ( int j = 0; j < 25; j++){
         if(mes_name->at(j) == "mu4"){
            for(int im = 0; im < probe_mesSA_mdtHitZ->front().size(); im++){
               if(probe_mesSA_mdtHitIsOutlier->at(j)[im] == 3){
                  mdtZRIs3->Draw("P,same");
               }
               if(probe_mesSA_mdtHitIsOutlier->at(j)[im] == 2){
                  mdtZRIs2->Draw("P,same");
               }
               if(probe_mesSA_mdtHitIsOutlier->at(j)[im] == 1){
                  mdtZRIs1->Draw("P,same");
               }
               if(probe_mesSA_mdtHitIsOutlier->at(j)[im] == 0){
                  mdtarc->DrawArc(mdtZ[im],mdtR[im],mdtSpace[im],0,360,"same");
               }
            }
         }
      }

      SPBIZR->Draw("P,same");
//      fSPBIslope->Draw("same");
      fSPBIinterslope->Draw("same");
      SPBMZR->Draw("P,same");
//      fSPBMslope->Draw("same");
      SPBOZR->Draw("P,same");
//      fSPBOslope->Draw("same");

/*
      SPEIZR->Draw("P,same");
      SPEMZR->Draw("P,same");
      SPEOZR->Draw("P,same");
      SPEEZR->Draw("P,same");
      SPCSCZR->Draw("P,same");
      SPBEEZR->Draw("P,same");
      SPBMEZR->Draw("P,same");
*/

//all station end //

      TCanvas *c2 = new TCanvas();

      TH1F *frame = c2->DrawFrame( -15000, 0, 15000, 15000);  
      frame->GetXaxis()->SetTitle("Z[mm]");
      frame->GetYaxis()->SetTitle("R[mm]");
      frame->SetTitle("Circle fit 1");

//      psegZR->Draw("P,same");      
//      for(int k = 0; k < psegn; k++){
//         fpsegslope[k]->Draw("same");
//      }

      SPBIZR->Draw("P,same");
//      fSPBIslope->Draw("same");
//      fSPBIinterslope->Draw("same");
      SPBMZR->Draw("P,same");
//      fSPBMslope->Draw("same");
      SPBOZR->Draw("P,same");
//      fSPBOslope->Draw("same");

      SParc1->DrawArc(Cz1[0],Cr1[0],Radius1[0],0,360,"same");
      SParc2->DrawArc(Cz2[0],Cr2[0],Radius2[0],0,360,"same");
//      SParc3->DrawArc(Cz3[0],Cr3[0],Radius3[0],0,360,"same");
      SParc4->DrawArc(Cz4[0],Cr4[0],Radius4[0],0,360,"same");

      }
   }
}


void NewEventDisplay::Circfit (int Nmeas,float *XI,float *YI,float *RI,float *WI,int *IG, float *A,float *B,float DAB[2][2], float *Chi2) {

  float XX[NMEAMX],YY[NMEAMX],Test,Toll,Xnor,Aold,Bold,Epsi;
  float SAA,SAB,SBB,Square;
  int j,Niter;

  Toll    = .1;
  Niter   = 0;
  //      *A      = 0.;
  //      *B      = 0.;
  //      SAA     = 0.;
  //      SAB     = 0.;
  //      SBB     = 0.;
  Square  = 0.;

  //    Many iterations ...

  do{
    Niter++;
    Xnor  = 1. / sqrt(1. + *A * *A);
    Aold  = *A;
    Bold  = *B;
      for(j=0;j<Nmeas;j++) {
        XX[j] = 0.;
        YY[j] = 0.;
      if(IG[j]==1) {
        XX[j] = XI[j];
        YY[j] = YI[j];
      }
      else if(IG[j]==2) {
        if(*A * XI[j] + *B - YI[j]>=0.) Epsi = 1.0;    // mod 961017
          else Epsi = -1.0;
          XX[j] = XI[j] - Epsi * Xnor * fabs(RI[j]) * *A;
          YY[j] = YI[j] + Epsi * Xnor * fabs(RI[j]);
        } else if(IG[j]==3) {
          XX[j] = XI[j] - Xnor * RI[j] * *A;
          YY[j] = YI[j] + Xnor * RI[j];
        }
      }

    Xline(XX,YY,WI,IG,Nmeas,A,B,&SAA,&SBB,&SAB,&Square);
    if(Square<=0.) break;
    Test = ((Aold-*A)*(Aold-*A))/ SAA + ((Bold-*B)*(Bold-*B))/ SBB;

  }

//  std::cout << "SAA, SAB, SBA, SBB = " << SAA << ", "<< SAB << ", " << SAB << ", "<< SBB << ", "<< std::endl;

  while(Test>=Toll&&Niter<=20);

  DAB[0][0] = SAA;
  DAB[0][1] = SAB;
  DAB[1][0] = SAB;
  DAB[1][1] = SBB;
  *Chi2     = Square;

  std::cout << "SAA, SAB, SBA, SBB, Chi2 = " << SAA << ", "<< SAB << ", " << SAB << ", "<< SBB << ", " << Square << std::endl;

}



void NewEventDisplay::Xline (float *X,float *Y,float *W,int *IG,int NP, float *A,float *B,float *SAA,float *SBB,float *SAB,float *Square) {

  int j;
  float S1,SX,SY,SXX,SXY,SYY,Deter,DY;

  *Square = -7.;
  S1      = 0.;
  SX      = 0.;
  SY      = 0.;
  SXX     = 0.;
  SXY     = 0.;
  SYY     = 0.;

  for(j=0;j<NP;j++) {
    if(IG[j]>=1) {
      S1  = S1  + W[j];
      SX  = SX  + W[j] * X[j];
      SY  = SY  + W[j] * Y[j];
      SXX = SXX + W[j] * X[j] * X[j];
      SXY = SXY + W[j] * X[j] * Y[j];
      SYY = SYY + W[j] * Y[j] * Y[j];
    }
  }
  std::cout << "S1, SX, SY, SXX, SXY, SYY = " << S1 << ", " << SX << ", "<< SY << ", "<< SXX << ", "<< SXY << ", "<< SYY << ", "<< std::endl;

  Deter  = S1 * SXX - SX * SX;

  if(Deter!=0.) {
    *A      = (S1 * SXY - SX * SY)  / Deter;
    *B      = (SY * SXX - SX * SXY) / Deter;
    *SAA    =   S1  / Deter;
    *SBB    =   SXX / Deter;
    *SAB    = - SX  / Deter;
    *Square = 0.;
    for(j=0;j<NP;j++) {
      if(IG[j]>=1) {
        DY =(Y[j] - *A * X[j] - *B)/sqrt(1 + *A * *A);

        //printf("Punto n.=%d , DY = %12.6f\n",j,DY);
        *Square = *Square + W[j] * DY * DY;
      }
    }
  }
}

// insert Circls() with 14 virticals //
// this code maybe makes SlopeCand, InterceptCand, and Chi2Cand, maybe//


void NewEventDisplay::Circles (int Nmeas,float *XI,float *YI,float *RI,float *WI,int *IG, float *A,float *B,float DAB[2][2],float *Chi2,float *Pchi2, float *SlopeCand, float *InterceptCand, float *Chi2Cand) {
   float RRi[NMEAMX],WIlim,CHbest,Abest,Bbest;
   float A0,B0,SAA,SBB,SAB,Square,Aj,Bj;
   int i,j,jj,Ntrue,Ngood,Nbest,Ntry,NOgo,Isig,Iflg,Ksign[4];
   int Ntrue2,Igg[NMEAMX];

   std::vector<float> st_chi2;
   std::vector<float> st_A;
   std::vector<float> st_B;

   for (int ii=0; ii<NCAND; ii++) {
     SlopeCand[ii]     = 0.;//99999;
     InterceptCand[ii] = 0.;//99999;
     Chi2Cand[ii]      = 0.;
   }

 //    Find the four "besrit" points (equispaced, wi.gt.1/2*wimax)

   *Chi2 = -1.;
   if (Nmeas<=2||Nmeas>=NMEAMX+1) return;
   Ntrue = 0;
   Ngood = 0;
   Nbest = 3;
   WIlim = 0.;

   Abest = 0.;
   Bbest = 0.;

   for (i=0;i<4;i++) Ksign[i] = -1; // i : 0,1,Nmeas-2,Nmeas-1

   for (j=0;j<Nmeas;j++) { //j : number of mdt hits
     if (IG[j]>=1) {
        Ntrue++;
        WIlim = (WIlim>=WI[j])? WIlim : WI[j];
        std::cout << "Ntrue = " << Ntrue << std::endl;
        std::cout << "WIlim = " << WIlim << std::endl;
     }
   }

   if (Ntrue<=2) return;

   WIlim = 0.1 * WIlim;

   for(j=0;j<Nmeas;j++) if(IG[j]>=1&&WI[j]>=WIlim) Ngood++;

   for (j=0;j<Nmeas;j++) {
     if (IG[j]>=1&&(WI[j]>=WIlim||Ngood<=3)) {

// select tube//

       if (Ksign[0]==-1.) {
         Ksign[0] = j;
       } else if (Ksign[1]==-1.) {
         Ksign[1] = j;
       } else if (Ksign[2]==-1.) {
         Ksign[2] = j;
       } else if(Ksign[3]==-1.) {
         Ksign[3] = j;
         Nbest    = 4;
       } else {
         Ksign[2] = Ksign[3];
         Ksign[3] = j;
       }
// select tube end//
     }
   }

   std::cout << "Ksign[0] = " << Ksign[0] << std::endl; 
   std::cout << "Ksign[1] = " << Ksign[1] << std::endl; 
   std::cout << "Ksign[2] = " << Ksign[2] << std::endl; 
   std::cout << "Ksign[3] = " << Ksign[3] << std::endl; 

   //    First attempt, try with a line through the centers

   Xline(XI,YI,WI,IG,Nmeas,&A0,&B0,&SAA,&SBB,&SAB,&Square);

 //    Then try 16 times trough the best four points
   st_A.clear(); st_B.clear(); st_chi2.clear();

   for (i=0;i<NMEAMX;i++) Igg[i] = -1;

   CHbest = 1.e25;
   Ntry   = (int)floor(pow(2.,Nbest)) - 1;         // 2**4 - 1 = 15

   for (j=0;j<=Ntry;j++) {

     NOgo = 0;

     for (jj=1;jj<=Nbest;jj++) {
       Isig = (j&(int)pow(2.,jj-1))? 1 : 0;
       //          Isig = ibits(&j,&jj1,&one);            // alternatively 0, 1
       Iflg = IG[Ksign[jj-1]];
       Igg[Ksign[jj-1]] = Iflg;
       RRi[Ksign[jj-1]] = RI[Ksign[jj-1]];

       if (Iflg==2) {

         Igg[Ksign[jj-1]] = 3;

        if (Isig==1) RRi[Ksign[jj-1]] = - RI[Ksign[jj-1]];

       } else if (Isig==1) {

         NOgo = 1;

       }
     }

//     std::cout << "NOgo = " << NOgo << std::endl;

     if (NOgo==0) {

       Aj = A0;
       Bj = B0;
       std::cout<< "---Circfit1---" << std::endl;

       Circfit(Nmeas,XI,YI,RRi,WI,Igg,&Aj,&Bj,DAB,Chi2);

       std::cout<< "---Circfit2---" << std::endl;

       Circfit(Nmeas,XI,YI,RI,WI,IG,&Aj,&Bj,DAB,Chi2);

       st_A.push_back(Aj); st_B.push_back(Bj); st_chi2.push_back(*Chi2);
 
       std::cout<< "---Circfit result---" << std::endl;
       std::cout << "Aj, Bj, *Chi2 = " << Aj << ", " << Bj << ", " << *Chi2<< std::endl;


       if (*Chi2>=0.0&&*Chi2<=CHbest) {
         Abest  = Aj;
         Bbest  = Bj;
         CHbest = *Chi2;
       }
     }
   }

   std::multimap<float, int>chi_map;
   chi_map.clear();
   std::vector<float> t_A;
   std::vector<float> t_B;
   std::vector<float> t_chi2;
   t_A.clear();
   t_B.clear();
   t_chi2.clear();

   for (int ir=0; ir<(int)st_chi2.size(); ir++) chi_map.insert(std::make_pair(st_chi2.at(ir), ir));

   for (std::multimap<float, int>::iterator jt = chi_map.begin(); jt != chi_map.end(); ++jt) {
     t_A.push_back(st_A.at(jt->second));
     t_B.push_back(st_B.at(jt->second));
     t_chi2.push_back(st_chi2.at(jt->second));
   }

   for (int nv=0; nv<6; nv++) {
     SlopeCand[nv]     = t_A[nv];
     InterceptCand[nv] = t_B[nv];
     Chi2Cand[nv]      = t_chi2[nv];
   }

    //    ... and finally with all the points
   *A = Abest;
   *B = Bbest;
   Circfit(Nmeas,XI,YI,RI,WI,IG,A,B,DAB,Chi2);

   std::cout << "finally A, B, *Chi2 = " << *A << ", " << *B << ", " << *Chi2<< std::endl;

   if (*Chi2>=0.0) {
     Ntrue2 = Ntrue - 2;
     *Pchi2 = TMath::Prob(*Chi2, Ntrue2);
  }

return;

}

void NewEventDisplay::Circle3point (double Z1, double R1, double Z2, double R2, double Z3, double R3, float *Cz, float *Cr, float *Radius) {

   float Zm12[1], Rm12[1], Slope12[1];
   float Zm23[1], Rm23[1], Slope23[1];

   for (int ii=0; ii<1; ii++) {

     Cz[ii]    = 0.;
     Cr[ii]    = 0.;
     Radius[ii] = 0.;

     Zm12[ii]    = 0.;
     Rm12[ii]    = 0.;
     Slope12[ii] = 0.;

     Zm23[ii]    = 0.;
     Rm23[ii]    = 0.;
     Slope23[ii] = 0.;

   }

   Bisector( Z1, R1, Z2, R2, Zm12, Rm12, Slope12 );
   Bisector( Z2, R2, Z3, R3, Zm23, Rm23, Slope23 );

   Cz[0] = (( Zm12[0] * Slope12[0] - Rm12[0])-( Zm23[0] * Slope23[0] -Rm23[0]) ) / (Slope12[0] - Slope23[0]);
   Cr[0] = ((-Zm23[0] * Slope23[0] + Rm23[0]) * Slope12[0] - (-Zm12[0] * Slope12[0] + Rm12[0]) * Slope23[0] ) / (Slope12[0] - Slope23[0]);
 
   Radius[0] = sqrt( ( Z1 - Cz[0] ) * ( Z1 - Cz[0] ) + ( R1 - Cr[0] ) * ( R1 - Cr[0] ) );
    
}

void NewEventDisplay::Circle2point1slope (double Z1, double R1, double Z2, double R2, double Slope2, float *Cz, float *Cr, float *Radius) {

   float Zm12[1], Rm12[1], Slope12[1];

   for (int ii=0; ii<1; ii++) {

     Cz[ii]    = 0.;
     Cr[ii]    = 0.;
     Radius[ii] = 0.;

     Zm12[ii]    = 0.;
     Rm12[ii]    = 0.;
     Slope12[ii] = 0.;

   }
   
   Bisector( Z1, R1, Z2, R2, Zm12, Rm12, Slope12 );

   Cz[0] = ( (Z2 / Slope2) + R2 + (Slope12[0] * Zm12[0]) - Rm12[0] ) / ( 1 / Slope2 + Slope12[0] );
   Cr[0] = ( ( (-Slope12[0] * Zm12[0]) + Rm12[0]) / Slope2 + Slope12[0] * ( Z2 / Slope2 + R2 ) ) / ( 1 / Slope2 + Slope12[0] );

   Radius[0] = sqrt( ( Z1 - Cz[0] ) * ( Z1 - Cz[0] ) + ( R1 - Cr[0] ) * ( R1 - Cr[0] ) );

}

void NewEventDisplay::Bisector(double Za, double Ra, double Zb, double Rb, float *Zm, float *Rm, float *Slope ){

  
   for (int ii=0; ii<1; ii++) {
     Zm[ii]    = 0.;
     Rm[ii]    = 0.;
     Slope[ii] = 0.;
   }

   std::cout << Za << ", " << Ra << ", " << Zb << ", " << Rb << std::endl; 

   Zm[0] = (Za + Zb) / 2;
   Rm[0] = (Ra + Rb) / 2;
   
   if (Za == Zb) {
      Slope[0] = 0;
   }
   else{
      Slope[0] = -1 / ( ( Ra - Rb ) / ( Za - Zb ) ); 
   }

   std::cout << Zm[0] << ", " << Rm[0] << ", "<< Slope[0] << std::endl; 

return; 

}

