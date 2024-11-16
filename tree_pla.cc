#include "TArtPlastic.hh"
#include "TArtPlasticPara.hh"
#include "TArtFSDBSDPla.hh"
#include "TArtHODPla.hh"

#include "TClonesArray.h"

#include "TArtBigRIPSParameters.hh"
#include "TArtSAMURAIParameters.hh"

#include "TArtCalibPlastic.hh"

#include "TArtCalibCoin.hh"
#include "TArtCalibFSDBSDPla.hh"
#include "TArtCalibHODPla.hh"

#include "TArtUtil.hh"
#include "TArtStoreManager.hh"
#include "TArtEventStore.hh"
#include "TArtEventInfo.hh"

#include <fstream>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <TSystem.h>
#include <TObjectTable.h>

using namespace std;


void tree_pla(Int_t nRun=372){

  TArtEventStore *estore = new TArtEventStore();
  estore->Open(Form("ridf/sdaq04/data%04d.ridf.gz",nRun));

  TArtBigRIPSParameters *brprm = TArtBigRIPSParameters::Instance();
  brprm->LoadParameter("db/SAMURAIPlastic.xml");
  
  // prepare tree
  ULong64_t neve=0;
  TTree *tree = new TTree("pla_tree","Plastic tree");
  
  // output file
  TFile *fout = new TFile(Form("rootfiles/tree/tree_pla%04d.root",nRun),"RECREATE");

  
  // Triger bit
  TArtCalibCoin *coin = new TArtCalibCoin();
  TClonesArray *triggerbit_array = (TClonesArray *)estore -> GetEventInfoArray();
  TArtEventInfo *triggerinfo = (TArtEventInfo *)triggerbit_array -> At(0);
  Int_t ftriggerbit[6] = {-1,-1,-1,-1,-1,-1};
  tree->Branch("ftriggerbit[6]",ftriggerbit,"ftriggerbit[6]/I");

  // Plastic
  TArtCalibPlastic *calibPlastic = new TArtCalibPlastic();
  TClonesArray* Plastic = (TClonesArray*)calibPlastic->GetPlasticArray();
  tree->Branch("Plastic",Plastic);

  // event
  tree->Branch("EventNum",&neve);
  
  // for p gate
  bool Bitp = false;
  tree->Branch("Bitp",&Bitp);
  
  // for d gate
  bool Bitd = false;
  tree->Branch("Bitd",&Bitd);
  
  // for 10Be gate
  bool Bit10Be = false;
  tree->Branch("Bit10Be",&Bit10Be);
   

  while(estore->GetNextEvent() //&& neve <100001
    ){

    //if (neve%1000==0) cout<<"\r  Running...     "<<neve<<" Events"<<flush;
   
    cout << "Event: " <<  neve << endl; 

    //trigger bit
    coin -> ClearData();
    coin -> LoadData();
    for(Int_t i=0;i<6;i++){
      ftriggerbit[i] = 0;
      ftriggerbit[i] = coin -> IsChTrue(i+1);
    }

    // Plastic
    calibPlastic->ClearData();
    calibPlastic->ReconstructData(); 
    
    TArtPlastic* PlasticF7 = (TArtPlastic*)TArtUtil::FindDataObject(Plastic,(char*)"F7pl");
    TArtPlastic* PlasticSBT1 = (TArtPlastic*)TArtUtil::FindDataObject(Plastic,(char*)"F13pl-1");
    TArtPlastic* PlasticSBT2 = (TArtPlastic*)TArtUtil::FindDataObject(Plastic,(char*)"F13pl-2");
    
    double F7_Time = PlasticF7->GetTime(); 
    double F7_Charge = PlasticF7->GetQAveRaw(); 

    double SBT1_Time = PlasticSBT1->GetTime(); 
    double SBT1_Charge = PlasticSBT1->GetQAveRaw(); 

    double SBT2_Time = PlasticSBT2->GetTime(); 
    double SBT2_Charge = PlasticSBT2->GetQAveRaw(); 
    
    double SBT_Time = (SBT1_Time+SBT2_Time)/2.; 
    double SBT_Charge = sqrt(SBT1_Charge*SBT2_Charge); 
    
    Bitp = false;
    if(SBT2_Time <= 73 && SBT2_Time >= 67 &&  SBT2_Charge <= 370 && SBT2_Charge >= 200){
      Bitp = true;
    }
    
    Bitd = false;
    if(SBT2_Time <= 73 && SBT2_Time >= 67 &&  SBT2_Charge <= 370 && SBT2_Charge >= 200){
      Bitd = true;
    }
    
    Bit10Be = false;
    if(SBT_Time-F7_Time<=-357 && SBT_Time-F7_Time>=-358 &&  F7_Charge<=680 && F7_Charge>=0){
      Bit10Be = true;
    }
    
    tree->Fill();
    
    estore->ClearData();
    neve ++;
  }
  
  cout << endl;
  tree->Write();
  fout->Close();

}//end of tree_pla__________________________________________________


