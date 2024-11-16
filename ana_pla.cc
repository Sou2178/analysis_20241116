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

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <TSystem.h>
#include <TObjectTable.h>

using namespace std;

void ana_pla(Int_t nRun=372) {

  // import file
  TFile *fpla = new TFile(Form("rootfiles/tree/tree_pla%04d.root",nRun));
 
  // import tree
  TTree *pla_tree = (TTree*) fpla->Get("pla_tree");
  
  // set branch
  Int_t ftriggerbit[6];
  bool Bitp=false;
  bool Bitd=false;
  bool Bit10Be=false;
  TClonesArray* Plastic = new TClonesArray("TArtPlastic");

  pla_tree->SetBranchAddress("ftriggerbit[6]",ftriggerbit);
  pla_tree->SetBranchAddress("Bitp",&Bitp);
  pla_tree->SetBranchAddress("Bitd",&Bitd);
  pla_tree->SetBranchAddress("Bit10Be",&Bit10Be);
  pla_tree->SetBranchAddress("Plastic",&Plastic);

  // output file
  TFile *fout = new TFile(Form("macros/chkpla/rootfile/ana_pla%04d.root",nRun),"RECREATE");

  // hist 
  TH1 *htrigger = new TH1F("htrigger","",6,0.5,6.5);
  TH1 *htrigger_d = new TH1F("htrigger_d","",6,0.5,6.5);
  TH1 *htrigger_be = new TH1F("htrigger_be","",6,0.5,6.5);
  
  TH1 *htdc3l = new TH1F("htdc3l","",1000,0,30000);
  TH1 *htdc5l = new TH1F("htdc5l","",1000,0,50000);
  TH1 *htdc7l = new TH1F("htdc7l","",1000,0,70000);
  TH1 *htdc1l = new TH1F("htdc1l","",1000,0,70000); // SBT1
  TH1 *htdc2l = new TH1F("htdc2l","",1000,0,70000); // SBT2
  TH1 *htdc3 = new TH1F("htdc3","",1000,0,30000);
  TH1 *htdc5 = new TH1F("htdc5","",1000,0,50000);
  TH1 *htdc7 = new TH1F("htdc7","",1000,0,70000);
  TH1 *htdc1 = new TH1F("htdc1","",1000,0,70000); // SBT1
  TH1 *htdc2 = new TH1F("htdc2","",1000,0,70000); // SBT2
  TH1 *htdc_sbt = new TH1F("htdc_sbt","",1000,0,70000); // SBT_Average
  TH1 *htime3 = new TH1F("htime3","",2000,0,2000);
  TH1 *htime5 = new TH1F("htime5","",2000,0,2000);
  TH1 *htime7 = new TH1F("htime7","",2000,0,2000);
  TH1 *htime1 = new TH1F("htime1","",2000,0,2000);
  TH1 *htime2 = new TH1F("htime2","",2000,0,2000);
  TH1 *htime_sbt = new TH1F("htime_sbt","",2000,0,2000);
  htime7->GetXaxis()->SetTitle("TDC[ns]");
  htime1->GetXaxis()->SetTitle("TDC[ns]");
  htime2->GetXaxis()->SetTitle("TDC[ns]");
  htime_sbt->GetXaxis()->SetTitle("TDC[ns]");
  TH1 *htime7_hodf = new TH1F("htime7_hodf","",2000,0,2000);
  TH1 *htime1_hodf = new TH1F("htime1_hodf","",2000,0,2000);
  TH1 *htime2_hodf = new TH1F("htime2_hodf","",2000,0,2000);
  TH1 *htime_sbt_hodf = new TH1F("htime_sbt_hodf","",2000,0,2000);
  TH1 *htime7_cp = new TH1F("htime7_cp","",2000,0,2000);
  TH1 *htime1_cp = new TH1F("htime1_cp","",2000,0,2000);
  TH1 *htime2_cp = new TH1F("htime2_cp","",2000,0,2000);
  TH1 *htime_sbt_cp = new TH1F("htime_sbt_cp","",2000,0,2000);
  TH1 *htime7_hodp = new TH1F("htime7_hodp","",2000,0,2000);
  TH1 *htime1_hodp = new TH1F("htime1_hodp","",2000,0,2000);
  TH1 *htime2_hodp = new TH1F("htime2_hodp","",2000,0,2000);
  TH1 *htime_sbt_hodp = new TH1F("htime_sbt_hodp","",2000,0,2000);
  TH1 *htime7_n = new TH1F("htime7_n","",2000,0,2000);
  TH1 *htime1_n = new TH1F("htime1_n","",2000,0,2000);
  TH1 *htime2_n = new TH1F("htime2_n","",2000,0,2000);
  TH1 *htime_sbt_n = new TH1F("htime_sbt_n","",2000,0,2000);
  TH1 *htime7_np = new TH1F("htime7_np","",2000,0,2000);
  TH1 *htime1_np = new TH1F("htime1_np","",2000,0,2000);
  TH1 *htime2_np = new TH1F("htime2_np","",2000,0,2000);
  TH1 *htime_sbt_np = new TH1F("htime_sbt_np","",2000,0,2000);
  TH1 *htime1_temp = new TH1F("htime1_temp","",10000,0,1000);
  TH1 *htime2_temp = new TH1F("htime2_temp","",10000,0,1000);
 
  TH1 *hqdc3l = new TH1F("hqdc3l","",1000,0,3000);
  TH1 *hqdc5l = new TH1F("hqdc5l","",1000,0,3000);
  TH1 *hqdc7l = new TH1F("hqdc7l","",1000,0,3000);
  TH1 *hqdc1l = new TH1F("hqdc1l","",1000,0,3000);
  TH1 *hqdc2l = new TH1F("hqdc2l","",1000,0,3000);
  TH1 *hq3 = new TH1F("hq3","",3000,0,3000);
  TH1 *hq5 = new TH1F("hq5","",3000,0,3000);
  TH1 *hq7 = new TH1F("hq7","",3000,0,3000);
  TH1 *hq1 = new TH1F("hq1","",3000,0,3000);
  TH1 *hq2 = new TH1F("hq2","",3000,0,3000);

  TH1 *hqt7 = new TH2F("hqt7","",1000,0,100,1500,0,1500);
  TH1 *hqt1 = new TH2F("hqt1","",1000,0,100,1500,0,1500);
  TH1 *hqt2 = new TH2F("hqt2","",1000,0,100,1500,0,1500);
  TH1 *hqt = new TH2F("hqt","",1000,0,100,1500,0,1500); // SBT1,2 average
  TH1 *hqt2_d = new TH2F("hqt2_d","",1000,0,100,1500,0,1500);
  TH1 *hqt_be = new TH2F("hqt_be","",1000,0,100,1500,0,1500);
   
  TH1 *htdc7lr = new TH2F("htdc7lr","",1000,0,30000,1000,0,30000); // SBT1 TDC L vs R
  TH1 *htdc1lr = new TH2F("htdc1lr","",1000,0,30000,1000,0,30000); // SBT1 TDC L vs R
  TH1 *htdc2lr = new TH2F("htdc2lr","",1000,0,30000,1000,0,30000); // SBT2 TDC L vs R
  TH1 *htime7lr = new TH2F("htime7lr","",1500,0,1500,1500,0,1500); // SBT1 time L vs R
  TH1 *htime1lr = new TH2F("htime1lr","",1500,0,1500,1500,0,1500); // SBT1 time L vs R
  TH1 *htime2lr = new TH2F("htime2lr","",1500,0,1500,1500,0,1500); // SBT2 time L vs R
  TH1 *hqdc7lr = new TH2F("hqdc7lr","",1000,0,5000,1000,0,5000); // SBT1 QDC L vs R
  TH1 *hqdc1lr = new TH2F("hqdc1lr","",1000,0,5000,1000,0,5000); // SBT1 QDC L vs R
  TH1 *hqdc2lr = new TH2F("hqdc2lr","",1000,0,5000,1000,0,5000); // SBT2 QDC L vs R
  
  TH1 *q1t12 = new TH2F("q1t12","",400,-10,10,1500,0,1500);
  TH1 *q2t12 = new TH2F("q2t12","",400,-10,10,1500,0,1500);
  TH1 *q7t71 = new TH2F("q7t71","",500,-370,-320,1500,0,1500);
  TH1 *q1t71 = new TH2F("q1t71","",500,-370,-320,1500,0,1500);
  TH1 *q7t72 = new TH2F("q7t72","",500,-370,-320,1500,0,1500);
  TH1 *q2t72 = new TH2F("q2t72","",500,-370,-320,1500,0,1500);
  TH1 *q7t7_sbt = new TH2F("q7t7_sbt","",500,-370,-320,1500,0,1500);
  q7t7_sbt->SetStats(0);
  q7t7_sbt->SetXTitle("TOF[ns]"); q7t7_sbt->SetYTitle("QDC[ch]");
  q7t7_sbt->GetYaxis()->SetRangeUser(0,1200);
  TH1 *q7t7_sbt_temp = new TH2F("q7t7_sbt_temp","",3000,-9999,20001,1500,0,1500);
  TH1 *q7t7_sbt_be = new TH2F("q7t7_sbt_be","",500,-370,-320,1500,0,1500);
  TH1 *qsbtt7_sbt = new TH2F("qsbtt7_sbt","",500,-370,-320,1500,0,1500);
 
  TH1 *corr_tdc71 = new TH2F("corr_tdc71","",1000,0,60000,1000,0,60000);  
  TH1 *corr_tdcl71 = new TH2F("corr_tdcl71","",1000,0,60000,1000,0,60000);  
  TH1 *corr_tdcr71 = new TH2F("corr_tdcr71","",1000,0,60000,1000,0,60000);  
  TH1 *corr_time71 = new TH2F("corr_time71","",1000,0,1000,1000,0,1000);  
 
  TH1 *q7q1 = new TH2F("q7q1","",1500,0,1500,1500,0,1500);
  TH1 *q7q2 = new TH2F("q7q2","",1500,0,1500,1500,0,1500);
  TH1 *q1q2 = new TH2F("q1q2","",1500,0,1500,1500,0,1500);
  TH1 *q7qsbt = new TH2F("q7qsbt","",1500,0,1500,1500,0,1500);

  double nentries = pla_tree->GetEntries();
  cout << nentries << endl;

  for(Long64_t ieve=0;ieve<nentries;++ieve){ // loop1
    //for(Long64_t ieve=0;ieve<10000;++ieve){ // loop1
    pla_tree->GetEntry(ieve);

    if (ieve%1000==0){
      cout<<"\r  Running...     "<<ieve<<" Events"<<flush;
    }

    for(int i=0;i<6;i++){ 
      if(ftriggerbit[i]==1) htrigger->Fill(i+1);
      if(Bitd && ftriggerbit[i]==1) htrigger_d->Fill(i+1);
      if(Bit10Be && ftriggerbit[i]==1) htrigger_be->Fill(i+1);
    }


    TArtPlastic *f3 = (TArtPlastic*)TArtUtil::FindDataObject(Plastic,(char*)"F3pl");
    TArtPlastic *f5 = (TArtPlastic*)TArtUtil::FindDataObject(Plastic,(char*)"F5pl");
    TArtPlastic *f7 = (TArtPlastic*)TArtUtil::FindDataObject(Plastic,(char*)"F7pl");
    TArtPlastic *sbt1 = (TArtPlastic*)TArtUtil::FindDataObject(Plastic,(char*)"F13pl-1");
    TArtPlastic *sbt2 = (TArtPlastic*)TArtUtil::FindDataObject(Plastic,(char*)"F13pl-2");

    double tdc3l = f3->GetTLRaw(), tdc3r = f3->GetTRRaw();
    double time3l = f3->GetTimeL(), time3r = f3->GetTimeR(), time3 = f3->GetTime();
    double qraw3l = f3->GetQLRaw(), qraw3r = f3->GetQRRaw(), qraw3 = f3->GetQAveRaw();
    // QAveRaw = sqrt(QL*QR)
    double tdc5l = f5->GetTLRaw(), tdc5r = f5->GetTRRaw();
    double time5l = f5->GetTimeL(), time5r = f5->GetTimeR(), time5 = f5->GetTime();
    double qraw5l = f5->GetQLRaw(), qraw5r = f5->GetQRRaw(), qraw5 = f5->GetQAveRaw();
    double tdc7l = f7->GetTLRaw(), tdc7r = f7->GetTRRaw();
    double time7l = f7->GetTimeL(), time7r = f7->GetTimeR(), time7 = f7->GetTime();
    double qraw7l = f7->GetQLRaw(), qraw7r = f7->GetQRRaw(), qraw7 = f7->GetQAveRaw();
    double tdc1l = sbt1->GetTLRaw(), tdc1r = sbt1->GetTRRaw();
    double time1l = sbt1->GetTimeL(), time1r = sbt1->GetTimeR(), time1 = sbt1->GetTime();
    double qraw1l = sbt1->GetQLRaw(), qraw1r = sbt1->GetQRRaw(), qraw1 = sbt1->GetQAveRaw();
    double tdc2l = sbt2->GetTLRaw(), tdc2r = sbt2->GetTRRaw();
    double time2l = sbt2->GetTimeL(), time2r = sbt2->GetTimeR(), time2 = sbt2->GetTime();
    double qraw2l = sbt2->GetQLRaw(), qraw2r = sbt2->GetQRRaw(), qraw2 = sbt2->GetQAveRaw();

    double tdc3 = -9999., tdc5 = -9999., tdc7 = -9999., tdc1 = -9999., tdc2 = -9999.;
    double tdc_sbt = -9999., time_sbt = -9999., qraw_sbt = -9999.;
    if(tdc3l>0 && tdc3r>0) tdc3 = (tdc3l+tdc3r)/2.; 
    if(tdc5l>0 && tdc5r>0) tdc5 = (tdc5l+tdc5r)/2.; 
    if(tdc7l>0 && tdc7r>0) tdc7 = (tdc7l+tdc7r)/2.; 
    if(tdc1l>0 && tdc1r>0) tdc1 = (tdc1l+tdc1r)/2.; 
    if(tdc2l>0 && tdc2r>0) tdc2 = (tdc2l+tdc2r)/2.; 
    if(tdc1>0 && tdc2>0) tdc_sbt = (tdc1+tdc2)/2.;
    if(time1>0 && time2>0) time_sbt = (time1+time2)/2.;
    if(qraw1>0 && qraw2>0) qraw_sbt = sqrt(qraw1*qraw2);

    if(ftriggerbit[0]==1){ // DSB trigger 2024/10/01   

      double tof12 =-9999.;
      if(time1>0 && time2>0) tof12 = time2 - time1;
      q1t12->Fill(tof12, qraw1);
      q2t12->Fill(tof12, qraw2);
      
      double tof7_1 =-9999.;
      if(time7>0 && time1>0) tof7_1 = time1 - time7;
      q7t71->Fill(tof7_1, qraw7);
      q1t71->Fill(tof7_1, qraw1);
      
      double tof7_2 =-9999.;
      if(time7>0 && time2>0) tof7_2 = time2 - time7;
      q7t72->Fill(tof7_2, qraw7);
      q2t72->Fill(tof7_2, qraw2);
      
      double tof7_sbt =-9999.;
      if(time7>0 && time_sbt>0) tof7_sbt = time_sbt - time7;
      q7t7_sbt->Fill(tof7_sbt, qraw7);
      q7t7_sbt_temp->Fill(tof7_sbt, qraw7);
      if(Bit10Be) q7t7_sbt_be->Fill(tof7_sbt, qraw7);
      qsbtt7_sbt->Fill(tof7_sbt, qraw_sbt);
      
      // left tdc[ch]
      htdc3l->Fill(tdc3l); htdc5l->Fill(tdc5l); htdc7l->Fill(tdc7l); htdc1l->Fill(tdc1l); htdc2l->Fill(tdc2l); 
      // tdc ave [ch]
      htdc3->Fill(tdc3); htdc5->Fill(tdc5); htdc7->Fill(tdc7); htdc1->Fill(tdc1); htdc2->Fill(tdc2); htdc_sbt->Fill(tdc_sbt);
      // tdc ave [ns]
      htime3->Fill(time3); htime5->Fill(time5); htime7->Fill(time7);
      htime1->Fill(time1);
      if((time1l>80 && time1r>80) || (time1l<80 && time1r<80)) htime1_temp->Fill(time1);
      htime2->Fill(time2);
      if((time2l>80 && time2r>80) || (time2l<80 && time2r<80)) htime2_temp->Fill(time2);
      htime_sbt->Fill(time_sbt);
      // left qdc [ch]
      hqdc3l->Fill(qraw3l); hqdc5l->Fill(qraw5l); hqdc7l->Fill(qraw7l); hqdc1l->Fill(qraw1l); hqdc2l->Fill(qraw2l);
      // qdc ave [ch]
      hq3->Fill(qraw3); hq5->Fill(qraw5); hq7->Fill(qraw7); hq1->Fill(qraw1); hq2->Fill(qraw2);
      // qdc[ch] vs tdc[ns]
      hqt7->Fill(time7, qraw7); hqt1->Fill(time1, qraw1); hqt2->Fill(time2, qraw2); hqt->Fill(time_sbt, qraw_sbt);
      
      if(Bitd) hqt2_d->Fill(time2, qraw2);
      if(Bit10Be) hqt_be->Fill(time_sbt, qraw_sbt);
     
      if(tdc7l>0 && tdc7r>0) htdc7lr->Fill(tdc7l,tdc7r);
      if(tdc1l>0 && tdc1r>0) htdc1lr->Fill(tdc1l,tdc1r);
      if(tdc2l>0 && tdc2r>0) htdc2lr->Fill(tdc2l,tdc2r);
      if(time7l>0 && time7r>0) htime7lr->Fill(time7l,time7r);
      if(time1l>0 && time1r>0) htime1lr->Fill(time1l,time1r);
      if(time2l>0 && time2r>0) htime2lr->Fill(time2l,time2r);
      if(qraw7l>0 && qraw7r>0) hqdc7lr->Fill(qraw7l,qraw7r);
      if(qraw1l>0 && qraw1r>0) hqdc1lr->Fill(qraw1l,qraw1r);
      if(qraw2l>0 && qraw2r>0) hqdc2lr->Fill(qraw2l,qraw2r);

      if(tdc7>0 && tdc1>0) corr_tdc71->Fill(tdc7,tdc1);
      if(tdc7l>0 && tdc1l>0) corr_tdcl71->Fill(tdc7l,tdc1l);
      if(tdc7r>0 && tdc1r>0) corr_tdcr71->Fill(tdc7r,tdc1r);
      if(time7>0 && time1>0) corr_time71->Fill(time7,time1);
      
      q7q1->Fill(qraw7,qraw1);
      q7q2->Fill(qraw7,qraw2);
      q1q2->Fill(qraw1,qraw2);
      q7qsbt->Fill(qraw7,qraw_sbt);
 
    }
    
    if(ftriggerbit[1]==1){ // B*Cg*HODF trigger 2024/10/15 
      htime7_hodf->Fill(time7);
      htime1_hodf->Fill(time1);
      htime2_hodf->Fill(time2);
      htime_sbt_hodf->Fill(time_sbt);
    }
    if(ftriggerbit[2]==1){ // B*CpLR trigger 2024/10/15 
      htime7_cp->Fill(time7);
      htime1_cp->Fill(time1);
      htime2_cp->Fill(time2);
      htime_sbt_cp->Fill(time_sbt);
    }
    if(ftriggerbit[3]==1){ // B*HODP trigger 2024/10/15 
      htime7_hodp->Fill(time7);
      htime1_hodp->Fill(time1);
      htime2_hodp->Fill(time2);
      htime_sbt_hodp->Fill(time_sbt);
    }
    if(ftriggerbit[4]==1){ // B*N trigger 2024/10/15 
      htime7_n->Fill(time7);
      htime1_n->Fill(time1);
      htime2_n->Fill(time2);
      htime_sbt_n->Fill(time_sbt);
    }
    if(ftriggerbit[5]==1){ // B*NP trigger 2024/10/15 
      htime7_np->Fill(time7);
      htime1_np->Fill(time1);
      htime2_np->Fill(time2);
      htime_sbt_np->Fill(time_sbt);
    }

  }// end of loop1
  cout << endl;
  fout->Write();
  //fout->Close();

  int bin_min = 45, bin_max = 80;
  double sum1_dsb = 0, sum1_hodf = 0, sum1_cp = 0, sum1_hodp = 0, sum1_n = 0, sum1_np = 0;
  double sum2_dsb = 0, sum2_hodf = 0, sum2_cp = 0, sum2_hodp = 0, sum2_n = 0, sum2_np = 0;
  double sum_sbt_dsb = 0, sum_sbt_hodf = 0, sum_sbt_cp = 0, sum_sbt_hodp = 0, sum_sbt_n = 0, sum_sbt_np = 0;
  for (int i = bin_min; i <= bin_max; i++) {
    sum1_dsb += htime1->GetBinContent(i);
    sum1_hodf += htime1_hodf->GetBinContent(i);
    sum1_cp += htime1_cp->GetBinContent(i);
    sum1_hodp += htime1_hodp->GetBinContent(i);
    sum1_n += htime1_n->GetBinContent(i);
    sum1_np += htime1_np->GetBinContent(i);
    
    sum2_dsb += htime2->GetBinContent(i);
    sum2_hodf += htime2_hodf->GetBinContent(i);
    sum2_cp += htime2_cp->GetBinContent(i);
    sum2_hodp += htime2_hodp->GetBinContent(i);
    sum2_n += htime2_n->GetBinContent(i);
    sum2_np += htime2_np->GetBinContent(i);
    
    sum_sbt_dsb += htime_sbt->GetBinContent(i);
    sum_sbt_hodf += htime_sbt_hodf->GetBinContent(i);
    sum_sbt_cp += htime_sbt_cp->GetBinContent(i);
    sum_sbt_hodp += htime_sbt_hodp->GetBinContent(i);
    sum_sbt_n += htime_sbt_n->GetBinContent(i);
    sum_sbt_np += htime_sbt_np->GetBinContent(i);
  } 
  cout << "sum1_dsb= "  << sum1_dsb << ", " << sum1_dsb/htime1->GetEntries() << endl;
  cout << "sum1_hodf= "  << sum1_hodf << ", " << sum1_hodf/htime1_hodf->GetEntries() << endl;
  cout << "sum1_cp= "  << sum1_cp << ", " << sum1_cp/htime1_cp->GetEntries() << endl;
  cout << "sum1_hodp= "  << sum1_hodp << ", " << sum1_hodp/htime1_hodp->GetEntries() << endl;
  cout << "sum1_n= "  << sum1_n << ", " << sum1_n/htime1_n->GetEntries() << endl;
  cout << "sum1_np= "  << sum1_np << ", " << sum1_np/htime1_np->GetEntries() << endl;
  
  cout << "sum2_dsb= "  << sum2_dsb << ", " << sum2_dsb/htime2->GetEntries() << endl;
  cout << "sum2_hodf= "  << sum2_hodf << ", " << sum2_hodf/htime2_hodf->GetEntries() << endl;
  cout << "sum2_cp= "  << sum2_cp << ", " << sum2_cp/htime2_cp->GetEntries() << endl;
  cout << "sum2_hodp= "  << sum2_hodp << ", " << sum2_hodp/htime2_hodp->GetEntries() << endl;
  cout << "sum2_n= "  << sum2_n << ", " << sum2_n/htime2_n->GetEntries() << endl;
  cout << "sum2_np= "  << sum2_np << ", " << sum2_np/htime2_np->GetEntries() << endl;

  cout << "sum_sbt_dsb= "  << sum_sbt_dsb << ", " << sum_sbt_dsb/htime_sbt->GetEntries() << endl;
  cout << "sum_sbt_hodf= "  << sum_sbt_hodf << ", " << sum_sbt_hodf/htime_sbt_hodf->GetEntries() << endl;
  cout << "sum_sbt_cp= "  << sum_sbt_cp << ", " << sum_sbt_cp/htime_sbt_cp->GetEntries() << endl;
  cout << "sum_sbt_hodp= "  << sum_sbt_hodp << ", " << sum_sbt_hodp/htime_sbt_hodp->GetEntries() << endl;
  cout << "sum_sbt_n= "  << sum_sbt_n << ", " << sum_sbt_n/htime_sbt_n->GetEntries() << endl;
  cout << "sum_sbt_np= "  << sum_sbt_np << ", " << sum_sbt_np/htime_sbt_np->GetEntries() << endl;

}//end of ana_pla__________________________________________________


