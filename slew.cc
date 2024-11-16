#include <iostream>
#include "TArtPlastic.hh"
#include "TArtPlasticPara.hh"
#include "TClonesArray.h"
#include "TArtBigRIPSParameters.hh"
#include "TArtCalibPlastic.hh"
#include "TArtCalibCoin.hh"
#include "TArtUtil.hh"
#include "TArtStoreManager.hh"
#include "TArtEventStore.hh"
#include "TArtEventInfo.hh"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

void slew(Int_t nRun=372){// run293: deuteron

  // import 
  TFile *fpla = new TFile(Form("macros/chkpla/rootfile/tree_pla%04d.root",nRun)); 
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
  TFile *fout = new TFile(Form("macros/chkpla/slew/slew%04d.root",nRun),"RECREATE");
  
  // hist 
  TH1 *htrigger = new TH1F("htrigger","r",6,0.5,6.5);
  TH1 *q1lt12 = new TH2F("q1lt12","",600,200,800,200,-5,5);
  TH1 *q2lt12 = new TH2F("q2lt12","",600,200,800,200,-5,5);
  TH1 *q1lt12_be = new TH2F("q1lt12_be","",600,200,800,200,-5,5);
  TH1 *q2lt12_be = new TH2F("q2lt12_be","",600,200,800,200,-5,5);
  TH1 *q1lt12_be_dsb = new TH2F("q1lt12_be_dsb","",600,200,800,200,-5,5);
  TH1 *q2lt12_be_dsb = new TH2F("q2lt12_be_dsb","",600,200,800,200,-5,5);
  q1lt12->GetXaxis()->SetTitle("QDC[ch]");
  q1lt12->GetYaxis()->SetTitle("TOF[ns]");
  q2lt12->GetXaxis()->SetTitle("QDC[ch]");
  q2lt12->GetYaxis()->SetTitle("TOF[ns]");
  q1lt12_be->GetXaxis()->SetTitle("QDC[ch]");
  q1lt12_be->GetYaxis()->SetTitle("TOF[ns]");
  q2lt12_be->GetXaxis()->SetTitle("QDC[ch]");
  q2lt12_be->GetYaxis()->SetTitle("TOF[ns]");
  q1lt12_be_dsb->GetXaxis()->SetTitle("QDC[ch]");
  q1lt12_be_dsb->GetYaxis()->SetTitle("TOF[ns]");
  q2lt12_be_dsb->GetXaxis()->SetTitle("QDC[ch]");
  q2lt12_be_dsb->GetYaxis()->SetTitle("TOF[ns]");

  // after slew
  TH1 *q1lt12_slew = new TH2F("q1lt12_slew","",600,200,800,200,-5,5);
  TH1 *q2lt12_slew = new TH2F("q2lt12_slew","",600,200,800,200,-5,5);
  q1lt12_slew->GetXaxis()->SetTitle("QDC[ch]");
  q1lt12_slew->GetYaxis()->SetTitle("TOF[ns]");
  q2lt12_slew->GetXaxis()->SetTitle("QDC[ch]");
  q2lt12_slew->GetYaxis()->SetTitle("TOF[ns]");
 
  // slew correction from
  ifstream ifs(Form("macros/chkpla/slew/slew%04d.txt",nRun));
  double pslew1a = 0,pslew1b = 0, pslew2a = 0, pslew2b = 0;
  ifs >> pslew1a >> pslew1b;
  ifs >> pslew2a >> pslew2b;
  cout << pslew1a << ", " << pslew1b << endl; 
  cout << pslew2a << ", " << pslew2b << endl; 

  //the number of entries
  double nentries = pla_tree->GetEntries();
  cout<< nentries << endl;
  
  //start of loop1 -------------------------
  for (Long64_t ieve=0;ieve<nentries;++ieve){
    pla_tree->GetEntry(ieve);
     
    if (ieve%1000==0){
      cout<<"\r  Running...     "<<ieve<<" Events"<<flush;
    }
    
    for(int i=0;i<6;i++){ 
      if(ftriggerbit[i]==1){ 
        htrigger->Fill(i+1);
      }
    }
    
    TArtPlastic *f3 = (TArtPlastic*)TArtUtil::FindDataObject(Plastic,(char*)"F3pl");
    TArtPlastic *f5 = (TArtPlastic*)TArtUtil::FindDataObject(Plastic,(char*)"F5pl");
    TArtPlastic *f7 = (TArtPlastic*)TArtUtil::FindDataObject(Plastic,(char*)"F7pl");
    TArtPlastic *sbt1 = (TArtPlastic*)TArtUtil::FindDataObject(Plastic,(char*)"F13pl-1");
    TArtPlastic *sbt2 = (TArtPlastic*)TArtUtil::FindDataObject(Plastic,(char*)"F13pl-2");

    double time3l = f3->GetTimeL(), time3r = f3->GetTimeR(), time3 = f3->GetTime();
    double qraw3l = f3->GetQLRaw(), qraw3r = f3->GetQRRaw(), qraw3 = f3->GetQAveRaw();
    double time5l = f5->GetTimeL(), time5r = f5->GetTimeR(), time5 = f5->GetTime();
    double qraw5l = f5->GetQLRaw(), qraw5r = f5->GetQRRaw(), qraw5 = f5->GetQAveRaw();
    double time7l = f7->GetTimeL(), time7r = f7->GetTimeR(), time7 = f7->GetTime();
    double qraw7l = f7->GetQLRaw(), qraw7r = f7->GetQRRaw(), qraw7 = f7->GetQAveRaw();
    double time1l = sbt1->GetTimeL(), time1r = sbt1->GetTimeR(), time1 = sbt1->GetTime();
    double qraw1l = sbt1->GetQLRaw(), qraw1r = sbt1->GetQRRaw(), qraw1 = sbt1->GetQAveRaw();
    double time2l = sbt2->GetTimeL(), time2r = sbt2->GetTimeR(), time2 = sbt2->GetTime();
    double qraw2l = sbt2->GetQLRaw(), qraw2r = sbt2->GetQRRaw(), qraw2 = sbt2->GetQAveRaw();

    double time_sbt = -9999., qraw_sbt = -9999.;
    if(time1>0 && time2>0) time_sbt = (time1+time2)/2.;
    if(qraw1>0 && qraw2>0) qraw_sbt = sqrt(qraw1*qraw2);
    
    double time1l_slew = -9999., time1r_slew = -9999., time1_slew = -9999.;
    if(time1l>0 && qraw1l>0 && time1r>0 && qraw1r>0) {
      time1l_slew = time1l - pslew1a/sqrt(qraw1l) - pslew1b;
      time1r_slew = time1r - pslew1a/sqrt(qraw1r) - pslew1b;
      time1_slew = (time1l_slew + time1r_slew)/2.;
    }
    double time2l_slew = -9999., time2r_slew = -9999., time2_slew = -9999.;
    if(time2l>0 && qraw2l>0 && time2r>0 && qraw2r>0) {
      time2l_slew = time2l - pslew2a/sqrt(qraw2l) - pslew2b;
      time2r_slew = time2r - pslew2a/sqrt(qraw2r) - pslew2b;
      time2_slew = (time2l_slew + time2r_slew)/2.;
    }

    double tof12 =-9999.;
    double tof12_slew =-9999.;
    if(time1>0 && time2>0) {
      tof12 = time2 - time1;
      //tof12_slew = tof12 - pslew2/sqrt(qraw2l);
    }
    if(time1_slew != -9999. && time2_slew != -9999.) tof12_slew = time2_slew - time1_slew;
    
    q1lt12->Fill(qraw1l, tof12);
    q2lt12->Fill(qraw2l, tof12);
    if(Bit10Be) {
      q1lt12_be->Fill(qraw1l, tof12);
      q2lt12_be->Fill(qraw2l, tof12);
      if(ftriggerbit[0]==1) {
        q1lt12_be_dsb->Fill(qraw1l, tof12);
        q2lt12_be_dsb->Fill(qraw2l, tof12);
        q1lt12_slew->Fill(qraw1l, tof12_slew);
        q2lt12_slew->Fill(qraw2l, tof12_slew);
      }
    }
    //tof12_slew = time2_slew - time1_slew;
           
  }// end of loop1----------------------------
  cout << endl;
  fout->Write();
  //fout->Close();
     
  TF1 *fslew = new TF1("fslew","[0]/sqrt(x)+[1]");
  
  fslew->SetParameter(0, 10);
  fslew->SetParameter(1, 10);

  q1lt12_be_dsb->Fit("fslew","","",200,800);
  double slew1a = fslew->GetParameter(0); 
  double slew1b = fslew->GetParameter(1); 
  q2lt12_be_dsb->Fit("fslew","","",200,800);
  double slew2a = fslew->GetParameter(0); 
  double slew2b = fslew->GetParameter(1); 
  
  ofstream ofs(Form("macros/chkpla/slew/slew%04d.txt",nRun));
  ofs.seekp(0,ios_base::cur);
  ofs  <<  slew1a << " " << slew1b << endl;
  ofs.seekp(0,ios_base::cur);
  ofs  <<  slew2a << " " << slew2b << endl;
   
   
} 
