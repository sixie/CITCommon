//root -l CITCommon/Examples/ReadZeeEventExample.C+'("/data/blue/sixie/LeptonScaleAndResolution/Electrons/ZeeNtuple.AllNtuple_HZZ4lNtuple_s12-zllm50-2-v9_noskim_0000.root.root")'
//root -l CITCommon/Examples/ReadZeeEventExample.C+'("/data/blue/sixie/LeptonScaleAndResolution/Electrons/Hgg/tree_zee_nosmear.root",2)'

//================================================================================================
//
// Example
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TH1F.h>                   // 
#include <TCanvas.h>                // 
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // class for 4-vector calculations
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

// structure for output ntuple
#include "CITCommon/CommonData/interface/ZeeEventTree.h"
#endif


//=== MAIN MACRO ================================================================================================= 

void ReadZeeEventExample(const string inputfile, Int_t Option = 0) {


  //*****************************************************************************************
  //Make a  Histogram
  //*****************************************************************************************
  TH1F *hist = new TH1F ("mass",";Mass [GeV/c^{2}];Number of Events",80,40,200);


  //*****************************************************************************************
  // Input
  //*****************************************************************************************
  citana::ZeeEventTree zeeTree;
  if (Option == 0) {
    cout << zeeTree.fTreeType << endl;
    zeeTree.LoadTree(inputfile.c_str());
    zeeTree.InitTree();
    cout << zeeTree.fTreeType << endl;
  }
  if (Option == 1) {
    zeeTree.LoadTree(inputfile.c_str(),citana::ZeeEventTree::kHggFutyanZeeDataEvent);
    zeeTree.InitTree(citana::ZeeEventTree::kHggFutyanZeeDataEvent);
  }
  if (Option == 2) {
    zeeTree.LoadTree(inputfile.c_str(),citana::ZeeEventTree::kHggFutyanZeeMCEvent);
    zeeTree.InitTree(citana::ZeeEventTree::kHggFutyanZeeMCEvent);
  }

  //*****************************************************************************************
  // Read Input File
  //*****************************************************************************************
  for (Long64_t ievt=0; ievt<zeeTree.tree_->GetEntries();ievt++) {
    
    if (ievt%10000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
    
    zeeTree.tree_->GetEntry(ievt);

    if (zeeTree.Ele1Pt() > 25 && zeeTree.Ele2Pt() > 25) {

      //For supercluster energy
      TLorentzVector ele1;
      TLorentzVector ele2;

      Bool_t useRegression = kFALSE;
      //Note: Et = E / Cosh(eta)

      if (useRegression) {
        ele1.SetPtEtaPhiM( zeeTree.fEle1EnergyRegression / TMath::CosH(zeeTree.Ele1Eta()) , zeeTree.Ele1Eta(), zeeTree.Ele1Phi(), 0.51099892e-3);
        ele2.SetPtEtaPhiM( zeeTree.fEle2EnergyRegression / TMath::CosH(zeeTree.Ele2Eta()) , zeeTree.Ele2Eta(), zeeTree.Ele2Phi(), 0.51099892e-3);
      } else {
        ele1.SetPtEtaPhiM( zeeTree.Ele1Pt() , zeeTree.Ele1Eta(), zeeTree.Ele1Phi(), 0.51099892e-3);
        ele2.SetPtEtaPhiM( zeeTree.Ele2Pt() , zeeTree.Ele2Eta(), zeeTree.Ele2Phi(), 0.51099892e-3);
      }
      
      double weight = 1; //zeeTree.fWeight  //if you want to use pileup reweighting      
      hist->Fill((ele1+ele2).M(), weight);

    }
   
  }
  

     
//--------------------------------------------------------------------------------------------------------------
// Plot
//==============================================================================================================
   
  TCanvas *cv = new TCanvas("cv","cv", 800, 600);
  hist->Draw("hist");
  cv->SaveAs("massExample.gif");

      
}
