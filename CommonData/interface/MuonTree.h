#ifndef MuonTree_H
#define MuonTree_H

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include <cmath>
#include "assert.h"

namespace citana
{
  class MuonTree {

    public:

      /// variables
      Float_t                 fWeight;
      UInt_t                  fRunNumber;
      UInt_t                  fLumiSectionNumber;
      UInt_t                  fEventNumber;
      Bool_t                  fMuEventNumberParity;
      Float_t                 fRho; 
      UInt_t                  fNVertices; 

      Float_t                 fMuPt; 
      Float_t                 fMuEta; 
      Float_t                 fMuPhi; 

      UInt_t  	            fMuTypeBits;
      Bool_t  		    fIsAllArbitrated;
      Float_t                 fMuTkNchi2; 
      Float_t                 fMuGlobalNchi2; 
      Float_t                 fMuNValidHits; 
      Float_t                 fMuNTrackerHits; 
      Float_t                 fMuNPixelHits; 
      Float_t                 fMuNMatches; 
      Float_t                 fMuD0; 
    
      //Additional Vars used in Likelihood
      Float_t                 fMuIP3d; 
      Float_t                 fMuIP3dSig; 
      Float_t                 fMuTrkKink; 
      Float_t                 fMuGlobalKink; 
      Float_t                 fMuSegmentCompatibility; 
      Float_t                 fMuCaloCompatibility; 
      Float_t                 fMuHadEnergy; 
      Float_t                 fMuHoEnergy; 
      Float_t                 fMuEmEnergy; 
      Float_t                 fMuHadS9Energy; 
      Float_t                 fMuHoS9Energy; 
      Float_t                 fMuEmS9Energy; 

      //Isolation Variables
      Float_t                 fMuTrkIso03; 
      Float_t                 fMuEMIso03; 
      Float_t                 fMuHadIso03; 
      Float_t                 fMuTrkIso05; 
      Float_t                 fMuEMIso05; 
      Float_t                 fMuHadIso05; 
      Float_t                 fMuPFIso04; 
      Float_t                 fChargedIso_DR0p0To0p1;
      Float_t                 fChargedIso_DR0p1To0p2;
      Float_t                 fChargedIso_DR0p2To0p3;
      Float_t                 fChargedIso_DR0p3To0p4;
      Float_t                 fChargedIso_DR0p4To0p5;
      Float_t                 fGammaIso_DR0p0To0p1;
      Float_t                 fGammaIso_DR0p1To0p2;
      Float_t                 fGammaIso_DR0p2To0p3;
      Float_t                 fGammaIso_DR0p3To0p4;
      Float_t                 fGammaIso_DR0p4To0p5;
      Float_t                 fNeutralHadronIso_DR0p0To0p1;
      Float_t                 fNeutralHadronIso_DR0p1To0p2;
      Float_t                 fNeutralHadronIso_DR0p2To0p3;
      Float_t                 fNeutralHadronIso_DR0p3To0p4;
      Float_t                 fNeutralHadronIso_DR0p4To0p5;

      Bool_t                  fMuPassTriggerDenominator;

    public:
      /// this is the main element
      TTree *tree_;
      TFile *f_;
  
      /// hold the names of variables to facilitate things (filled during Init)
      std::vector<std::string> variables_;

      /// default constructor  
      MuonTree()  {};
      /// default destructor
      ~MuonTree(){ 
        if (f_) f_->Close();  
      };
    
      /// initialize varibles and fill list of available variables
      void InitVariables() {
        fWeight			   = 0.0;
        fRunNumber		   = 0.0;
        fLumiSectionNumber	   = 0.0;
        fEventNumber		   = 0.0;
        fMuEventNumberParity	   = 0.0;
        fRho 			   = 0.0;
        fNVertices 		   = 0.0;
        fMuPt  			   = 0.0;
        fMuEta 			   = 0.0;
        fMuPhi 			   = 0.0;
        fMuTypeBits		   = 0.0;
        fIsAllArbitrated		   = 0.0;
        fMuTkNchi2 		   = 0.0;
        fMuGlobalNchi2 		   = 0.0;
        fMuNValidHits  		   = 0.0;
        fMuNTrackerHits 		   = 0.0;
        fMuNPixelHits  		   = 0.0;
        fMuNMatches 		   = 0.0;
        fMuD0  			   = 0.0;
        fMuIP3d 			   = 0.0;
        fMuIP3dSig 		   = 0.0;
        fMuTrkKink 		   = 0.0;
        fMuGlobalKink  		   = 0.0;
        fMuSegmentCompatibility 	   = 0.0;
        fMuCaloCompatibility 	   = 0.0;
        fMuHadEnergy 		   = 0.0;
        fMuHoEnergy 		   = 0.0;
        fMuEmEnergy 		   = 0.0;
        fMuHadS9Energy 		   = 0.0;
        fMuHoS9Energy  		   = 0.0;
        fMuEmS9Energy  		   = 0.0;
        fMuTrkIso03 		   = 0.0;
        fMuEMIso03 		   = 0.0;
        fMuHadIso03 		   = 0.0;
        fMuTrkIso05 		   = 0.0;
        fMuEMIso05 		   = 0.0;
        fMuHadIso05 		   = 0.0;
        fMuPFIso04 		   = 0.0;
        fChargedIso_DR0p0To0p1 	   = 0.0;
        fChargedIso_DR0p1To0p2 	   = 0.0;
        fChargedIso_DR0p2To0p3 	   = 0.0;
        fChargedIso_DR0p3To0p4 	   = 0.0;
        fChargedIso_DR0p4To0p5 	   = 0.0;
        fGammaIso_DR0p0To0p1	   = 0.0;
        fGammaIso_DR0p1To0p2	   = 0.0;
        fGammaIso_DR0p2To0p3	   = 0.0;
        fGammaIso_DR0p3To0p4	   = 0.0;
        fGammaIso_DR0p4To0p5	   = 0.0;
        fNeutralHadronIso_DR0p0To0p1 = 0.0;
        fNeutralHadronIso_DR0p1To0p2 = 0.0;
        fNeutralHadronIso_DR0p2To0p3 = 0.0;
        fNeutralHadronIso_DR0p3To0p4 = 0.0;
        fNeutralHadronIso_DR0p4To0p5 = 0.0;
        fMuPassTriggerDenominator = kFALSE;
      }
    
      /// load a MuonTree
      void LoadTree(const char* file){
        f_ = TFile::Open(file);
        assert(f_);
        tree_ = dynamic_cast<TTree*>(f_->Get("Muons"));
        assert(tree_);
      }
    
      /// create a MuonTree
      void CreateTree(){
        tree_ = new TTree("Muons","Muons");
        f_ = 0;

        //book the branches
        tree_->Branch("weight",&fWeight,"weight/F");
        tree_->Branch("run",&fRunNumber,"run/i");
        tree_->Branch("lumi",&fLumiSectionNumber,"lumi/i");
        tree_->Branch("event",&fEventNumber,"event/i");
        tree_->Branch("EventNumberParity",&fMuEventNumberParity,"EventNumberParity/O"); 
        tree_->Branch("Rho",&fRho,"Rho/F"); 
        tree_->Branch("NVertices",&fNVertices,"NVertices/F"); 
        tree_->Branch("pt",&fMuPt,"pt/F"); 
        tree_->Branch("eta",&fMuEta,"eta/F"); 
        tree_->Branch("phi",&fMuPhi,"phi/F"); 
        tree_->Branch("typeBits",&fMuTypeBits,"typeBits/i"); 
        tree_->Branch("isAllArbitrated",&fIsAllArbitrated,"isAllArbitrated/i"); 
        tree_->Branch("TkNchi2",&fMuTkNchi2,"TkNchi2/F"); 
        tree_->Branch("GlobalNchi2",&fMuGlobalNchi2,"GlobalNchi2/F"); 
        tree_->Branch("NValidHits",&fMuNValidHits,"NValidHits/F"); 
        tree_->Branch("NTrackerHits",&fMuNTrackerHits,"NTrackerHits/F"); 
        tree_->Branch("NPixelHits",&fMuNPixelHits,"NPixelHits/F"); 
        tree_->Branch("NMatches",&fMuNMatches,"NMatches/F"); 
        tree_->Branch("D0",&fMuD0,"D0/F"); 
        tree_->Branch("IP3d",&fMuIP3d,"IP3d/F"); 
        tree_->Branch("IP3dSig",&fMuIP3dSig,"IP3dSig/F"); 
        tree_->Branch("TrkKink",&fMuTrkKink,"TrkKink/F"); 
        tree_->Branch("GlobalKink",&fMuGlobalKink,"GlobalKink/F"); 
        tree_->Branch("SegmentCompatibility",&fMuSegmentCompatibility,"SegmentCompatibility/F"); 
        tree_->Branch("CaloCompatibility",&fMuCaloCompatibility,"CaloCompatibility/F"); 
        tree_->Branch("HadEnergy",&fMuHadEnergy,"HadEnergy/F"); 
        tree_->Branch("HoEnergy",&fMuHoEnergy,"HoEnergy/F"); 
        tree_->Branch("EmEnergy",&fMuEmEnergy,"EmEnergy/F"); 
        tree_->Branch("HadS9Energy",&fMuHadS9Energy,"HadS9Energy/F"); 
        tree_->Branch("HoS9Energy",&fMuHoS9Energy,"HoS9Energy/F"); 
        tree_->Branch("EmS9Energy",&fMuEmS9Energy,"EmS9Energy/F"); 
        tree_->Branch("TrkIso03",&fMuTrkIso03,"TrkIso03/F"); 
        tree_->Branch("EMIso03",&fMuEMIso03,"EMIso03/F"); 
        tree_->Branch("HadIso03",&fMuHadIso03,"HadIso03/F"); 
        tree_->Branch("TrkIso05",&fMuTrkIso05,"TrkIso05/F"); 
        tree_->Branch("EMIso05",&fMuEMIso05,"EMIso05/F"); 
        tree_->Branch("HadIso05",&fMuHadIso05,"HadIso05/F"); 
        tree_->Branch("PFIso04",&fMuPFIso04,"PFIso04/F"); 
        tree_->Branch("ChargedIso_DR0p0To0p1",&fChargedIso_DR0p0To0p1,"ChargedIso_DR0p0To0p1/F");
        tree_->Branch("ChargedIso_DR0p1To0p2",&fChargedIso_DR0p1To0p2,"ChargedIso_DR0p1To0p2/F");
        tree_->Branch("ChargedIso_DR0p2To0p3",&fChargedIso_DR0p2To0p3,"ChargedIso_DR0p2To0p3/F");
        tree_->Branch("ChargedIso_DR0p3To0p4",&fChargedIso_DR0p3To0p4,"ChargedIso_DR0p3To0p4/F");
        tree_->Branch("ChargedIso_DR0p4To0p5",&fChargedIso_DR0p4To0p5,"ChargedIso_DR0p4To0p5/F");
        tree_->Branch("GammaIso_DR0p0To0p1",&fGammaIso_DR0p0To0p1,"GammaIso_DR0p0To0p1/F");
        tree_->Branch("GammaIso_DR0p1To0p2",&fGammaIso_DR0p1To0p2,"GammaIso_DR0p1To0p2/F");
        tree_->Branch("GammaIso_DR0p2To0p3",&fGammaIso_DR0p2To0p3,"GammaIso_DR0p2To0p3/F");
        tree_->Branch("GammaIso_DR0p3To0p4",&fGammaIso_DR0p3To0p4,"GammaIso_DR0p3To0p4/F");
        tree_->Branch("GammaIso_DR0p4To0p5",&fGammaIso_DR0p4To0p5,"GammaIso_DR0p4To0p5/F");
        tree_->Branch("NeutralHadronIso_DR0p0To0p1",&fNeutralHadronIso_DR0p0To0p1,"NeutralHadronIso_DR0p0To0p1/F");
        tree_->Branch("NeutralHadronIso_DR0p1To0p2",&fNeutralHadronIso_DR0p1To0p2,"NeutralHadronIso_DR0p1To0p2/F");
        tree_->Branch("NeutralHadronIso_DR0p2To0p3",&fNeutralHadronIso_DR0p2To0p3,"NeutralHadronIso_DR0p2To0p3/F");
        tree_->Branch("NeutralHadronIso_DR0p3To0p4",&fNeutralHadronIso_DR0p3To0p4,"NeutralHadronIso_DR0p3To0p4/F");
        tree_->Branch("NeutralHadronIso_DR0p4To0p5",&fNeutralHadronIso_DR0p4To0p5,"NeutralHadronIso_DR0p4To0p5/F");
        tree_->Branch("PassTriggerDenominator",&fMuPassTriggerDenominator,"PassTriggerDenominator/O");

      } 

      // initialze a MuonTree
      void InitTree(){
        assert(tree_);
        // don't forget to set pointers to zero before you set address
        // or you will fully appreciate that "ROOT sucks" :)
        InitVariables();
        //Set branch address
        Int_t currentState = gErrorIgnoreLevel;
        // gErrorIgnoreLevel = kError;
        gErrorIgnoreLevel = kBreak;

        tree_->SetBranchAddress("weight",&fWeight);
        tree_->SetBranchAddress("run",&fRunNumber);
        tree_->SetBranchAddress("lumi",&fLumiSectionNumber);
        tree_->SetBranchAddress("event",&fEventNumber);
        tree_->SetBranchAddress("EventNumberParity",&fMuEventNumberParity);
        tree_->SetBranchAddress("Rho",&fRho);
        tree_->SetBranchAddress("NVertices",&fNVertices);
        tree_->SetBranchAddress("pt",&fMuPt);
        tree_->SetBranchAddress("eta",&fMuEta);
        tree_->SetBranchAddress("phi",&fMuPhi);
        tree_->SetBranchAddress("typeBits",&fMuTypeBits);
        tree_->SetBranchAddress("isAllArbitrated",&fIsAllArbitrated);
        tree_->SetBranchAddress("TkNchi2",&fMuTkNchi2);
        tree_->SetBranchAddress("GlobalNchi2",&fMuGlobalNchi2);
        tree_->SetBranchAddress("NValidHits",&fMuNValidHits);
        tree_->SetBranchAddress("NTrackerHits",&fMuNTrackerHits);
        tree_->SetBranchAddress("NPixelHits",&fMuNPixelHits);
        tree_->SetBranchAddress("NMatches",&fMuNMatches);
        tree_->SetBranchAddress("D0",&fMuD0);
        tree_->SetBranchAddress("IP3d",&fMuIP3d);
        tree_->SetBranchAddress("IP3dSig",&fMuIP3dSig);
        tree_->SetBranchAddress("TrkKink",&fMuTrkKink);
        tree_->SetBranchAddress("GlobalKink",&fMuGlobalKink);
        tree_->SetBranchAddress("SegmentCompatibility",&fMuSegmentCompatibility);
        tree_->SetBranchAddress("CaloCompatibility",&fMuCaloCompatibility);
        tree_->SetBranchAddress("HadEnergy",&fMuHadEnergy);
        tree_->SetBranchAddress("HoEnergy",&fMuHoEnergy);
        tree_->SetBranchAddress("EmEnergy",&fMuEmEnergy);
        tree_->SetBranchAddress("HadS9Energy",&fMuHadS9Energy);
        tree_->SetBranchAddress("HoS9Energy",&fMuHoS9Energy);
        tree_->SetBranchAddress("EmS9Energy",&fMuEmS9Energy);
        tree_->SetBranchAddress("TrkIso03",&fMuTrkIso03);
        tree_->SetBranchAddress("EMIso03",&fMuEMIso03);
        tree_->SetBranchAddress("HadIso03",&fMuHadIso03);
        tree_->SetBranchAddress("TrkIso05",&fMuTrkIso05);
        tree_->SetBranchAddress("EMIso05",&fMuEMIso05);
        tree_->SetBranchAddress("HadIso05",&fMuHadIso05);
        tree_->SetBranchAddress("PFIso04",&fMuPFIso04);
        tree_->SetBranchAddress("ChargedIso_DR0p0To0p1",&fChargedIso_DR0p0To0p1);
        tree_->SetBranchAddress("ChargedIso_DR0p1To0p2",&fChargedIso_DR0p1To0p2);
        tree_->SetBranchAddress("ChargedIso_DR0p2To0p3",&fChargedIso_DR0p2To0p3);
        tree_->SetBranchAddress("ChargedIso_DR0p3To0p4",&fChargedIso_DR0p3To0p4);
        tree_->SetBranchAddress("ChargedIso_DR0p4To0p5",&fChargedIso_DR0p4To0p5);
        tree_->SetBranchAddress("GammaIso_DR0p0To0p1",&fGammaIso_DR0p0To0p1);
        tree_->SetBranchAddress("GammaIso_DR0p1To0p2",&fGammaIso_DR0p1To0p2);
        tree_->SetBranchAddress("GammaIso_DR0p2To0p3",&fGammaIso_DR0p2To0p3);
        tree_->SetBranchAddress("GammaIso_DR0p3To0p4",&fGammaIso_DR0p3To0p4);
        tree_->SetBranchAddress("GammaIso_DR0p4To0p5",&fGammaIso_DR0p4To0p5);
        tree_->SetBranchAddress("NeutralHadronIso_DR0p0To0p1",&fNeutralHadronIso_DR0p0To0p1);
        tree_->SetBranchAddress("NeutralHadronIso_DR0p1To0p2",&fNeutralHadronIso_DR0p1To0p2);
        tree_->SetBranchAddress("NeutralHadronIso_DR0p2To0p3",&fNeutralHadronIso_DR0p2To0p3);
        tree_->SetBranchAddress("NeutralHadronIso_DR0p3To0p4",&fNeutralHadronIso_DR0p3To0p4);
        tree_->SetBranchAddress("NeutralHadronIso_DR0p4To0p5",&fNeutralHadronIso_DR0p4To0p5);
        tree_->SetBranchAddress("PassTriggerDenominator",&fMuPassTriggerDenominator);

        gErrorIgnoreLevel = currentState;
      }

  }; 

}

#endif
