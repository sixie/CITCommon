#ifndef ZeeEventTree_H
#define ZeeEventTree_H

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include <cmath>
#include "assert.h"

namespace citana
{
  class ZeeEventTree {

    public:

      enum ZeeEventTreeType { kCITZeeEvent                                     = 0, 
                              kHggFutyanZeeDataEvent, 
                              kHggFutyanZeeMCEvent,
      };

      enum DiElecTrigBits { kEle17SC8Trigger                                   = 1UL<<0, 
                            kEle20SC4Trigger                                   = 1UL<<1, 
                            kEle32SC17Trigger                                  = 1UL<<2, 
                            kEle17Ele8Loose                                    = 1UL<<3,
                            kEle17Ele8Tight                                    = 1UL<<4,
                            kSingleEleTight                                    = 1UL<<5
      };
      
      /// variables
      UInt_t                  fTreeType;
      Float_t                 fWeight;
      UInt_t                  fRunNumber;
      UInt_t                  fLumiSectionNumber;
      UInt_t                  fEventNumber;
      UInt_t                  fNPU;
      Float_t                 fRho; 
      UInt_t                  fNVertices; 
      UInt_t                  fEventTriggerBits;

      Float_t                 fMass;
      Float_t                 fMassRegression;


      Float_t                 fEle1Pt; 
      Float_t                 fEle1Eta; 
      Float_t                 fEle1Phi; 
      Float_t                 fEle1SCEt; 
      Float_t                 fEle1SCEta; 
      Float_t                 fEle1SCPhi; 
      Float_t                 fEle1GenPt; 
      Float_t                 fEle1EnergyCorrAndSmeared; 
      Float_t                 fEle1Energy; 
      Float_t                 fEle1EnergyRegression; 
      Float_t                 fEle1EnergyRegressionV0; 
      Float_t                 fEle1EnergyRegressionV1; 
      Float_t                 fEle1EnergyRegressionV2;
      Float_t                 fEle1EnergyRegressionErrorV0; 
      Float_t                 fEle1EnergyRegressionErrorV1; 
      Float_t                 fEle1EnergyRegressionErrorV2;
      Int_t                   fEle1Charge;
      Float_t                 fEle1HZZICHEP2012IDMVA; 
      Float_t                 fEle1PFIso04;
      Float_t                 fEle1R9;
      Bool_t                  fEle1PassLooseSimpleCuts; 
      Bool_t                  fEle1PassMediumSimpleCuts; 
      Bool_t                  fEle1PassTightSimpleCuts; 
      Bool_t                  fEle1PassHZZICHEP2012; 


      Float_t                 fEle2Pt; 
      Float_t                 fEle2Eta; 
      Float_t                 fEle2Phi; 
      Float_t                 fEle2SCEt; 
      Float_t                 fEle2SCEta; 
      Float_t                 fEle2SCPhi; 
      Float_t                 fEle2GenPt; 
      Float_t                 fEle2EnergyCorrAndSmeared; 
      Float_t                 fEle2Energy; 
      Float_t                 fEle2EnergyRegression; 
      Float_t                 fEle2EnergyRegressionV0; 
      Float_t                 fEle2EnergyRegressionV1; 
      Float_t                 fEle2EnergyRegressionV2; 
      Float_t                 fEle2EnergyRegressionErrorV0; 
      Float_t                 fEle2EnergyRegressionErrorV1; 
      Float_t                 fEle2EnergyRegressionErrorV2;
      Int_t                   fEle2Charge;
      Float_t                 fEle2HZZICHEP2012IDMVA; 
      Float_t                 fEle2PFIso04;
      Float_t                 fEle2R9;
      Bool_t                  fEle2PassLooseSimpleCuts; 
      Bool_t                  fEle2PassMediumSimpleCuts; 
      Bool_t                  fEle2PassTightSimpleCuts; 
      Bool_t                  fEle2PassHZZICHEP2012; 


      Int_t                   fRunNumber_int;
      Int_t                   fLumiSectionNumber_int;
      Int_t                   fEventNumber_int;
      Int_t                   fNVertices_int; 
      Double_t                dEle1Pt;
      Double_t                dEle1Eta;
      Double_t                dEle1Phi;
      Double_t                dEle1EnergyCorrAndSmeared;
      Double_t                dEle2Pt;
      Double_t                dEle2Eta;
      Double_t                dEle2Phi;
      Double_t                dEle2EnergyCorrAndSmeared;



    public:
      /// this is the main element
      TTree *tree_;
      TFile *f_;
  
      /// hold the names of variables to facilitate things (filled during Init)
      std::vector<std::string> variables_;
      
      /// default constructor  
      ZeeEventTree()  {
        fTreeType = kCITZeeEvent;
      };
      /// default destructor
      ~ZeeEventTree(){ 
        if (f_) f_->Close();  
      };
    
      /// initialize varibles and fill list of available variables
      void InitVariables() {
        fWeight			       = 0.0;
        fRunNumber		       = 0.0;
        fLumiSectionNumber	       = 0.0;
        fEventNumber		       = 0.0;
        fNPU  			       = 0.0;
        fRho  			       = 0.0;
        fNVertices 		       = 0.0;
        fEventTriggerBits              = 0;
        fMass                            = 0.0;
        fEle1Pt                          = 0.0; 
        fEle1Eta                         = 0.0; 
        fEle1Phi                         = 0.0; 
        fEle1SCEt                        = 0.0; 
        fEle1SCEta                       = 0.0; 
        fEle1SCPhi                       = 0.0; 
        fEle1GenPt                       = 0.0; 
        fEle1EnergyCorrAndSmeared        = 0.0; 
        fEle1Energy                      = 0.0; 
        fEle1EnergyRegression            = 0.0; 
        fEle1EnergyRegressionV0          = 0.0; 
        fEle1EnergyRegressionV1          = 0.0; 
        fEle1EnergyRegressionV2          = 0.0; 
        fEle1EnergyRegressionErrorV0     = 0.0; 
        fEle1EnergyRegressionErrorV1     = 0.0; 
        fEle1EnergyRegressionErrorV2     = 0.0; 
        fEle1Charge                      = 0; 
        fEle1HZZICHEP2012IDMVA           = 0.0; 
        fEle1PFIso04                     = 0.0;
        fEle1R9                          = 0.0;
        fEle1PassLooseSimpleCuts         = kFALSE; 
        fEle1PassMediumSimpleCuts         = kFALSE; 
        fEle1PassTightSimpleCuts         = kFALSE; 
        fEle1PassHZZICHEP2012            = kFALSE; 
        fEle2Pt                          = 0.0; 
        fEle2Eta                         = 0.0; 
        fEle2Phi                         = 0.0; 
        fEle2SCEt                        = 0.0; 
        fEle2SCEta                       = 0.0; 
        fEle2SCPhi                       = 0.0; 
        fEle2GenPt                       = 0.0; 
        fEle2EnergyCorrAndSmeared        = 0.0; 
        fEle2Energy                      = 0.0; 
        fEle2EnergyRegression            = 0.0; 
        fEle2EnergyRegressionV0          = 0.0; 
        fEle2EnergyRegressionV1          = 0.0; 
        fEle2EnergyRegressionV2          = 0.0; 
        fEle2EnergyRegressionErrorV0     = 0.0; 
        fEle2EnergyRegressionErrorV1     = 0.0; 
        fEle2EnergyRegressionErrorV2     = 0.0; 
        fEle2Charge                      = 0; 
        fEle2HZZICHEP2012IDMVA           = 0.0; 
        fEle2PFIso04                     = 0.0;
        fEle2R9                          = 0.0;
        fEle2PassLooseSimpleCuts         = kFALSE; 
        fEle2PassMediumSimpleCuts        = kFALSE; 
        fEle2PassTightSimpleCuts         = kFALSE; 
        fEle2PassHZZICHEP2012            = kFALSE; 

        fRunNumber_int                   = 0;
        fLumiSectionNumber_int           = 0;
        fEventNumber_int                 = 0;
        fNVertices_int                   = 0; 
        dEle1Pt                          = 0;
        dEle1Eta                         = 0;
        dEle1Phi                         = 0;
        dEle1EnergyCorrAndSmeared        = 0;
        dEle2Pt                          = 0;
        dEle2Eta                         = 0;
        dEle2Phi                         = 0;
        dEle2EnergyCorrAndSmeared        = 0;

      }
    
      /// load a ZeeEventTree
      void LoadTree(const char* file, UInt_t Type = kCITZeeEvent){
        fTreeType = Type;
        f_ = TFile::Open(file);
        assert(f_);
        if (Type == kCITZeeEvent) {
          tree_ = dynamic_cast<TTree*>(f_->Get("ZeeEvent"));
        } else if (Type == kHggFutyanZeeDataEvent) {
          tree_ = dynamic_cast<TTree*>(f_->Get("Data"));
        } else if (Type == kHggFutyanZeeMCEvent) {
          tree_ = dynamic_cast<TTree*>(f_->Get("DYJetsToLL"));
        } else {
          cout << "Warning: Type " <<  Type << " is not supported \n";
        }
        InitTree(Type);
        assert(tree_);
      }
    
      /// create a ZeeEventTree
      void CreateTree(){
        tree_ = new TTree("ZeeEvent","ZeeEvent");
        f_ = 0;

        //book the branches
        tree_->Branch("weight",&fWeight,"weight/F");
        tree_->Branch("run",&fRunNumber,"run/i");
        tree_->Branch("lumi",&fLumiSectionNumber,"lumi/i");
        tree_->Branch("event",&fEventNumber,"event/i");
        tree_->Branch("npu",&fNPU,"npu/i"); 
        tree_->Branch("rho",&fRho,"rho/F"); 
        tree_->Branch("vertices",&fNVertices,"vertices/i"); 
        tree_->Branch("triggerbits",&fEventTriggerBits,"triggerbits/i"); 
        tree_->Branch("mass",&fMass,"mass/F"); 

        tree_->Branch("Ele1Pt",&fEle1Pt,"Ele1Pt/F"); 
        tree_->Branch("Ele1Eta",&fEle1Eta,"Ele1Eta/F"); 
        tree_->Branch("Ele1Phi",&fEle1Phi,"Ele1Phi/F"); 
        tree_->Branch("Ele1SCEt",&fEle1SCEt,"Ele1SCEt/F"); 
        tree_->Branch("Ele1SCEta",&fEle1SCEta,"Ele1SCEta/F"); 
        tree_->Branch("Ele1SCPhi",&fEle1SCPhi,"Ele1SCPhi/F"); 
        tree_->Branch("Ele1GenPt",&fEle1GenPt,"Ele1GenPt/F"); 
        tree_->Branch("Ele1EnergyCorrAndSmeared",&fEle1EnergyCorrAndSmeared,"Ele1EnergyCorrAndSmeared/F"); 
        tree_->Branch("Ele1Energy",&fEle1Energy,"Ele1Energy/F"); 
        tree_->Branch("Ele1EnergyRegression",&fEle1EnergyRegression,"Ele1EnergyRegression/F"); 
        tree_->Branch("Ele1EnergyRegressionV0",&fEle1EnergyRegressionV0,"Ele1EnergyRegressionV0/F"); 
        tree_->Branch("Ele1EnergyRegressionV1",&fEle1EnergyRegressionV1,"Ele1EnergyRegressionV1/F"); 
        tree_->Branch("Ele1EnergyRegressionV2",&fEle1EnergyRegressionV2,"Ele1EnergyRegressionV2/F");
        tree_->Branch("Ele1EnergyRegressionErrorV0",&fEle1EnergyRegressionErrorV0,"Ele1EnergyRegressionErrorV0/F"); 
        tree_->Branch("Ele1EnergyRegressionErrorV1",&fEle1EnergyRegressionErrorV1,"Ele1EnergyRegressionErrorV1/F"); 
        tree_->Branch("Ele1EnergyRegressionErrorV2",&fEle1EnergyRegressionErrorV2,"Ele1EnergyRegressionErrorV2/F");
        tree_->Branch("Ele1Charge",&fEle1Charge,"Ele1Charge/I"); 
        tree_->Branch("Ele1HZZICHEP2012IDMVA",&fEle1HZZICHEP2012IDMVA,"Ele1HZZICHEP2012IDMVA/F"); 
        tree_->Branch("Ele1PFIso04",&fEle1PFIso04,"Ele1PFIso04/F"); 
        tree_->Branch("Ele1R9",&fEle1R9,"Ele1R9/F"); 
        tree_->Branch("Ele1PassLooseSimpleCuts",&fEle1PassLooseSimpleCuts,"Ele1PassLooseSimpleCuts/O"); 
        tree_->Branch("Ele1PassMediumSimpleCuts",&fEle1PassMediumSimpleCuts,"Ele1PassMediumSimpleCuts/O"); 
        tree_->Branch("Ele1PassTightSimpleCuts",&fEle1PassTightSimpleCuts,"Ele1PassTightSimpleCuts/O"); 
        tree_->Branch("Ele1PassHZZICHEP2012",&fEle1PassHZZICHEP2012,"Ele1PassHZZICHEP2012/O");
        tree_->Branch("Ele2Pt",&fEle2Pt,"Ele2Pt/F"); 
        tree_->Branch("Ele2Eta",&fEle2Eta,"Ele2Eta/F"); 
        tree_->Branch("Ele2Phi",&fEle2Phi,"Ele2Phi/F"); 
        tree_->Branch("Ele2SCEt",&fEle2SCEt,"Ele2SCEt/F"); 
        tree_->Branch("Ele2SCEta",&fEle2SCEta,"Ele2SCEta/F"); 
        tree_->Branch("Ele2SCPhi",&fEle2SCPhi,"Ele2SCPhi/F"); 
        tree_->Branch("Ele2GenPt",&fEle2GenPt,"Ele2GenPt/F"); 
        tree_->Branch("Ele2EnergyCorrAndSmeared",&fEle2EnergyCorrAndSmeared,"Ele2EnergyCorrAndSmeared/F"); 
        tree_->Branch("Ele2Energy",&fEle2Energy,"Ele2Energy/F"); 
        tree_->Branch("Ele2EnergyRegression",&fEle2EnergyRegression,"Ele2EnergyRegression/F"); 
        tree_->Branch("Ele2EnergyRegressionV0",&fEle2EnergyRegressionV0,"Ele2EnergyRegressionV0/F"); 
        tree_->Branch("Ele2EnergyRegressionV1",&fEle2EnergyRegressionV1,"Ele2EnergyRegressionV1/F"); 
        tree_->Branch("Ele2EnergyRegressionV2",&fEle2EnergyRegressionV2,"Ele2EnergyRegressionV2/F");
        tree_->Branch("Ele2EnergyRegressionErrorV0",&fEle2EnergyRegressionErrorV0,"Ele2EnergyRegressionErrorV0/F"); 
        tree_->Branch("Ele2EnergyRegressionErrorV1",&fEle2EnergyRegressionErrorV1,"Ele2EnergyRegressionErrorV1/F"); 
        tree_->Branch("Ele2EnergyRegressionErrorV2",&fEle2EnergyRegressionErrorV2,"Ele2EnergyRegressionErrorV2/F");
        tree_->Branch("Ele2Charge",&fEle2Charge,"Ele2Charge/I"); 
        tree_->Branch("Ele2HZZICHEP2012IDMVA",&fEle2HZZICHEP2012IDMVA,"Ele2HZZICHEP2012IDMVA/F"); 
        tree_->Branch("Ele2PFIso04",&fEle2PFIso04,"Ele2PFIso04/F"); 
        tree_->Branch("Ele2R9",&fEle2R9,"Ele2R9/F"); 
        tree_->Branch("Ele2PassLooseSimpleCuts",&fEle2PassLooseSimpleCuts,"Ele2PassLooseSimpleCuts/O"); 
        tree_->Branch("Ele2PassMediumSimpleCuts",&fEle2PassMediumSimpleCuts,"Ele2PassMediumSimpleCuts/O"); 
        tree_->Branch("Ele2PassTightSimpleCuts",&fEle2PassTightSimpleCuts,"Ele2PassTightSimpleCuts/O"); 
        tree_->Branch("Ele2PassHZZICHEP2012",&fEle2PassHZZICHEP2012,"Ele2PassHZZICHEP2012/O"); 

      }

      // initialze a ZeeEventTree
      void InitTree(UInt_t Type = kCITZeeEvent){
        assert(tree_);
        // don't forget to set pointers to zero before you set address
        // or you will fully appreciate that "ROOT sucks" :)
        InitVariables();
        //Set branch address
        Int_t currentState = gErrorIgnoreLevel;

        if (Type == kCITZeeEvent) {
          cout << "CITZeeEvent\n";
          tree_->SetBranchAddress("weight",&fWeight);
          tree_->SetBranchAddress("run",&fRunNumber);
          tree_->SetBranchAddress("lumi",&fLumiSectionNumber);
          tree_->SetBranchAddress("event",&fEventNumber);
          tree_->SetBranchAddress("npu",&fNPU);
          tree_->SetBranchAddress("rho",&fRho);
          tree_->SetBranchAddress("vertices",&fNVertices);
          tree_->SetBranchAddress("triggerbits",&fEventTriggerBits); 
          tree_->SetBranchAddress("mass",&fMass);
          tree_->SetBranchAddress("Ele1Pt",&fEle1Pt);
          tree_->SetBranchAddress("Ele1Eta",&fEle1Eta);
          tree_->SetBranchAddress("Ele1Phi",&fEle1Phi);
          tree_->SetBranchAddress("Ele1SCEt",&fEle1SCEt);
          tree_->SetBranchAddress("Ele1SCEta",&fEle1SCEta);
          tree_->SetBranchAddress("Ele1SCPhi",&fEle1SCPhi);
          tree_->SetBranchAddress("Ele1GenPt",&fEle1GenPt);
          tree_->SetBranchAddress("Ele1EnergyCorrAndSmeared",&fEle1EnergyCorrAndSmeared);
          tree_->SetBranchAddress("Ele1Energy",&fEle1Energy);
          tree_->SetBranchAddress("Ele1EnergyRegression",&fEle1EnergyRegression);
          tree_->SetBranchAddress("Ele1EnergyRegressionV0",&fEle1EnergyRegressionV0);
          tree_->SetBranchAddress("Ele1EnergyRegressionV1",&fEle1EnergyRegressionV1);
          tree_->SetBranchAddress("Ele1EnergyRegressionV2",&fEle1EnergyRegressionV2);
          tree_->SetBranchAddress("Ele1EnergyRegressionErrorV0",&fEle1EnergyRegressionErrorV0);
          tree_->SetBranchAddress("Ele1EnergyRegressionErrorV1",&fEle1EnergyRegressionErrorV1);
          tree_->SetBranchAddress("Ele1EnergyRegressionErrorV2",&fEle1EnergyRegressionErrorV2);
          tree_->SetBranchAddress("Ele1Charge",&fEle1Charge); 
          tree_->SetBranchAddress("Ele1HZZICHEP2012IDMVA",&fEle1HZZICHEP2012IDMVA);
          tree_->SetBranchAddress("Ele1PFIso04",&fEle1PFIso04);
          tree_->SetBranchAddress("Ele1R9",&fEle1R9);
          tree_->SetBranchAddress("Ele1PassLooseSimpleCuts",&fEle1PassLooseSimpleCuts);
          tree_->SetBranchAddress("Ele1PassMediumSimpleCuts",&fEle1PassMediumSimpleCuts);
          tree_->SetBranchAddress("Ele1PassTightSimpleCuts",&fEle1PassTightSimpleCuts);
          tree_->SetBranchAddress("Ele1PassHZZICHEP2012",&fEle1PassHZZICHEP2012);
          tree_->SetBranchAddress("Ele2Pt",&fEle2Pt);
          tree_->SetBranchAddress("Ele2Eta",&fEle2Eta);
          tree_->SetBranchAddress("Ele2Phi",&fEle2Phi);
          tree_->SetBranchAddress("Ele2SCEt",&fEle2SCEt);
          tree_->SetBranchAddress("Ele2SCEta",&fEle2SCEta);
          tree_->SetBranchAddress("Ele2SCPhi",&fEle2SCPhi);
          tree_->SetBranchAddress("Ele2GenPt",&fEle2GenPt);
          tree_->SetBranchAddress("Ele2EnergyCorrAndSmeared",&fEle2EnergyCorrAndSmeared);
          tree_->SetBranchAddress("Ele2Energy",&fEle2Energy);
          tree_->SetBranchAddress("Ele2EnergyRegression",&fEle2EnergyRegression);
          tree_->SetBranchAddress("Ele2EnergyRegressionV0",&fEle2EnergyRegressionV0);
          tree_->SetBranchAddress("Ele2EnergyRegressionV1",&fEle2EnergyRegressionV1);
          tree_->SetBranchAddress("Ele2EnergyRegressionV2",&fEle2EnergyRegressionV2);
          tree_->SetBranchAddress("Ele2EnergyRegressionErrorV0",&fEle2EnergyRegressionErrorV0);
          tree_->SetBranchAddress("Ele2EnergyRegressionErrorV1",&fEle2EnergyRegressionErrorV1);
          tree_->SetBranchAddress("Ele2EnergyRegressionErrorV2",&fEle2EnergyRegressionErrorV2);
          tree_->SetBranchAddress("Ele2Charge",&fEle2Charge); 
          tree_->SetBranchAddress("Ele2HZZICHEP2012IDMVA",&fEle2HZZICHEP2012IDMVA);
          tree_->SetBranchAddress("Ele2PFIso04",&fEle2PFIso04);
          tree_->SetBranchAddress("Ele2R9",&fEle2R9);
          tree_->SetBranchAddress("Ele2PassLooseSimpleCuts",&fEle2PassLooseSimpleCuts);
          tree_->SetBranchAddress("Ele2PassMediumSimpleCuts",&fEle2PassMediumSimpleCuts);
          tree_->SetBranchAddress("Ele2PassTightSimpleCuts",&fEle2PassTightSimpleCuts);
          tree_->SetBranchAddress("Ele2PassHZZICHEP2012",&fEle2PassHZZICHEP2012);
        } else if (Type == kHggFutyanZeeDataEvent || Type == kHggFutyanZeeMCEvent  ) {
          cout << "other type\n";
          //************************************************************
          //The commented out branches don't exist in this version
          //************************************************************

          tree_->SetBranchAddress("weight",&fWeight);
          tree_->SetBranchAddress("run",&fRunNumber_int);
          tree_->SetBranchAddress("lumi",&fLumiSectionNumber_int);
          tree_->SetBranchAddress("event",&fEventNumber_int);
          //tree_->SetBranchAddress("rho",&fRho); 
          tree_->SetBranchAddress("nvtx",&fNVertices_int);
          //tree_->SetBranchAddress("triggerbits",&fEventTriggerBits); 
          tree_->SetBranchAddress("dipho_mass",&fMass);
          tree_->SetBranchAddress("pho1_pt",&dEle1Pt);
          tree_->SetBranchAddress("pho1_eta",&dEle1Eta);
          tree_->SetBranchAddress("pho1_phi",&dEle1Phi);
          //tree_->SetBranchAddress("Ele1SCEt",&fEle1SCEt);
          tree_->SetBranchAddress("pho1_sceta",&fEle1SCEta);
          tree_->SetBranchAddress("pho1_scphi",&fEle1SCPhi);
          //tree_->SetBranchAddress("Ele1GenPt",&fEle1GenPt);
          tree_->SetBranchAddress("pho1_energy",&dEle1EnergyCorrAndSmeared);
          tree_->SetBranchAddress("pho1_energy_noregr",&fEle1Energy);
          tree_->SetBranchAddress("pho1_energy_regr",&fEle1EnergyRegression);
          //tree_->SetBranchAddress("Ele1Charge",&fEle1Charge); 
          //tree_->SetBranchAddress("Ele1HZZICHEP2012IDMVA",&fEle1HZZICHEP2012IDMVA);
          //tree_->SetBranchAddress("Ele1PFIso04",&fEle1PFIso04);
          tree_->SetBranchAddress("pho1_r9",&fEle1R9);
          //tree_->SetBranchAddress("Ele1PassLooseSimpleCuts",&fEle1PassLooseSimpleCuts);
          //tree_->SetBranchAddress("Ele1PassMediumSimpleCuts",&fEle1PassMediumSimpleCuts);
          //tree_->SetBranchAddress("Ele1PassTightSimpleCuts",&fEle1PassTightSimpleCuts);
          //tree_->SetBranchAddress("Ele1PassHZZICHEP2012",&fEle1PassHZZICHEP2012);

          tree_->SetBranchAddress("pho2_pt",&dEle2Pt);
          tree_->SetBranchAddress("pho2_eta",&dEle2Eta);
          tree_->SetBranchAddress("pho2_phi",&dEle2Phi);
          //tree_->SetBranchAddress("Ele2SCEt",&fEle2SCEt);
          tree_->SetBranchAddress("pho2_sceta",&fEle2SCEta);
          tree_->SetBranchAddress("pho2_scphi",&fEle2SCPhi);
          //tree_->SetBranchAddress("Ele2GenPt",&fEle2GenPt);
          tree_->SetBranchAddress("pho2_energy",&dEle2EnergyCorrAndSmeared);
          tree_->SetBranchAddress("pho2_energy_noregr",&fEle2Energy);
          tree_->SetBranchAddress("pho2_energy_regr",&fEle2EnergyRegression);
          //tree_->SetBranchAddress("Ele2Charge",&fEle2Charge); 
          //tree_->SetBranchAddress("Ele2HZZICHEP2012IDMVA",&fEle2HZZICHEP2012IDMVA);
          //tree_->SetBranchAddress("Ele2PFIso04",&fEle2PFIso04);
          tree_->SetBranchAddress("pho2_r9",&fEle2R9);
          //tree_->SetBranchAddress("Ele2PassLooseSimpleCuts",&fEle2PassLooseSimpleCuts);
          //tree_->SetBranchAddress("Ele2PassMediumSimpleCuts",&fEle2PassMediumSimpleCuts);
          //tree_->SetBranchAddress("Ele2PassTightSimpleCuts",&fEle2PassTightSimpleCuts);
          //tree_->SetBranchAddress("Ele2PassHZZICHEP2012",&fEle2PassHZZICHEP2012);
        } else {
          cout << "Warning: Type " <<  Type << " is not supported \n";
        }

        gErrorIgnoreLevel = currentState;
      }

      Float_t Ele1Pt() {
        if (fTreeType == kHggFutyanZeeDataEvent || fTreeType == kHggFutyanZeeMCEvent) return dEle1Pt;
        else return fEle1Pt;
      }
      Float_t Ele1Eta() {
        if (fTreeType == kHggFutyanZeeDataEvent || fTreeType == kHggFutyanZeeMCEvent) return dEle1Eta;
        else return fEle1Eta;
      }
      Float_t Ele1Phi() {
        if (fTreeType == kHggFutyanZeeDataEvent || fTreeType == kHggFutyanZeeMCEvent) return dEle1Phi;
        else return fEle1Phi;
      }
      Float_t Ele1EnergyCorrAndSmeared() {
        if (fTreeType == kHggFutyanZeeDataEvent || fTreeType == kHggFutyanZeeMCEvent) return dEle1EnergyCorrAndSmeared;
        else return fEle1EnergyCorrAndSmeared;
      }
      Float_t Ele2Pt() {
        if (fTreeType == kHggFutyanZeeDataEvent || fTreeType == kHggFutyanZeeMCEvent) return dEle2Pt;
        else return fEle2Pt;
      }
      Float_t Ele2Eta() {
        if (fTreeType == kHggFutyanZeeDataEvent || fTreeType == kHggFutyanZeeMCEvent) return dEle2Eta;
        else return fEle2Eta;
      }
      Float_t Ele2Phi() {
        if (fTreeType == kHggFutyanZeeDataEvent || fTreeType == kHggFutyanZeeMCEvent) return dEle2Phi;
        else return fEle2Phi;
      }
      Float_t Ele2EnergyCorrAndSmeared() {
        if (fTreeType == kHggFutyanZeeDataEvent || fTreeType == kHggFutyanZeeMCEvent) return dEle2EnergyCorrAndSmeared;
        else return fEle2EnergyCorrAndSmeared;
      }
      UInt_t RunNumber() {
        if (fTreeType == kHggFutyanZeeDataEvent || fTreeType == kHggFutyanZeeMCEvent) return fRunNumber_int;
        else return fRunNumber;
      }
      UInt_t LumiSectionNumber() {
        if (fTreeType == kHggFutyanZeeDataEvent || fTreeType == kHggFutyanZeeMCEvent) return fLumiSectionNumber_int;
        else return fLumiSectionNumber;
      }
      UInt_t EventNumber() {
        if (fTreeType == kHggFutyanZeeDataEvent || fTreeType == kHggFutyanZeeMCEvent) return fEventNumber_int;
        else return fEventNumber;
      }
      UInt_t NVertices() {
        if (fTreeType == kHggFutyanZeeDataEvent || fTreeType == kHggFutyanZeeMCEvent) return fNVertices_int;
        else return fNVertices;
      }


  }; 
}


#endif
