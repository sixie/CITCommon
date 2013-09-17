#ifndef ZmumuGammaEventTree_H
#define ZmumuGammaEventTree_H

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include <cmath>
#include "assert.h"

namespace citana
{
  class ZmumuGammaEventTree {

    public:

      enum ZmumuGammaEventTreeType { kCITZmumuGammaEvent    = 0
      };

      
      /// setup variables
      UInt_t                  fTreeType;

      //Tree variables
      Float_t                 fWeight;
      UInt_t                  fRunNumber;
      UInt_t                  fLumiSectionNumber;
      UInt_t                  fEventNumber;
      UInt_t                  fNPU;
      Float_t                 fRho; 
      UInt_t                  fNVertices; 
      Float_t                 fMass;
      Float_t                 fDileptonMass;
      Float_t                 fPhotonPt; 
      Float_t                 fPhotonEta; 
      Float_t                 fPhotonPhi; 
      Float_t                 fPhotonR9; 
      Float_t                 fPhotonIsEB; 
      Float_t                 fPhotonHoE; 
      Float_t                 fSCEt; 
      Float_t                 fSCE; 
      Float_t                 fSCRawE; 
      Float_t                 fSCEta; 
      Float_t                 fSCEtaWidth; 
      Float_t                 fSCPhiWidth; 
      Float_t                 fPhoToTrackDeltaR; 
      Float_t                 fPhoPassElectronVeto; 
      Float_t                 fPhoHasMatchedConversion; 
      Float_t                 fPhoHasPixelMatch; 
      Float_t                 fPhoSigmaIEtaIEta; 
      Float_t                 fPhoTrackIso; 
      Float_t                 fPhoEcalIso; 
      Float_t                 fPhoHcalIso; 
      Float_t                 fPhoPdgId; 
      Float_t                 fPhoCrackCorr; 
      Float_t                 fMu1Pt; 
      Float_t                 fMu1Eta; 
      Float_t                 fMu1Phi; 
      Float_t                 fMu1Charge;
      Float_t                 fMu1TrackChi2; 
      Float_t                 fMu1TrackNormalizedChi2; 
      Float_t                 fMu1DeltaR; 
      Float_t                 fMu1CalEnergyEm; 
      Float_t                 fMu1CalEnergyEmMax; 
      Float_t                 fMu1CalEnergyEmHad; 
      Float_t                 fMu2Pt; 
      Float_t                 fMu2Eta; 
      Float_t                 fMu2Phi; 
      Float_t                 fMu2Charge;
      Float_t                 fMu2TrackChi2; 
      Float_t                 fMu2TrackNormalizedChi2; 
      Float_t                 fMu2DeltaR; 
      Float_t                 fMu2CalEnergyEm; 
      Float_t                 fMu2CalEnergyEmMax; 
      Float_t                 fMu2CalEnergyEmHad; 
      Float_t                 fMinDeltaEta;
      Float_t                 fMinDeltaR;
      Float_t                 fMinDeltaPhi;
      Float_t                 fkRatio;
      Float_t                 fPreshowerE;

      Float_t                 fGenMu1Pt; 
      Float_t                 fGenMu1Eta; 
      Float_t                 fGenMu1Phi; 
      Float_t                 fGenMu2Pt; 
      Float_t                 fGenMu2Eta; 
      Float_t                 fGenMu2Phi;
      Float_t                 fGenPhoE;  
      Float_t                 fGenPhoEt; 
      Float_t                 fGenPhoEta; 
      Float_t                 fGenPhoPhi; 
      Float_t                 fGenPhoMotherPdgId; 
      Bool_t                  fIsFSR;
      Bool_t                  fIsISR;

      Int_t                   fPhoIEtaX; 
      Int_t                   fPhoIPhiY; 
      Int_t                   fMuNearIEtaX; 
      Int_t                   fMuNearIPhiY; 
      Bool_t                  fMuNearIsEB; 
      Int_t                   fMuNearIndex;




    public:
      /// this is the main element
      TTree *tree_;
      TFile *f_;
  
      /// hold the names of variables to facilitate things (filled during Init)
      std::vector<std::string> variables_;
      
      /// default constructor  
      ZmumuGammaEventTree()  {
        fTreeType = kCITZmumuGammaEvent;
      };
      /// default destructor
      ~ZmumuGammaEventTree(){ 
        if (f_) f_->Close();  
      };
    
      /// initialize varibles and fill list of available variables
      void InitVariables() {

        fWeight = 0.0;
        fRunNumber = 0.0;
        fLumiSectionNumber = 0.0;
        fEventNumber = 0.0;
        fNPU = 0.0;
        fRho = 0.0;
        fNVertices = 0.0;
        fMass = 0.0;
        fDileptonMass = 0.0;
        fPhotonPt = 0.0;
        fPhotonEta = 0.0;
        fPhotonPhi = 0.0;
        fPhotonR9 = 0.0;
        fPhotonIsEB = 0.0;
        fPhotonHoE = 0.0;
        fSCEt = 0.0;
        fSCE = 0.0;
        fSCRawE = 0.0;
        fSCEta = 0.0;
        fSCEtaWidth = 0.0;
        fSCPhiWidth = 0.0;
        fPhoToTrackDeltaR = 0.0;
        fPhoPassElectronVeto = 0.0;
        fPhoHasMatchedConversion = 0.0;
        fPhoHasPixelMatch = 0.0;
        fPhoSigmaIEtaIEta = 0.0;
        fPhoTrackIso = 0.0;
        fPhoEcalIso = 0.0;
        fPhoHcalIso = 0.0;
        fPhoPdgId = 0.0;
        fPhoCrackCorr = 0.0;
        fMu1Pt = 0.0;
        fMu1Eta = 0.0;
        fMu1Phi = 0.0;
        fMu1Charge = 0;
        fMu1TrackChi2 = 0.0;
        fMu1TrackNormalizedChi2 = 0.0;
        fMu1DeltaR = 0.0;
        fMu1CalEnergyEm = 0.0;
        fMu1CalEnergyEmMax = 0.0;
        fMu1CalEnergyEmHad = 0.0;
        fMu2Pt = 0.0;
        fMu2Eta = 0.0;
        fMu2Phi = 0.0;
        fMu2Charge = 0;
        fMu2TrackChi2 = 0.0;
        fMu2TrackNormalizedChi2 = 0.0;
        fMu2DeltaR = 0.0;
        fMu2CalEnergyEm = 0.0;
        fMu2CalEnergyEmMax = 0.0;
        fMu2CalEnergyEmHad = 0.0;
        fMinDeltaEta = 0.0;
        fMinDeltaR = 0.0;
        fMinDeltaPhi = 0.0;
        fkRatio = 0.0;
        fPreshowerE = 0.0;
        fGenMu1Pt = 0.0;
        fGenMu1Eta = 0.0;
        fGenMu1Phi = 0.0;
        fGenMu2Pt = 0.0;
        fGenMu2Eta = 0.0;
        fGenMu2Phi = 0.0;
        fGenPhoE = 0.0;
        fGenPhoEt = 0.0;
        fGenPhoEta = 0.0;
        fGenPhoPhi = 0.0;
        fGenPhoMotherPdgId = 0.0;
        fIsFSR = false;
        fIsISR = false;
        fPhoIEtaX = 0.0;
        fPhoIPhiY = 0.0;
        fMuNearIEtaX = 0.0;
        fMuNearIPhiY = 0.0;
        fMuNearIsEB = false;
        fMuNearIndex = 0.0;

      }
    
      /// load a ZmumuGammaEventTree
      void LoadTree(const char* file, UInt_t Type = kCITZmumuGammaEvent){
        fTreeType = Type;
        f_ = TFile::Open(file);
        assert(f_);
        if (Type == kCITZmumuGammaEvent) {
          tree_ = dynamic_cast<TTree*>(f_->Get("ZmumuGammaEvent"));
        } else {
          cout << "Warning: Type " <<  Type << " is not supported \n";
        }
        InitTree(Type);
        assert(tree_);
      }
    
      /// create a ZmumuGammaEventTree
      void CreateTree(){
        tree_ = new TTree("ZmumuGammaEvent","ZmumuGammaEvent");
        f_ = 0;

        //book the branches
        tree_->Branch("Weight",&fWeight,"Weight/F");
        tree_->Branch("RunNumber",&fRunNumber,"RunNumber/i");
        tree_->Branch("LumiSectionNumber",&fLumiSectionNumber,"LumiSectionNumber/i");
        tree_->Branch("EventNumber",&fEventNumber,"EventNumber/i");
        tree_->Branch("NPU",&fNPU,"NPU/i");
        tree_->Branch("Rho",&fRho,"Rho/F");
        tree_->Branch("NVertices",&fNVertices,"NVertices/i");
        tree_->Branch("Mass",&fMass,"Mass/F");
        tree_->Branch("DileptonMass",&fDileptonMass,"DileptonMass/F");
        tree_->Branch("PhotonPt",&fPhotonPt,"PhotonPt/F");
        tree_->Branch("PhotonEta",&fPhotonEta,"PhotonEta/F");
        tree_->Branch("PhotonPhi",&fPhotonPhi,"PhotonPhi/F");
        tree_->Branch("PhotonR9",&fPhotonR9,"PhotonR9/F");
        tree_->Branch("PhotonIsEB",&fPhotonIsEB,"PhotonIsEB/F");
        tree_->Branch("PhotonHoE",&fPhotonHoE,"PhotonHoE/F");
        tree_->Branch("SCEt",&fSCEt,"SCEt/F");
        tree_->Branch("SCE",&fSCE,"SCE/F");
        tree_->Branch("SCRawE",&fSCRawE,"SCRawE/F");
        tree_->Branch("SCEta",&fSCEta,"SCEta/F");
        tree_->Branch("SCEtaWidth",&fSCEtaWidth,"SCEtaWidth/F");
        tree_->Branch("SCPhiWidth",&fSCPhiWidth,"SCPhiWidth/F");
        tree_->Branch("PhoToTrackDeltaR",&fPhoToTrackDeltaR,"PhoToTrackDeltaR/F");
        tree_->Branch("PhoPassElectronVeto",&fPhoPassElectronVeto,"PhoPassElectronVeto/F");
        tree_->Branch("PhoHasMatchedConversion",&fPhoHasMatchedConversion,"PhoHasMatchedConversion/F");
        tree_->Branch("PhoHasPixelMatch",&fPhoHasPixelMatch,"PhoHasPixelMatch/F");
        tree_->Branch("PhoSigmaIEtaIEta",&fPhoSigmaIEtaIEta,"PhoSigmaIEtaIEta/F");
        tree_->Branch("PhoTrackIso",&fPhoTrackIso,"PhoTrackIso/F");
        tree_->Branch("PhoEcalIso",&fPhoEcalIso,"PhoEcalIso/F");
        tree_->Branch("PhoHcalIso",&fPhoHcalIso,"PhoHcalIso/F");
        tree_->Branch("PhoPdgId",&fPhoPdgId,"PhoPdgId/F");
        tree_->Branch("PhoCrackCorr",&fPhoCrackCorr,"PhoCrackCorr/F");
        tree_->Branch("Mu1Pt",&fMu1Pt,"Mu1Pt/F");
        tree_->Branch("Mu1Eta",&fMu1Eta,"Mu1Eta/F");
        tree_->Branch("Mu1Phi",&fMu1Phi,"Mu1Phi/F");
        tree_->Branch("Mu1Charge",&fMu1Charge,"Mu1Charge/F");
        tree_->Branch("Mu1TrackChi2",&fMu1TrackChi2,"Mu1TrackChi2/F");
        tree_->Branch("Mu1TrackNormalizedChi2",&fMu1TrackNormalizedChi2,"Mu1TrackNormalizedChi2/F");
        tree_->Branch("Mu1DeltaR",&fMu1DeltaR,"Mu1DeltaR/F");
        tree_->Branch("Mu1CalEnergyEm",&fMu1CalEnergyEm,"Mu1CalEnergyEm/F");
        tree_->Branch("Mu1CalEnergyEmMax",&fMu1CalEnergyEmMax,"Mu1CalEnergyEmMax/F");
        tree_->Branch("Mu1CalEnergyEmHad",&fMu1CalEnergyEmHad,"Mu1CalEnergyEmHad/F");
        tree_->Branch("Mu2Pt",&fMu2Pt,"Mu2Pt/F");
        tree_->Branch("Mu2Eta",&fMu2Eta,"Mu2Eta/F");
        tree_->Branch("Mu2Phi",&fMu2Phi,"Mu2Phi/F");
        tree_->Branch("Mu2Charge",&fMu2Charge,"Mu2Charge/F");
        tree_->Branch("Mu2TrackChi2",&fMu2TrackChi2,"Mu2TrackChi2/F");
        tree_->Branch("Mu2TrackNormalizedChi2",&fMu2TrackNormalizedChi2,"Mu2TrackNormalizedChi2/F");
        tree_->Branch("Mu2DeltaR",&fMu2DeltaR,"Mu2DeltaR/F");
        tree_->Branch("Mu2CalEnergyEm",&fMu2CalEnergyEm,"Mu2CalEnergyEm/F");
        tree_->Branch("Mu2CalEnergyEmMax",&fMu2CalEnergyEmMax,"Mu2CalEnergyEmMax/F");
        tree_->Branch("Mu2CalEnergyEmHad",&fMu2CalEnergyEmHad,"Mu2CalEnergyEmHad/F");
        tree_->Branch("MinDeltaEta",&fMinDeltaEta,"MinDeltaEta/F");
        tree_->Branch("MinDeltaR",&fMinDeltaR,"MinDeltaR/F");
        tree_->Branch("MinDeltaPhi",&fMinDeltaPhi,"MinDeltaPhi/F");
        tree_->Branch("kRatio",&fkRatio,"kRatio/F");
        tree_->Branch("PreshowerE",&fPreshowerE,"PreshowerE/F");
        tree_->Branch("GenMu1Pt",&fGenMu1Pt,"GenMu1Pt/F");
        tree_->Branch("GenMu1Eta",&fGenMu1Eta,"GenMu1Eta/F");
        tree_->Branch("GenMu1Phi",&fGenMu1Phi,"GenMu1Phi/F");
        tree_->Branch("GenMu2Pt",&fGenMu2Pt,"GenMu2Pt/F");
        tree_->Branch("GenMu2Eta",&fGenMu2Eta,"GenMu2Eta/F");
        tree_->Branch("GenMu2Phi",&fGenMu2Phi,"GenMu2Phi/F");
        tree_->Branch("GenPhoE",&fGenPhoE,"GenPhoE/F");
        tree_->Branch("GenPhoEt",&fGenPhoEt,"GenPhoEt/F");
        tree_->Branch("GenPhoEta",&fGenPhoEta,"GenPhoEta/F");
        tree_->Branch("GenPhoPhi",&fGenPhoPhi,"GenPhoPhi/F");
        tree_->Branch("GenPhoMotherPdgId",&fGenPhoMotherPdgId,"GenPhoMotherPdgId/F");
        tree_->Branch("IsFSR",&fIsFSR,"IsFSR/O");
        tree_->Branch("IsISR",&fIsISR,"IsISR/O");
        tree_->Branch("PhoIEtaX",&fPhoIEtaX,"PhoIEtaX/I");
        tree_->Branch("PhoIPhiY",&fPhoIPhiY,"PhoIPhiY/I");
        tree_->Branch("MuNearIEtaX",&fMuNearIEtaX,"MuNearIEtaX/I");
        tree_->Branch("MuNearIPhiY",&fMuNearIPhiY,"MuNearIPhiY/I");
        tree_->Branch("MuNearIsEB",&fMuNearIsEB,"MuNearIsEB/O");
        tree_->Branch("MuNearIndex",&fMuNearIndex,"MuNearIndex/I");
      }

      // initialze a ZmumuGammaEventTree
      void InitTree(UInt_t Type = kCITZmumuGammaEvent){
        assert(tree_);
        // don't forget to set pointers to zero before you set address
        // or you will fully appreciate that "ROOT sucks" :)
        InitVariables();
        //Set branch address
        Int_t currentState = gErrorIgnoreLevel;

        if (Type == kCITZmumuGammaEvent) {
          cout << "CITZmumuGammaEvent\n";
          tree_->SetBranchAddress("Weight",&fWeight);
          tree_->SetBranchAddress("RunNumber",&fRunNumber);
          tree_->SetBranchAddress("LumiSectionNumber",&fLumiSectionNumber);
          tree_->SetBranchAddress("EventNumber",&fEventNumber);
          tree_->SetBranchAddress("NPU",&fNPU);
          tree_->SetBranchAddress("Rho",&fRho);
          tree_->SetBranchAddress("NVertices",&fNVertices);
          tree_->SetBranchAddress("Mass",&fMass);
          tree_->SetBranchAddress("DileptonMass",&fDileptonMass);
          tree_->SetBranchAddress("PhotonPt",&fPhotonPt);
          tree_->SetBranchAddress("PhotonEta",&fPhotonEta);
          tree_->SetBranchAddress("PhotonPhi",&fPhotonPhi);
          tree_->SetBranchAddress("PhotonR9",&fPhotonR9);
          tree_->SetBranchAddress("PhotonIsEB",&fPhotonIsEB);
          tree_->SetBranchAddress("PhotonHoE",&fPhotonHoE);
          tree_->SetBranchAddress("SCEt",&fSCEt);
          tree_->SetBranchAddress("SCE",&fSCE);
          tree_->SetBranchAddress("SCRawE",&fSCRawE);
          tree_->SetBranchAddress("SCEta",&fSCEta);
          tree_->SetBranchAddress("SCEtaWidth",&fSCEtaWidth);
          tree_->SetBranchAddress("SCPhiWidth",&fSCPhiWidth);
          tree_->SetBranchAddress("PhoToTrackDeltaR",&fPhoToTrackDeltaR);
          tree_->SetBranchAddress("PhoPassElectronVeto",&fPhoPassElectronVeto);
          tree_->SetBranchAddress("PhoHasMatchedConversion",&fPhoHasMatchedConversion);
          tree_->SetBranchAddress("PhoHasPixelMatch",&fPhoHasPixelMatch);
          tree_->SetBranchAddress("PhoSigmaIEtaIEta",&fPhoSigmaIEtaIEta);
          tree_->SetBranchAddress("PhoTrackIso",&fPhoTrackIso);
          tree_->SetBranchAddress("PhoEcalIso",&fPhoEcalIso);
          tree_->SetBranchAddress("PhoHcalIso",&fPhoHcalIso);
          tree_->SetBranchAddress("PhoPdgId",&fPhoPdgId);
          tree_->SetBranchAddress("PhoCrackCorr",&fPhoCrackCorr);
          tree_->SetBranchAddress("Mu1Pt",&fMu1Pt);
          tree_->SetBranchAddress("Mu1Eta",&fMu1Eta);
          tree_->SetBranchAddress("Mu1Phi",&fMu1Phi);
          tree_->SetBranchAddress("Mu1Charge",&fMu1Charge);
          tree_->SetBranchAddress("Mu1TrackChi2",&fMu1TrackChi2);
          tree_->SetBranchAddress("Mu1TrackNormalizedChi2",&fMu1TrackNormalizedChi2);
          tree_->SetBranchAddress("Mu1DeltaR",&fMu1DeltaR);
          tree_->SetBranchAddress("Mu1CalEnergyEm",&fMu1CalEnergyEm);
          tree_->SetBranchAddress("Mu1CalEnergyEmMax",&fMu1CalEnergyEmMax);
          tree_->SetBranchAddress("Mu1CalEnergyEmHad",&fMu1CalEnergyEmHad);
          tree_->SetBranchAddress("Mu2Pt",&fMu2Pt);
          tree_->SetBranchAddress("Mu2Eta",&fMu2Eta);
          tree_->SetBranchAddress("Mu2Phi",&fMu2Phi);
          tree_->SetBranchAddress("Mu2Charge",&fMu2Charge);
          tree_->SetBranchAddress("Mu2TrackChi2",&fMu2TrackChi2);
          tree_->SetBranchAddress("Mu2TrackNormalizedChi2",&fMu2TrackNormalizedChi2);
          tree_->SetBranchAddress("Mu2DeltaR",&fMu2DeltaR);
          tree_->SetBranchAddress("Mu2CalEnergyEm",&fMu2CalEnergyEm);
          tree_->SetBranchAddress("Mu2CalEnergyEmMax",&fMu2CalEnergyEmMax);
          tree_->SetBranchAddress("Mu2CalEnergyEmHad",&fMu2CalEnergyEmHad);
          tree_->SetBranchAddress("MinDeltaEta",&fMinDeltaEta);
          tree_->SetBranchAddress("MinDeltaR",&fMinDeltaR);
          tree_->SetBranchAddress("MinDeltaPhi",&fMinDeltaPhi);
          tree_->SetBranchAddress("kRatio",&fkRatio);
          tree_->SetBranchAddress("PreshowerE",&fPreshowerE);
          tree_->SetBranchAddress("GenMu1Pt",&fGenMu1Pt);
          tree_->SetBranchAddress("GenMu1Eta",&fGenMu1Eta);
          tree_->SetBranchAddress("GenMu1Phi",&fGenMu1Phi);
          tree_->SetBranchAddress("GenMu2Pt",&fGenMu2Pt);
          tree_->SetBranchAddress("GenMu2Eta",&fGenMu2Eta);
          tree_->SetBranchAddress("GenMu2Phi",&fGenMu2Phi);
          tree_->SetBranchAddress("GenPhoE",&fGenPhoE);
          tree_->SetBranchAddress("GenPhoEt",&fGenPhoEt);
          tree_->SetBranchAddress("GenPhoEta",&fGenPhoEta);
          tree_->SetBranchAddress("GenPhoPhi",&fGenPhoPhi);
          tree_->SetBranchAddress("GenPhoMotherPdgId",&fGenPhoMotherPdgId);
          tree_->SetBranchAddress("IsFSR",&fIsFSR);
          tree_->SetBranchAddress("IsISR",&fIsISR);
          tree_->SetBranchAddress("PhoIEtaX",&fPhoIEtaX);
          tree_->SetBranchAddress("PhoIPhiY",&fPhoIPhiY);
          tree_->SetBranchAddress("MuNearIEtaX",&fMuNearIEtaX);
          tree_->SetBranchAddress("MuNearIPhiY",&fMuNearIPhiY);
          tree_->SetBranchAddress("MuNearIsEB",&fMuNearIsEB);
          tree_->SetBranchAddress("MuNearIndex",&fMuNearIndex);
        } else {
          cout << "Warning: Type " <<  Type << " is not supported \n";
        }

        gErrorIgnoreLevel = currentState;
      }

  }; 
}


#endif
