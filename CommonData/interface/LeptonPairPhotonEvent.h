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

      enum ZmumuGammaEventTreeType { kCITLeptonPairPhotonEvent  = 0
      };

      
      /// setup variables
      UInt_t                  fTreeType;

      //Tree variables
      Float_t Weight;
      UInt_t  RunNumber;
      UInt_t  LumiSectionNumber;
      UInt_t  EventNumber;

      Float_t electronZmass;
      Float_t mllg;// 3-body mass
      Bool_t  muonZgVeto;
      Float_t muonZmass;
      Int_t   year;   

      Float_t ele1MVA;// Electron 1 MVA Value
      Float_t ele1charge;// Electron 1 Charge
      Float_t ele1energy;// Electron 1 Energy
      Float_t ele1px;// Electron 1 Px
      Float_t ele1py;// Electron 1 Py
      Float_t ele1pz;// Electron 1 Pz
      Float_t ele1pt;// Electron 1 Pt
      Float_t ele1eta;// Electron 1 Eta
      Float_t ele1mass;// Electron 1 Mass
      Float_t ele1phi;// Electron 1 Phi
      Float_t ele1dEtaIn;
      Float_t ele1dPhiIn;
      Float_t ele1sigmaIEtaIEta;
      Float_t ele1HadOverEm;
      Float_t ele1D0;
      Float_t ele1DZ;
      Float_t ele1OneOverEMinusOneOverP;
      Float_t ele1PFIsoOverPt;
      Bool_t  ele1Conversion;
      Float_t ele1missinghits;

      Float_t ele2MVA;// Electron 2 MVA Value
      Float_t ele2charge;// Electron 2 Charge
      Float_t ele2energy;// Electron 2 Energy
      Float_t ele2px;// Electron 2 Px
      Float_t ele2py;// Electron 2 Py
      Float_t ele2pz;// Electron 2 Pz
      Float_t ele2pt;// Electron 2 Pt
      Float_t ele2eta;// Electron 2 Eta
      Float_t ele2mass;// Electron 2 Mass
      Float_t ele2phi;// Electron 2 Phi
      Float_t ele2dEtaIn;
      Float_t ele2dPhiIn;
      Float_t ele2sigmaIEtaIEta;
      Float_t ele2HadOverEm;
      Float_t ele2D0;
      Float_t ele2DZ;
      Float_t ele2OneOverEMinusOneOverP;
      Float_t ele2PFIsoOverPt;
      Bool_t  ele2Conversion;
      Float_t ele2missinghits;

      Float_t chargediso_ele1;
      Float_t gammaiso_ele1;
      Float_t neutraliso_ele1;
      Float_t rho_ele1;
      Float_t effectivearea_ele1;
      Float_t chargediso_ele2;
      Float_t gammaiso_ele2;
      Float_t neutraliso_ele2;
      Float_t rho_ele2;
      Float_t effectivearea_ele2;

      Float_t photonidmva;// Photon MVA Value
      Float_t photonenergy;// Photon Energy
      Float_t photonpx;// Photon Px
      Float_t photonpy;// Photon Py
      Float_t photonpz;// Photon Pz
      Float_t photonpt;// Photon Pt
      Float_t photoneta;// Photon Eta
      Float_t photonmass;// Photon Mass??? photon->Mass() 
      Float_t photonphi;// Photon Phi
      Float_t photonr9;// Photon R9

      Float_t m1E;// Muon 1 Energy
      Float_t m1Pt;// Muon 1 Pt
      Float_t m1Mass;// Muon 1 Mass
      Float_t m1Px;// Muon 1 Px
      Float_t m1Py;// Muon 1 Py
      Float_t m1Pz;// Muon 1 Pz
      Float_t m1Eta;// Muon 1 Eta
      Float_t m1Phi;// Muon 1 Phi
      Float_t m1Charge;// Muon 1 Charge
      Float_t m2E;// Muon 2 Energy
      Float_t m2Pt;// Muon 2 Pt
      Float_t m2Mass;// Muon 2 Mass
      Float_t m2Px;// Muon 2 Px
      Float_t m2Py;// Muon 2 Py
      Float_t m2Pz;// Muon 2 Pz
      Float_t m2Eta;// Muon 2 Eta
      Float_t m2Phi;// Muon 2 Phi
      Float_t m2Charge;// Muon 2 Charge

      Float_t NPu; //Number of Pileup events
      Float_t NPuPlus;//Number of Pileup events in next signal readout
      Float_t NPuMinus;//Number of Pileup events in previous signal readout
      Float_t photonmatchmc;   

      Float_t costheta_lm_electrons;
      Float_t costheta_lp_electrons;
      Float_t phi_electrons;
      Float_t cosTheta_electrons;
      Float_t cosThetaG_electrons;
      Float_t costheta_lm_muons;
      Float_t costheta_lp_muons;
      Float_t phi_muons;
      Float_t cosTheta_muons;
      Float_t cosThetaG_muons;


    public:
      /// this is the main element
      TTree *tree_;
      TFile *f_;
  
      /// hold the names of variables to facilitate things (filled during Init)
      std::vector<std::string> variables_;
      
      /// default constructor  
      ZmumuGammaEventTree()  {
        fTreeType = kCITLeptonPairPhotonEvent;
      };
      /// default destructor
      ~ZmumuGammaEventTree(){ 
        if (f_) f_->Close();  
      };
    
      /// initialize varibles and fill list of available variables
      void InitVariables() {

        Weight = 0.0;
        RunNumber = 0.0;
        LumiSectionNumber = 0.0;
        EventNumber = 0.0;

        electronZmass = 0;
        mllg = 0;
        muonZgVeto = 0;
        muonZmass = 0;
        year = 0;

        ele1MVA = 0;
        ele1charge = 0;
        ele1energy = 0;
        ele1px = 0;
        ele1py = 0;
        ele1pz = 0;
        ele1pt = 0;
        ele1eta = 0;
        ele1mass = 0;
        ele1phi = 0;
        ele1dEtaIn = 0;
        ele1dPhiIn = 0;
        ele1sigmaIEtaIEta = 0;
        ele1HadOverEm = 0;
        ele1D0 = 0;
        ele1DZ = 0;
        ele1OneOverEMinusOneOverP = 0;
        ele1PFIsoOverPt = 0;
        ele1Conversion = 0;
        ele1missinghits = 0;

        ele2MVA = 0;
        ele2charge = 0;
        ele2energy = 0;
        ele2px = 0;
        ele2py = 0;
        ele2pz = 0;
        ele2pt = 0;
        ele2eta = 0;
        ele2mass = 0;
        ele2phi = 0;
        ele2dEtaIn = 0;
        ele2dPhiIn = 0;
        ele2sigmaIEtaIEta = 0;
        ele2HadOverEm = 0;
        ele2D0 = 0;
        ele2DZ = 0;
        ele2OneOverEMinusOneOverP = 0;
        ele2PFIsoOverPt = 0;
        ele2Conversion = 0;
        ele2missinghits = 0;

        chargediso_ele1 = 0;
        gammaiso_ele1 = 0;
        neutraliso_ele1 = 0;
        rho_ele1 = 0;
        effectivearea_ele1 = 0;
        chargediso_ele2 = 0;
        gammaiso_ele2 = 0;
        neutraliso_ele2 = 0;
        rho_ele2 = 0;
        effectivearea_ele2 = 0;

        photonidmva = 0;
        photonenergy = 0;
        photonpx = 0;
        photonpy = 0;
        photonpz = 0;
        photonpt = 0;
        photoneta = 0;
        photonmass = 0;
        photonphi = 0;
        photonr9 = 0;

        m1E = 0;
        m1Pt = 0;
        m1Mass = 0;
        m1Px = 0;
        m1Py = 0;
        m1Pz = 0;
        m1Eta = 0;
        m1Phi = 0;
        m1Charge = 0;
        m2E = 0;
        m2Pt = 0;
        m2Mass = 0;
        m2Px = 0;
        m2Py = 0;
        m2Pz = 0;
        m2Eta = 0;
        m2Phi = 0;
        m2Charge = 0;

        NPu = 0;
        NPuPlusNumber = 0;
        NPuMinusNumber = 0;
        photonmatchmc = 0;

        costheta_lm_electrons = 0;
        costheta_lp_electrons = 0;
        phi_electrons = 0;
        cosTheta_electrons = 0;
        cosThetaG_electrons = 0;
        costheta_lm_muons = 0;
        costheta_lp_muons = 0;
        phi_muons = 0;
        cosTheta_muons = 0;
        cosThetaG_muons = 0;


      }
    
      /// load a ZmumuGammaEventTree
      void LoadTree(const char* file, UInt_t Type = kCITLeptonPairPhotonEvent){
        fTreeType = Type;
        f_ = TFile::Open(file);
        assert(f_);
        if (Type == kCITLeptonPairPhotonEvent) {
          tree_ = dynamic_cast<TTree*>(f_->Get("h2LepPhotonTree"));
        } else {
          cout << "Warning: Type " <<  Type << " is not supported \n";
        }
        InitTree(Type);
        assert(tree_);
      }
    
      /// create a ZmumuGammaEventTree
      void CreateTree(){
        tree_ = new TTree("ZmumuGammaEvent","ZmumuGammaEvent");
        tree->SetAutoSave(300e9);

        f_ = 0;

        tree_->Branch("Weight",&Weight,"Weight/F");
        tree_->Branch("RunNumber",&RunNumber,"RunNumber/i");
        tree_->Branch("LumiSectionNumber",&LumiSectionNumber,"LumiSectionNumber/i");
        tree_->Branch("EventNumber",&EventNumber,"EventNumber/i");

        tree->Branch("electronZmass",&electronZmass,"electronZmass/F");
        tree->Branch("mllg",&mllg,"mllg/F");
        tree->Branch("ele1MVA",&ele1MVA,"ele1MVA/F");
        tree->Branch("ele2MVA",&ele2MVA,"ele2MVA/F");
  
        tree->Branch("ele1charge",&ele1charge,"ele1charge/F");
        tree->Branch("ele1energy",&ele1energy,"ele1energy/F");
        tree->Branch("ele1px",&ele1px,"ele1px/F");
        tree->Branch("ele1py",&ele1py,"ele1py/F");
        tree->Branch("ele1pz",&ele1pz,"ele1pz/F");
        tree->Branch("ele1pt",&ele1pt,"ele1pt/F");
        tree->Branch("ele1eta",&ele1eta,"ele1eta/F");
        tree->Branch("ele1mass",&ele1mass,"ele1mass/F");
        tree->Branch("ele1phi",&ele1phi,"ele1phi/F");

        tree->Branch("ele2charge",&ele2charge,"ele2charge/F");
        tree->Branch("ele2energy",&ele2energy,"ele2energy/F");
        tree->Branch("ele2px",&ele2px,"ele2px/F");
        tree->Branch("ele2py",&ele2py,"ele2py/F");
        tree->Branch("ele2pz",&ele2pz,"ele2pz/F");
        tree->Branch("ele2pt",&ele2pt,"ele2pt/F");
        tree->Branch("ele2eta",&ele2eta,"ele2eta/F");
        tree->Branch("ele2mass",&ele2mass,"ele2mass/F");
        tree->Branch("ele2phi",&ele2phi,"ele2phi/F");

        tree->Branch("ele1dEtaIn",&ele1dEtaIn,"ele1dEtaIn/F");
        tree->Branch("ele1dPhiIn",&ele1dPhiIn,"ele1dPhiIn/F");
        tree->Branch("ele1sigmaIEtaIEta",&ele1sigmaIEtaIEta, "ele1sigmaIEtaIEta/F");
        tree->Branch("ele1HadOverEm",&ele1HadOverEm,"ele1HadOverEm/F");
        tree->Branch("ele1D0",&ele1D0,"ele1D0/F");
        tree->Branch("ele1DZ",&ele1DZ,"ele1DZ/F");
        tree->Branch("ele1OneOverEMinusOneOverP",&ele1OneOverEMinusOneOverP,"ele1OneOverEMinusOneOverP/F");
        tree->Branch("ele1PFIsoOverPt",&ele1PFIsoOverPt,"ele1PFIsoOverPt/F");
        tree->Branch("ele1Conversion",&ele1Conversion,"ele1Conversion/O");
        tree->Branch("ele1missinghits",&ele1missinghits,"ele1missinghits/F");

        tree->Branch("ele2dEtaIn",&ele2dEtaIn,"ele2dEtaIn/F");
        tree->Branch("ele2dPhiIn",&ele2dPhiIn,"ele2dPhiIn/F");
        tree->Branch("ele2sigmaIEtaIEta",&ele2sigmaIEtaIEta, "ele2sigmaIEtaIEta/F");
        tree->Branch("ele2HadOverEm",&ele2HadOverEm,"ele2HadOverEm/F");
        tree->Branch("ele2D0",&ele2D0,"ele2D0/F");
        tree->Branch("ele2DZ",&ele2DZ,"ele2DZ/F");
        tree->Branch("ele2OneOverEMinusOneOverP",&ele2OneOverEMinusOneOverP,"ele2OneOverEMinusOneOverP/F");
        tree->Branch("ele2PFIsoOverPt",&ele2PFIsoOverPt,"ele2PFIsoOverPt/F");
        tree->Branch("ele2Conversion",&ele2Conversion,"ele2Conversion/O");
        tree->Branch("ele2missinghits",&ele2missinghits,"ele2missinghits/F");

        tree->Branch("chargediso_ele1",&chargediso_ele1,"chargediso_ele1/F");
        tree->Branch("gammaiso_ele1",&gammaiso_ele1,"gammaiso_ele1/F");
        tree->Branch("neutraliso_ele1",&neutraliso_ele1,"neutraliso_ele1/F");
        tree->Branch("rho_ele1",&rho_ele1,"rho_ele1/F");
        tree->Branch("effectivearea_ele1",&effectivearea_ele1,"effectivearea_ele1/F");
        tree->Branch("chargediso_ele2",&chargediso_ele2,"chargediso_ele2/F");
        tree->Branch("gammaiso_ele2",&gammaiso_ele2,"gammaiso_ele2/F");
        tree->Branch("neutraliso_ele2",&neutraliso_ele2,"neutraliso_ele2/F");
        tree->Branch("rho_ele2",&rho_ele2,"rho_ele2/F");
        tree->Branch("effectivearea_ele2",&effectivearea_ele2,"effectivearea_ele2/F");

        tree->Branch("costheta_lm_electrons",&costheta_lm_electrons,"costheta_lm_electrons/F");
        tree->Branch("costheta_lp_electrons",&costheta_lp_electrons,"costheta_lp_electrons/F");
        tree->Branch("phi_electrons",&phi_electrons,"phi_electrons/F");
        tree->Branch("cosTheta_electrons",&cosTheta_electrons,"cosTheta_electrons/F");
        tree->Branch("cosThetaG_electrons",&cosThetaG_electrons,"cosThetaG_electrons");
        tree->Branch("costheta_lm_muons",&costheta_lm_muons,"costheta_lm_muons/F");
        tree->Branch("costheta_lp_muons",&costheta_lp_muons,"costheta_lp_muons/F");
        tree->Branch("phi_muons",&phi_muons,"phi_muons/F");
        tree->Branch("cosTheta_muons",&cosTheta_muons,"cosTheta_muons/F");
        tree->Branch("cosThetaG_muons",&cosThetaG_muons,"cosThetaG_muons/F");

        tree->Branch("muonZgVeto",&muonZgVeto,"muonZgVeto/O");
        tree->Branch("muonZmass",&muonZmass,"muonZmass/F");
        tree->Branch("m1E",&m1E,"m1E/F");
        tree->Branch("m1Pt",&m1Pt,"m1Pt/F");
        tree->Branch("m1Mass",&m1Mass,"m1Mass/F");
        tree->Branch("m1Px",&m1Px,"m1Px/F");
        tree->Branch("m1Py",&m1Py,"m1Py/F");
        tree->Branch("m1Pz",&m1Pz,"m1Pz/F");
        tree->Branch("m1Eta",&m1Eta,"m1Eta/F");
        tree->Branch("m1Phi",&m1Phi,"m1Phi/F");
        tree->Branch("m1Charge",&m1Charge,"m1Charge/F");
        tree->Branch("m2E",&m2E,"m2E/F");
        tree->Branch("m2Pt",&m2Pt,"m2Pt/F");
        tree->Branch("m2Mass",&m2Mass,"m2Mass/F");
        tree->Branch("m2Px",&m2Px,"m2Px/F");
        tree->Branch("m2Py",&m2Py,"m2Py/F");
        tree->Branch("m2Pz",&m2Pz,"m2Pz/F");
        tree->Branch("m2Eta",&m2Eta,"m2Eta/F");
        tree->Branch("m2Phi",&m2Phi,"m2Phi/F");
        tree->Branch("m2Charge",&m2Charge,"m2Charge/F");

        tree->Branch("photonidmva",&photonidmva,"photonidmva/F");
        tree->Branch("photonr9",&photonr9,"photonr9/F");
        tree->Branch("photonenergy",&photonenergy,"photonenergy/F");
        tree->Branch("photonpx",&photonpx,"photonpx/F");
        tree->Branch("photonpy",&photonpy,"photonpy/F");
        tree->Branch("photonpz",&photonpz,"photonpz/F");
        tree->Branch("photonpt",&photonpt,"photonpt/F");
        tree->Branch("photoneta",&photoneta,"photoneta/F");
        tree->Branch("photonmass",&photonmass,"photonmass/F");
        tree->Branch("photonphi",&photonphi,"photonphi/F");
 
        tree->Branch("NPu",&NPu,"NPu/F");
        tree->Branch("NPuPlus",&NPuPlus,"NPuPlus/F");
        tree->Branch("NPuMinus",&NPuMinus,"NPuMinus/F");
 
        tree->Branch("photonmatchmc",&photonmatchmc,"photonmatchmc/F");

      }

      // initialze a ZmumuGammaEventTree
      void InitTree(UInt_t Type = kCITLeptonPairPhotonEvent){
        assert(tree_);
        // don't forget to set pointers to zero before you set address
        // or you will fully appreciate that "ROOT sucks" :)
        InitVariables();
        //Set branch address
        Int_t currentState = gErrorIgnoreLevel;

        if (Type == kCITLeptonPairPhotonEvent) {
          cout << "kCITLeptonPairPhotonEvent\n";
//           tree_->SetBranchAddress("weight",&fWeight);


        tree_->SetBranchAddress("Weight",&Weight);
        tree_->SetBranchAddress("RunNumber",&RunNumber);
        tree_->SetBranchAddress("LumiSectionNumber",&LumiSectionNumber);
        tree_->SetBranchAddress("EventNumber",&EventNumber);

        tree->SetBranchAddress("electronZmass",&electronZmass);
        tree->SetBranchAddress("mllg",&mllg);
        tree->SetBranchAddress("ele1MVA",&ele1MVA);
        tree->SetBranchAddress("ele2MVA",&ele2MVA);
  
        tree->SetBranchAddress("ele1charge",&ele1charge);
        tree->SetBranchAddress("ele1energy",&ele1energy);
        tree->SetBranchAddress("ele1px",&ele1px);
        tree->SetBranchAddress("ele1py",&ele1py);
        tree->SetBranchAddress("ele1pz",&ele1pz);
        tree->SetBranchAddress("ele1pt",&ele1pt);
        tree->SetBranchAddress("ele1eta",&ele1eta);
        tree->SetBranchAddress("ele1mass",&ele1mass);
        tree->SetBranchAddress("ele1phi",&ele1phi);

        tree->SetBranchAddress("ele2charge",&ele2charge);
        tree->SetBranchAddress("ele2energy",&ele2energy);
        tree->SetBranchAddress("ele2px",&ele2px);
        tree->SetBranchAddress("ele2py",&ele2py);
        tree->SetBranchAddress("ele2pz",&ele2pz);
        tree->SetBranchAddress("ele2pt",&ele2pt);
        tree->SetBranchAddress("ele2eta",&ele2eta);
        tree->SetBranchAddress("ele2mass",&ele2mass);
        tree->SetBranchAddress("ele2phi",&ele2phi);

        tree->SetBranchAddress("ele1dEtaIn",&ele1dEtaIn);
        tree->SetBranchAddress("ele1dPhiIn",&ele1dPhiIn);
        tree->SetBranchAddress("ele1sigmaIEtaIEta",&ele1sigmaIEtaIEta, "ele1sigmaIEtaIEta/F");
        tree->SetBranchAddress("ele1HadOverEm",&ele1HadOverEm);
        tree->SetBranchAddress("ele1D0",&ele1D0);
        tree->SetBranchAddress("ele1DZ",&ele1DZ);
        tree->SetBranchAddress("ele1OneOverEMinusOneOverP",&ele1OneOverEMinusOneOverP);
        tree->SetBranchAddress("ele1PFIsoOverPt",&ele1PFIsoOverPt);
        tree->SetBranchAddress("ele1Conversion",&ele1Conversion);
        tree->SetBranchAddress("ele1missinghits",&ele1missinghits);

        tree->SetBranchAddress("ele2dEtaIn",&ele2dEtaIn);
        tree->SetBranchAddress("ele2dPhiIn",&ele2dPhiIn);
        tree->SetBranchAddress("ele2sigmaIEtaIEta",&ele2sigmaIEtaIEta, "ele2sigmaIEtaIEta/F");
        tree->SetBranchAddress("ele2HadOverEm",&ele2HadOverEm);
        tree->SetBranchAddress("ele2D0",&ele2D0);
        tree->SetBranchAddress("ele2DZ",&ele2DZ);
        tree->SetBranchAddress("ele2OneOverEMinusOneOverP",&ele2OneOverEMinusOneOverP);
        tree->SetBranchAddress("ele2PFIsoOverPt",&ele2PFIsoOverPt);
        tree->SetBranchAddress("ele2Conversion",&ele2Conversion);
        tree->SetBranchAddress("ele2missinghits",&ele2missinghits);

        tree->SetBranchAddress("chargediso_ele1",&chargediso_ele1);
        tree->SetBranchAddress("gammaiso_ele1",&gammaiso_ele1);
        tree->SetBranchAddress("neutraliso_ele1",&neutraliso_ele1);
        tree->SetBranchAddress("rho_ele1",&rho_ele1);
        tree->SetBranchAddress("effectivearea_ele1",&effectivearea_ele1);
        tree->SetBranchAddress("chargediso_ele2",&chargediso_ele2);
        tree->SetBranchAddress("gammaiso_ele2",&gammaiso_ele2);
        tree->SetBranchAddress("neutraliso_ele2",&neutraliso_ele2);
        tree->SetBranchAddress("rho_ele2",&rho_ele2);
        tree->SetBranchAddress("effectivearea_ele2",&effectivearea_ele2);

        tree->SetBranchAddress("costheta_lm_electrons",&costheta_lm_electrons);
        tree->SetBranchAddress("costheta_lp_electrons",&costheta_lp_electrons);
        tree->SetBranchAddress("phi_electrons",&phi_electrons);
        tree->SetBranchAddress("cosTheta_electrons",&cosTheta_electrons);
        tree->SetBranchAddress("cosThetaG_electrons",&cosThetaG_electrons);
        tree->SetBranchAddress("costheta_lm_muons",&costheta_lm_muons);
        tree->SetBranchAddress("costheta_lp_muons",&costheta_lp_muons);
        tree->SetBranchAddress("phi_muons",&phi_muons);
        tree->SetBranchAddress("cosTheta_muons",&cosTheta_muons);
        tree->SetBranchAddress("cosThetaG_muons",&cosThetaG_muons);

        tree->SetBranchAddress("muonZgVeto",&muonZgVeto);
        tree->SetBranchAddress("muonZmass",&muonZmass);
        tree->SetBranchAddress("m1E",&m1E);
        tree->SetBranchAddress("m1Pt",&m1Pt);
        tree->SetBranchAddress("m1Mass",&m1Mass);
        tree->SetBranchAddress("m1Px",&m1Px);
        tree->SetBranchAddress("m1Py",&m1Py);
        tree->SetBranchAddress("m1Pz",&m1Pz);
        tree->SetBranchAddress("m1Eta",&m1Eta);
        tree->SetBranchAddress("m1Phi",&m1Phi);
        tree->SetBranchAddress("m1Charge",&m1Charge);
        tree->SetBranchAddress("m2E",&m2E);
        tree->SetBranchAddress("m2Pt",&m2Pt);
        tree->SetBranchAddress("m2Mass",&m2Mass);
        tree->SetBranchAddress("m2Px",&m2Px);
        tree->SetBranchAddress("m2Py",&m2Py);
        tree->SetBranchAddress("m2Pz",&m2Pz);
        tree->SetBranchAddress("m2Eta",&m2Eta);
        tree->SetBranchAddress("m2Phi",&m2Phi);
        tree->SetBranchAddress("m2Charge",&m2Charge);

        tree->SetBranchAddress("photonidmva",&photonidmva);
        tree->SetBranchAddress("photonr9",&photonr9);
        tree->SetBranchAddress("photonenergy",&photonenergy);
        tree->SetBranchAddress("photonpx",&photonpx);
        tree->SetBranchAddress("photonpy",&photonpy);
        tree->SetBranchAddress("photonpz",&photonpz);
        tree->SetBranchAddress("photonpt",&photonpt);
        tree->SetBranchAddress("photoneta",&photoneta);
        tree->SetBranchAddress("photonmass",&photonmass);
        tree->SetBranchAddress("photonphi",&photonphi);
 
        tree->SetBranchAddress("NPu",&NPu);
        tree->SetBranchAddress("NPuPlus",&NPuPlus);
        tree->SetBranchAddress("NPuMinus",&NPuMinus);
 
        tree->SetBranchAddress("photonmatchmc",&photonmatchmc);


        } else {
          cout << "Warning: Type " <<  Type << " is not supported \n";
        }

        gErrorIgnoreLevel = currentState;
      }

  }; 
}


#endif
