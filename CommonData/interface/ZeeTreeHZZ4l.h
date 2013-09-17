//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri May 17 11:37:45 2013 by ROOT version 5.32/03
// from TTree probe_tree/probe_tree
// found on file: zEE_lineshape_Data_2012.root
//////////////////////////////////////////////////////////

#ifndef ZeeTreeHZZ4l_h
#define ZeeTreeHZZ4l_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

namespace citana
{
  class ZeeTreeHZZ4l {
  public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    // Declaration of leaf types
    Float_t         intimeSimVertices;
    Float_t         l1bdtID;
    Float_t         l1bdtIso;
    Float_t         l1classification;
    Float_t         l1combDetIso;
    Float_t         l1deta;
    Float_t         l1dphi;
    Float_t         l1ecalIso;
    Float_t         l1esOverSC;
    Float_t         l1eta;
    Float_t         l1etawidth;
    Float_t         l1gsfp;
    Float_t         l1gsfpmode;
    Float_t         l1hOverE;
    Float_t         l1hcalIso;
    Float_t         l1p;
    Float_t         l1pdgId;
    Float_t         l1pfIsoChHad04;
    Float_t         l1pfIsoComb04EACorr;
    Float_t         l1pfIsoEAtot;
    Float_t         l1pfIsoNHad04;
    Float_t         l1pfIsoPhoton04;
    Float_t         l1phi;
    Float_t         l1phiwidth;
    Float_t         l1pt;
    Float_t         l1r9;
    Float_t         l1scenergy;
    Float_t         l1sceta;
    Float_t         l1scphi;
    Float_t         l1sietaieta;
    Float_t         l1sip;
    Float_t         l1tkIso;
    Float_t         l2bdtID;
    Float_t         l2bdtIso;
    Float_t         l2classification;
    Float_t         l2combDetIso;
    Float_t         l2deta;
    Float_t         l2dphi;
    Float_t         l2ecalIso;
    Float_t         l2esOverSC;
    Float_t         l2eta;
    Float_t         l2etawidth;
    Float_t         l2gsfp;
    Float_t         l2gsfpmode;
    Float_t         l2hOverE;
    Float_t         l2hcalIso;
    Float_t         l2p;
    Float_t         l2pdgId;
    Float_t         l2pfIsoChHad04;
    Float_t         l2pfIsoComb04EACorr;
    Float_t         l2pfIsoEAtot;
    Float_t         l2pfIsoNHad04;
    Float_t         l2pfIsoPhoton04;
    Float_t         l2phi;
    Float_t         l2phiwidth;
    Float_t         l2pt;
    Float_t         l2r9;
    Float_t         l2scenergy;
    Float_t         l2sceta;
    Float_t         l2scphi;
    Float_t         l2sietaieta;
    Float_t         l2sip;
    Float_t         l2tkIso;
    Float_t         massErr;
    Float_t         numTrueInteractions;
    Float_t         nvtx;
    Float_t         pairCombDetIso;
    Float_t         phoeta;
    Float_t         phophi;
    Float_t         phopt;
    Float_t         rho;
    Float_t         rhoAA;
    Float_t         rhoMuCentral;
    Float_t         zeta;
    Float_t         zmass;
    Float_t         zmll;
    Float_t         zphi;
    Float_t         zpt;
    Int_t           fsr;
    Int_t           l1ConvR;
    Int_t           l1idMVA;
    Int_t           l1idNew;
    Int_t           l1idPRL;
    Int_t           l1isoMVA;
    Int_t           l2ConvR;
    Int_t           l2idMVA;
    Int_t           l2idNew;
    Int_t           l2idPRL;
    Int_t           l2isoMVA;
    UInt_t          run;
    UInt_t          lumi;
    UInt_t          event;

    // List of branches
    TBranch        *b_intimeSimVertices;   //!
    TBranch        *b_l1bdtID;   //!
    TBranch        *b_l1bdtIso;   //!
    TBranch        *b_l1classification;   //!
    TBranch        *b_l1combDetIso;   //!
    TBranch        *b_l1deta;   //!
    TBranch        *b_l1dphi;   //!
    TBranch        *b_l1ecalIso;   //!
    TBranch        *b_l1esOverSC;   //!
    TBranch        *b_l1eta;   //!
    TBranch        *b_l1etawidth;   //!
    TBranch        *b_l1gsfp;   //!
    TBranch        *b_l1gsfpmode;   //!
    TBranch        *b_l1hOverE;   //!
    TBranch        *b_l1hcalIso;   //!
    TBranch        *b_l1p;   //!
    TBranch        *b_l1pdgId;   //!
    TBranch        *b_l1pfIsoChHad04;   //!
    TBranch        *b_l1pfIsoComb04EACorr;   //!
    TBranch        *b_l1pfIsoEAtot;   //!
    TBranch        *b_l1pfIsoNHad04;   //!
    TBranch        *b_l1pfIsoPhoton04;   //!
    TBranch        *b_l1phi;   //!
    TBranch        *b_l1phiwidth;   //!
    TBranch        *b_l1pt;   //!
    TBranch        *b_l1r9;   //!
    TBranch        *b_l1scenergy;   //!
    TBranch        *b_l1sceta;   //!
    TBranch        *b_l1scphi;   //!
    TBranch        *b_l1sietaieta;   //!
    TBranch        *b_l1sip;   //!
    TBranch        *b_l1tkIso;   //!
    TBranch        *b_l2bdtID;   //!
    TBranch        *b_l2bdtIso;   //!
    TBranch        *b_l2classification;   //!
    TBranch        *b_l2combDetIso;   //!
    TBranch        *b_l2deta;   //!
    TBranch        *b_l2dphi;   //!
    TBranch        *b_l2ecalIso;   //!
    TBranch        *b_l2esOverSC;   //!
    TBranch        *b_l2eta;   //!
    TBranch        *b_l2etawidth;   //!
    TBranch        *b_l2gsfp;   //!
    TBranch        *b_l2gsfpmode;   //!
    TBranch        *b_l2hOverE;   //!
    TBranch        *b_l2hcalIso;   //!
    TBranch        *b_l2p;   //!
    TBranch        *b_l2pdgId;   //!
    TBranch        *b_l2pfIsoChHad04;   //!
    TBranch        *b_l2pfIsoComb04EACorr;   //!
    TBranch        *b_l2pfIsoEAtot;   //!
    TBranch        *b_l2pfIsoNHad04;   //!
    TBranch        *b_l2pfIsoPhoton04;   //!
    TBranch        *b_l2phi;   //!
    TBranch        *b_l2phiwidth;   //!
    TBranch        *b_l2pt;   //!
    TBranch        *b_l2r9;   //!
    TBranch        *b_l2scenergy;   //!
    TBranch        *b_l2sceta;   //!
    TBranch        *b_l2scphi;   //!
    TBranch        *b_l2sietaieta;   //!
    TBranch        *b_l2sip;   //!
    TBranch        *b_l2tkIso;   //!
    TBranch        *b_massErr;   //!
    TBranch        *b_numTrueInteractions;   //!
    TBranch        *b_nvtx;   //!
    TBranch        *b_pairCombDetIso;   //!
    TBranch        *b_phoeta;   //!
    TBranch        *b_phophi;   //!
    TBranch        *b_phopt;   //!
    TBranch        *b_rho;   //!
    TBranch        *b_rhoAA;   //!
    TBranch        *b_rhoMuCentral;   //!
    TBranch        *b_zeta;   //!
    TBranch        *b_zmass;   //!
    TBranch        *b_zmll;   //!
    TBranch        *b_zphi;   //!
    TBranch        *b_zpt;   //!
    TBranch        *b_fsr;   //!
    TBranch        *b_l1ConvR;   //!
    TBranch        *b_l1idMVA;   //!
    TBranch        *b_l1idNew;   //!
    TBranch        *b_l1idPRL;   //!
    TBranch        *b_l1isoMVA;   //!
    TBranch        *b_l2ConvR;   //!
    TBranch        *b_l2idMVA;   //!
    TBranch        *b_l2idNew;   //!
    TBranch        *b_l2idPRL;   //!
    TBranch        *b_l2isoMVA;   //!
    TBranch        *b_run;   //!
    TBranch        *b_lumi;   //!
    TBranch        *b_event;   //!

  ZeeTreeHZZ4l(TTree *tree) : fChain(0) 
    {
      // if parameter tree is not specified (or zero), connect the file
      // used to generate this class and read the Tree.
      if (tree == 0) {
	TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("zEE_lineshape_Data_2012.root");
	if (!f || !f->IsOpen()) {
	  f = new TFile("zEE_lineshape_Data_2012.root");
	}
	TDirectory * dir = (TDirectory*)f->Get("zEE_lineshape_Data_2012.root:/zeetree");
	dir->GetObject("probe_tree",tree);

      }
      Init(tree);
    }

  ~ZeeTreeHZZ4l()
    {
      if (!fChain) return;
      delete fChain->GetCurrentFile();
    }

  Int_t GetEntry(Long64_t entry)
    {
      // Read contents of entry.
      if (!fChain) return 0;
      return fChain->GetEntry(entry);
    }
  Long64_t LoadTree(Long64_t entry)
    {
      // Set the environment to read one entry
      if (!fChain) return -5;
      Long64_t centry = fChain->LoadTree(entry);
      if (centry < 0) return centry;
      if (fChain->GetTreeNumber() != fCurrent) {
	fCurrent = fChain->GetTreeNumber();
	Notify();
      }
      return centry;
    }

  void Init(TTree *tree)
  {
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("intimeSimVertices", &intimeSimVertices, &b_intimeSimVertices);
    fChain->SetBranchAddress("l1bdtID", &l1bdtID, &b_l1bdtID);
    fChain->SetBranchAddress("l1bdtIso", &l1bdtIso, &b_l1bdtIso);
    fChain->SetBranchAddress("l1classification", &l1classification, &b_l1classification);
    fChain->SetBranchAddress("l1combDetIso", &l1combDetIso, &b_l1combDetIso);
    fChain->SetBranchAddress("l1deta", &l1deta, &b_l1deta);
    fChain->SetBranchAddress("l1dphi", &l1dphi, &b_l1dphi);
    fChain->SetBranchAddress("l1ecalIso", &l1ecalIso, &b_l1ecalIso);
    fChain->SetBranchAddress("l1esOverSC", &l1esOverSC, &b_l1esOverSC);
    fChain->SetBranchAddress("l1eta", &l1eta, &b_l1eta);
    fChain->SetBranchAddress("l1etawidth", &l1etawidth, &b_l1etawidth);
    fChain->SetBranchAddress("l1gsfp", &l1gsfp, &b_l1gsfp);
    fChain->SetBranchAddress("l1gsfpmode", &l1gsfpmode, &b_l1gsfpmode);
    fChain->SetBranchAddress("l1hOverE", &l1hOverE, &b_l1hOverE);
    fChain->SetBranchAddress("l1hcalIso", &l1hcalIso, &b_l1hcalIso);
    fChain->SetBranchAddress("l1p", &l1p, &b_l1p);
    fChain->SetBranchAddress("l1pdgId", &l1pdgId, &b_l1pdgId);
    fChain->SetBranchAddress("l1pfIsoChHad04", &l1pfIsoChHad04, &b_l1pfIsoChHad04);
    fChain->SetBranchAddress("l1pfIsoComb04EACorr", &l1pfIsoComb04EACorr, &b_l1pfIsoComb04EACorr);
    fChain->SetBranchAddress("l1pfIsoEAtot", &l1pfIsoEAtot, &b_l1pfIsoEAtot);
    fChain->SetBranchAddress("l1pfIsoNHad04", &l1pfIsoNHad04, &b_l1pfIsoNHad04);
    fChain->SetBranchAddress("l1pfIsoPhoton04", &l1pfIsoPhoton04, &b_l1pfIsoPhoton04);
    fChain->SetBranchAddress("l1phi", &l1phi, &b_l1phi);
    fChain->SetBranchAddress("l1phiwidth", &l1phiwidth, &b_l1phiwidth);
    fChain->SetBranchAddress("l1pt", &l1pt, &b_l1pt);
    fChain->SetBranchAddress("l1r9", &l1r9, &b_l1r9);
    fChain->SetBranchAddress("l1scenergy", &l1scenergy, &b_l1scenergy);
    fChain->SetBranchAddress("l1sceta", &l1sceta, &b_l1sceta);
    fChain->SetBranchAddress("l1scphi", &l1scphi, &b_l1scphi);
    fChain->SetBranchAddress("l1sietaieta", &l1sietaieta, &b_l1sietaieta);
    fChain->SetBranchAddress("l1sip", &l1sip, &b_l1sip);
    fChain->SetBranchAddress("l1tkIso", &l1tkIso, &b_l1tkIso);
    fChain->SetBranchAddress("l2bdtID", &l2bdtID, &b_l2bdtID);
    fChain->SetBranchAddress("l2bdtIso", &l2bdtIso, &b_l2bdtIso);
    fChain->SetBranchAddress("l2classification", &l2classification, &b_l2classification);
    fChain->SetBranchAddress("l2combDetIso", &l2combDetIso, &b_l2combDetIso);
    fChain->SetBranchAddress("l2deta", &l2deta, &b_l2deta);
    fChain->SetBranchAddress("l2dphi", &l2dphi, &b_l2dphi);
    fChain->SetBranchAddress("l2ecalIso", &l2ecalIso, &b_l2ecalIso);
    fChain->SetBranchAddress("l2esOverSC", &l2esOverSC, &b_l2esOverSC);
    fChain->SetBranchAddress("l2eta", &l2eta, &b_l2eta);
    fChain->SetBranchAddress("l2etawidth", &l2etawidth, &b_l2etawidth);
    fChain->SetBranchAddress("l2gsfp", &l2gsfp, &b_l2gsfp);
    fChain->SetBranchAddress("l2gsfpmode", &l2gsfpmode, &b_l2gsfpmode);
    fChain->SetBranchAddress("l2hOverE", &l2hOverE, &b_l2hOverE);
    fChain->SetBranchAddress("l2hcalIso", &l2hcalIso, &b_l2hcalIso);
    fChain->SetBranchAddress("l2p", &l2p, &b_l2p);
    fChain->SetBranchAddress("l2pdgId", &l2pdgId, &b_l2pdgId);
    fChain->SetBranchAddress("l2pfIsoChHad04", &l2pfIsoChHad04, &b_l2pfIsoChHad04);
    fChain->SetBranchAddress("l2pfIsoComb04EACorr", &l2pfIsoComb04EACorr, &b_l2pfIsoComb04EACorr);
    fChain->SetBranchAddress("l2pfIsoEAtot", &l2pfIsoEAtot, &b_l2pfIsoEAtot);
    fChain->SetBranchAddress("l2pfIsoNHad04", &l2pfIsoNHad04, &b_l2pfIsoNHad04);
    fChain->SetBranchAddress("l2pfIsoPhoton04", &l2pfIsoPhoton04, &b_l2pfIsoPhoton04);
    fChain->SetBranchAddress("l2phi", &l2phi, &b_l2phi);
    fChain->SetBranchAddress("l2phiwidth", &l2phiwidth, &b_l2phiwidth);
    fChain->SetBranchAddress("l2pt", &l2pt, &b_l2pt);
    fChain->SetBranchAddress("l2r9", &l2r9, &b_l2r9);
    fChain->SetBranchAddress("l2scenergy", &l2scenergy, &b_l2scenergy);
    fChain->SetBranchAddress("l2sceta", &l2sceta, &b_l2sceta);
    fChain->SetBranchAddress("l2scphi", &l2scphi, &b_l2scphi);
    fChain->SetBranchAddress("l2sietaieta", &l2sietaieta, &b_l2sietaieta);
    fChain->SetBranchAddress("l2sip", &l2sip, &b_l2sip);
    fChain->SetBranchAddress("l2tkIso", &l2tkIso, &b_l2tkIso);
    fChain->SetBranchAddress("massErr", &massErr, &b_massErr);
    fChain->SetBranchAddress("numTrueInteractions", &numTrueInteractions, &b_numTrueInteractions);
    fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
    fChain->SetBranchAddress("pairCombDetIso", &pairCombDetIso, &b_pairCombDetIso);
    fChain->SetBranchAddress("phoeta", &phoeta, &b_phoeta);
    fChain->SetBranchAddress("phophi", &phophi, &b_phophi);
    fChain->SetBranchAddress("phopt", &phopt, &b_phopt);
    fChain->SetBranchAddress("rho", &rho, &b_rho);
    fChain->SetBranchAddress("rhoAA", &rhoAA, &b_rhoAA);
    fChain->SetBranchAddress("rhoMuCentral", &rhoMuCentral, &b_rhoMuCentral);
    fChain->SetBranchAddress("zeta", &zeta, &b_zeta);
    fChain->SetBranchAddress("zmass", &zmass, &b_zmass);
    fChain->SetBranchAddress("zmll", &zmll, &b_zmll);
    fChain->SetBranchAddress("zphi", &zphi, &b_zphi);
    fChain->SetBranchAddress("zpt", &zpt, &b_zpt);
    fChain->SetBranchAddress("fsr", &fsr, &b_fsr);
    fChain->SetBranchAddress("l1ConvR", &l1ConvR, &b_l1ConvR);
    fChain->SetBranchAddress("l1idMVA", &l1idMVA, &b_l1idMVA);
    fChain->SetBranchAddress("l1idNew", &l1idNew, &b_l1idNew);
    fChain->SetBranchAddress("l1idPRL", &l1idPRL, &b_l1idPRL);
    fChain->SetBranchAddress("l1isoMVA", &l1isoMVA, &b_l1isoMVA);
    fChain->SetBranchAddress("l2ConvR", &l2ConvR, &b_l2ConvR);
    fChain->SetBranchAddress("l2idMVA", &l2idMVA, &b_l2idMVA);
    fChain->SetBranchAddress("l2idNew", &l2idNew, &b_l2idNew);
    fChain->SetBranchAddress("l2idPRL", &l2idPRL, &b_l2idPRL);
    fChain->SetBranchAddress("l2isoMVA", &l2isoMVA, &b_l2isoMVA);
    fChain->SetBranchAddress("run", &run, &b_run);
    fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
    fChain->SetBranchAddress("event", &event, &b_event);
    Notify();
  }

  Bool_t Notify()
    {
      // The Notify() function is called when a new file is opened. This
      // can be either for a new TTree in a TChain or when when a new TTree
      // is started when using PROOF. It is normally not necessary to make changes
      // to the generated code, but the routine can be extended by the
      // user if needed. The return value is currently not used.

      return kTRUE;
    }

  void Show(Long64_t entry)
  {
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
  }
  Int_t Cut(Long64_t entry)
    {
      // This function may be called from Loop.
      // returns  1 if entry is accepted.
      // returns -1 otherwise.
      return 1;
    }

  };
}
#endif

