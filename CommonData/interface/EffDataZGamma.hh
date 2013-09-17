#ifndef EFF_DATAZGAMMA_HH
#define EFF_DATAZGAMMA_HH

struct EffDataZGamma
{
    Float_t mass, pt, eta, phi, weight;
    Int_t q;
    UInt_t npv, npu, pass;
    UInt_t runNum, lumiSec, evtNum;
    Float_t rho;
    Float_t phosigiEtaiEta, phoR9;
    Float_t massll, masstagpho, massprobepho, drprobepho, drtagpho, phoet, phoeta;
    Float_t ptprobepho;
    Float_t dphitagtoprobepho;
    Float_t ptllpho;        
    Bool_t phopasspixelveto;
    Bool_t phoisreal, tagisreal, probeisreal;

};

// "mass/F:pt:eta:phi:weight:q/I:npv/i:npu:pass:runNum:lumiSec:evtNum:rho/F"
// "mass/F:pt:eta:phi:weight:q/I:npv/i:npu:pass:runNum:lumiSec:evtNum:rho/F:phosigiEtaiEta:phoR9:massll:masstagpho:massprobepho:drprobepho:drtagpho:phoet:phoeta:ptprobepho:dphitagtoprobepho:ptllpho:phopasspixelveto/O:phoisreal:tagisreal:probeisreal"
#endif
