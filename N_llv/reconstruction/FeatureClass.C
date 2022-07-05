// Define the class of Features, to be stored and used in BDT code
// specific to the llv final states

#include "Rtypes.h"
class Features {
public:
    Int_t iEvt;  // event index

    // lepton means the one decayed from HNL
    Float_t ptLep;   // pt of lepton
    Float_t etaLep;  // eta of lepton
    Float_t phiLep;  // phi of lepton
    Float_t ELep;    // energy of lepton
    Float_t pxLep;
    Float_t pyLep;
    Float_t pzLep;

    // typeLep: is the one from HNL
    Int_t typeLep;  // type of lepton (11: electron, 13: muon)
    // typeLep2 is the one from W bosn (but seems like not so useful?)
    Int_t typeLepW;  // type of lepton (11: electron, 13: muon)

    Float_t ptLepW;   // pt of of the lepton from W bosn (with the nu which gives heavier HNL mass)
    Float_t etaLepW;  // eta of the lepton from the W bosonn (with the nu which gives heavier HNL mass)
    Float_t phiLepW;  // phi of the lepton from the W bosonn (with the nu which gives heavier HNL mass)
    Float_t ELepW;    // E of the lepton from the W boson (with the nu which gives heavier HNL mass)

    Float_t DeltaPhillv1;  // Delta phi of lepton and lv (from W) (v which gives heavier HNL mass)
    Float_t DeltaRlv1;     // Delta R between lepton and nuetrino from the W bossn (v which gives heavier HNL mass)
    Float_t DeltaRllv1;    // Delta R between lepton and lv (from W) (v which gives heavier HNL mass)
    Float_t ptlv1;         // pt of W (v which gives heavier HNL mass)
    Float_t etalv1;        // eta of W (v which gives heavier HNL mass)
    Float_t philv1;        // phi of W (v which gives heavier HNL mass)
    Float_t Elv1;          // mass of W (v which gives heavier HNL mass)

    Float_t DeltaPhillv2;  // Delta phi of lepton and lv (from W)
    Float_t DeltaRlv2;     // Delta R between lepton and nuetrino from the W bossn
    Float_t DeltaRllv2;    // Delta R between lepton and lv (from W)
    Float_t ptlv2;         // pt of W
    Float_t etalv2;        // eta of W
    Float_t philv2;        // phi of W
    Float_t Elv2;          // mass of W

    Float_t mN11;  // mass of N (pair with the larger pT lepton and the solution that give heavier HNL mass)
    Float_t mN12;  // mass of N (pair with the larger pT lepton and the solution that give lighter HNL mass)
    Float_t mN21;  // mass of N (pair with the smaller pT lepton and the solution that give heavier HNL mass)
    Float_t mN22;  // mass of N (pair with the smaller pT lepton and the solution that give ligther HNL mass)

    // following N would represent only pairing up with the larger pt lepton
    Float_t ptN1;   // pt of N (heavier N mass solution)
    Float_t etaN1;  // eta of N (heavier N mass solution)
    Float_t phiN1;  // phi of N (heavier N mass solution)
    Float_t pxN1;   // pz of N (heavier N mass solution)
    Float_t pyN1;   // px of N (heavier N mass solution)
    Float_t pzN1;   // py of N (heavier N mass solution)

    Float_t ptN2;   // pt of N (ligther N mass solution)
    Float_t etaN2;  // eta of N (ligther N mass solution)
    Float_t phiN2;  // phi of N (ligther N mass solution)
    Float_t pxN2;   // pz of N (ligther N mass solution)
    Float_t pyN2;   // px of N (ligther N mass solution)
    Float_t pzN2;   // py of N (ligther N mass solution)

    Int_t chargeLep;  // charge of lepton (lepton means the one decayed from HNL)

    Float_t MET;
    Float_t pTheta;  // proxy of W polarisation from 2008.04318

    Float_t MinPtLep;
    Float_t DeltaPhiNMET;

    Float_t ptLepTrue;
    Float_t etaLepTrue;
    Float_t phiLepTrue;
    Float_t ELepTrue;
    Float_t ptJet1True;
    Float_t etaJet1True;
    Float_t phiJet1True;
    Float_t EJet1True;
    Float_t ptJet2True;
    Float_t etaJet2True;
    Float_t phiJet2True;
    Float_t EJet2True;
    Float_t DeltaRjjTrue;
    Float_t DeltaRjjlTrue;
    Float_t ptJJTrue;
    Float_t etaJJTrue;
    Float_t phiJJTrue;
    Float_t mJJTrue;
    Float_t mNTrue;
    Float_t ptNTrue;
    Float_t etaNTrue;
    Float_t phiNTrue;
    Float_t ENTrue;
    Float_t pzNTrue;
    Int_t chargeLepTrue;
    Int_t typeLepTrue;

    Float_t pxLepTrue;
    Float_t pyLepTrue;
    Float_t pzLepTrue;
    Float_t pxNTrue;
    Float_t pyNTrue;

    Float_t pxNu11;
    Float_t pyNu11;
    Float_t pzNu11;
    Float_t pxNu12;
    Float_t pyNu12;
    Float_t pzNu12;

    Float_t pxNu21;
    Float_t pyNu21;
    Float_t pzNu21;
    Float_t pxNu22;
    Float_t pyNu22;
    Float_t pzNu22;

    Float_t pxNuTrue;
    Float_t pyNuTrue;
    Float_t pzNuTrue;

    Float_t imbalance11;
    Float_t imbalance12;
    Float_t imbalance21;
    Float_t imbalance22;

    Int_t bkgPIDV1 = 99999;
    Int_t bkgPIDV2 = 99999;
};
