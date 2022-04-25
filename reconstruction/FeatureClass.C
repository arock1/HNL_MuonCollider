
#include "Rtypes.h"
class Features {
public:
    Int_t iEvt;           // event index
    Float_t ptLep;        // pt of lepton
    Float_t etaLep;       // eta of lepton
    Float_t phiLep;       // phi of lepton
    Float_t ELep;         // energy of lepton
    Int_t typeLep;        // type of lepton (11: electron, 13: muon)
    Int_t typeLep2;       // type of lepton (11: electron, 13: muon)
    Float_t ptJet1;       // pt of Jet1, if any
    Float_t etaJet1;      // eta of Jet1, if any
    Float_t phiJet1;      // phi of Jet1, if any
    Float_t EJet1;        // E of Jet1, if any
    Float_t ptJet2;       // pt of Jet2, if any
    Float_t etaJet2;      // eta of Jet2, if any
    Float_t phiJet2;      // phi of Jet2, if any
    Float_t EJet2;        // E of Jet2, if any
    Float_t DeltaPhijjl;  // Delta phi of lepton and JJ (from W)
    Float_t DeltaRjj;     // Delta R between 2 Jets, if any
    Float_t DeltaRjjl;    // Delta R between lepton and JJ (Jet(s) from W)
    Int_t nJets;          // number of Jets
    Float_t ptJJ;         // pt of W
    Float_t etaJJ;        // eta of W
    Float_t phiJJ;        // phi of W
    Float_t mJJ;          // mass of W
    Float_t mN;           // mass of N
    Float_t ptN;          // pt of N
    Float_t etaN;         // eta of N
    Float_t phiN;         // phi of N
    Float_t pzN;          // pz of N
    Int_t chargeLep;      // charge of lepton
    Float_t MET;

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

    Int_t bkgPIDV1 = 99999;
    Int_t bkgPIDV2 = 99999;
};
