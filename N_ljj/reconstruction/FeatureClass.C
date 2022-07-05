// Define the class of Features, to be stored and used in BDT code

#include "Rtypes.h"
class Features {
public:
    Int_t iEvt;      // event index
    Float_t ptLep;   // pt of lepton
    Float_t etaLep;  // eta of lepton
    Float_t phiLep;  // phi of lepton
    Float_t ELep;    // energy of lepton
    Float_t pxLep;
    Float_t pyLep;
    Float_t pzLep;
    Int_t chargeLep;    // charge of lepton
    Int_t lepisMu = 0;  // type of lepton (1: muon; 0: eletron)

    Float_t mJJ;    // mass of W
    Float_t ptJJ;   // pt of W
    Float_t etaJJ;  // eta of W
    Float_t phiJJ;  // phi of W

    Float_t mN;    // mass of N
    Float_t ptN;   // pt of N
    Float_t etaN;  // eta of N
    Float_t phiN;  // phi of N
    Float_t pxN;   // pz of N
    Float_t pyN;   // px of N
    Float_t pzN;   // py of N

    Float_t DeltaRjj;         // Delta R between 2 Jets, if any
    Float_t EJet1 = 99999;    // E of Jet1, if any
    Float_t ptJet1 = 99999;   // pt of Jet1, if any
    Float_t etaJet1 = 99999;  // eta of Jet1, if any
    Float_t phiJet1 = 99999;  // phi of Jet1, if any
    Float_t EJet2 = 99999;    // E of Jet2, if any
    Float_t ptJet2 = 99999;   // pt of Jet2, if any
    Float_t etaJet2 = 99999;  // eta of Jet2, if any
    Float_t phiJet2 = 99999;  // phi of Jet2, if any
    Float_t DeltaRjjl;        // Delta R between lepton and JJ (Jet(s) from W)
    Float_t DeltaPhijjl;      // Delta phi of lepton and JJ (from W)

    Float_t ptFwMu = 99999;  // pt of the forward muon
    // =============================================================
    // extra, just for trial
    Float_t MET;
    Float_t DeltaPhiNMET;
};
