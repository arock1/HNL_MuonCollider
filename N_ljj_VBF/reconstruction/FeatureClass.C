// Define the class of Features, to be stored and used in BDT code

#include "Rtypes.h"
class Features {
public:
    Int_t iEvt;      // event index

    Float_t ptLep1;   // pt of lepton
    Float_t etaLep1;  // eta of lepton
    Float_t phiLep1;  // phi of lepton
    Float_t ELep1;    // energy of lepton
    Float_t pxLep1;
    Float_t pyLep1;
    Float_t pzLep1;
    Int_t chargeLep1;     // charge of lepton
    Int_t lep1isEle = 0;  // type of lepton (0: muon; 1: eletron; 0: tau)
    Int_t lep1isMu = 0;   // type of lepton (1: muon; 0: eletron; 0:tau)

    Float_t ptLep2;   // pt of lepton
    Float_t etaLep2;  // eta of lepton
    Float_t phiLep2;  // phi of lepton
    Float_t ELep2;    // energy of lepton
    Float_t pxLep2;
    Float_t pyLep2;
    Float_t pzLep2;
    Int_t chargeLep2;     // charge of lepton
    Int_t lep2isEle = 0;  // type of lepton (0: muon; 1: eletron; 0: tau)
    Int_t lep2isMu = 0;   // type of lepton (1: muon; 0: eletron; 0:tau)

    Float_t mJJ;    // mass of W
    Float_t EJJ;    // energy of W
    Float_t PJJ;    // energy of W
    Float_t ptJJ;   // pt of W
    Float_t etaJJ;  // eta of W
    Float_t phiJJ;  // phi of W
    Float_t pxJJ;   // phi of W
    Float_t pyJJ;   // phi of W
    Float_t pzJJ;   // phi of W

    Float_t mN1;    // mass of N
    Float_t EN1;    // mass of N
    Float_t PN1;    // mass of N
    Float_t ptN1;   // pt of N
    Float_t etaN1;  // eta of N
    Float_t phiN1;  // phi of N
    Float_t pxN1;   // pz of N
    Float_t pyN1;   // px of N
    Float_t pzN1;   // py of N

    Float_t mN2;    // mass of N
    Float_t EN2;    // mass of N
    Float_t PN2;    // mass of N
    Float_t ptN2;   // pt of N
    Float_t etaN2;  // eta of N
    Float_t phiN2;  // phi of N
    Float_t pxN2;   // pz of N
    Float_t pyN2;   // px of N
    Float_t pzN2;   // py of N

    Float_t mJJLep1Lep2;  

    Float_t DeltaRjj;         // Delta R between 2 Jets, if any
    Float_t EJet1 = 99999;    // E of Jet1, if any
    Float_t ptJet1 = 99999;   // pt of Jet1, if any
    Float_t etaJet1 = 99999;  // eta of Jet1, if any
    Float_t phiJet1 = 99999;  // phi of Jet1, if any
    Float_t EJet2 = 99999;    // E of Jet2, if any
    Float_t ptJet2 = 99999;   // pt of Jet2, if any
    Float_t etaJet2 = 99999;  // eta of Jet2, if any
    Float_t phiJet2 = 99999;  // phi of Jet2, if any
    Float_t DeltaRjjl1;        // Delta R between lepton and JJ (Jet(s) from W)
    Float_t DeltaPhijjl1;      // Delta phi of lepton and JJ (from W)
    Float_t DeltaRjjl2;        // Delta R between lepton and JJ (Jet(s) from W)
    Float_t DeltaPhijjl2;      // Delta phi of lepton and JJ (from W)

    Float_t ptFwMu = 99999;  // pt of the forward muon

    Float_t tau1 = 99999;
    Float_t tau2 = -99999;


    // =============================================================
    // Truth level
    Float_t ptNTrue;
    Float_t ENTrue;
    Float_t etaNTrue;

    Float_t ptJJTrue;
    Float_t EJJTrue;
    Float_t PJJTrue;
    Float_t etaJJTrue;
    Float_t phiJJTrue;

    Float_t ptLepTrue;
    Float_t ELepTrue;
    Float_t etaLepTrue;
    Float_t phiLepTrue;
    Float_t pxLepTrue;
    Float_t pyLepTrue;
    Float_t pzLepTrue;

    Float_t ptLep2True;
    // =============================================================
    // extra, just for trial
    Float_t MET;
    Float_t DeltaPhiNMET;

};
