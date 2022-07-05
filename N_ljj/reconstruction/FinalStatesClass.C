// Define the class of the targeted final states, and used later in reconstruction

#pragma once

#include <vector>
using namespace std;
#include "Rtypes.h"
#include "TLorentzVector.h"

struct iFinalStates {
    Int_t foundAll = 0;
    // Int_t iW;
    vector<TLorentzVector> iLeps;
    vector<Int_t> iElectronIndeces;
    vector<Int_t> iMuonIndeces;
    vector<Int_t> iLepCharges;
    vector<TLorentzVector> i1Jets;
    vector<TLorentzVector> i2Jets;
};
