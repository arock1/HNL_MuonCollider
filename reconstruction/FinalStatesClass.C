#pragma once

#include <vector>
using namespace std;
#include "Rtypes.h"

struct iFinalStates {
    Int_t foundAll = 0;
    Int_t iW;
    vector<Int_t> iLeps;
    vector<Float_t> mLeps;
    vector<Int_t> iJets;
};
