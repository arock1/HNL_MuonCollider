#pragma once

#include <iostream>

using namespace std;

#include <vector>

#include "FinalStatesClass.C"
#include "Rtypes.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "classes/DelphesClasses.h"

iFinalStates FindFinalStatesIndex(TClonesArray* branchTrack, TClonesArray* branchVLC1Jet, TClonesArray* branchVLC2Jet) {
    iFinalStates iFinalStatesIndexes;
    Int_t nTracks = branchTrack->GetEntries();
    Int_t nVLC1Jets = branchVLC1Jet->GetEntries();
    Int_t nVLC2Jets = branchVLC2Jet->GetEntries();

    // require there be 1 or 2 jets, rom W>qq
    if (nVLC1Jets != 1 && nVLC2Jets != 2) return iFinalStatesIndexes;

    Track* track1;
    Float_t mLep;
    for (Int_t it = 0; it < nTracks; it++) {
        track1 = (Track*)branchTrack->At(it);
        if (abs(track1->PID) == 13) {  // muon
            iFinalStatesIndexes.mLeps.push_back(0.1056583745);
            iFinalStatesIndexes.iLeps.push_back(it);
        } else if (abs(track1->PID) == 11) {  // electron
            iFinalStatesIndexes.mLeps.push_back(0.00051099894615109989461);
            iFinalStatesIndexes.iLeps.push_back(it);
        }
    }

    Jet* jet1;
    Jet* jet2;
    // required 2 jets event
    if (nVLC1Jets == 1) {
        iFinalStatesIndexes.iJets.push_back(0);  // save jet index
    } else if (nVLC2Jets == 2) {
        jet1 = (Jet*)branchVLC2Jet->At(0);
        jet2 = (Jet*)branchVLC2Jet->At(1);
        if (jet1->PT > jet2->PT) {
            iFinalStatesIndexes.iJets.push_back(0);  // save jet index
            iFinalStatesIndexes.iJets.push_back(1);  // save jet index
        } else {
            iFinalStatesIndexes.iJets.push_back(1);  // save jet index
            iFinalStatesIndexes.iJets.push_back(0);  // save jet index
        }
    } else {
        return iFinalStatesIndexes;
    }

    iFinalStatesIndexes.foundAll = 1;  // flag to process

    return iFinalStatesIndexes;
}
