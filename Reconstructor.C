#pragma once

#include <iostream>

using namespace std;

#include <vector>

#include "FinalStatesClass.C"
#include "Rtypes.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "classes/DelphesClasses.h"

iFinalStates FindFinalStatesIndex(TClonesArray* branchTrack, TClonesArray* branchJet) {
    iFinalStates iFinalStatesIndexes;
    Int_t nTracks = branchTrack->GetEntries();
    Int_t nJets = branchJet->GetEntries();

    Track* track1;
    Jet* jet1;
    Jet* jet2;

    if (nJets != 2) return iFinalStatesIndexes;

    Int_t iLep = 99999;
    Int_t nEle = 0;
    Int_t nMu = 0;
    for (Int_t it = 0; it < nTracks; it++) {
        track1 = (Track*)branchTrack->At(it);
        if (abs(track1->PID) == 13) {
            nMu += 1;
            iLep = it;
        } else if (abs(track1->PID) == 11) {
            nEle += 1;
            iLep = it;
        }
    }
    if (nMu + nEle != 1) return iFinalStatesIndexes;

    TLorentzVector j1, j2, lep, N;
    jet1 = (Jet*)branchJet->At(0);
    jet2 = (Jet*)branchJet->At(1);
    track1 = (Track*)branchTrack->At(iLep);
    j1.SetPtEtaPhiM(jet1->PT, jet1->Eta, jet1->Phi, jet1->Mass);
    j2.SetPtEtaPhiM(jet2->PT, jet2->Eta, jet2->Phi, jet2->Mass);

    Float_t mLep;
    if (nMu == 1) mLep = 0.1056583745;
    if (nEle == 1) mLep = 0.51099894615109989461;

    lep.SetPtEtaPhiM(track1->PT, track1->Eta, track1->Phi, mLep);

    if (not(lep.Px() * j1.Px() + lep.Py() * j1.Py() + lep.Pz() * j1.Pz() >= 0)) return iFinalStatesIndexes;
    if (not(lep.Px() * j2.Px() + lep.Py() * j2.Py() + lep.Pz() * j2.Pz() >= 0)) return iFinalStatesIndexes;
    if (not(j2.Px() * j1.Px() + j2.Py() * j1.Py() + j2.Pz() * j1.Pz() >= 0)) return iFinalStatesIndexes;

    iFinalStatesIndexes.iLep = iLep;
    iFinalStatesIndexes.iJet1 = 0;
    iFinalStatesIndexes.iJet2 = 0;
    iFinalStatesIndexes.mLep = mLep;
    iFinalStatesIndexes.foundAll = 1;

    return iFinalStatesIndexes;
}
