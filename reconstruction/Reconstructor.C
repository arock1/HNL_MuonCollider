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

    Track* track1;
    Jet* jet1;
    Jet* jet2;

    if (nVLC1Jets != 1 && nVLC2Jets != 2) return iFinalStatesIndexes;

    Int_t iLep = 99999;
    Int_t nEle = 0;
    Int_t nMu = 0;
    Float_t ptMax = 0;
    for (Int_t it = 0; it < nTracks; it++) {
        track1 = (Track*)branchTrack->At(it);
        // muon
        if (abs(track1->PID) == 13) {
            nMu += 1;
            if (ptMax < track1->PT) {
                ptMax = track1->PT;
                iLep = it;
            }
            // electron
        } else if (abs(track1->PID) == 11) {
            nEle += 1;
            if (ptMax < track1->PT) {
                ptMax = track1->PT;
                iLep = it;
            }
        }
    }
    //  required at least lepton, either muon or electron
    if (nMu + nEle < 1) return iFinalStatesIndexes;
    track1 = (Track*)branchTrack->At(iLep);

    // required 2 jets event
    if (nVLC1Jets == 1) {
        iFinalStatesIndexes.nJets = nVLC1Jets;
        iFinalStatesIndexes.iJet1 = 0;  // save jet index
    } else if (nVLC2Jets == 2) {
        iFinalStatesIndexes.nJets = nVLC2Jets;
        jet1 = (Jet*)branchVLC2Jet->At(0);
        jet2 = (Jet*)branchVLC2Jet->At(1);
        if (jet1->PT > jet2->PT) {
            iFinalStatesIndexes.iJet1 = 0;  // save jet index
            iFinalStatesIndexes.iJet2 = 1;  // save jet index
        } else {
            iFinalStatesIndexes.iJet1 = 1;  // save jet index
            iFinalStatesIndexes.iJet2 = 0;  // save jet index
        }
    } else {
        iFinalStatesIndexes.nJets = 0;
        return iFinalStatesIndexes;
    }

    // set the lepton mass, check if it is muon or electron
    Float_t mLep;
    if (nMu == 1) mLep = 0.1056583745;
    if (nEle == 1) mLep = 0.00051099894615109989461;

    // TLorentzVector j1, j2, lep;
    // j1.SetPtEtaPhiM(jet1->PT, jet1->Eta, jet1->Phi, jet1->Mass);
    // j2.SetPtEtaPhiM(jet2->PT, jet2->Eta, jet2->Phi, jet2->Mass);
    // lep.SetPtEtaPhiM(track1->PT, track1->Eta, track1->Phi, mLep);

    iFinalStatesIndexes.iLep = iLep;   // save lepton index
    iFinalStatesIndexes.mLep = mLep;   // save lepton mass
    iFinalStatesIndexes.foundAll = 1;  // flag to process

    return iFinalStatesIndexes;
}
