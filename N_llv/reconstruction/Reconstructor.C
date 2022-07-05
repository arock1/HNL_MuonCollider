// To find the targeted final states

#pragma once

#include <iostream>

using namespace std;

#include <vector>

#include "FinalStatesClass.C"
#include "Rtypes.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "classes/DelphesClasses.h"

TParticlePDG* electronPDG = TDatabasePDG::Instance()->GetParticle(11);
extern Float_t mElectronPDG;
Float_t mElectronPDG = electronPDG->Mass();

TParticlePDG* muonPDG = TDatabasePDG::Instance()->GetParticle(13);
extern Float_t mMuonPDG;
Float_t mMuonPDG = muonPDG->Mass();

TParticlePDG* WPDG = TDatabasePDG::Instance()->GetParticle(24);
extern Float_t mWPDG, widthWPDG;
Float_t mWPDG = WPDG->Mass();
Float_t widthWPDG = WPDG->Width();

iFinalStates FindFinalStatesIndex(TClonesArray* branchElectron, TClonesArray* branchMuon) {
    iFinalStates iFinalStatesIndexes;
    Int_t nElectrons = branchElectron->GetEntries();
    Int_t nMuons = branchMuon->GetEntries();

    // require at least two leptons
    if (nElectrons + nMuons < 2) return iFinalStatesIndexes;

    TLorentzVector lep;
    Electron* electron1;
    for (Int_t iel = 0; iel < nElectrons; iel++) {
        electron1 = (Electron*)branchElectron->At(iel);
        lep.SetPtEtaPhiM(electron1->PT, electron1->Eta, electron1->Phi, 0);  // assumed lepton to be massless
        iFinalStatesIndexes.iLeps.push_back(lep);
        iFinalStatesIndexes.iElectronIndeces.push_back(iel);
        iFinalStatesIndexes.iLepCharges.push_back(electron1->Charge);
    }
    Muon* muon1;
    for (Int_t imu = 0; imu < nMuons; imu++) {
        muon1 = (Muon*)branchMuon->At(imu);
        lep.SetPtEtaPhiM(muon1->PT, muon1->Eta, muon1->Phi, 0);  // assumed lepton to be massless
        iFinalStatesIndexes.iLeps.push_back(lep);
        iFinalStatesIndexes.iMuonIndeces.push_back(imu);
        iFinalStatesIndexes.iLepCharges.push_back(muon1->Charge);
    }

    iFinalStatesIndexes.foundAll = 1;  // flag to process

    return iFinalStatesIndexes;
}
