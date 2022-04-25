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

iFinalStates FindFinalStatesIndex(TClonesArray* branchElectron, TClonesArray* branchMuon, TClonesArray* branchVLC1Jet, TClonesArray* branchVLC2Jet) {
    iFinalStates iFinalStatesIndexes;
    Int_t nElectrons = branchElectron->GetEntries();
    Int_t nMuons = branchMuon->GetEntries();
    Int_t nVLC1Jets = branchVLC1Jet->GetEntries();
    Int_t nVLC2Jets = branchVLC2Jet->GetEntries();

    // require at least one lepton
    if (nElectrons + nMuons < 1) return iFinalStatesIndexes;
    // require there be 1 or 2 jets, from W>qq
    if (nVLC1Jets != 1 && nVLC2Jets != 2) return iFinalStatesIndexes;

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

    TLorentzVector jet;
    Jet* jet1;
    Jet* jet2;
    if (nVLC1Jets == 1) {
        jet1 = (Jet*)branchVLC1Jet->At(0);
        jet.SetPtEtaPhiM(jet1->PT, jet1->Eta, jet1->Phi, jet1->Mass);
        iFinalStatesIndexes.iJets.push_back(jet);  // save jet
    } else if (nVLC2Jets == 2) {
        jet1 = (Jet*)branchVLC2Jet->At(0);
        jet2 = (Jet*)branchVLC2Jet->At(1);
        if (jet1->PT > jet2->PT) {
            jet.SetPtEtaPhiM(jet1->PT, jet1->Eta, jet1->Phi, jet1->Mass);
            iFinalStatesIndexes.iJets.push_back(jet);  // save jet
            jet.SetPtEtaPhiM(jet2->PT, jet2->Eta, jet2->Phi, jet2->Mass);
            iFinalStatesIndexes.iJets.push_back(jet);  // save jet
        } else {
            jet.SetPtEtaPhiM(jet2->PT, jet2->Eta, jet2->Phi, jet2->Mass);
            iFinalStatesIndexes.iJets.push_back(jet);  // save jet
            jet.SetPtEtaPhiM(jet1->PT, jet1->Eta, jet1->Phi, jet1->Mass);
            iFinalStatesIndexes.iJets.push_back(jet);  // save jet
        }
    } else {
        return iFinalStatesIndexes;
    }

    iFinalStatesIndexes.foundAll = 1;  // flag to process

    return iFinalStatesIndexes;
}
