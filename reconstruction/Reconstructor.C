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

    // =============================================================
    //  trying new method: use the case when the JJ mass is closer to the W mass
    Int_t storeN = 0;
    TLorentzVector jet1;
    TLorentzVector jet21, jet22, jet2;
    if (nVLC1Jets == 1) {
        Jet* jetVLC1 = (Jet*)branchVLC1Jet->At(0);
        jet1.SetPtEtaPhiM(jetVLC1->PT, jetVLC1->Eta, jetVLC1->Phi, jetVLC1->Mass);
    }
    if (nVLC2Jets == 2) {
        Jet* jetVLC21 = (Jet*)branchVLC2Jet->At(0);
        Jet* jetVLC22 = (Jet*)branchVLC2Jet->At(1);
        jet21.SetPtEtaPhiM(jetVLC21->PT, jetVLC21->Eta, jetVLC21->Phi, jetVLC21->Mass);
        jet22.SetPtEtaPhiM(jetVLC22->PT, jetVLC22->Eta, jetVLC22->Phi, jetVLC22->Mass);
        jet2 = jet21 + jet22;
    }

    if (nVLC1Jets == 1 && nVLC2Jets == 2) {  // if there are both, then pick the one with mass closest to W boson
        if (abs(mWPDG - jet1.M()) < abs(mWPDG - jet2.M())) {
            storeN = 1;
        } else {
            storeN = 2;
        }
    } else if (nVLC1Jets == 1 && nVLC2Jets == 0) {  // if there is only 1 fat jet
        storeN = 1;
    } else if (nVLC2Jets == 2 && nVLC1Jets == 0) {  // if there is ony 2 jets
        storeN = 2;
    }

    if (storeN == 1) {
        iFinalStatesIndexes.iJets.push_back(jet1);  // save jet
    } else if (storeN == 2) {
        // if there are two jets, let the first one be the one with higher pT
        if (jet21.Pt() > jet22.Pt()) {
            iFinalStatesIndexes.iJets.push_back(jet21);  // save jet
            iFinalStatesIndexes.iJets.push_back(jet22);  // save jet
        } else {
            iFinalStatesIndexes.iJets.push_back(jet22);  // save jet
            iFinalStatesIndexes.iJets.push_back(jet21);  // save jet
        }
    } else {
        return iFinalStatesIndexes;
    }

    // =============================================================
    // TLorentzVector jet;
    // Jet* jet1;
    // Jet* jet2;
    // if (nVLC1Jets == 1) {
    // jet1 = (Jet*)branchVLC1Jet->At(0);
    // jet.SetPtEtaPhiM(jet1->PT, jet1->Eta, jet1->Phi, jet1->Mass);
    // iFinalStatesIndexes.iJets.push_back(jet);  // save jet
    //} else if (nVLC2Jets == 2) {
    // jet1 = (Jet*)branchVLC2Jet->At(0);
    // jet2 = (Jet*)branchVLC2Jet->At(1);
    // if (jet1->PT > jet2->PT) {
    // jet.SetPtEtaPhiM(jet1->PT, jet1->Eta, jet1->Phi, jet1->Mass);
    // iFinalStatesIndexes.iJets.push_back(jet);  // save jet
    // jet.SetPtEtaPhiM(jet2->PT, jet2->Eta, jet2->Phi, jet2->Mass);
    // iFinalStatesIndexes.iJets.push_back(jet);  // save jet
    //} else {
    // jet.SetPtEtaPhiM(jet2->PT, jet2->Eta, jet2->Phi, jet2->Mass);
    // iFinalStatesIndexes.iJets.push_back(jet);  // save jet
    // jet.SetPtEtaPhiM(jet1->PT, jet1->Eta, jet1->Phi, jet1->Mass);
    // iFinalStatesIndexes.iJets.push_back(jet);  // save jet
    //}
    //} else {
    // return iFinalStatesIndexes;
    //}

    iFinalStatesIndexes.foundAll = 1;  // flag to process

    return iFinalStatesIndexes;
}
