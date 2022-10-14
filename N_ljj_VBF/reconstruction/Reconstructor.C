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

TParticlePDG* tauPDG = TDatabasePDG::Instance()->GetParticle(15);
extern Float_t mTauPDG;
Float_t mTauPDG = tauPDG->Mass();

TParticlePDG* WPDG = TDatabasePDG::Instance()->GetParticle(24);
extern Float_t mWPDG, widthWPDG;
Float_t mWPDG = WPDG->Mass();
Float_t widthWPDG = WPDG->Width();

/*
iFinalStates FindFinalStatesIndex(TClonesArray* branchElectron, TClonesArray* branchMuon, TClonesArray* branchVLC1Jet, TClonesArray* branchVLC2Jet) {
    iFinalStates iFinalStatesIndexes;
    Int_t nElectrons = branchElectron->GetEntries();
    Int_t nMuons = branchMuon->GetEntries();
    Int_t nVLC1Jets = branchVLC1Jet->GetEntries();
    Int_t nVLC2Jets = branchVLC2Jet->GetEntries();

    //  require at least two leptons
    // if (nElectrons + nMuons < 1 && tauTagged == 0) return iFinalStatesIndexes;
    if (nElectrons + nMuons < 2) return iFinalStatesIndexes;
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
        iFinalStatesIndexes.typeLeps.push_back(11);
    }
    Muon* muon1;
    for (Int_t imu = 0; imu < nMuons; imu++) {
        muon1 = (Muon*)branchMuon->At(imu);
        lep.SetPtEtaPhiM(muon1->PT, muon1->Eta, muon1->Phi, 0);  // assumed lepton to be massless
        iFinalStatesIndexes.iLeps.push_back(lep);
        iFinalStatesIndexes.iMuonIndeces.push_back(imu);
        iFinalStatesIndexes.iLepCharges.push_back(muon1->Charge);
        iFinalStatesIndexes.typeLeps.push_back(13);
    }

    // =============================================================
    TLorentzVector jet1;
    TLorentzVector jet21, jet22, jet2;
    if (nVLC1Jets == 1) {
        Jet* jetVLC1 = (Jet*)branchVLC1Jet->At(0);
        jet1.SetPtEtaPhiM(jetVLC1->PT, jetVLC1->Eta, jetVLC1->Phi, jetVLC1->Mass);
        iFinalStatesIndexes.i1Jets.push_back(jet1);  // save jet
        iFinalStatesIndexes.tau1 = jetVLC1->Tau[0];
        iFinalStatesIndexes.tau2 = jetVLC1->Tau[1];
    }
    if (nVLC2Jets == 2) {
        Jet* jetVLC21 = (Jet*)branchVLC2Jet->At(0);
        Jet* jetVLC22 = (Jet*)branchVLC2Jet->At(1);
        jet21.SetPtEtaPhiM(jetVLC21->PT, jetVLC21->Eta, jetVLC21->Phi, jetVLC21->Mass);
        jet22.SetPtEtaPhiM(jetVLC22->PT, jetVLC22->Eta, jetVLC22->Phi, jetVLC22->Mass);
        if (jet21.Pt() > jet22.Pt()) {
            iFinalStatesIndexes.i2Jets.push_back(jet21);  // save jet
            iFinalStatesIndexes.i2Jets.push_back(jet22);  // save jet
        } else {
            iFinalStatesIndexes.i2Jets.push_back(jet22);  // save jet
            iFinalStatesIndexes.i2Jets.push_back(jet21);  // save jet
        }
    }

    iFinalStatesIndexes.foundAll = 1;  // flag to process

    return iFinalStatesIndexes;
}
*/


iFinalStates FindFinalStatesLepsIndex(TClonesArray* branchElectron, TClonesArray* branchMuon, TClonesArray* branchVLC1Jet, TClonesArray* branchVLC2Jet) {
    iFinalStates iFinalStatesIndexes;
    Int_t nElectrons = branchElectron->GetEntries();
    Int_t nMuons = branchMuon->GetEntries();

    //  require at least two leptons
    // if (nElectrons + nMuons < 1 && tauTagged == 0) return iFinalStatesIndexes;
    if (nElectrons + nMuons < 2) return iFinalStatesIndexes;

    TLorentzVector lep;
    Electron* electron1;
    for (Int_t iel = 0; iel < nElectrons; iel++) {
        electron1 = (Electron*)branchElectron->At(iel);
        lep.SetPtEtaPhiM(electron1->PT, electron1->Eta, electron1->Phi, 0);  // assumed lepton to be massless
        iFinalStatesIndexes.iLeps.push_back(lep);
        iFinalStatesIndexes.iElectronIndeces.push_back(iel);
        iFinalStatesIndexes.iLepCharges.push_back(electron1->Charge);
        iFinalStatesIndexes.typeLeps.push_back(11);
    }
    Muon* muon1;
    for (Int_t imu = 0; imu < nMuons; imu++) {
        muon1 = (Muon*)branchMuon->At(imu);
        lep.SetPtEtaPhiM(muon1->PT, muon1->Eta, muon1->Phi, 0);  // assumed lepton to be massless
        iFinalStatesIndexes.iLeps.push_back(lep);
        iFinalStatesIndexes.iMuonIndeces.push_back(imu);
        iFinalStatesIndexes.iLepCharges.push_back(muon1->Charge);
        iFinalStatesIndexes.typeLeps.push_back(13);
    }

    iFinalStatesIndexes.foundAll = 1;  // flag to process

    return iFinalStatesIndexes;
}

iFinalStates FindFinalStatesJetsIndex(TClonesArray* branchElectron, TClonesArray* branchMuon, TClonesArray* branchVLC1Jet, TClonesArray* branchVLC2Jet) {
    iFinalStates iFinalStatesIndexes;
    Int_t nVLC1Jets = branchVLC1Jet->GetEntries();
    Int_t nVLC2Jets = branchVLC2Jet->GetEntries();
    if (nVLC1Jets != 1 && nVLC2Jets != 2) return iFinalStatesIndexes;

    TLorentzVector jet1;
    TLorentzVector jet21, jet22, jet2;
    if (nVLC1Jets == 1) {
        Jet* jetVLC1 = (Jet*)branchVLC1Jet->At(0);
        jet1.SetPtEtaPhiM(jetVLC1->PT, jetVLC1->Eta, jetVLC1->Phi, jetVLC1->Mass);
        iFinalStatesIndexes.i1Jets.push_back(jet1);  // save jet
        iFinalStatesIndexes.tau1 = jetVLC1->Tau[0];
        iFinalStatesIndexes.tau2 = jetVLC1->Tau[1];
    }
    if (nVLC2Jets == 2) {
        Jet* jetVLC21 = (Jet*)branchVLC2Jet->At(0);
        Jet* jetVLC22 = (Jet*)branchVLC2Jet->At(1);
        jet21.SetPtEtaPhiM(jetVLC21->PT, jetVLC21->Eta, jetVLC21->Phi, jetVLC21->Mass);
        jet22.SetPtEtaPhiM(jetVLC22->PT, jetVLC22->Eta, jetVLC22->Phi, jetVLC22->Mass);
        if (jet21.Pt() > jet22.Pt()) {
            iFinalStatesIndexes.i2Jets.push_back(jet21);  // save jet
            iFinalStatesIndexes.i2Jets.push_back(jet22);  // save jet
        } else {
            iFinalStatesIndexes.i2Jets.push_back(jet22);  // save jet
            iFinalStatesIndexes.i2Jets.push_back(jet21);  // save jet
        }
    }

    iFinalStatesIndexes.foundAll = 1;  // flag to process

    return iFinalStatesIndexes;
}
TLorentzVector getWJet(iFinalStates iFS) {
    // if there are both one jet and two jets, then pick the one with mass closest to W boson
    // if there is only one(two) jet, then assign that as an W jet
    TLorentzVector jet1, jet21, jet22, jet2, jj;
    if (iFS.i1Jets.size() == 1 && iFS.i2Jets.size() == 2) {
        jet1 = iFS.i1Jets[0];
        jet21 = iFS.i2Jets[0];
        jet22 = iFS.i2Jets[1];
        jet2 = jet21 + jet22;
        if (abs(jet1.M() - mWPDG) < abs(jet2.M() - mWPDG)) {
            jj = jet1;
        } else {
            jj = jet2;
        }
    } else if (iFS.i1Jets.size() == 1 && iFS.i2Jets.size() == 0) {
        jj = iFS.i1Jets[0];
    } else if (iFS.i1Jets.size() == 0 && iFS.i2Jets.size() == 2) {
        jet21 = iFS.i2Jets[0];
        jet22 = iFS.i2Jets[1];
        jj = jet21 + jet22;
    }

    return jj;
}

TLorentzVector deduceTau(TLorentzVector tauJet, TLorentzVector jj) {
    TLorentzVector lepSum = tauJet + jj;
    TLorentzVector lepi = tauJet;

    // solving quadratic eqn (solve for pz(nu))
    Float_t A = lepSum.Px() * lepSum.Px() + lepSum.Py() * lepSum.Py();
    Float_t B = -(lepi.Px() * lepSum.Px() + lepi.Py() * lepSum.Py()) + mTauPDG * mTauPDG / 2;

    Float_t a = (lepi.Pz() * lepi.Pz() - lepi.E() * lepi.E());
    Float_t b = 2 * B * lepi.Pz();
    Float_t c = B * B - A * lepi.E() * lepi.E();
    Float_t pz1 = (-b + pow(b * b - 4 * a * c, 0.5)) / (2 * a);
    Float_t pz2 = (-b - pow(b * b - 4 * a * c, 0.5)) / (2 * a);

    // if no solution, just use -b/2a
    if (isnan(pz1) && isnan(pz2)) {
        pz1 = -b / (2 * a);
        pz2 = -b / (2 * a);
    }

    // defining the 4 momentum
    TLorentzVector nu1_, nu2_, tau1_, tau2_;  // two possible solution of the nu, and the corresponding taus
    Float_t E1 = pow(lepSum.Px() * lepSum.Px() + lepSum.Py() * lepSum.Py() + pz1 * pz1, 0.5);
    Float_t E2 = pow(lepSum.Px() * lepSum.Px() + lepSum.Py() * lepSum.Py() + pz2 * pz2, 0.5);
    nu1_.SetPxPyPzE(-lepSum.Px(), -lepSum.Py(), pz1, E1);
    nu2_.SetPxPyPzE(-lepSum.Px(), -lepSum.Py(), pz2, E2);

    tau1_ = nu1_ + lepi;
    tau2_ = nu2_ + lepi;

    // if (pz1 != pz2) cout << "/////////////////\n";

    TLorentzVector tau_ = (tau1_ * 0.5 + tau2_ * 0.5);

    return tau_;
}

// since in VBF, there are at least two leptons
// pick the one with smallest P as the one decayed from HNL
void getLep(iFinalStates iFS, TLorentzVector jj, pair<Int_t, Int_t>* typeLeps, pair<Int_t, Int_t>* chargeLeps, pair<TLorentzVector, TLorentzVector>* leps) {
    // loop over all the leptons
    pair<Int_t, Int_t> chargeLeps_;
    pair<Int_t, Int_t> typeLeps_;
    pair<TLorentzVector, TLorentzVector> leps_;
    
    vector<TLorentzVector> twoLeps = {}; 
    vector<Int_t> twoCharges = {};
    vector<Int_t> twoTypes = {};
    if (iFS.iLepCharges.size() == 2) { // if there are exactly two leptons, then just append them
        twoLeps.push_back(iFS.iLeps[0]); twoLeps.push_back(iFS.iLeps[1]);
        twoCharges.push_back(iFS.iLepCharges[0]); twoCharges.push_back(iFS.iLepCharges[1]);
        twoTypes.push_back(iFS.typeLeps[0]); twoTypes.push_back(iFS.typeLeps[1]);
        
    } else  {   // if more than two leptons, use the two with min |eta|
        Int_t indexMinFirst = 99999; // index for the smallest eta
        for (Int_t _ = 0; _ < 2; _++) { // looking foe two leptons
            // init the index and min eta for this round
            Int_t indexMin = 99999;
            Float_t etaMin = 99999;
            for (Int_t il = 0; il < iFS.iLepCharges.size(); il++) {
                if (il == indexMinFirst) continue; // skip the smallest eta one 
                TLorentzVector lepI = iFS.iLeps[il];
                if (abs(lepI.Eta()) < etaMin) { // find the smallest 
                    etaMin = abs(lepI.Eta());
                    indexMin = il;
                }
            }
            if (_ == 0) indexMinFirst = indexMin; 
            TLorentzVector lepI = iFS.iLeps[indexMin];
            twoLeps.push_back(iFS.iLeps[indexMin]);
            twoCharges.push_back(iFS.iLepCharges[indexMin]);
            twoTypes.push_back(iFS.typeLeps[indexMin]);
        }
    }


    // store the lepton according to their pt 
    // larger pT be the first 
    TLorentzVector lep1, lep2;
    lep1 = twoLeps[0]; lep2 = twoLeps[1];
    if (lep1.Pt() > lep2.Pt()) {
        leps_.first = lep1; leps_.second = lep2; 
        typeLeps_.first = twoTypes[0]; typeLeps_.second= twoTypes[1]; 
        chargeLeps_.first = twoCharges[0]; chargeLeps_.second= twoCharges[1]; 
    } else {
        leps_.first = lep2; leps_.second = lep1; 
        typeLeps_.first = twoTypes[1]; typeLeps_.second= twoTypes[0]; 
        chargeLeps_.first = twoCharges[1]; chargeLeps_.second= twoCharges[0]; 
    }

    *typeLeps = typeLeps_;
    *chargeLeps = chargeLeps_;
    *leps = leps_;
}

