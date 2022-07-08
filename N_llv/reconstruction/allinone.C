//
//
//
// ===========================================================================
// || This Code Does:                                                       ||
// || 1. Classifiy the event type (Signal / What Background type)           ||
// ||    and store the corresponding truth level info.                      ||
// || 2. Find the targeted final state in detector level:                   ||
// ||       a. Two leptons (electron or muon)                               ||
// || 3. Reconstruct W jet, by deducing the nutrino info (N>(W>l nu) l)     ||
// ||       a. Since N general has pT~0, assume pT(nu) = -pT(2l)            ||
// ||       b. Then solve the quadratic eqn of (W>l nu)                     ||
// ||          Find what is the pz(nu) s.t. (p(l) + p(nu))^2 = m_W^2        ||
// || 6. Store features for later BDT analysis                              ||
// ===========================================================================
//

#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>
using namespace std;

#include "FeatureClass.C"
#include "FinalStatesClass.C"
#include "Geometry.C"
#include "Reconstructor.C"
#include "Rtypes.h"
#include "SignalBackgroundClassifier.C"
#include "TChain.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

R__ADD_LIBRARY_PATH($DELPHES)
R__LOAD_LIBRARY(libDelphes)

void allinone(
    const string type,
    const Bool_t save = false,
    Int_t num_test = 0,
    Int_t print_detail = 1) {
    // formmating the input output files
    string inputFile_st;
    string outputFile_st;

    Int_t foundFiles;
    // foundFiles = getFileNames(type, &inputFile_st, &outputFile_st);

    // inputFile_st = "../data/detector/signal_tchannel_E-3TeV_N-1TeV_llv.root";
    // outputFile_st = "../data/features/signal_tchannel_E-3TeV_N-1TeV_reco_llv.root";
    //   outputFile_st = "../data/features/signal_tchannel_E-3TeV_N-1TeV_reco_llv_cheat.root";

    inputFile_st = "../data/detector/signal_tchannel_E-3TeV_N-1TeV_llv_Dirac.root";
    outputFile_st = "../data/features/signal_tchannel_E-3TeV_N-1TeV_reco_llv_Dirac.root";

    const char* inputFile = inputFile_st.c_str();
    const char* outputFile = outputFile_st.c_str();
    cout << "\nReading: " << inputFile << "\n";
    cout << "Expected outputFile: " << outputFile << "\n\n";
    if (gSystem->AccessPathName(inputFile)) {
        cout << "inputFile: " << inputFile << " does not exist" << endl;
        return;
    }

    if (not save) outputFile = "../data/_dummy.root";

    // Load lib, and read data
    gSystem->Load("libDelphes");
    TChain chain("Delphes");
    chain.Add(inputFile);
    ExRootTreeReader* treeReader = new ExRootTreeReader(&chain);
    Int_t numberOfEntries = treeReader->GetEntries();

    // clone the branches
    TClonesArray* branchParticle = treeReader->UseBranch("Particle");
    TClonesArray* branchElectron = treeReader->UseBranch("Electron");
    TClonesArray* branchMuon = treeReader->UseBranch("Muon");
    TClonesArray* branchVLC1Jet = treeReader->UseBranch("VLCjetR12N1");
    TClonesArray* branchVLC2Jet = treeReader->UseBranch("VLCjetR02N2");
    TClonesArray* branchMET = treeReader->UseBranch("MissingET");

    // book feature storing tree and file
    TFile fea(outputFile, "recreate");
    TTree tr("t", "features");
    Features* features = new Features;
    tr.Branch("features", &features);

    if (num_test == 0) num_test = numberOfEntries;

    // tracing the number of events are each cut
    Int_t nEv = 0;      // # of events in truth level
    Int_t nFS = 0;      // # of events having identified final states
    Int_t nLepEta = 0;  // # of events after lep eta cut
    Int_t nLepPt = 0;   // # of events after lep pt cut
    Int_t nJJPt = 0;    // # of events after jj pt cut
    Int_t nWM = 0;      // # of events after W reco

    Float_t lepEtaCut = 2.5;  // lep eta cut
    Float_t lepPtCut = 100;   // lep pt cut

    // loop the events
    for (Int_t i_en = 0; i_en < num_test; i_en++) {
        // progress
        if ((i_en % 1000) == 0) cout << "Reconstruction Progress: " << i_en << "/" << numberOfEntries << "\r";
        cout.flush();
        // cout << "\nEvent: " << i_en << endl;

        treeReader->ReadEntry(i_en);  // reading the entry

        //========================================================================
        //=======================   Classifiy Event Type   =======================
        //========================================================================
        iFinalStates iFSTrue;                     // indeces of the true level final states
        TLorentzVector lepTrue, lepTrue2, NTrue;  // truth level lepton, jets, and HNL
        Int_t typeLepTrue;
        Int_t passing = 0;
        if (type.at(0) == 'i' || type.at(0) == 's' || type.at(0) == 't') {  // search for signals
            if (type.at(2) == 'M') {                                        // Majorana
                passing = ClassifySingal(branchParticle, &iFSTrue, &NTrue, 9900012);
            } else if (type.at(2) == 'D') {  // Dirac
                passing = ClassifySingal(branchParticle, &iFSTrue, &NTrue, 9990012);
            }
        }

        if (passing == 0) continue;
        nEv += 1;

        //========================================================================
        //=======================      Reconstruction      =======================
        //========================================================================
        iFinalStates iFS;
        // finding final states: 2 leptons
        iFS = FindFinalStatesIndex(branchElectron, branchMuon);

        if (iFS.foundAll == 0) continue;
        nFS += 1;

        TLorentzVector lep1, lep2, lepSum;  // two tagged leptons, the sum of them (used to deduce the pT(nu))
        lep1 = iFS.iLeps[0];
        lep2 = iFS.iLeps[1];
        lepSum = lep1 + lep2;
        // lepSum = lep1 + lep2 - NTrue;

        // lepton eta cut
        if ((abs(lep1.Eta()) >= lepEtaCut)) continue;
        if ((abs(lep2.Eta()) >= lepEtaCut)) continue;
        nLepEta += 1;

        // lepton pt cut
        if (not(lep1.Pt() >= lepPtCut)) continue;
        if (not(lep2.Pt() >= lepPtCut)) continue;
        nLepPt += 1;

        Int_t haveWboson = 0;               // flag for reconstructing a W boson
        vector<TLorentzVector> sortedLeps;  // sorted the lepton by their pT
        if (iFS.iLeps[0].Pt() > iFS.iLeps[1].Pt()) {
            // if (iFS.iLeps[0].E() > iFS.iLeps[1].E()) {
            sortedLeps.push_back(iFS.iLeps[0]);
            sortedLeps.push_back(iFS.iLeps[1]);
        } else {
            sortedLeps.push_back(iFS.iLeps[1]);
            sortedLeps.push_back(iFS.iLeps[0]);
        }

        // cout << " lep1:    " << sortedLeps[0].Px() / sortedLeps[0].Pt() << "; " << sortedLeps[0].Py() / sortedLeps[0].Pt() << "\n";
        // cout << " lep2:    " << sortedLeps[1].Px() / sortedLeps[1].Pt() << "; " << sortedLeps[1].Py() / sortedLeps[1].Pt() << "\n";
        // cout << " NTrue:   " << NTrue.Px() / NTrue.Pt() << "; " << NTrue.Py() / NTrue.Pt() << "\n";
        // cout << endl;

        TLorentzVector W1, W2;
        for (Int_t iLeps = 0; iLeps < 2; iLeps++) {  // loop over the tagged leptons
            TLorentzVector lepi = sortedLeps[iLeps];

            // solving quadratic eqn (solve for pz(nu))
            Float_t A = lepSum.Px() * lepSum.Px() + lepSum.Py() * lepSum.Py();
            Float_t B = -(lepi.Px() * lepSum.Px() + lepi.Py() * lepSum.Py()) + mWPDG * mWPDG / 2;

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
            haveWboson = 1;

            // defining the 4 momentum
            TLorentzVector nu1_, nu2_, W1_, W2_;  // two possible solution of the nu, and the corresponding W bosons
            Float_t E1 = pow(lepSum.Px() * lepSum.Px() + lepSum.Py() * lepSum.Py() + pz1 * pz1, 0.5);
            Float_t E2 = pow(lepSum.Px() * lepSum.Px() + lepSum.Py() * lepSum.Py() + pz2 * pz2, 0.5);
            nu1_.SetPxPyPzE(-lepSum.Px(), -lepSum.Py(), pz1, E1);
            nu2_.SetPxPyPzE(-lepSum.Px(), -lepSum.Py(), pz2, E2);

            W1_ = nu1_ + lepi;
            W2_ = nu2_ + lepi;

            // reconstructing the HNL
            TLorentzVector N1_, N2_;
            N1_ = nu1_ + lep1 + lep2;
            N2_ = nu2_ + lep1 + lep2;

            // define m1 as the larger mN solution, and m2 as the smaller mN solution
            Float_t m1, m2;

            Float_t pxnu1, pynu1, pznu1;
            Float_t pxnu2, pynu2, pznu2;
            TLorentzVector nu1, nu2, N1, N2, W1, W2;
            if (N1_.M() > N2_.M()) {
                N1 = N1_;
                nu1 = nu1_;
                W1 = W1_;
                N2 = N2_;
                nu2 = nu2_;
                W2 = W2_;
            } else {
                N1 = N2_;
                nu1 = nu2_;
                W1 = W2_;
                N2 = N1_;
                nu2 = nu1_;
                W2 = W1_;
            }
            if (iLeps == 0) {  // if the current lepton is the larger pT one
                features->mN11 = N1_.M();
                features->pxNu11 = nu1_.Px();
                features->pyNu11 = nu1_.Py();
                features->pzNu11 = nu1_.Pz();
                features->imbalance11 = abs(nu1_.E() - lepi.E()) / W1_.P();
                features->mN12 = N2_.M();
                features->pxNu12 = nu2_.Px();
                features->pyNu12 = nu2_.Py();
                features->pzNu12 = nu2_.Pz();
                features->imbalance12 = abs(nu2_.E() - lepi.E()) / W2_.P();
                W1 = W1_;
                W2 = W2_;
            } else if (iLeps == 1) {  // if the current lepton is the smaller pT one
                features->mN21 = N1_.M();
                features->pxNu21 = nu1_.Px();
                features->pyNu21 = nu1_.Py();
                features->pzNu21 = nu1_.Pz();
                features->imbalance21 = abs(nu1_.E() - lepi.E()) / W1_.P();
                features->mN22 = N2_.M();
                features->pxNu22 = nu2_.Px();
                features->pyNu22 = nu2_.Py();
                features->pzNu22 = nu2_.Pz();
                features->imbalance22 = abs(nu2_.E() - lepi.E()) / W2_.P();
            }
        }
        nWM += 1;
        TLorentzVector lep, lepW;
        lep = sortedLeps[0];
        lepW = sortedLeps[1];

        //========================================================================
        //=======================         Features         =======================
        //========================================================================

        // Float_t DeltaRjj, DeltaRjjl;
        // if (iFS.i2Jets.size() == 2) {
        // Float_t jjEtaDiff = jet1.Eta() - jet2.Eta();
        // Float_t jjPhiDiff = deltaPhi(jet1.Phi(), jet2.Phi());
        // DeltaRjj = pow(jjEtaDiff * jjEtaDiff + jjPhiDiff * jjPhiDiff, 0.5);
        //}

        // Float_t jjlEtaDiff = jj.Eta() - lep.Eta();
        // Float_t jjlPhiDiff = deltaPhi(jj.Phi(), lep.Phi());
        // DeltaRjjl = pow(jjlEtaDiff * jjlEtaDiff + jjlPhiDiff * jjlPhiDiff, 0.5);

        // Float_t DeltaRjjTrue, DeltaRjjlTrue;
        // Float_t jjEtaDiffTrue = jet1True.Eta() - jet2True.Eta();
        // Float_t jjPhiDiffTrue = deltaPhi(jet1True.Phi(), jet2True.Phi());
        // DeltaRjjTrue = pow(jjEtaDiffTrue * jjEtaDiffTrue + jjPhiDiffTrue * jjPhiDiffTrue, 0.5);
        // Float_t jjlEtaDiffTrue = jjTrue.Eta() - lepTrue.Eta();
        // Float_t jjlPhiDiffTrue = deltaPhi(jjTrue.Phi(), lepTrue.Phi());
        // DeltaRjjlTrue = pow(jjlEtaDiffTrue * jjlEtaDiffTrue + jjlPhiDiffTrue * jjlPhiDiffTrue, 0.5);

        features->iEvt = i_en;
        features->pxNuTrue = iFSTrue.iNu[0].Px();
        features->pyNuTrue = iFSTrue.iNu[0].Py();
        features->pzNuTrue = iFSTrue.iNu[0].Pz();

        features->ptLep = lep.Pt();
        features->etaLep = lep.Eta();
        features->phiLep = lep.Phi();
        features->ELep = lep.E();
        features->pxLep = lep.Px();
        features->pyLep = lep.Py();
        features->pzLep = lep.Pz();

        // features->DeltaPhijjl = jjlPhiDiff;
        // features->DeltaRjj = DeltaRjj;
        // features->DeltaRjjl = DeltaRjjl;
        // if (iFS.i2Jets.size() == 2) features->DeltaRjj = DeltaRjj;
        // if (iFS.i2Jets.size() == 2) features->DeltaRjjl = DeltaRjjl;

        // if (iFS.i2Jets.size() == 2) {
        // features->pTheta = abs(jet21.E() - jet22.E()) / jj.P();
        //}

        //// if (iFS.i2Jets.size() != 2) cout << " .....................: ";
        // if (iFS.i2Jets.size() == 2) {
        // features->EJet1 = jet21.E();
        // features->EJet2 = jet22.E();
        //}
        //// features->nJets = iFS.iJets.size();
        // features->mJJ = jj.M();
        // features->ptJJ = jj.Pt();
        // features->etaJJ = jj.Eta();
        // features->phiJJ = jj.Phi();
        // features->mN = N.M();
        // features->ptN = N.Pt();
        // features->etaN = N.Eta();
        // features->phiN = N.Phi();
        // features->pxN = N.Px();
        // features->pyN = N.Py();
        // features->pzN = N.Pz();

        // features->ptLepTrue = lepTrue.Pt();
        // features->etaLepTrue = lepTrue.Eta();
        // features->phiLepTrue = lepTrue.Phi();
        // features->ELepTrue = lepTrue.E();

        // features->ptJet1True = jet1True.Pt();
        // features->etaJet1True = jet1True.Eta();
        // features->phiJet1True = jet1True.Phi();
        // features->EJet1True = jet1True.E();

        // features->ptJet2True = jet2True.Pt();
        // features->etaJet2True = jet2True.Eta();
        // features->phiJet2True = jet2True.Phi();
        // features->EJet2True = jet2True.E();

        // features->DeltaRjjTrue = DeltaRjjTrue;
        // features->DeltaRjjlTrue = DeltaRjjlTrue;
        // features->mJJTrue = jjTrue.M();
        // features->ptJJTrue = jjTrue.Pt();
        // features->etaJJTrue = jjTrue.Eta();
        // features->phiJJTrue = jjTrue.Phi();
        // features->mNTrue = NTrue.M();
        // features->ptNTrue = NTrue.Pt();
        // features->etaNTrue = NTrue.Eta();
        // features->phiNTrue = NTrue.Phi();
        // features->ENTrue = NTrue.E();
        // features->pzNTrue = NTrue.Pz();

        // features->chargeLep = chargeLep;
        // if (type.at(0) != 'b') {
        // features->chargeLepTrue = iFSTrue.iLepCharges[0];
        //}

        // features->typeLep = typeLep;
        // features->typeLepTrue = typeLepTrue;

        // Preliminary: Try to find the second lepton for the 2VBF case, for LNV
        // finding the muon that give cloeset mass to W?
        // TLorentzVector lep2;
        // Float_t mDiffWMin = 99999;
        // Int_t chargeLep2;
        // Int_t lep2Index = 99999;
        // Int_t typeLep2;
        // TLorentzVector lN;  // second lepton + N
        // for (Int_t il = 0; il < iFS.iLepCharges.size(); il++) {
        // if (il == lepIndex) continue;
        // TLorentzVector lepI = iFS.iLeps[il];
        //// eta cut
        // if (not(abs(lepI.Eta()) <= lepEtaCut)) continue;

        // lN = lep2 + N;
        // if (abs(mWPDG - lN.M()) < mDiffWMin) {
        // chargeLep2 = iFS.iLepCharges[il];
        // mDiffWMin = abs(mWPDG - lN.M());
        // if (not(lN.M() >= WMLowCut && lN.M() <= WMHighCut)) continue;
        // lep2 = lepI;
        //}
        //}
        // features->typeLep2 = typeLep2;

        // Int_t nMET = branchMET->GetEntries();
        // if (nMET != 0) {
        // MissingET* met = (MissingET*)branchMET->At(0);
        // features->MET = met->MET;
        // features->DeltaPhiNMET = deltaPhi(met->Phi, N.Phi());
        //}

        // Float_t minptlep = 99999;
        // for (Int_t il = 0; il < iFS.iLepCharges.size(); il++) {
        // TLorentzVector lepI = iFS.iLeps[il];
        // if (lepI.Pt() < minptlep) minptlep = lepI.Pt();
        //}
        // features->MinPtLep = minptlep;

        // passing = ClassifyiBackground(branchParticle, &bkgTypesReco);
        tr.Fill();

        // cout << "M: " << lepTrue.M() << "; " << lep.M() << endl;
        // cout << "E: " << lepTrue.E() << "; " << lep.E() << endl;
        // cout << "P: " << lepTrue.P() << "; " << lep.P() << endl;
        // cout << "PT: " << lepTrue.Pt() << "; " << lep.Pt() << endl;
        // cout << endl;
    }

    cout << "Reconstruction Progress: " << numberOfEntries << "/" << numberOfEntries << "\n\n";

    if (print_detail) {
        cout << "\t---------------------------------------------------------------------------------" << endl;
        cout << "\t==============================       Cut eff.      ==============================\n";
        cout << "\t---------------------------------------------------------------------------------" << endl;

        cout << "\t# of targeted events in truth level:\t\t\t" << nEv << endl;
        cout << "\t# after identified final states:\t\t\t" << nFS << "\t(" << 100 * float(nFS) / float(nEv) << "%)" << endl;
        cout << "\t# after lepton eta cut \t(eta(l) <= " << float(lepEtaCut) << "):\t\t" << nLepEta << "\t(" << 100 * float(nLepEta) / float(nFS) << "%)" << endl;
        cout << "\t# after lepton pt cut \t(pT(l) >= " << float(lepPtCut) << "):\t\t\t" << nLepPt << "\t(" << 100 * float(nLepPt) / float(nLepEta) << "%)" << endl;
        // cout << "\t# after jj pt cut \t(pT(jj) >= " << float(jjPtCut) << "):\t\t" << nJJPt << "\t(" << 100 * float(nJJPt) / float(nLepPt) << "%)" << endl;
        cout << "\t# after W reco:  \t\t\t\t\t" << nWM << "\t(" << 100 * float(nWM) / float(nLepPt) << "%)" << endl;
        //  cout << "# after jet pt cut (pT(j1, j2) >= " << float(jetPtCut) << "):\t\t" << nJetPt << "\t(" << 100 * float(nJetPt) / float(nLepPt) << "%)" << endl;
        cout << "\t---------------------------------------------------------------------------------" << endl;
        cout << "\t\t\t\t\t\tTotal eff.:\t" << nWM << " / " << nEv << "\t(" << 100 * float(nWM) / float(nEv) << " %) " << endl;
        cout << "\t---------------------------------------------------------------------------------" << endl;
        cout << "\n\n\n";
    }

    tr.Write();
    fea.Close();
}
