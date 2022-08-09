//
// ===========================================================================
// || This Code Does:                                                       ||
// || 1. Classifiy the event type (Signal / What Background type)           ||
// ||    and store the corresponding truth level info.                      ||
// || 2. Find the targeted final state in detector level:                   ||
// ||       a. At least two leptons (electron or muon)                      ||
// ||       b. Either one 1VLC jet or two 2VLC jets found                   ||
// || 3. Reconstruct W jet                                                  ||
// ||       a. If there is only (1VLC case found), use it as W jet          ||
// ||       b. If there is only (2VLC case found), use their sum as W jet   ||
// || 4. Get HNL (Need to have |eta|<threshold to simulate detector)        ||
// ||       a. Since there are at least two leptons,                        ||
// ||          store the two with smallest eta                              ||
// ||       b. Reconstruct two possible HNL with each lepton,               || 
// ||          and let BDT decide how to use the info                       ||
// || 5. Apply cuts                                                         ||
// ||       a. each lepton pT                                               ||
// ||       b. W jet pT                                                     ||
// ||       c. W jet mass                                                   ||
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

// Formatting input and output file name, based on the input "type" variable {{{
Int_t getFileNames(string type, string* inputFile_st_, string* outputFile_st_) {
    string inputFile_st;
    string outputFile_st;

    // expect input like: "s_M_3_1000" means signal, Majorana, sqrt{s}=3TeV, m_N=1000GeV
    if (type.at(0) == 's') {  // signals
        cout << "Processing Signal data" << endl;

        vector<string> words{};
        for (Int_t i = 0; i < 4; i++) {
            Int_t pos = type.find("_");
            words.push_back(type.substr(0, pos));
            // cout << " type: " << type.substr(0, pos) << "\n";
            type.erase(0, pos + 1);
        }

        // expect: "../data/VBF_Maj_E-3_m-1000.root"
        // inputFile_st = "../data/detector/VBF_";
        inputFile_st = "../data/detector/sig_";
        string fermionType;
        if (words[1] == 'M') {
            fermionType = "Maj";
        } else if (words[1] == 'D') {
            fermionType = "Dir";
        } else {
            cout << " Wrong input\n";
            return 0;
        }
        inputFile_st += fermionType;
        inputFile_st += "_E-";
        inputFile_st += words[2];
        inputFile_st += "_m-";
        inputFile_st += words[3];
        // inputFile_st += "_tau";
        // inputFile_st += "_l";
        inputFile_st += ".root";

        // expect: "../features_ISR/sig_Maj_E-3_m-1000_reco.root"
        outputFile_st = "../data/features/sig_";
        outputFile_st += fermionType;
        outputFile_st += "_E-";
        outputFile_st += words[2];
        outputFile_st += "_m-";
        outputFile_st += words[3];
        outputFile_st += "_reco.root";

        *inputFile_st_ = inputFile_st;
        *outputFile_st_ = outputFile_st;
        return 1;

        // expect input like: "b_MuMu_qqll_3" means, mu-mu interaction at sqrt{s}=3TeV
        // expect input like: "b_aMu_qql_3" means, photon-mu interaction at sqrt{s}=3TeV
    } else if (type.at(0) == 'b') {
        cout << "Processing Background data" << endl;

        vector<string> words{};
        for (Int_t i = 0; i < 4; i++) {
            Int_t pos = type.find("_");
            words.push_back(type.substr(0, pos));
            type.erase(0, pos + 1);
        }

        // expect: "../data/bg_MuMu_qqll_E-3.root"
        inputFile_st = "../data/detector/bg_";
        inputFile_st += words[1];
        inputFile_st += "_";
        inputFile_st += words[2];
        inputFile_st += "_E-";
        inputFile_st += words[3];
        inputFile_st += ".root";

        // expect: "../data/bg_2W_E-3_reco.root"
        outputFile_st = "../data/features/bg_";
        outputFile_st += words[1];
        outputFile_st += "_";
        outputFile_st += words[2];
        outputFile_st += "_E-";
        outputFile_st += words[3];
        outputFile_st += "_reco.root";

        *inputFile_st_ = inputFile_st;
        *outputFile_st_ = outputFile_st;
        return 1;

    } else {
        cout << "Wrong input type" << endl;
        return 0;
    }
}
///}}}

void allinone_ISR(
    const string type,
    const Bool_t save = false,
    Int_t num_test = 0,  // 0: run all data samples
    Int_t print_detail = 1,
    Int_t all_on = 1) {  // turn on all the cuts
    // formmating the input output files
    string inputFile_st;
    string outputFile_st;

    Int_t foundFiles;
    foundFiles = getFileNames(type, &inputFile_st, &outputFile_st);

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
    // TClonesArray* branchVLC1Jet = treeReader->UseBranch("VLCjetR17N1");
    TClonesArray* branchVLC2Jet = treeReader->UseBranch("VLCjetR02N2");
    // TClonesArray* branchVLC2Jet = treeReader->UseBranch("VLCjetR05N2");
    TClonesArray* branchVLC3Jet = treeReader->UseBranch("VLCjetR02N3");
    TClonesArray* branchMET = treeReader->UseBranch("MissingET");
    TClonesArray* branchFwMu = treeReader->UseBranch("ForwardMuon");
    TClonesArray* branchTrack = treeReader->UseBranch("Track");

    // book feature storing tree and file
    TFile fea(outputFile, "recreate");
    TTree tr("t", "features");
    Features* features = new Features;
    tr.Branch("features", &features);

    if (num_test == 0) num_test = numberOfEntries;

    // tracing the number of events are each cut
    Int_t nEv = 0;          // # of events in truth level
    Int_t nFS = 0;          // # of events having identified final states
    Int_t nLepEta = 0;      // # of events after lep eta cut
    Int_t nLepPt = 0;       // # of events after lep pt cut
    Int_t nJJPt = 0;        // # of events after jj pt cut
    Int_t nJJM = 0;         // # of events after jj mas cut
    BkgTypes bkgTypes;      // Count the type of bkg
    BkgTypes bkgTypesReco;  // count the type of bkg after reconstruction

    Float_t lepEtaCut = 2.5;                    // lep eta cut
    Float_t lepPtCut = 50;                     // lep pt cut
    Float_t jjPtCut = 50;                      // jj pt cut
    Float_t WMLowCut = mWPDG - 5 * widthWPDG;   // jj mass cut
    Float_t WMHighCut = mWPDG + 5 * widthWPDG;  // jj mass cut

    // loop the events
    for (Int_t i_en = 0; i_en < num_test; i_en++) {
        // progress
        if ((i_en % 1000) == 0) cout << "Reconstruction Progress: " << i_en << "/" << numberOfEntries << "\r";
        //cout << i_en << endl;
        cout.flush();
        // cout << "\nEvent: " << i_en << endl;
        // cout << "............................................." << endl;

        treeReader->ReadEntry(i_en);  // reading the entry

        //========================================================================
        //=======================   Classifiy Event Type   =======================
        //========================================================================
        iFinalStates iFSTrue;                                 // indeces of the true level final states
        TLorentzVector lepTrue, jet1True_, jet2True_, NTrue;  // truth level lepton, jets, and HNL
        Int_t typeLepTrue;
        Int_t passing = 0;
        if (type.at(0) == 's') {      // search for signals
            if (type.at(2) == 'M') {  // Majorana
                passing = ClassifySingal(branchParticle, &iFSTrue, &NTrue, 9900012);
            } else if (type.at(2) == 'D') {  // Dirac
                passing = ClassifySingal(branchParticle, &iFSTrue, &NTrue, 9990012);
            }
        } else {
            passing = ClassifyiBackground(branchParticle, &bkgTypes);
        }
        if (passing == 0) continue;

        nEv += 1;
        if (type.at(0) == 's') {  // search for signals
            lepTrue = iFSTrue.iLeps[0];
            jet1True_ = iFSTrue.i2Jets[0];
            jet2True_ = iFSTrue.i2Jets[1];
            if (iFSTrue.iElectronIndeces.size() == 1) typeLepTrue = 11;
            if (iFSTrue.iMuonIndeces.size() == 1) typeLepTrue = 13;
            if (iFSTrue.iTauIndeces.size() == 1) typeLepTrue = 15;
        }

        TLorentzVector jjTrue;
        jjTrue = jet1True_ + jet2True_;

        //========================================================================
        //=======================      Reconstruction      =======================
        //========================================================================
        iFinalStates iFS;
        // finding final states: at least two leptons + 1 (or 2) jet
        iFS = FindFinalStatesIndex(branchElectron, branchMuon, branchVLC1Jet, branchVLC2Jet);
        if (all_on && iFS.foundAll == 0) continue;
        nFS += 1;  // found the targeted final states

        TLorentzVector jet1, jet21, jet22, jet2;
        TLorentzVector lep1, lep2, jj, N1, N2, JJLep1Lep2;
        pair<Int_t, Int_t> chargeLeps;
        pair<Int_t, Int_t> typeLeps;
        pair<TLorentzVector, TLorentzVector> leps;

        jj = getWJet(iFS);

        getLep(iFS, jj, &typeLeps, &chargeLeps, &leps);
        lep1 = leps.first; lep2 = leps.second;
        nLepEta += 1;

        // reconstruct N
        N1 = jj + lep1;
        N2 = jj + lep2;
        JJLep1Lep2 = jj + lep1 + lep2;

        //========================================================================
        //=======================           Cuts           =======================
        //========================================================================

        if (all_on && not(lep1.Pt() >= lepPtCut && lep2.Pt() >= lepPtCut)) continue;
        nLepPt += 1;

        if (all_on && not(jj.Pt() >= jjPtCut)) continue;
        nJJPt += 1;

        if (all_on && not(jj.M() >= WMLowCut && jj.M() <= WMHighCut)) continue;
        nJJM += 1;
        // cout << " .................................................................................. " <<endl;

        //========================================================================
        //=======================         Features         =======================
        //========================================================================
        // event index
        features->iEvt = i_en;

        // lepton 4-momentum info
        features->ptLep1 = lep1.Pt();
        features->etaLep1 = lep1.Eta();
        features->phiLep1 = lep1.Phi();
        features->ELep1 = lep1.E();
        features->pxLep1 = lep1.Px();
        features->pyLep1 = lep1.Py();
        features->pzLep1 = lep1.Pz();
        // lepton 4-momentum info
        features->ptLep2 = lep2.Pt();
        features->etaLep2 = lep2.Eta();
        features->phiLep2 = lep2.Phi();
        features->ELep2 = lep2.E();
        features->pxLep2 = lep2.Px();
        features->pyLep2 = lep2.Py();
        features->pzLep2 = lep2.Pz();

        // lepton info
        features->chargeLep1 = chargeLeps.first;
        Int_t lep1isMu_ = 0, lep1isEle_ = 0;
        if (typeLeps.first == 13) {
            lep1isMu_ = 1; lep1isEle_ = 0;
        } else if (typeLeps.first == 11) {
            lep1isMu_ = 0; lep1isEle_ = 1;
        }
        features->lep1isMu = lep1isMu_;
        features->lep1isEle = lep1isEle_;

        // lepton info
        features->chargeLep2 = chargeLeps.second;
        Int_t lep2isMu_ = 0, lep2isEle_ = 0;
        if (typeLeps.second == 13) {
            lep2isMu_ = 1; lep2isEle_ = 0;
        } else if (typeLeps.second == 11) {
            lep2isMu_ = 0; lep2isEle_ = 1;
        }
        features->lep2isMu = lep2isMu_;
        features->lep2isEle = lep2isEle_;

        // W boson 4-momentum info
        features->mJJ = jj.M();
        features->EJJ = jj.E();
        features->PJJ = jj.P();
        features->ptJJ = jj.Pt();
        features->etaJJ = jj.Eta();
        features->phiJJ = jj.Phi();
        features->pxJJ = jj.Px();
        features->pyJJ = jj.Py();
        features->pzJJ = jj.Pz();

        // HNL 4-momentum info
        features->mN1 = N1.M();
        features->EN1 = N1.E();
        features->PN1 = N1.P();
        features->ptN1 = N1.Pt();
        features->etaN1 = N1.Eta();
        features->phiN1 = N1.Phi();
        features->pxN1 = N1.Px();
        features->pyN1 = N1.Py();
        features->pzN1 = N1.Pz();

        // HNL 4-momentum info
        features->mN2 = N2.M();
        features->EN2 = N2.E();
        features->PN2 = N2.P();
        features->ptN2 = N2.Pt();
        features->etaN2 = N2.Eta();
        features->phiN2 = N2.Phi();
        features->pxN2 = N2.Px();
        features->pyN2 = N2.Py();
        features->pzN2 = N2.Pz();

        features->mJJLep1Lep2 = JJLep1Lep2.M();

        // angular distance between the (two jets) and (W and lepton)
        Float_t DeltaRjj = 99999;
        // angular distance between the two jets, only if there are two jets in VLC branch
        Float_t EJ1 = 99999, ptJ1 = 99999, etaJ1 = 99999, phiJ1 = 99999;
        Float_t EJ2 = 99999, ptJ2 = 99999, etaJ2 = 99999, phiJ2 = 99999;
        if (iFS.i2Jets.size() == 2) {
            Float_t jjEtaDiff = jet1.Eta() - jet2.Eta();
            Float_t jjPhiDiff = deltaPhi(jet1.Phi(), jet2.Phi());
            DeltaRjj = pow(jjEtaDiff * jjEtaDiff + jjPhiDiff * jjPhiDiff, 0.5);
            // the energy of the jets, used as info to measure the imbalance
            // (since the W boson in t-ch is mostly longitudnal, so generally small imbalance)
            EJ1 = jet21.E();
            ptJ1 = jet21.Pt();
            etaJ1 = jet21.Pt();
            phiJ1 = jet21.Phi();

            EJ2 = jet22.E();
            ptJ2 = jet22.Pt();
            etaJ2 = jet22.Pt();
            phiJ2 = jet22.Phi();
        }
        features->DeltaRjj = DeltaRjj;
        features->EJet1 = EJ1;
        features->ptJet1 = ptJ1;
        features->etaJet1 = etaJ1;
        features->phiJet1 = phiJ1;
        features->EJet2 = EJ2;
        features->ptJet2 = ptJ2;
        features->etaJet2 = etaJ2;
        features->phiJet2 = phiJ2;

        Float_t tau1_ = 99999, tau2_ = 99999;
        if (iFS.i1Jets.size() == 1) {
            tau1_ = iFS.tau1;
            tau2_ = iFS.tau2;
        }
        features->tau1 = tau1_;
        features->tau2 = tau2_;

        // angular distance between the W and lepton
        Float_t DeltaRjjl1 = 99999, DeltaRjjl2 = 99999;
        Float_t jjl1EtaDiff = jj.Eta() - lep1.Eta();
        Float_t jjl1PhiDiff = deltaPhi(jj.Phi(), lep1.Phi());
        DeltaRjjl1 = pow(jjl1EtaDiff * jjl1EtaDiff + jjl1PhiDiff * jjl1PhiDiff, 0.5);
        features->DeltaRjjl1 = DeltaRjjl1;
        features->DeltaPhijjl1 = jjl1PhiDiff;  // the phi difference between the W and lepton
        Float_t jjl2EtaDiff = jj.Eta() - lep2.Eta();
        Float_t jjl2PhiDiff = deltaPhi(jj.Phi(), lep2.Phi());
        DeltaRjjl2 = pow(jjl2EtaDiff * jjl2EtaDiff + jjl2PhiDiff * jjl2PhiDiff, 0.5);
        features->DeltaRjjl2 = DeltaRjjl2;
        features->DeltaPhijjl2 = jjl2PhiDiff;  // the phi difference between the W and lepton

        // extra, just for trial
        // the forward muon info, for distinguishing the ISR mu-photon interaction
        Int_t nFwMu = branchFwMu->GetEntries();
        Float_t ptFwmu = 99999;
        if (nFwMu != 0) {
            Muon* fwmu = (Muon*)branchFwMu->At(0);
            ptFwmu = fwmu->PT;
            cout << " fwmu: " << fwmu->PT << "\n";
        }
        features->ptFwMu = ptFwmu;

        features->ptNTrue = NTrue.Pt();
        features->ENTrue = NTrue.E();
        features->etaNTrue = NTrue.Eta();

        features->ptJJTrue = jjTrue.Pt();
        features->EJJTrue = jjTrue.E();
        features->PJJTrue = jjTrue.P();
        features->etaJJTrue = jjTrue.Eta();
        features->phiJJTrue = jjTrue.Phi();

        features->etaLepTrue = lepTrue.Eta();
        features->phiLepTrue = lepTrue.Phi();
        features->ptLepTrue = lepTrue.Pt();
        features->ELepTrue = lepTrue.E();
        features->pxLepTrue = lepTrue.Px();
        features->pyLepTrue = lepTrue.Py();
        features->pzLepTrue = lepTrue.Pz();

        Int_t nJetcon;
        Jet* wjet1 = (Jet*)branchVLC1Jet->At(0);

        Int_t nTracks;
        nTracks = branchTrack->GetEntries();




        tr.Fill();
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
        cout << "\t# after jj pt cut \t(pT(jj) >= " << float(jjPtCut) << "):\t\t\t" << nJJPt << "\t(" << 100 * float(nJJPt) / float(nLepPt) << "%)" << endl;
        cout << "\t# after jj mass cut \t(" << WMLowCut << " <= m(jj) <= " << float(WMHighCut) << "):\t" << nJJM << "\t(" << 100 * float(nJJM) / float(nJJPt) << "%)" << endl;
        cout << "\t---------------------------------------------------------------------------------" << endl;
        cout << "\t\t\t\t\t\tTotal eff.:\t" << nJJM << " / " << nEv << "\t(" << 100 * float(nJJM) / float(nEv) << " %) " << endl;
        cout << "\t---------------------------------------------------------------------------------" << endl;
        cout << "\n\n\n";
    }

    tr.Write();
    fea.Close();
    cout << "Writing to: " << outputFile << "\n\n";
    cout << "STORE: " << outputFile << "; e=" << float(nJJM) / float(nEv) << " END STORE";
}
