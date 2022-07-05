//
// ===========================================================================
// || This Code Does:                                                       ||
// || 1. Classifiy the event type (Signal / What Background type)           ||
// ||    and store the corresponding truth level info.                      ||
// || 2. Find the targeted final state in detector level:                   ||
// ||       a. At least one lepton (electron or muon)                       ||
// ||       b. Either one 1VLC jet or two 2VLC jets found                   ||
// || 3. Reconstruct W jet                                                  ||
// ||       a. If there is only (1VLC case found), use it as W jet          ||
// ||       b. If there is only (2VLC case found), use their sum as W jet   ||
// || 4. Select lepton (Need to have |eta|<threshold to simulate detector)  ||
// ||       a. If there is only one, use it                                 ||
// ||       b. If there are multiple, use the one with highest pT           ||
// || 5. Apply cuts                                                         ||
// ||       a. lepton pT                                                    ||
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

        // expect: "../data_ISR/sig_Maj_E-3_m-1000.root"
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

        // expect: "../data_ISR/bg_MuMu_qqll_E-3.root"
        inputFile_st = "../data/detector/bg_";
        inputFile_st += words[1];
        inputFile_st += "_";
        inputFile_st += words[2];
        inputFile_st += "_E-";
        inputFile_st += words[3];
        inputFile_st += ".root";

        // expect: "../features_ISR/bg_2W_E-3_reco.root"
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
    Int_t num_test = 0,
    Int_t print_detail = 1) {
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
    TClonesArray* branchVLC2Jet = treeReader->UseBranch("VLCjetR02N2");
    TClonesArray* branchMET = treeReader->UseBranch("MissingET");
    TClonesArray* branchFwMu = treeReader->UseBranch("ForwardMuon");

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

    Float_t lepEtaCut = 2.5;  // lep eta cut
    Float_t lepPtCut = 100;   // lep pt cut
    Float_t jjPtCut = 100;    // jj pt cut

    Float_t WMLowCut = mWPDG - 5 * widthWPDG;   // jj mass cut
    Float_t WMHighCut = mWPDG + 5 * widthWPDG;  // jj mass cut

    // loop the events
    for (Int_t i_en = 0; i_en < num_test; i_en++) {
        // progress
        if ((i_en % 1000) == 0) cout << "Reconstruction Progress: " << i_en << "/" << numberOfEntries << "\r";
        cout.flush();
        // cout << "\nEvent: " << i_en << endl;

        treeReader->ReadEntry(i_en);  // reading the entry
        nEv += 1;

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
            lepTrue = iFSTrue.iLeps[0];
            jet1True_ = iFSTrue.i2Jets[0];
            jet2True_ = iFSTrue.i2Jets[1];
            if (iFSTrue.iElectronIndeces.size() == 1) typeLepTrue = 11;
            if (iFSTrue.iMuonIndeces.size() == 1) typeLepTrue = 13;
        } else {
            passing = ClassifyiBackground(branchParticle, &bkgTypes);
        }

        if (passing == 0) continue;

        TLorentzVector jjTrue;
        jjTrue = jet1True_ + jet2True_;
        //========================================================================
        //=======================      Reconstruction      =======================
        //========================================================================
        iFinalStates iFS;
        // finding final states: lepton + 1 (or 2) jet
        iFS = FindFinalStatesIndex(branchElectron, branchMuon, branchVLC1Jet, branchVLC2Jet);

        if (iFS.foundAll == 0) continue;
        nFS += 1;  // found the targeted final states

        TLorentzVector jet1, jet21, jet22, jet2;
        TLorentzVector lep, jj, N;

        // if there are both one jet and two jets, then pick the one with mass closest to W boson
        // if there is only one(two) jet, then assign that as an W jet
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

        // loop over all the leptons
        Float_t ptLepMax = -99999;
        Int_t chargeLep;
        Int_t lepIndex;
        Int_t typeLep;
        Int_t foundLep = 0;
        for (Int_t il = 0; il < iFS.iLepCharges.size(); il++) {
            TLorentzVector lepI = iFS.iLeps[il];
            // eta cut
            if (not(abs(lepI.Eta()) <= lepEtaCut)) continue;
            // store the lepton with max. pT among all
            if (lepI.Pt() > ptLepMax) {
                chargeLep = iFS.iLepCharges[il];
                foundLep = 1;
                ptLepMax = lepI.Pt();
                lep = lepI;
                lepIndex = il;
                if (il < iFS.iElectronIndeces.size()) typeLep = 11;
                if (il >= iFS.iElectronIndeces.size()) typeLep = 13;
            }
        }
        // if no such lepton, then pass
        if (foundLep == 0) continue;
        nLepEta += 1;

        // reconstruct N
        N = jj + lep;

        //========================================================================
        //=======================           Cuts           =======================
        //========================================================================

        if (not(lep.Pt() >= lepPtCut)) continue;
        nLepPt += 1;

        if (not(jj.Pt() >= jjPtCut)) continue;
        nJJPt += 1;

        if (not(jj.M() >= WMLowCut && jj.M() <= WMHighCut)) continue;
        nJJM += 1;

        //========================================================================
        //=======================         Features         =======================
        //========================================================================
        // event index
        features->iEvt = i_en;

        // lepton 4-momentum info
        features->ptLep = lep.Pt();
        features->etaLep = lep.Eta();
        features->phiLep = lep.Phi();
        features->ELep = lep.E();
        features->pxLep = lep.Px();
        features->pyLep = lep.Py();
        features->pzLep = lep.Pz();
        // lepton info for distinguishing Maj/Dir
        features->chargeLep = chargeLep;
        if (typeLep == 13) features->lepisMu = 1;
        if (typeLep == 11) features->lepisMu = 0;

        // W boson 4-momentum info
        features->mJJ = jj.M();
        features->ptJJ = jj.Pt();
        features->etaJJ = jj.Eta();
        features->phiJJ = jj.Phi();

        // HNL 4-momentum info
        features->mN = N.M();
        features->ptN = N.Pt();
        features->etaN = N.Eta();
        features->phiN = N.Phi();
        features->pxN = N.Px();
        features->pyN = N.Py();
        features->pzN = N.Pz();

        // angular distance between the (two jets) and (W and lepton)
        Float_t DeltaRjj, DeltaRjjl;
        // angular distance between the two jets, only if there are two jets in VLC branch
        if (iFS.i2Jets.size() == 2) {
            Float_t jjEtaDiff = jet1.Eta() - jet2.Eta();
            Float_t jjPhiDiff = deltaPhi(jet1.Phi(), jet2.Phi());
            DeltaRjj = pow(jjEtaDiff * jjEtaDiff + jjPhiDiff * jjPhiDiff, 0.5);
            features->DeltaRjj = DeltaRjj;
            // the energy of the jets, used as info to measure the imbalance
            // (since the W boson in t-ch is mostly longitudnal, so generally small imbalance)
            features->EJet1 = jet21.E();
            features->ptJet1 = jet21.Pt();
            features->etaJet1 = jet21.Pt();
            features->phiJet1 = jet21.Phi();

            features->EJet2 = jet22.E();
            features->ptJet2 = jet22.Pt();
            features->etaJet2 = jet22.Pt();
            features->phiJet2 = jet22.Phi();
        }

        // angular distance between the W and lepton
        Float_t jjlEtaDiff = jj.Eta() - lep.Eta();
        Float_t jjlPhiDiff = deltaPhi(jj.Phi(), lep.Phi());
        DeltaRjjl = pow(jjlEtaDiff * jjlEtaDiff + jjlPhiDiff * jjlPhiDiff, 0.5);
        features->DeltaRjjl = DeltaRjjl;
        features->DeltaPhijjl = jjlPhiDiff;  // the phi difference between the W and lepton

        // extra, just for trial
        Int_t nMET = branchMET->GetEntries();
        if (nMET != 0) {
            MissingET* met = (MissingET*)branchMET->At(0);
            features->MET = met->MET;
            features->DeltaPhiNMET = deltaPhi(met->Phi, N.Phi());
        }

        // extra, just for trial
        // the forward muon info, for distinguishing the ISR mu-photon interaction
        Int_t nFwMu = branchFwMu->GetEntries();
        if (nFwMu != 0) {
            Muon* fwmu = (Muon*)branchFwMu->At(0);
            features->ptFwMu = fwmu->PT;
        }

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
        cout << "\t# after jj pt cut \t(pT(jj) >= " << float(jjPtCut) << "):\t\t" << nJJPt << "\t(" << 100 * float(nJJPt) / float(nLepPt) << "%)" << endl;
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
