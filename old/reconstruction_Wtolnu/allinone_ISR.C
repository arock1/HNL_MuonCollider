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
        inputFile_st = "../data_ISR/sig_";
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
        outputFile_st = "../features_ISR/sig_";
        outputFile_st += fermionType;
        outputFile_st += "_E-";
        outputFile_st += words[2];
        outputFile_st += "_m-";
        outputFile_st += words[3];
        outputFile_st += "_reco.root";

        *inputFile_st_ = inputFile_st;
        *outputFile_st_ = outputFile_st;
        return 1;

        // expect input like: "b_sZ_3" means background, s-channel Z boson subtype, sqrt{s}=3TeV
    } else if (type.at(0) == 'b') {
        cout << "Processing Background data" << endl;

        vector<string> words{};
        for (Int_t i = 0; i < 3; i++) {
            Int_t pos = type.find("_");
            words.push_back(type.substr(0, pos));
            type.erase(0, pos + 1);
        }

        // expect: "../data_ISR/bg_2W_E-3.root"
        inputFile_st = "../data_ISR/bg_";
        inputFile_st += words[1];
        inputFile_st += "_E-";
        inputFile_st += words[2];
        inputFile_st += ".root";

        // expect: "../features_ISR/bg_2W_E-3_reco.root"
        outputFile_st = "../features_ISR/bg_";
        outputFile_st += words[1];
        outputFile_st += "_E-";
        outputFile_st += words[2];
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
            if (passing == 1) {
                cout << " passing: " << passing << "\n";
                lepTrue = iFSTrue.iLeps[0];
                jet1True_ = iFSTrue.i2Jets[0];
                jet2True_ = iFSTrue.i2Jets[1];
                if (iFSTrue.iElectronIndeces.size() == 1) typeLepTrue = 11;
                if (iFSTrue.iMuonIndeces.size() == 1) typeLepTrue = 13;
            }
        } else {
            passing = ClassifyiBackground(branchParticle, &bkgTypes);
        }

        // if (passing == 0) continue;
        nEv += 1;

        // TLorentzVector jjTrue;
        // jjTrue = jet1True_ + jet2True_;

        //========================================================================
        //=======================      Reconstruction      =======================
        //========================================================================
        iFinalStates iFS;
        // finding final states: lepton + 1 (or 2) jet
        iFS = FindFinalStatesIndex(branchElectron, branchMuon, branchVLC1Jet, branchVLC2Jet);

        if (iFS.foundAll == 0) continue;
        nFS += 1;

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
        // Tagging truth level
        // TLorentzVector jet1True, jet2True;
        // if (iFS.i2Jets.size() == 2) {
        //// if there are two jets, match them to the reconstcuted jet according to the eta
        // if (abs(jet1.Eta() - jet1True_.Eta()) < abs(jet1.Eta() - jet2True_.Eta())) {
        // jet1True = jet1True_;
        // jet2True = jet2True_;
        //} else {
        // jet1True = jet2True_;
        // jet2True = jet1True_;
        //}

        //} else if (iFS.i1Jets.size() == 1) {
        // jet1True = jet1True_;
        // jet2True = jet2True_;
        //}

        Float_t DeltaRjj, DeltaRjjl;
        if (iFS.i2Jets.size() == 2) {
            Float_t jjEtaDiff = jet1.Eta() - jet2.Eta();
            Float_t jjPhiDiff = deltaPhi(jet1.Phi(), jet2.Phi());
            DeltaRjj = pow(jjEtaDiff * jjEtaDiff + jjPhiDiff * jjPhiDiff, 0.5);
        }

        Float_t jjlEtaDiff = jj.Eta() - lep.Eta();
        Float_t jjlPhiDiff = deltaPhi(jj.Phi(), lep.Phi());
        DeltaRjjl = pow(jjlEtaDiff * jjlEtaDiff + jjlPhiDiff * jjlPhiDiff, 0.5);

        // Float_t DeltaRjjTrue, DeltaRjjlTrue;
        // Float_t jjEtaDiffTrue = jet1True.Eta() - jet2True.Eta();
        // Float_t jjPhiDiffTrue = deltaPhi(jet1True.Phi(), jet2True.Phi());
        // DeltaRjjTrue = pow(jjEtaDiffTrue * jjEtaDiffTrue + jjPhiDiffTrue * jjPhiDiffTrue, 0.5);
        // Float_t jjlEtaDiffTrue = jjTrue.Eta() - lepTrue.Eta();
        // Float_t jjlPhiDiffTrue = deltaPhi(jjTrue.Phi(), lepTrue.Phi());
        // DeltaRjjlTrue = pow(jjlEtaDiffTrue * jjlEtaDiffTrue + jjlPhiDiffTrue * jjlPhiDiffTrue, 0.5);

        features->iEvt = i_en;
        features->ptLep = lep.Pt();
        features->etaLep = lep.Eta();
        features->phiLep = lep.Phi();
        features->ELep = lep.E();
        features->pxLep = lep.Px();
        features->pyLep = lep.Py();
        features->pzLep = lep.Pz();

        features->DeltaPhijjl = jjlPhiDiff;
        features->DeltaRjj = DeltaRjj;
        features->DeltaRjjl = DeltaRjjl;
        if (iFS.i2Jets.size() == 2) features->DeltaRjj = DeltaRjj;
        if (iFS.i2Jets.size() == 2) features->DeltaRjjl = DeltaRjjl;

        if (iFS.i2Jets.size() == 2) {
            features->pTheta = abs(jet21.E() - jet22.E()) / jj.P();
        }

        // if (iFS.i2Jets.size() != 2) cout << " .....................: ";
        if (iFS.i2Jets.size() == 2) {
            features->EJet1 = jet21.E();
            features->EJet2 = jet22.E();
        }
        // features->nJets = iFS.iJets.size();
        features->mJJ = jj.M();
        features->ptJJ = jj.Pt();
        features->etaJJ = jj.Eta();
        features->phiJJ = jj.Phi();
        features->mN = N.M();
        features->ptN = N.Pt();
        features->etaN = N.Eta();
        features->phiN = N.Phi();
        features->pxN = N.Px();
        features->pyN = N.Py();
        features->pzN = N.Pz();

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

        features->chargeLep = chargeLep;
        // if (type.at(0) != 'b') {
        // features->chargeLepTrue = iFSTrue.iLepCharges[0];
        //}

        features->typeLep = typeLep;
        // features->typeLepTrue = typeLepTrue;

        // Preliminary: Try to find the second lepton for the 2VBF case, for LNV
        // finding the muon that give cloeset mass to W?
        TLorentzVector lep2;
        Float_t mDiffWMin = 99999;
        Int_t chargeLep2;
        Int_t lep2Index = 99999;
        Int_t typeLep2;
        TLorentzVector lN;  // second lepton + N
        for (Int_t il = 0; il < iFS.iLepCharges.size(); il++) {
            if (il == lepIndex) continue;
            TLorentzVector lepI = iFS.iLeps[il];
            // eta cut
            if (not(abs(lepI.Eta()) <= lepEtaCut)) continue;

            lN = lep2 + N;
            if (abs(mWPDG - lN.M()) < mDiffWMin) {
                chargeLep2 = iFS.iLepCharges[il];
                mDiffWMin = abs(mWPDG - lN.M());
                if (not(lN.M() >= WMLowCut && lN.M() <= WMHighCut)) continue;
                lep2 = lepI;
            }
        }
        features->typeLep2 = typeLep2;

        Int_t nMET = branchMET->GetEntries();
        if (nMET != 0) {
            MissingET* met = (MissingET*)branchMET->At(0);
            features->MET = met->MET;
            features->DeltaPhiNMET = deltaPhi(met->Phi, N.Phi());
        }

        Float_t minptlep = 99999;
        for (Int_t il = 0; il < iFS.iLepCharges.size(); il++) {
            TLorentzVector lepI = iFS.iLeps[il];
            if (lepI.Pt() < minptlep) minptlep = lepI.Pt();
        }
        features->MinPtLep = minptlep;

        passing = ClassifyiBackground(branchParticle, &bkgTypesReco);
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
        cout << "\t# after jj pt cut \t(pT(jj) >= " << float(jjPtCut) << "):\t\t" << nJJPt << "\t(" << 100 * float(nJJPt) / float(nLepPt) << "%)" << endl;
        cout << "\t# after jj mass cut \t(" << WMLowCut << " <= m(jj) <= " << float(WMHighCut) << "):\t" << nJJM << "\t(" << 100 * float(nJJM) / float(nJJPt) << "%)" << endl;
        // cout << "# after jet pt cut (pT(j1, j2) >= " << float(jetPtCut) << "):\t\t" << nJetPt << "\t(" << 100 * float(nJetPt) / float(nLepPt) << "%)" << endl;
        cout << "\t---------------------------------------------------------------------------------" << endl;
        cout << "\t\t\t\t\t\tTotal eff.:\t" << nJJM << " / " << nEv << "\t(" << 100 * float(nJJM) / float(nEv) << " %) " << endl;
        cout << "\t---------------------------------------------------------------------------------" << endl;
        cout << "\n\n\n";

        // if (type.at(0) == 'b') {
        // cout << "\t---------------------------------------------------------------------------------" << endl;
        // cout << "\tBkg. Types\t\tNumber in truth level \t\tReconstructed Numebr" << endl;
        //// cout << "\t(of same 'mother')" << endl;
        // cout << "\t---------------------------------------------------------------------------------" << endl;
        // cout << "\t# of W-W:\t\t" << bkgTypes.WW << "\t(" << 100 * float(bkgTypes.WW) / float(nEv) << "%)"
        //<< "\t\t" << bkgTypesReco.WW << "\t(" << 100 * float(bkgTypesReco.WW) / float(nJJM) << "%)" << endl;
        // cout << "\t# of photon-photon:\t" << bkgTypes.phph << "\t(" << 100 * float(bkgTypes.phph) / float(nEv) << "%)"
        //<< "\t\t\t" << bkgTypesReco.phph << "\t(" << 100 * float(bkgTypesReco.phph) / float(nJJM) << "%)" << endl;
        // cout << "\t# of Z-photon:\t\t" << bkgTypes.Zph << "\t(" << 100 * float(bkgTypes.Zph) / float(nEv) << "%)"
        //<< "\t\t\t" << bkgTypesReco.Zph << "\t(" << 100 * float(bkgTypesReco.Zph) / float(nJJM) << "%)" << endl;
        // cout << "\t# of Z-Z:\t\t" << bkgTypes.ZZ << "\t(" << 100 * float(bkgTypes.ZZ) / float(nEv) << "%)"
        //<< "\t\t\t" << bkgTypesReco.ZZ << "\t(" << 100 * float(bkgTypesReco.ZZ) / float(nJJM) << "%)" << endl;
        // cout << "\t# of W-Z:\t\t" << bkgTypes.WZ << "\t(" << 100 * float(bkgTypes.WZ) / float(nEv) << "%)"
        //<< "\t\t\t" << bkgTypesReco.WZ << "\t(" << 100 * float(bkgTypesReco.WZ) / float(nJJM) << "%)" << endl;
        // cout << "\t# of W-photon:\t\t" << bkgTypes.Wph << "\t(" << 100 * float(bkgTypes.Wph) / float(nEv) << "%)"
        //<< "\t\t" << bkgTypesReco.Wph << "\t(" << 100 * float(bkgTypesReco.Wph) / float(nJJM) << "%)" << endl;
        // cout << endl;

        // cout << "\t# of lepFromPhys:\t" << bkgTypes.lepFromPhys << "\t(" << 100 * float(bkgTypes.lepFromPhys) / float(nEv) << "%)"
        //<< "\t\t\t" << bkgTypesReco.lepFromPhys << "\t(" << 100 * float(bkgTypesReco.lepFromPhys) / float(nJJM) << "%)" << endl;
        // cout << "\t\t lep=e: \t  [" << bkgTypes.eFromPhys << "\t  (" << 100 * float(bkgTypes.eFromPhys) / float(bkgTypes.lepFromPhys) << "%)]"
        //<< "\t\t  [" << bkgTypesReco.eFromPhys << "\t  (" << 100 * float(bkgTypesReco.eFromPhys) / float(bkgTypesReco.lepFromPhys) << "%)]" << endl;
        // cout << "\t\t lep=mu:\t  [" << bkgTypes.muFromPhys << "\t  (" << 100 * float(bkgTypes.muFromPhys) / float(bkgTypes.lepFromPhys) << "%)]"
        //<< "\t\t  [" << bkgTypesReco.muFromPhys << "\t  (" << 100 * float(bkgTypesReco.muFromPhys) / float(bkgTypesReco.lepFromPhys) << "%)]" << endl;
        // cout << endl;

        // cout << "\t# of others:\t\t" << bkgTypes.others << "\t(" << 100 * float(bkgTypes.others) / float(nEv) << "%)"
        //<< "\t\t" << bkgTypesReco.others << "\t(" << 100 * float(bkgTypesReco.others) / float(nJJM) << "%)" << endl;
        // cout << "\t\t Have electron:\t  [" << bkgTypes.others_haveEle << "\t  (" << 100 * float(bkgTypes.others_haveEle) / float(bkgTypes.others) << "%)]"
        //<< "\t\t  [" << bkgTypesReco.others_haveEle << "\t  (" << 100 * float(bkgTypesReco.others_haveEle) / float(bkgTypesReco.others) << "%)]" << endl;
        // cout << "\t\t No electron:  \t  [" << bkgTypes.others - bkgTypes.others_haveEle << "\t  (" << 100 * float(bkgTypes.others - bkgTypes.others_haveEle) / float(bkgTypes.others) << "%)]"
        //<< "\t  [" << bkgTypesReco.others - bkgTypesReco.others_haveEle << "\t  (" << 100 * float(bkgTypesReco.others - bkgTypesReco.others_haveEle) / float(bkgTypesReco.others) << "%)]" << endl;
        // cout << "\t---------------------------------------------------------------------------------" << endl;
        //}
    }

    tr.Write();
    fea.Close();
    cout << "Writing to: " << outputFile << "\n\n";
    cout << "STORE: " << outputFile << "; e=" << float(nJJM) / float(nEv) << " END STORE";
    // Int_t pos = type.find("_");
    // if (type[0] == 'b') {
    // cout << "STORE: t=b; E=" << type.substr(1, -1) << "; m=" << NAN << "; e=" << float(nJJM) / float(nEv) << " END STORE";
    //} else {
    // cout << "STORE: t=" << type.substr(0, 2) << "; E=" << type.substr(2, pos - 2) << "; m=" << type.substr(pos + 1, -1) << "; e=" << float(nJJM) / float(nEv) << " END STORE";
    //}
}
