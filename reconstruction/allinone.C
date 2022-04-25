
// ==========================================================================
// || This Code
// ==========================================================================

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

    if (type.at(0) == 'b') {  // background (e.g. "b3": background with E=3TeV)
        inputFile_st = "../data/background_inclusive_E-";
        inputFile_st += type.substr(1, -1);  // energy (sqrt s)
        inputFile_st += "TeV_CustomizedJet.root";

        outputFile_st = "../features/background_reco_E-";
        outputFile_st += type.substr(1, -1);
        outputFile_st += "TeV.root";

        *inputFile_st_ = inputFile_st;
        *outputFile_st_ = outputFile_st;
        return 1;

    } else if (type.at(0) == 'i' || type.at(0) == 's' || type.at(0) == 't') {  // signals
        cout << "Processing Signal data" << endl;
        string fermionType;
        Int_t pos = type.find("_");
        if (type.at(0) == 'i') {  // inclusive signal (e.g. "iM3_01": inclusive signal with Majorana, E=3TeV, mN=0.1TeV)
            inputFile_st = "../data/signal_E-";
        } else if (type.at(0) == 's') {  // s-channel
            inputFile_st = "../data/signal_schannel_E-";
        } else if (type.at(0) == 't') {  // t-channel
            inputFile_st = "../data/signal_tchannel_E-";
        }
        inputFile_st += type.substr(2, pos - 2);
        inputFile_st += "TeV_N-";
        inputFile_st += type.substr(pos + 1, -1);
        inputFile_st += "TeV";
        if (type.at(1) == 'D') {
            fermionType = "_Dirac";
        } else if (type.at(1) == 'M') {
            fermionType = "";
        } else {
            cout << "Wrong input type" << endl;
            return 0;
        }
        inputFile_st += fermionType;
        inputFile_st += "_CustomizedJet.root";

        if (type.at(0) == 'i') {
            outputFile_st = "../features/signal_reco_E-";
        } else if (type.at(0) == 's') {
            outputFile_st = "../features/signal_schannel_reco_E-";
        } else if (type.at(0) == 't') {
            outputFile_st = "../features/signal_tchannel_reco_E-";
        }
        outputFile_st += type.substr(2, pos - 2);
        outputFile_st += "TeV_N-";
        outputFile_st += type.substr(pos + 1, -1);
        outputFile_st += "TeV";
        outputFile_st += fermionType;
        outputFile_st += ".root";
        *inputFile_st_ += inputFile_st;
        *outputFile_st_ = outputFile_st;
        return 1;

    } else {
        cout << "Wrong input type" << endl;
        return 0;
    }
}
///}}}

void allinone(
    const string type,
    const Bool_t save = false,
    Int_t num_test = 0) {
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
    Int_t nEv = 0;      // # of events in truth level
    Int_t nFS = 0;      // # of events having identified final states
    Int_t nLepEta = 0;  // # of events after lep eta cut
    Int_t nLepPt = 0;   // # of events after lep pt cut
    Int_t nJJPt = 0;    // # of events after jj pt cut
    Int_t nJJM = 0;     // # of events after jj mas cut
    BkgTypes bkgTypes;
    BkgTypes bkgTypesReco;

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
        if (type.at(0) == 'i' || type.at(0) == 's' || type.at(0) == 't') {  // search for signals
            if (type.at(1) == 'M') {                                        // Majorana
                passing = ClassifySingal(branchParticle, &iFSTrue, &NTrue, 9900012);
            } else if (type.at(1) == 'D') {
                passing = ClassifySingal(branchParticle, &iFSTrue, &NTrue, 9990012);
            }
            lepTrue = iFSTrue.iLeps[0];
            jet1True_ = iFSTrue.iJets[0];
            jet2True_ = iFSTrue.iJets[1];
            // cout << lepTrue.M() << endl;
            // cout << iFSTrue.iElectronIndeces.size() << "; " << iFSTrue.iMuonIndeces.size() << endl;
            if (iFSTrue.iElectronIndeces.size() == 1) typeLepTrue = 11;
            if (iFSTrue.iMuonIndeces.size() == 1) typeLepTrue = 13;
        } else {
            passing = ClassifyiBackground(branchParticle, &bkgTypes);
        }

        if (passing == 0) continue;
        nEv += 1;

        TLorentzVector jjTrue;
        jjTrue = jet1True_ + jet2True_;

        //========================================================================
        //=======================      Reconstruction      =======================
        //========================================================================
        iFinalStates iFS;
        // finding final states: lepton + 1 (or 2) jet
        iFS = FindFinalStatesIndex(branchElectron, branchMuon, branchVLC1Jet, branchVLC2Jet);

        if (iFS.foundAll == 0) continue;
        nFS += 1;

        TLorentzVector lep, jet1, jet2, jj, N;

        // if there is only one jet, then it is corresponding to the first index
        // and the W is this only jet
        // if there are two jets, then corresponding to the first and second index
        // and the W is their sum
        if (iFS.iJets.size() == 1) {
            jet1 = iFS.iJets[0];
            jj = jet1;
        } else if (iFS.iJets.size() == 2) {
            jet1 = iFS.iJets[0];
            jet2 = iFS.iJets[1];
            jj = jet1 + jet2;
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
        TLorentzVector jet1True, jet2True;
        Float_t DeltaRjj, DeltaRjjl;
        if (iFS.iJets.size() == 2) {
            // if there are two jets, match them to the reconstcuted jet according to the eta
            if (abs(jet1.Eta() - jet1True_.Eta()) < abs(jet1.Eta() - jet2True_.Eta())) {
                jet1True = jet1True_;
                jet2True = jet2True_;
            } else {
                jet1True = jet2True_;
                jet2True = jet1True_;
            }
            Float_t jjEtaDiff = jet1.Eta() - jet2.Eta();
            Float_t jjPhiDiff = deltaPhi(jet1.Phi(), jet2.Phi());
            DeltaRjj = pow(jjEtaDiff * jjEtaDiff + jjPhiDiff * jjPhiDiff, 0.5);

        } else if (iFS.iJets.size() == 1) {
            jet1True = jet1True_;
            jet2True = jet2True_;
        }

        Float_t jjlEtaDiff = jj.Eta() - lep.Eta();
        Float_t jjlPhiDiff = deltaPhi(jj.Phi(), lep.Phi());
        DeltaRjjl = pow(jjlEtaDiff * jjlEtaDiff + jjlPhiDiff * jjlPhiDiff, 0.5);

        Float_t DeltaRjjTrue, DeltaRjjlTrue;
        Float_t jjEtaDiffTrue = jet1True.Eta() - jet2True.Eta();
        Float_t jjPhiDiffTrue = deltaPhi(jet1True.Phi(), jet2True.Phi());
        DeltaRjjTrue = pow(jjEtaDiffTrue * jjEtaDiffTrue + jjPhiDiffTrue * jjPhiDiffTrue, 0.5);
        Float_t jjlEtaDiffTrue = jjTrue.Eta() - lepTrue.Eta();
        Float_t jjlPhiDiffTrue = deltaPhi(jjTrue.Phi(), lepTrue.Phi());
        DeltaRjjlTrue = pow(jjlEtaDiffTrue * jjlEtaDiffTrue + jjlPhiDiffTrue * jjlPhiDiffTrue, 0.5);

        features->iEvt = i_en;
        features->ptLep = lep.Pt();
        features->etaLep = lep.Eta();
        features->phiLep = lep.Phi();
        features->ELep = lep.E();

        features->ptJet1 = jet1.Pt();
        features->etaJet1 = jet1.Eta();
        features->phiJet1 = jet1.Phi();
        features->EJet1 = jet1.E();

        features->ptJet2 = jet2.Pt();
        features->etaJet2 = jet2.Eta();
        features->phiJet2 = jet2.Phi();
        features->EJet2 = jet2.E();

        features->DeltaPhijjl = jjlPhiDiff;
        features->DeltaRjj = DeltaRjj;
        features->DeltaRjjl = DeltaRjjl;
        if (iFS.iJets.size() == 2) features->DeltaRjj = DeltaRjj;
        if (iFS.iJets.size() == 2) features->DeltaRjjl = DeltaRjjl;

        features->nJets = iFS.iJets.size();
        features->mJJ = jj.M();
        features->ptJJ = jj.Pt();
        features->etaJJ = jj.Eta();
        features->phiJJ = jj.Phi();
        features->mN = N.M();
        features->ptN = N.Pt();
        features->etaN = N.Eta();
        features->phiN = N.Phi();
        features->pzN = N.Pz();

        features->ptLepTrue = lepTrue.Pt();
        features->etaLepTrue = lepTrue.Eta();
        features->phiLepTrue = lepTrue.Phi();
        features->ELepTrue = lepTrue.E();

        features->ptJet1True = jet1True.Pt();
        features->etaJet1True = jet1True.Eta();
        features->phiJet1True = jet1True.Phi();
        features->EJet1True = jet1True.E();

        features->ptJet2True = jet2True.Pt();
        features->etaJet2True = jet2True.Eta();
        features->phiJet2True = jet2True.Phi();
        features->EJet2True = jet2True.E();

        features->DeltaRjjTrue = DeltaRjjTrue;
        features->DeltaRjjlTrue = DeltaRjjlTrue;
        features->mJJTrue = jjTrue.M();
        features->ptJJTrue = jjTrue.Pt();
        features->etaJJTrue = jjTrue.Eta();
        features->phiJJTrue = jjTrue.Phi();
        features->mNTrue = NTrue.M();
        features->ptNTrue = NTrue.Pt();
        features->etaNTrue = NTrue.Eta();
        features->phiNTrue = NTrue.Phi();
        features->ENTrue = NTrue.E();
        features->pzNTrue = NTrue.Pz();

        features->chargeLep = chargeLep;
        if (type.at(0) != 'b') {
            features->chargeLepTrue = iFSTrue.iLepCharges[0];
        }

        features->typeLep = typeLep;
        features->typeLepTrue = typeLepTrue;

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
        }

        passing = ClassifyiBackground(branchParticle, &bkgTypesReco);
        tr.Fill();

        // cout << "M: " << lepTrue.M() << "; " << lep.M() << endl;
        // cout << "E: " << lepTrue.E() << "; " << lep.E() << endl;
        // cout << "P: " << lepTrue.P() << "; " << lep.P() << endl;
        // cout << "PT: " << lepTrue.Pt() << "; " << lep.Pt() << endl;
        // cout << endl;
    }
    cout << "Reconstruction Progress: " << numberOfEntries << "/" << numberOfEntries << "\n\n";

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

    if (type.at(0) == 'b') {
        cout << "\t---------------------------------------------------------------------------------" << endl;
        cout << "\tBkg. Types\t\tNumber in truth level \t\tReconstructed Numebr" << endl;
        // cout << "\t(of same 'mother')" << endl;
        cout << "\t---------------------------------------------------------------------------------" << endl;
        cout << "\t# of W-W:\t\t" << bkgTypes.WW << "\t(" << 100 * float(bkgTypes.WW) / float(nEv) << "%)"
             << "\t\t" << bkgTypesReco.WW << "\t(" << 100 * float(bkgTypesReco.WW) / float(nJJM) << "%)" << endl;
        cout << "\t# of photon-photon:\t" << bkgTypes.phph << "\t(" << 100 * float(bkgTypes.phph) / float(nEv) << "%)"
             << "\t\t\t" << bkgTypesReco.phph << "\t(" << 100 * float(bkgTypesReco.phph) / float(nJJM) << "%)" << endl;
        cout << "\t# of Z-photon:\t\t" << bkgTypes.Zph << "\t(" << 100 * float(bkgTypes.Zph) / float(nEv) << "%)"
             << "\t\t\t" << bkgTypesReco.Zph << "\t(" << 100 * float(bkgTypesReco.Zph) / float(nJJM) << "%)" << endl;
        cout << "\t# of Z-Z:\t\t" << bkgTypes.ZZ << "\t(" << 100 * float(bkgTypes.ZZ) / float(nEv) << "%)"
             << "\t\t\t" << bkgTypesReco.ZZ << "\t(" << 100 * float(bkgTypesReco.ZZ) / float(nJJM) << "%)" << endl;
        cout << "\t# of W-Z:\t\t" << bkgTypes.WZ << "\t(" << 100 * float(bkgTypes.WZ) / float(nEv) << "%)"
             << "\t\t\t" << bkgTypesReco.WZ << "\t(" << 100 * float(bkgTypesReco.WZ) / float(nJJM) << "%)" << endl;
        cout << "\t# of W-photon:\t\t" << bkgTypes.Wph << "\t(" << 100 * float(bkgTypes.Wph) / float(nEv) << "%)"
             << "\t\t" << bkgTypesReco.Wph << "\t(" << 100 * float(bkgTypesReco.Wph) / float(nJJM) << "%)" << endl;
        cout << endl;

        cout << "\t# of lepFromPhys:\t" << bkgTypes.lepFromPhys << "\t(" << 100 * float(bkgTypes.lepFromPhys) / float(nEv) << "%)"
             << "\t\t\t" << bkgTypesReco.lepFromPhys << "\t(" << 100 * float(bkgTypesReco.lepFromPhys) / float(nJJM) << "%)" << endl;
        cout << "\t\t lep=e: \t  [" << bkgTypes.eFromPhys << "\t  (" << 100 * float(bkgTypes.eFromPhys) / float(bkgTypes.lepFromPhys) << "%)]"
             << "\t\t  [" << bkgTypesReco.eFromPhys << "\t  (" << 100 * float(bkgTypesReco.eFromPhys) / float(bkgTypesReco.lepFromPhys) << "%)]" << endl;
        cout << "\t\t lep=mu:\t  [" << bkgTypes.muFromPhys << "\t  (" << 100 * float(bkgTypes.muFromPhys) / float(bkgTypes.lepFromPhys) << "%)]"
             << "\t\t  [" << bkgTypesReco.muFromPhys << "\t  (" << 100 * float(bkgTypesReco.muFromPhys) / float(bkgTypesReco.lepFromPhys) << "%)]" << endl;
        cout << endl;

        cout << "\t# of others:\t\t" << bkgTypes.others << "\t(" << 100 * float(bkgTypes.others) / float(nEv) << "%)"
             << "\t\t" << bkgTypesReco.others << "\t(" << 100 * float(bkgTypesReco.others) / float(nJJM) << "%)" << endl;
        cout << "\t\t Have electron:\t  [" << bkgTypes.others_haveEle << "\t  (" << 100 * float(bkgTypes.others_haveEle) / float(bkgTypes.others) << "%)]"
             << "\t\t  [" << bkgTypesReco.others_haveEle << "\t  (" << 100 * float(bkgTypesReco.others_haveEle) / float(bkgTypesReco.others) << "%)]" << endl;
        cout << "\t\t No electron:  \t  [" << bkgTypes.others - bkgTypes.others_haveEle << "\t  (" << 100 * float(bkgTypes.others - bkgTypes.others_haveEle) / float(bkgTypes.others) << "%)]"
             << "\t  [" << bkgTypesReco.others - bkgTypesReco.others_haveEle << "\t  (" << 100 * float(bkgTypesReco.others - bkgTypesReco.others_haveEle) / float(bkgTypesReco.others) << "%)]" << endl;
        cout << "\t---------------------------------------------------------------------------------" << endl;
    }

    tr.Write();
    fea.Close();
    cout << "Writing to: " << outputFile << "\n\n";
}
