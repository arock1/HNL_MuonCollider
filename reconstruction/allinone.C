
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>
using namespace std;

#include "FeatureClass.C"
#include "FinalStatesClass.C"
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
        Int_t pos = type.find("_");
        if (type.at(0) == 'i') {  // inclusive signal (e.g. "i3_01": inclusive signal with E=3TeV, mN=0.1TeV)
            inputFile_st = "../data/signal_E-";
        } else if (type.at(0) == 's') {  // s-channel
            inputFile_st = "../data/signal_schannel_E-";
        } else if (type.at(0) == 't') {  // t-channel
            inputFile_st = "../data/signal_tchannel_E-";
        }
        inputFile_st += type.substr(1, pos - 1);
        inputFile_st += "TeV_N-";
        inputFile_st += type.substr(pos + 1, -1);
        inputFile_st += "TeV_CustomizedJet.root";

        if (type.at(0) == 'i') {
            outputFile_st = "../features/signal_reco_E-";
        } else if (type.at(0) == 's') {
            outputFile_st = "../features/signal_schannel_reco_E-";
        } else if (type.at(0) == 't') {
            outputFile_st = "../features/signal_tchannel_reco_E-";
        }
        outputFile_st += type.substr(1, pos - 1);
        outputFile_st += "TeV_N-";
        outputFile_st += type.substr(pos + 1, -1);
        outputFile_st += "TeV.root";
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
    if (not save) outputFile = "../dummy.root";

    cout << "\nReading: " << inputFile << "\n\n";
    if (gSystem->AccessPathName(inputFile)) {
        cout << "inputFile: " << inputFile << " does not exist" << endl;
        return;
    }

    // Load lib, and read data
    gSystem->Load("libDelphes");
    TChain chain("Delphes");
    chain.Add(inputFile);
    ExRootTreeReader* treeReader = new ExRootTreeReader(&chain);
    Int_t numberOfEntries = treeReader->GetEntries();

    // clone the branches
    TClonesArray* branchParticle = treeReader->UseBranch("Particle");
    TClonesArray* branchTrack = treeReader->UseBranch("Track");
    TClonesArray* branchKTJet = treeReader->UseBranch("KTjet");
    TClonesArray* branchVLC1Jet = treeReader->UseBranch("VLCjetR12N1");
    TClonesArray* branchVLC2Jet = treeReader->UseBranch("VLCjetR02N2");
    TClonesArray* branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
    TClonesArray* branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");

    GenParticle* particle;
    Track* track;
    Jet* jet;
    Tower* eflowphoton;
    Tower* eflowneutralhadron;

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

    Float_t lepEtaCut = 2.5;  // lep eta cut
    Float_t lepPtCut = 100;   // lep pt cut
    Float_t jjPtCut = 100;    // jj pt cut

    Float_t mW = 80.379;
    Float_t mWWidth = 2.085;

    Float_t jjMLowCut = mW - 5 * mWWidth;   // jj mass cut
    Float_t jjMHighCut = mW + 5 * mWWidth;  // jj mass cut

    for (Int_t i_en = 0; i_en < num_test; i_en++) {
        // progress
        if ((i_en % 1000) == 0) cout << "Reconstruction Progress: " << i_en << "/" << numberOfEntries << "\r";
        cout.flush();

        treeReader->ReadEntry(i_en);  // reading the entry

        //========================================================================
        //=======================   Classifiy Event Type   =======================
        //========================================================================
        iFinalStates iFSTrue;
        TLorentzVector lepTrue, jet1True_, jet2True_, NTrue;
        Int_t passing = 0;
        if (type.at(0) == 'i' || type.at(0) == 's' || type.at(0) == 't') {
            passing = ClassifySingal(branchParticle, &iFSTrue, &lepTrue, &jet1True_, &jet2True_, &NTrue);
        } else {
            passing = ClassifyiBackground(branchParticle);
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
        iFS = FindFinalStatesIndex(branchTrack, branchVLC1Jet, branchVLC2Jet);

        if (iFS.foundAll == 0) continue;
        nFS += 1;

        TLorentzVector lep, jet1, jet2, jj, N;

        // if there is only one jet, then it is corresponding to the first index
        // and the W is this only jet
        // if there are two jets, then corresponding to the first and second index
        // and the W is their sum
        if (iFS.iJets.size() == 1) {
            Jet* jet1Jet = (Jet*)branchVLC1Jet->At(iFS.iJets[0]);
            jet1.SetPtEtaPhiM(jet1Jet->PT, jet1Jet->Eta, jet1Jet->Phi, jet1Jet->Mass);
            jj = jet1;
        } else if (iFS.iJets.size() == 2) {
            Jet* jet1Jet = (Jet*)branchVLC2Jet->At(iFS.iJets[0]);
            Jet* jet2Jet = (Jet*)branchVLC2Jet->At(iFS.iJets[1]);
            jet1.SetPtEtaPhiM(jet1Jet->PT, jet1Jet->Eta, jet1Jet->Phi, jet1Jet->Mass);
            jet2.SetPtEtaPhiM(jet2Jet->PT, jet2Jet->Eta, jet2Jet->Phi, jet2Jet->Mass);
            jj = jet1 + jet2;
        }

        // loop over all the leptons
        Float_t ptLepMax = -99999;
        Int_t chargeLep;
        Int_t foundLep = 0;
        for (Int_t iLep = 0; iLep < iFS.iLeps.size(); iLep++) {
            Track* lepTrack = (Track*)branchTrack->At(iFS.iLeps[iLep]);
            // eta cut
            if (not(abs(lepTrack->Eta) <= lepEtaCut)) continue;

            // exclude the lepton, if it is inside the jet cone
            // Float_t DeltaRjjl;
            // Float_t jjlEtaDiff = jj.Eta() - lepTrack->Eta;
            // Float_t jjlPhiDiff = jj.Phi() - lepTrack->Phi;
            // DeltaRjjl = pow(jjlEtaDiff * jjlEtaDiff + jjlPhiDiff * jjlPhiDiff, 0.5);
            // if (iFS.iJets.size() == 1 && DeltaRjjl < 1.2) continue;
            // if (iFS.iJets.size() == 2 && DeltaRjjl < 0.2) continue;

            // Int_t nTracks = branchTrack->GetEntries();
            // Float_t closestDeltaR = 99999;
            // for (Int_t itr = 0; itr < nTracks; itr++) {
            // if (itr == iFS.iLeps[iLep]) continue;
            // Float_t DeltaRanyl;
            // Track* anyTrack = (Track*)branchTrack->At(itr);
            // Float_t anylEtaDiff = anyTrack->Eta - lepTrack->Eta;
            // Float_t anylPhiDiff = anyTrack->Phi - lepTrack->Phi;
            // DeltaRanyl = pow(anylEtaDiff * anylEtaDiff + anylPhiDiff * anylPhiDiff, 0.5);
            // if (DeltaRanyl < closestDeltaR) closestDeltaR = DeltaRanyl;
            //}
            //// cout << closestDeltaR << endl;
            //// cout << endl;
            // if (closestDeltaR < 1) continue;

            // store the lepton with max. pT among all
            if (lepTrack->PT > ptLepMax) {
                chargeLep = lepTrack->Charge;
                foundLep = 1;
                ptLepMax = lepTrack->PT;
                lep.SetPtEtaPhiM(lepTrack->PT, lepTrack->Eta, lepTrack->Phi, iFS.mLeps[iLep]);
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

        if (not(jj.M() >= jjMLowCut && jj.M() <= jjMHighCut)) continue;
        nJJM += 1;

        //========================================================================
        //=======================         Features         =======================
        //========================================================================

        TLorentzVector jet1True,
            jet2True;
        Float_t DeltaRjj, DeltaRjjl;
        if (iFS.iJets.size() == 2) {
            if (abs(jet1.Eta() - jet1True_.Eta()) < abs(jet1.Eta() - jet2True_.Eta())) {
                jet1True = jet1True_;
                jet2True = jet2True_;
            } else {
                jet1True = jet2True_;
                jet2True = jet1True_;
            }
            Float_t jjEtaDiff = jet1.Eta() - jet2.Eta();
            Float_t jjPhiDiff = jet1.Phi() - jet2.Phi();
            DeltaRjj = pow(jjEtaDiff * jjEtaDiff + jjPhiDiff * jjPhiDiff, 0.5);

        } else if (iFS.iJets.size() == 1) {
            jet1True = jet1True_;
            jet2True = jet2True_;
        }

        Float_t jjlEtaDiff = jj.Eta() - lep.Eta();
        Float_t jjlPhiDiff = jj.Phi() - lep.Phi();
        DeltaRjjl = pow(jjlEtaDiff * jjlEtaDiff + jjlPhiDiff * jjlPhiDiff, 0.5);

        Float_t DeltaRjjTrue, DeltaRjjlTrue;
        Float_t jjEtaDiffTrue = jet1True.Eta() - jet2True.Eta();
        Float_t jjPhiDiffTrue = jet1True.Phi() - jet2True.Phi();
        DeltaRjjTrue = pow(jjEtaDiffTrue * jjEtaDiffTrue + jjPhiDiffTrue * jjPhiDiffTrue, 0.5);
        Float_t jjlEtaDiffTrue = jjTrue.Eta() - lepTrue.Eta();
        Float_t jjlPhiDiffTrue = jjTrue.Phi() - lepTrue.Phi();
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

        features->chargeLep = chargeLep;
        if (type.at(0) != 'b') {
            GenParticle* lTrue = (GenParticle*)branchParticle->At(iFSTrue.iLeps[0]);
            features->chargeLepTrue = lTrue->Charge;
        }

        tr.Fill();
    }
    cout << "Reconstruction Progress: " << numberOfEntries << "/" << numberOfEntries << "\r";

    cout << "---------------------------------------------------------------------------------" << endl;
    cout << "==============================       Cut eff.      ==============================\n";
    cout << "---------------------------------------------------------------------------------" << endl;

    cout << "# of targeted events in truth level:\t\t\t" << nEv << endl;
    cout << "# after identified final states:\t\t\t" << nFS << "\t(" << 100 * float(nFS) / float(nEv) << "%)" << endl;
    cout << "# after lepton eta cut \t(eta(l) <= " << float(lepEtaCut) << "):\t\t" << nLepEta << "\t(" << 100 * float(nLepEta) / float(nFS) << "%)" << endl;
    cout << "# after lepton pt cut \t(pT(l) >= " << float(lepPtCut) << "):\t\t\t" << nLepPt << "\t(" << 100 * float(nLepPt) / float(nLepEta) << "%)" << endl;
    cout << "# after jj pt cut \t(pT(jj) >= " << float(jjPtCut) << "):\t\t" << nJJPt << "\t(" << 100 * float(nJJPt) / float(nLepPt) << "%)" << endl;
    cout << "# after jj mass cut \t(" << jjMLowCut << " <= m(jj) <= " << float(jjMHighCut) << "):\t" << nJJM << "\t(" << 100 * float(nJJM) / float(nJJPt) << "%)" << endl;
    // cout << "# after jet pt cut (pT(j1, j2) >= " << float(jetPtCut) << "):\t\t" << nJetPt << "\t(" << 100 * float(nJetPt) / float(nLepPt) << "%)" << endl;
    cout << "---------------------------------------------------------------------------------" << endl;
    cout << "\t\t\t\t\tTotal eff.:\t" << nJJM << " / " << nEv << "\t(" << 100 * float(nJJM) / float(nEv) << " %) " << endl;
    cout << "---------------------------------------------------------------------------------" << endl;
    cout << "\n\n\n";

    tr.Write();
    fea.Close();
    cout << "\nWriting: " << outputFile << "\n\n";
}
