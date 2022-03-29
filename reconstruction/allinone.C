
#include <cmath>
#include <iostream>
#include <random>
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

void allinone(
    const string type,
    const Bool_t save = false,
    Int_t num_test = 0) {
    const char* inputFile;
    const char* outputFile;

    outputFile = "../features/dummy.root";
    if (type == "s3_1") {
        // inputFile = "../data/signal_E-3TeV_N-1TeV.root";
        inputFile = "../data/signal_E-3TeV_N-1TeV_testJet.root";
        outputFile = "../features/signal_reco_E-3TeV_N-1TeV.root";
    } else if (type == "s10_1") {
        // inputFile = "../data/signal_E-10TeV_N-1TeV.root";
        inputFile = "../data/signal_E-10TeV_N-1TeV_testJet.root";
        outputFile = "../features/signal_reco_E-10TeV_N-1TeV.root";
    } else if (type == "b3") {
        // inputFile = "../data/background_inclusive_E-3TeV.root";
        inputFile = "../data/background_inclusive_E-3TeV_testJet.root";
        outputFile = "../features/background_reco_E-3TeV.root";
    } else if (type == "b10") {
        // inputFile = "../data/background_inclusive_E-10TeV.root";
        inputFile = "../data/background_inclusive_E-10TeV_testJet.root";
        outputFile = "../features/background_reco_E-10TeV.root";
    } else {
        cout << "Wrong input" << endl;
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
    // TClonesArray* branchJet = treeReader->UseBranch("VLCjetR12N2");
    // TClonesArray* branchJet = treeReader->UseBranch("GenJet");
    TClonesArray* branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
    TClonesArray* branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");

    GenParticle* particle;
    Track* track;
    Jet* jet;
    Tower* eflowphoton;
    Tower* eflowneutralhadron;

    TFile fea(outputFile, "recreate");
    TTree tr("t", "features");
    Features* features = new Features;
    tr.Branch("features", &features);

    if (num_test == 0) num_test = numberOfEntries;
    Int_t nEv = 0;
    Int_t nFS = 0;
    Int_t nLepEta = 0;
    Int_t nLepPt = 0;
    Int_t nJetPt = 0;

    Float_t lepEtaCut = 2.5;
    Float_t lepPtCut = 200;
    Float_t jetPtCut = 50;

    for (Int_t i_en = 0; i_en < num_test; i_en++) {
        // progress
        if ((i_en % 1000) == 0) cout << "Reconstruction Progress: " << i_en << "/" << numberOfEntries << endl;

        treeReader->ReadEntry(i_en);  // reading the entry

        //========================================================================
        //=======================   Classifiy Event Type   =======================
        //========================================================================
        iFinalStates iFSTrue;
        TLorentzVector lepTrue, jet1True_, jet2True_, NTrue;
        Int_t passing = 0;
        if (type.at(0) == 's') {
            passing = ClassifySingal(branchParticle, &iFSTrue, &lepTrue, &jet1True_, &jet2True_, &NTrue);
        } else {
            passing = 1;
        }

        if (passing == 0) continue;
        nEv += 1;

        TLorentzVector jjTrue;
        jjTrue = jet1True_ + jet2True_;

        //========================================================================
        //=======================      Reconstruction      =======================
        //========================================================================
        iFinalStates iFS;
        // finding final states: lepton + 2jets
        // iFS = FindFinalStatesIndex(branchTrack, branchJet);
        iFS = FindFinalStatesIndex(branchTrack, branchVLC1Jet, branchVLC2Jet);
        if (iFS.foundAll == 0) continue;
        nFS += 1;
        // if (iFS.nJets != 2) continue;

        TLorentzVector lep, jet1, jet2, jj, N;

        Track* lepTrack = (Track*)branchTrack->At(iFS.iLep);
        lep.SetPtEtaPhiM(lepTrack->PT, lepTrack->Eta, lepTrack->Phi, iFS.mLep);

        if (iFS.nJets == 1) {
            Jet* jet1Jet = (Jet*)branchVLC1Jet->At(iFS.iJet1);
            jet1.SetPtEtaPhiM(jet1Jet->PT, jet1Jet->Eta, jet1Jet->Phi, jet1Jet->Mass);
            jj = jet1;
        } else if (iFS.nJets == 2) {
            Jet* jet1Jet = (Jet*)branchVLC2Jet->At(iFS.iJet1);
            Jet* jet2Jet = (Jet*)branchVLC2Jet->At(iFS.iJet2);
            jet1.SetPtEtaPhiM(jet1Jet->PT, jet1Jet->Eta, jet1Jet->Phi, jet1Jet->Mass);
            jet2.SetPtEtaPhiM(jet2Jet->PT, jet2Jet->Eta, jet2Jet->Phi, jet2Jet->Mass);
            jj = jet1 + jet2;
        }
        N = jj + lep;

        //========================================================================
        //=======================           Cuts           =======================
        //========================================================================
        // if (not(abs(lep.Eta()) <= lepEtaCut)) continue;
        nLepEta += 1;

        // if (not(lep.Pt() >= lepPtCut)) continue;
        nLepPt += 1;

        // if (not(abs(jj.M() - 80.379) < 10)) continue;

        // if (not(jet1.Pt() >= jetPtCut && jet2.Pt() >= jetPtCut)) continue;
        nJetPt += 1;

        //========================================================================
        //=======================         Features         =======================
        //========================================================================

        TLorentzVector jet1True, jet2True;
        Float_t DeltaRjj, DeltaRjjl;
        if (iFS.nJets == 2) {
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

        } else if (iFS.nJets == 1) {
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
        if (iFS.nJets == 2) features->DeltaRjj = DeltaRjj;
        if (iFS.nJets == 2) features->DeltaRjjl = DeltaRjjl;

        features->nJets = iFS.nJets;
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

        tr.Fill();
    }

    cout << "\n\n===== Done! =====\n\n";

    cout << "# of targeted events in truth level:\t\t" << nEv << endl;
    cout << "# after identified final states:\t\t" << nFS << "\t(" << 100 * float(nFS) / float(nEv) << "%)" << endl;
    cout << "# after lepton eta cut (eta(l) <= " << float(lepEtaCut) << "):\t\t" << nLepEta << "\t(" << 100 * float(nLepEta) / float(nFS) << "%)" << endl;
    cout << "# after lepton pt cut (pT(l)) >= " << float(lepPtCut) << "):\t\t" << nLepPt << "\t(" << 100 * float(nLepPt) / float(nLepEta) << "%)" << endl;
    cout << "# after jet pt cut (pT(j1, j2) >= " << float(jetPtCut) << "):\t\t" << nJetPt << "\t(" << 100 * float(nJetPt) / float(nLepPt) << "%)" << endl;
    cout << endl;
    cout << "Total eff.:\t\t\t\t\t" << nJetPt << "/" << nEv << "\t(" << 100 * float(nJetPt) / float(nEv) << "%)" << endl;
    cout << endl;

    tr.Write();
    fea.Close();
}
