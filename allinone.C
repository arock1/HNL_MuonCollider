
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

void inv_mass_jjl(
    const string type,
    const Bool_t save = false,
    Int_t num_test = 0) {
    const char* inputFile;
    const char* outputFile;

    if (type == "s1") {
        inputFile = "./signal_pythia8_events_1TeV.root";
        outputFile = "./signal_reco_1TeV.root";
    } else if (type == "b1") {
        inputFile = "./background_pythia8_events_1TeV.root";
        outputFile = "./background_reco_1TeV.root";
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
    // TClonesArray* branchJet = treeReader->UseBranch("KTjet");
    TClonesArray* branchJet = treeReader->UseBranch("VLCjetR05N2");
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
        if (type == "s1") {
            passing = ClassifySingal(branchParticle, &iFSTrue, &lepTrue, &jet1True_, &jet2True_);
        } else if (type == "b1") {
            passing = 1;
        }

        if (passing == 0) continue;
        nEv += 1;

        NTrue = lepTrue + jet1True_ + jet2True_;
        TLorentzVector jjTrue;
        jjTrue = jet1True_ + jet2True_;

        //========================================================================
        //=======================      Reconstruction      =======================
        //========================================================================
        iFinalStates iFS;
        // finding final states: lepton + 2jets
        iFS = FindFinalStatesIndex(branchTrack, branchJet);
        if (iFS.foundAll == 0) continue;
        nFS += 1;

        Track* lepTrack = (Track*)branchTrack->At(iFS.iLep);
        Jet* jet1Jet = (Jet*)branchJet->At(iFS.iJet1);
        Jet* jet2Jet = (Jet*)branchJet->At(iFS.iJet2);

        TLorentzVector lep, jet1, jet2, N;

        lep.SetPtEtaPhiM(lepTrack->PT, lepTrack->Eta, lepTrack->Phi, iFS.mLep);
        jet1.SetPtEtaPhiM(jet1Jet->PT, jet1Jet->Eta, jet1Jet->Phi, jet1Jet->Mass);
        jet2.SetPtEtaPhiM(jet2Jet->PT, jet2Jet->Eta, jet2Jet->Phi, jet2Jet->Mass);

        N = jet1 + jet2 + lep;
        TLorentzVector jj;
        jj = jet1 + jet2;

        //========================================================================
        //=======================           Cuts           =======================
        //========================================================================
        // require in same direction
        // if (not(lep.Px() * jet1.Px() + lep.Py() * jet1.Py() + lep.Pz() * jet1.Pz() >= 0)) continue;
        // if (not(lep.Px() * jet2.Px() + lep.Py() * jet2.Py() + lep.Pz() * jet2.Pz() >= 0)) continue;
        // if (not(jet2.Px() * jet1.Px() + jet2.Py() * jet1.Py() + jet2.Pz() * jet1.Pz() >= 0)) continue;
        if (not(abs(lep.Eta()) <= 2.5)) continue;
        nLepEta += 1;
        if (not(lep.Pt() >= 200)) continue;
        nLepPt += 1;
        if (not(jet1.Pt() >= 100 || jet2.Pt() >= 100)) continue;
        if (not(jet1.Pt() >= 50 && jet2.Pt() >= 50)) continue;
        nJetPt += 1;

        /*
        cout << jet1.Eta() << "; " << jet1True.Eta() << endl;
        cout << jet2.Eta() << "; " << jet2True.Eta() << endl;
        cout << jet1.E() * jet1.E() - jet1.P() * jet1.P() << "; " << jet2.E() * jet2.E() - jet2.P() * jet2.P() << endl;
        cout << jet1True.E() * jet1True.E() - jet1True.P() * jet1True.P() << "; " << jet2True.E() * jet2True.E() - jet2True.P() * jet2True.P() << endl;
        cout << endl;
        cout << jj.Eta() << "; " << jjTrue.Eta() << endl;
        cout << endl;
        cout << N.M() << "; " << NTrue.M() << endl;

        cout << endl;
        cout << endl;
        cout << endl;
        */
        //========================================================================
        //=======================         Features         =======================
        //========================================================================
        TLorentzVector jet1True, jet2True;
        if (abs(jet1.Eta() - jet1True_.Eta()) < abs(jet1.Eta() - jet2True_.Eta())) {
            jet1True = jet1True_;
            jet2True = jet2True_;
        } else {
            jet1True = jet2True_;
            jet2True = jet1True_;
        }

        Float_t DeltaRjj, DeltaRjjl;
        Float_t DeltaRjjTrue, DeltaRjjlTrue;

        Float_t jjEtaDiff = jet1.Eta() - jet2.Eta();
        Float_t jjPhiDiff = jet1.Phi() - jet2.Phi();
        DeltaRjj = pow(jjEtaDiff * jjEtaDiff + jjPhiDiff * jjPhiDiff, 0.5);
        Float_t jjlEtaDiff = jj.Eta() - lep.Eta();
        Float_t jjlPhiDiff = jj.Phi() - lep.Phi();
        DeltaRjjl = pow(jjlEtaDiff * jjlEtaDiff + jjlPhiDiff * jjlPhiDiff, 0.5);

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
        features->ptJet1 = jet1.Pt();
        features->etaJet1 = jet1.Eta();
        features->phiJet1 = jet1.Phi();
        features->ptJet2 = jet2.Pt();
        features->etaJet2 = jet2.Eta();
        features->phiJet2 = jet2.Phi();
        features->DeltaRjj = DeltaRjj;
        features->DeltaRjjl = DeltaRjjl;
        features->mN = N.M();

        features->ptLepTrue = lepTrue.Pt();
        features->etaLepTrue = lepTrue.Eta();
        features->phiLepTrue = lepTrue.Phi();
        features->ptJet1True = jet1True.Pt();
        features->etaJet1True = jet1True.Eta();
        features->phiJet1True = jet1True.Phi();
        features->ptJet2True = jet2True.Pt();
        features->etaJet2True = jet2True.Eta();
        features->phiJet2True = jet2True.Phi();
        features->DeltaRjjTrue = DeltaRjjTrue;
        features->DeltaRjjlTrue = DeltaRjjlTrue;

        tr.Fill();
    }

    cout << "\n\n===== Done! =====\n\n";

    cout << "# of targeted events in truth level: " << nEv << endl;
    cout << "# after identified final states:     " << nFS << "\t(" << 100 * float(nFS) / float(nEv) << "%)" << endl;
    cout << "# after lepton eta cut:              " << nLepEta << "\t(" << 100 * float(nLepEta) / float(nFS) << "%)" << endl;
    cout << "# after lepton pt cut:               " << nLepPt << "\t(" << 100 * float(nLepPt) / float(nLepEta) << "%)" << endl;
    cout << "# after jet pt cut:                  " << nJetPt << "\t(" << 100 * float(nJetPt) / float(nLepPt) << "%)" << endl;
    tr.Write();
    fea.Close();
}
