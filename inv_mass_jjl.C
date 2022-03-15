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

void inv_mass_jjl(Int_t num_test = 0) {
    // const char* inputFile = "./background_pythia8_events.root";
    const char* inputFile = "./signal_pythia8_events.root";

    const char* outputFile = "./dummy.root";

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
    TClonesArray* branchJet = treeReader->UseBranch("VLCjetR12N2");
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
    for (Int_t i_en = 0; i_en < num_test; i_en++) {
        // progress
        cout << endl;
        if ((i_en % 1) == 0) cout << "Reconstruction Progress: " << i_en << "/" << numberOfEntries << endl;

        treeReader->ReadEntry(i_en);  // reading the entry

        iFinalStates iFSTrue;
        Int_t passing = 0;
        passing = ClassifySingal(branchParticle, &iFSTrue);
        if (passing == 0) continue;

        GenParticle* lepParticle = (GenParticle*)branchParticle->At(iFSTrue.iLep);
        GenParticle* jet1Particle = (GenParticle*)branchParticle->At(iFSTrue.iJet1);
        GenParticle* jet2Particle = (GenParticle*)branchParticle->At(iFSTrue.iJet2);

        TLorentzVector lepTrue, jet1True, jet2True, NTrue;
        lepTrue.SetPtEtaPhiE(lepParticle->PT, lepParticle->Eta, lepParticle->Phi, lepParticle->E);
        jet1True.SetPtEtaPhiE(jet1Particle->PT, jet1Particle->Eta, jet1Particle->Phi, jet1Particle->E);
        jet2True.SetPtEtaPhiE(jet2Particle->PT, jet2Particle->Eta, jet2Particle->Phi, jet2Particle->E);

        NTrue = lepTrue + jet1True + jet2True;
        /*
        cout << "=== Truth ===" << endl;
        cout << "N1 Mass: " << NTrue.M() << endl;
        cout << "Lep p: " << lepTrue.P() << "; lep eta: " << lepTrue.Eta() << endl;
        cout << " ============" << endl;
        cout << endl;
        */
        iFinalStates iFS;
        iFS = FindFinalStatesIndex(branchTrack, branchJet);
        if (iFS.foundAll == 0) continue;

        Track* lepTrack = (Track*)branchTrack->At(iFS.iLep);
        Jet* jet1Jet = (Jet*)branchJet->At(iFS.iJet1);
        Jet* jet2Jet = (Jet*)branchJet->At(iFS.iJet2);

        TLorentzVector lep, jet1, jet2, N;

        lep.SetPtEtaPhiM(lepTrack->PT, lepTrack->Eta, lepTrack->Phi, iFS.mLep);
        jet1.SetPtEtaPhiM(jet1Jet->PT, jet1Jet->Eta, jet1Jet->Phi, jet1Jet->Mass);
        jet2.SetPtEtaPhiM(jet2Jet->PT, jet2Jet->Eta, jet2Jet->Phi, jet2Jet->Mass);

        cout << lepTrue.Pt() << "; " << lep.Pt() << endl;
    }
    cout << "\n\n===== Done! =====\n\n";
}
