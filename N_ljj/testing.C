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

#include "Rtypes.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "reconstruction/FeatureClass.C"
#include "reconstruction/FinalStatesClass.C"
#include "reconstruction/Geometry.C"
#include "reconstruction/Reconstructor.C"
#include "reconstruction/SignalBackgroundClassifier.C"

R__ADD_LIBRARY_PATH($DELPHES)
R__LOAD_LIBRARY(libDelphes)

void testing(
    const string type,
    const Bool_t save = false,
    Int_t num_test = 0,
    Int_t print_detail = 1,
    Int_t all_on = 1) {
    // formmating the input output files
    string inputFile_st;
    string outputFile_st;

    Int_t foundFiles;
    // foundFiles = getFileNames(type, &inputFile_st, &outputFile_st);
    inputFile_st = "./data/detector/bg_MuMu_qqll_E-3.root";
    outputFile_st = "./testing_qqll.root";

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

        TLorentzVector l1_, l1, l2_, l2, j1, j2, jj;
        Int_t nl = 0;
        Int_t nj = 0;

        Int_t nParticles = branchParticle->GetEntries();
        for (Int_t ip = 0; ip < nParticles; ip++) {
            GenParticle* particle = (GenParticle*)branchParticle->At(ip);
            if (particle->M1 == -1) continue;
            GenParticle* particleM = (GenParticle*)branchParticle->At(particle->M1);

            if (not(abs(particle->PID) == 13 || abs(particle->PID) == 11 || abs(particle->PID) / 10 == 0)) continue;
            if (abs(particleM->PID) != 13) continue;

            if (abs(particle->PID) == 13 || abs(particle->PID) == 11) {
                if (nl == 0) l1_.SetPtEtaPhiE(particle->PT, particle->Eta, particle->Phi, particle->E);
                if (nl == 1) l2_.SetPtEtaPhiE(particle->PT, particle->Eta, particle->Phi, particle->E);
                nl += 1;
            }
            if (abs(particle->PID) / 10 == 0) {
                if (nj == 0) j1.SetPtEtaPhiE(particle->PT, particle->Eta, particle->Phi, particle->E);
                if (nj == 1) j2.SetPtEtaPhiE(particle->PT, particle->Eta, particle->Phi, particle->E);
                nj += 1;
            }
        }
        if (l1_.E() > l2_.E()) {
            l1 = l1_;
            l2 = l2_;
        } else {
            l1 = l2_;
            l2 = l1_;
        }

        jj = j1 + j2;

        TLorentzVector jjl = jj + l1;
        // cout << " jjl: " << jjl.M() << "\n";
        // cout << endl;
        TLorentzVector jjl2 = jj + l2;

        features->mN = jjl.M();
        features->mJJl2 = jjl2.M();
        features->ptLep = l1.Pt();
        features->ptl2 = l2.Pt();
        features->pzl2 = l2.Pz();
        features->El2 = l2.E();
        features->etaLep = l1.Eta();
        features->pzLep = l1.Pz();
        features->ptJJ = jj.Pt();
        features->mJJ = jj.M();
        features->etaJJ = jj.Eta();
        Float_t jjlEtaDiff = jj.Eta() - l1.Eta();
        Float_t jjlPhiDiff = deltaPhi(jj.Phi(), l1.Phi());
        Float_t DeltaRjjl = pow(jjlEtaDiff * jjlEtaDiff + jjlPhiDiff * jjlPhiDiff, 0.5);
        features->DeltaRjjl = DeltaRjjl;
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
