

// ==========================================================================
// || This Code Does:
// || Store the truth level lepton info, without reconstruction
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

    if (type.at(0) == 'i' || type.at(0) == 's' || type.at(0) == 't') {  // signals
        cout << "Processing Signal data" << endl;
        string fermionType;
        Int_t pos = type.find("_");
        if (type.at(0) == 't') {  // t-channel
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

        if (type.at(0) == 't') {
            outputFile_st = "../features/truthLevelLepE_E-";
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

void truthLevelLeptonEnergy(
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

    // book feature storing tree and file
    TFile fea(outputFile, "recreate");
    TTree tr("t", "features");
    Features* features = new Features;
    tr.Branch("features", &features);

    if (num_test == 0) num_test = numberOfEntries;
    // loop the events
    for (Int_t i_en = 0; i_en < num_test; i_en++) {
        // progress
        if ((i_en % 1000) == 0) cout << "Reconstruction Progress: " << i_en << "/" << numberOfEntries << "\r";
        cout.flush();

        treeReader->ReadEntry(i_en);  // reading the entry
        //========================================================================
        //=======================   Classifiy Event Type   =======================
        //========================================================================
        iFinalStates iFSTrue;                                 // indeces of the true level final states
        TLorentzVector lepTrue, jet1True_, jet2True_, NTrue;  // truth level lepton, jets, and HNL
        Int_t typeLepTrue;
        Int_t passing = 0;
        if (type.at(1) == 'M') {  // Majorana
            passing = ClassifySingal(branchParticle, &iFSTrue, &NTrue, 9900012);
        } else if (type.at(1) == 'D') {  // Dirac
            passing = ClassifySingal(branchParticle, &iFSTrue, &NTrue, 9990012);
        }
        lepTrue = iFSTrue.iLeps[0];
        if (passing == 0) continue;

        features->ELepTrue = lepTrue.E();
        features->pxLepTrue = lepTrue.Px();
        features->pyLepTrue = lepTrue.Py();
        features->pzLepTrue = lepTrue.Pz();
        features->pxNTrue = NTrue.Px();
        features->pyNTrue = NTrue.Py();
        features->pzNTrue = NTrue.Pz();
        tr.Fill();
        // cout << " lepton energy : " << lepTrue.E() << "\n";
    }
    tr.Write();
    fea.Close();
    cout << "Reconstruction Progress: " << numberOfEntries << "/" << numberOfEntries << "\n\n";
}
