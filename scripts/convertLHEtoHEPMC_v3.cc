// main11.cc is a part of the PYTHIA event generator.
// Copyright (C) 2021 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: basic usage; LHE file;

// This is a simple test program.
// It illustrates how Les Houches Event File input can be used in Pythia8.
// It uses the ttsample.lhe(.gz) input file, the latter only with 100 events.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
using namespace Pythia8;

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "1. LHEF" << endl
             << "2. HEPMC" << endl;
    }

    // Specify file where HepMC events will be stored.
    HepMC::Pythia8ToHepMC ToHepMC;

    // Specify file where HepMC events will be stored.
    HepMC::IO_GenEvent ascii_io(argv[2], std::ios::out);

    // Generator. We here stick with default values, but changes
    // could be inserted with readString or readFile.
    Pythia pythia;

    // Read input LHE File
    string iPath = string(argv[1]);

    // Initialize Les Houches Event File run. List initialization information.
    pythia.readString("Beams:frameType = 4");
    pythia.readString("LesHouches:matchInOut = off");
    // pythia.readString("PartonLevel:ISR = off");
    pythia.readString("ProcessLevel:all = off");
    // pythia.readString("PartonLevel:all = off");
    // pythia.readString("HadronLevel:all = on");
    pythia.settings.word("Beams:LHEF", iPath);
    pythia.init();

    // Begin event loop; generate until none left in input file.
    for (int iEvent = 0;; ++iEvent) {
        // Generate events, and check whether generation failed.
        if (!pythia.next()) {
            // If failure because reached end of file then exit event loop.
            if (pythia.info.atEndOfFile()) break;
        }

        HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
        ToHepMC.fill_next_event(pythia, hepmcevt);
        ascii_io << hepmcevt;
        delete hepmcevt;
    }

    // Give statistics. Print histogram.
    // pythia.stat();

    // Done.
    return 0;
}
