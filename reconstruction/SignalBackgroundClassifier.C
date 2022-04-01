#include <iostream>
using namespace std;

#include "FinalStatesClass.C"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "classes/DelphesClasses.h"

//{{{
Int_t ClassifySingal(TClonesArray* branchParticle, iFinalStates* iFSTrue, TLorentzVector* lepTrue, TLorentzVector* Q1True, TLorentzVector* Q2True, TLorentzVector* NTrue) {
    GenParticle* particle;
    GenParticle* particle0;
    GenParticle* particle0M;
    GenParticle* particle0MM;

    GenParticle* particleNM1;
    GenParticle* particleNM2;

    Int_t iNTrue = 99999;
    Int_t foundN = 0;
    Int_t iWTrue = 99999;
    Int_t foundW = 0;
    Int_t iLepTrue = 99999;
    Int_t foundLep = 0;
    Int_t iQ1True = 99999;  // down type quark
    Int_t foundQ1 = 0;
    Int_t iQ2True = 99999;  // up type quark
    Int_t foundQ2 = 0;

    Int_t nParticles = branchParticle->GetEntries();
    // cout << "\n\n\n Num of p: " << nParticles << endl;

    for (Int_t iN = 0; iN < nParticles; iN++) {
        foundN = 0;
        particle = (GenParticle*)branchParticle->At(iN);
        // cout << particle->PID << endl;

        if (abs(particle->PID) == 9900012) {
            GenParticle* particle_mu1 = (GenParticle*)branchParticle->At(particle->M1);
            GenParticle* particle_mu2 = (GenParticle*)branchParticle->At(particle->M2);
            // cout << particle_mu1->Eta << "; " << particle_mu1->PID << endl;
            // cout << particle_mu2->Eta << "; " << particle_mu2->PID << endl;
            // cout << endl;
            particleNM1 = (GenParticle*)branchParticle->At(particle->M1);
            particleNM2 = (GenParticle*)branchParticle->At(particle->M2);
            // cout << "\nN1 Mother1: " << particleNM1->PID << "\nN1 Mother2: " << particleNM2->PID << "\n";

            // cout << "found N1" << endl;
            foundN = 1;
            iNTrue = iN;

            for (Int_t i = 0; i < nParticles; i++) {
                particle0 = (GenParticle*)branchParticle->At(i);
                if ((abs(particle0->PID) == 11 || abs(particle0->PID) == 13) && particle0->M1 == iNTrue) {
                    // cout << particle0->PID << endl;
                    //  cout << "Found Lepton: " << particle0->PID << endl;
                    foundLep = 1;
                    iLepTrue = i;
                    break;
                }
            }

            for (Int_t ii = 0; ii < nParticles; ii++) {
                particle0 = (GenParticle*)branchParticle->At(ii);
                if (particle0->M1 == -1) continue;
                Int_t iM = particle0->M1;
                particle0M = (GenParticle*)branchParticle->At(iM);
                if (abs(particle0->PID) <= 8) {
                    while (true) {
                        if (abs(particle0M->PID) == 9900012) {
                            if (abs(particle0->PID) % 2 == 1) {  // odd: down type quark
                                // cout << "down: " << particle0->PID << endl;
                                iQ1True = ii;
                                foundQ1 = 1;
                            }
                            if (abs(particle0->PID) % 2 == 0) {  // even: up type quark
                                // cout << "up  : " << particle0->PID << endl;
                                iQ2True = ii;
                                foundQ2 = 1;
                            }
                        }
                        if (foundQ1 == 1 && foundQ2 == 1) break;

                        particle0M = (GenParticle*)branchParticle->At(iM);
                        if (abs(particle0M->PID) != 24 && abs(particle0M->PID) != 9900012) break;
                        iM = particle0M->M1;
                    }
                }
                if (foundQ1 == 1 && foundQ2 == 1) break;
            }
            if (foundN == 1 && foundQ1 == 1 && foundQ2 == 1) {
                // cout << "............" << endl;
                iFinalStates iFSTrue_;
                iFSTrue_.iLeps.push_back(iLepTrue);
                iFSTrue_.iJets.push_back(iQ1True);
                iFSTrue_.iJets.push_back(iQ2True);
                *iFSTrue = iFSTrue_;

                TLorentzVector lepTrue_, Q1True_, Q2True_, NTrue_;
                GenParticle* lepParticle;
                GenParticle* Q1Particle;
                GenParticle* Q2Particle;
                GenParticle* NParticle;
                lepParticle = (GenParticle*)branchParticle->At(iLepTrue);
                Q1Particle = (GenParticle*)branchParticle->At(iQ1True);
                Q2Particle = (GenParticle*)branchParticle->At(iQ2True);
                NParticle = (GenParticle*)branchParticle->At(iNTrue);

                lepTrue_.SetPtEtaPhiE(lepParticle->PT, lepParticle->Eta, lepParticle->Phi, lepParticle->E);
                Q1True_.SetPtEtaPhiE(Q1Particle->PT, Q1Particle->Eta, Q1Particle->Phi, Q1Particle->E);
                Q2True_.SetPtEtaPhiE(Q2Particle->PT, Q2Particle->Eta, Q2Particle->Phi, Q2Particle->E);
                NTrue_.SetPtEtaPhiE(NParticle->PT, NParticle->Eta, NParticle->Phi, NParticle->E);

                *lepTrue = lepTrue_;
                *Q1True = Q1True_;
                *Q2True = Q2True_;
                *NTrue = NTrue_;

                return 1;
            }
        }
    }
    // cout << foundN1 << "; " << foundQ1 << "; " << foundQ2 << endl;

    return 0;
}
//}}}

Int_t ClassifyiBackground(TClonesArray* branchParticle) {
    GenParticle* particle;
    GenParticle* particleM;
    GenParticle* particleMM;

    Int_t nParticles = branchParticle->GetEntries();

    Int_t found2W = 0;

    Int_t iW1;
    Int_t iWMother;
    for (Int_t ip1 = 0; ip1 < nParticles; ip1++) {
        particle = (GenParticle*)branchParticle->At(ip1);
        if (particle->M1 < 0) continue;
        if (abs(particle->PID) == 13 || abs(particle->PID) == 11) {      // lepton
            particleM = (GenParticle*)branchParticle->At(particle->M1);  // targeting W boson
            if (abs(particleM->PID) != 24) continue;
            if (particleM->M1 < 0) continue;

            particleMM = (GenParticle*)branchParticle->At(particleM->M1);  // W boson's mother
            iW1 = ip1;
            iWMother = particleM->M1;
        }
    }

    Int_t iW2;
    for (Int_t ip2 = 0; ip2 < nParticles; ip2++) {
        particle = (GenParticle*)branchParticle->At(ip2);
        if (particle->M1 < 0) continue;
        if (abs(particle->PID) >= 1 && abs(particle->PID) <= 6) {        // quark
            particleM = (GenParticle*)branchParticle->At(particle->M1);  // targeting W boson
            if (abs(particleM->PID) != 24) continue;
            if (particle->M1 == iW1) continue;  // exclude the tagged W boson
            if (particleM->M1 < 0) continue;

            particleMM = (GenParticle*)branchParticle->At(particleM->M1);
            if (particleM->M1 != iWMother) continue;
            iW2 = ip2;

            found2W = 1;
        }
    }

    if (found2W == 0) return 1;
    particle = (GenParticle*)branchParticle->At(iW1);
    // cout << endl;
    // cout << particle->PID << endl;
    particle = (GenParticle*)branchParticle->At(iW2);
    // cout << particle->PID << endl;
    // cout << endl;

    return 1;
}
