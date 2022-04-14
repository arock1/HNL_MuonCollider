#include <iostream>
using namespace std;

#include "FinalStatesClass.C"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "classes/DelphesClasses.h"

//{{{
Int_t ClassifySingal(TClonesArray* branchParticle, iFinalStates* iFSTrue, TLorentzVector* NTrue) {
    GenParticle* particle;
    GenParticle* particle0;
    GenParticle* particle0M;
    GenParticle* particle0MM;

    GenParticle* particleNM1;
    GenParticle* particleNM2;

    Int_t iNTrue = 99999;  // index for HNL
    Int_t foundN = 0;
    Int_t iWTrue = 99999;  // index for W from N
    Int_t foundW = 0;
    Int_t iLepTrue = 99999;  // index for lepton from N
    Int_t foundLep = 0;
    Int_t iQ1True = 99999;  // down type quark
    Int_t foundQ1 = 0;
    Int_t iQ2True = 99999;  // up type quark
    Int_t foundQ2 = 0;

    Int_t typeLep;
    Int_t lepCharge;

    Int_t nParticles = branchParticle->GetEntries();
    // cout << "\n\n\n Num of p: " << nParticles << endl;

    for (Int_t iN = 0; iN < nParticles; iN++) {
        foundN = 0;
        particle = (GenParticle*)branchParticle->At(iN);
        // cout << particle->PID << endl;

        if (abs(particle->PID) == 9900012) {                               // find HNL
            particleNM1 = (GenParticle*)branchParticle->At(particle->M1);  // mother of HNL
            particleNM2 = (GenParticle*)branchParticle->At(particle->M2);  // mother of HNL
            // cout << "\nN1 Mother1: " << particleNM1->PID << "\nN1 Mother2: " << particleNM2->PID << "\n";
            // cout << particleNM1->Charge << "; " << particleNM1->Eta << endl;
            // cout << endl;

            // cout << "found N1" << endl;
            foundN = 1;
            iNTrue = iN;

            for (Int_t i = 0; i < nParticles; i++) {  // find lepton (electron or muon only)
                particle0 = (GenParticle*)branchParticle->At(i);
                if ((abs(particle0->PID) == 11 || abs(particle0->PID) == 13) && particle0->M1 == iNTrue) {
                    // cout << particle0->PID << endl;
                    //  cout << "Found Lepton: " << particle0->PID << endl;
                    foundLep = 1;
                    iLepTrue = i;
                    typeLep = abs(particle0->PID);
                    lepCharge = particle0->Charge;
                    break;
                }
            }

            for (Int_t ii = 0; ii < nParticles; ii++) {            // find qq<W<N
                particle0 = (GenParticle*)branchParticle->At(ii);  // quark
                if (particle0->M1 == -1) continue;                 // need to find quark's mother (W boson), if no mother, then skip
                Int_t iM = particle0->M1;
                particle0M = (GenParticle*)branchParticle->At(iM);

                if (abs(particle0->PID) <= 8) {  // if it is quark
                    while (true) {               // loop to trace its history
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

                GenParticle* lepParticle = (GenParticle*)branchParticle->At(iLepTrue);
                GenParticle* Q1Particle = (GenParticle*)branchParticle->At(iQ1True);
                GenParticle* Q2Particle = (GenParticle*)branchParticle->At(iQ2True);
                GenParticle* NParticle = (GenParticle*)branchParticle->At(iNTrue);

                TLorentzVector lepTrue_, Q1True_, Q2True_, NTrue_;
                lepTrue_.SetPtEtaPhiE(lepParticle->PT, lepParticle->Eta, lepParticle->Phi, lepParticle->E);
                Q1True_.SetPtEtaPhiE(Q1Particle->PT, Q1Particle->Eta, Q1Particle->Phi, Q1Particle->E);
                Q2True_.SetPtEtaPhiE(Q2Particle->PT, Q2Particle->Eta, Q2Particle->Phi, Q2Particle->E);
                NTrue_.SetPtEtaPhiE(NParticle->PT, NParticle->Eta, NParticle->Phi, NParticle->E);

                iFSTrue_.iLeps.push_back(lepTrue_);
                iFSTrue_.iJets.push_back(Q1True_);
                iFSTrue_.iJets.push_back(Q2True_);
                iFSTrue_.iLepCharges.push_back(lepCharge);
                *iFSTrue = iFSTrue_;

                *NTrue = NTrue_;

                return 1;
            }
        }
    }
    // cout << foundN1 << "; " << foundQ1 << "; " << foundQ2 << endl;
    return 0;
}
//}}}

struct BkgTypes {
    Int_t WW = 0;
    Int_t phph = 0;
    Int_t Zph = 0;
    Int_t ZZ = 0;
    Int_t WZ = 0;
    Int_t Wph = 0;
    Int_t lepFromPhys = 0;  // lepton from physics decay
    Int_t eFromPhys = 0;    // e from physics decay
    Int_t muFromPhys = 0;   // mu from physics decay
    Int_t others = 0;
    Int_t others_haveEle = 0;
    Int_t motherMu = 0;
};

Int_t ClassifyiBackground(TClonesArray* branchParticle, BkgTypes* bkgTypes) {
    GenParticle* particle1;
    GenParticle* particle2;
    GenParticle* particle1M;
    GenParticle* particle2M;
    GenParticle* particle2MM;

    Int_t nParticles = branchParticle->GetEntries();

    Int_t found1V = 0;
    Int_t found2V = 0;

    Int_t iV1, iV2;
    Int_t PIDV1 = 99999;
    Int_t PIDV2 = 99999;
    Int_t lep = 0;
    Int_t haveEle = 0;

    Int_t maxPt = -99999;
    for (Int_t ip1 = 0; ip1 < nParticles; ip1++) {
        found1V = 0;
        particle1 = (GenParticle*)branchParticle->At(ip1);

        if (abs(particle1->PID) == 13 || abs(particle1->PID) == 11) {  // lepton
            if (particle1->M1 < 0) continue;
            // only focus on where the max pT lepton from
            if (particle1->PT < maxPt) continue;
            maxPt = particle1->PT;

            if (abs(particle1->PID) == 11) haveEle = 1;
            particle1M = (GenParticle*)branchParticle->At(particle1->M1);                                  // lepton mother
            if (abs(particle1M->PID) == 22 || abs(particle1M->PID) == 23 || abs(particle1M->PID) == 24) {  // vector bosons
                iV1 = particle1->M1;
                PIDV1 = particle1M->PID;
                lep = particle1->PID;
                // break;
            }
            if ((abs(particle1->PID) == 13 && abs(particle1M->PID) != 13) || (abs(particle1->PID) == 11)) {  // if not vector boson, check if it is from phys decay
                lep = particle1->PID;
                PIDV1 = particle1M->PID;
                // DO NOT break. It might be prossible that in later loop, there is another lepton from V
            }
            // if (abs(particle1->PID) == 13 && abs(particle1M->PID) == 13) {
            //  cout << ">????" << endl;
            // PIDV1 = particle1M->PID;
            //}
        }
    }

    for (Int_t ip2 = 0; ip2 < nParticles; ip2++) {
        found2V = 0;
        particle2 = (GenParticle*)branchParticle->At(ip2);
        if (abs(particle2->PID) >= 1 && abs(particle2->PID) <= 6) {  // quark
            if (particle2->M1 < 0) continue;
            particle2M = (GenParticle*)branchParticle->At(particle2->M1);  // targeting Vector boson
            if (particle2->M1 == iV1) continue;                            // exclude the tagged Vector boson
            if (abs(particle2M->PID) == 22 || abs(particle2M->PID) == 23 || abs(particle2M->PID) == 24) {
                iV2 = particle2->M1;
                PIDV2 = particle2M->PID;
                break;
            }
            if (abs(PIDV2) != 22 && abs(PIDV2) != 23 && abs(PIDV2) != 24) {
                PIDV2 = particle2M->PID;
            }
            // if (abs(particle2M->PID) == 21) {
            // cout << " ..." << endl;
            //}
            found2V = 1;
            PIDV2 = particle2M->PID;
        }
    }

    if (abs(PIDV1) == 24 && abs(PIDV2) == 24) {
        bkgTypes->WW += 1;
    } else if ((abs(PIDV1) == 22 && abs(PIDV2) == 22)) {
        bkgTypes->phph += 1;
    } else if ((abs(PIDV1) == 22 && abs(PIDV2) == 23) || (abs(PIDV1) == 23 && abs(PIDV2) == 22)) {
        bkgTypes->Zph += 1;
    } else if ((abs(PIDV1) == 23 && abs(PIDV2) == 23)) {
        bkgTypes->ZZ += 1;
    } else if ((abs(PIDV1) == 24 && abs(PIDV2) == 23) || (abs(PIDV1) == 23 && abs(PIDV2) == 24)) {
        bkgTypes->WZ += 1;
    } else if ((abs(PIDV1) == 24 && abs(PIDV2) == 22) || (abs(PIDV1) == 22 && abs(PIDV2) == 24)) {
        bkgTypes->Wph += 1;
    } else if (abs(PIDV1) != 99999 && abs(PIDV1) != 22 && abs(PIDV1) != 23 && abs(PIDV1) != 24) {
        //} else if (abs(PIDV1) != 13 && abs(PIDV1) != 99999 && abs(PIDV1) != 22 && abs(PIDV1) != 23 && abs(PIDV1) != 24) {
        if (abs(lep) == 13) bkgTypes->muFromPhys += 1;
        if (abs(lep) == 11) bkgTypes->eFromPhys += 1;
        bkgTypes->lepFromPhys += 1;
        // if (abs(lep) == 11) {
        // cout << PIDV1 << endl;
        //}

        //} else if (abs(PIDV1) == 13) {
        // cout << "...." << endl;
        // bkgTypes->motherMu += 1;
    } else {
        // cout << particle1M->PID << "; " << particle2M->PID << endl;
        if (haveEle == 1) bkgTypes->others_haveEle += 1;
        bkgTypes->others += 1;
    }

    return 1;
}
