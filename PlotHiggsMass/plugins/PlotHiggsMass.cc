// -*- C++ -*-
//
// Package:    PlotHiggsMass
// Class:      PlotHiggsMass
//
/**\class PlotHiggsMass PlotHiggsMass.cc GenLevelTools/PlotHiggsMass/plugins/PlotHiggsMass.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  David Sperka
//         Created:  Fri, 30 Jan 2015 09:47:24 GMT
// $Id$
//
//


// system include files
#include <memory>
#include <vector>
#include <iostream>
using namespace std;

// user include files
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

#include "computeAngles.h"

//
// class declaration
//

class PlotHiggsMass : public edm::EDAnalyzer {
   public:
      explicit PlotHiggsMass(const edm::ParameterSet&);
      ~PlotHiggsMass();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      static bool sortByPt( const reco::GenParticle &p1, const reco::GenParticle &p2 ){ return (p1.pt() > p2.pt()); };
      static bool sortByM( const reco::GenParticle &p1, const reco::GenParticle &p2 ){ return (p1.mass() > p2.mass()); };
      static bool sortJetsByPt( const reco::GenJet &p1, const reco::GenJet &p2 ){ return (p1.pt() > p2.pt()); };
      static bool sortPointerByPz( const reco::Candidate *&p1, const reco::Candidate *&p2 ){ return (p1->pz() > p2->pz()); };

   private:

      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      void bookPassedEventTree(TString treeName, TTree *tree);

      edm::ConsumesCollector iC;
      edm::EDGetTokenT<reco::GenParticleCollection> particleCollectionToken_;
      edm::EDGetTokenT<reco::GenJetCollection> jetCollectionToken_;

      TH1F* h_nleptons;

      TTree *GenEventsTree;

      ULong64_t Run, Event, LumiSect;

      int nleptons, njets;
      double pT_H, mH, pT_jet0;

      TClonesArray *p4_lep; TClonesArray *p4_lepS1; TClonesArray *p4_jet;
      TClonesArray *p4_4l; TClonesArray *p4_4lS1; TClonesArray *p4_Z; TClonesArray *p4_ZZ;

      std::vector<int> id_lep; std::vector<int> status_lep; std::vector<int> motherid_lep;

      std::vector<double> *AssociatedParticlePt;
      std::vector<double> *AssociatedParticleEta;
      std::vector<double> *AssociatedParticlePhi;
      std::vector<double> *AssociatedParticleMass;
      std::vector<int> *AssociatedParticleId;

      std::vector<double> *GenAssociatedParticlePt;
      std::vector<double> *GenAssociatedParticleEta;
      std::vector<double> *GenAssociatedParticlePhi;
      std::vector<double> *GenAssociatedParticleMass;
      std::vector<int> *GenAssociatedParticleId;

      std::vector<double> *GenMotherPz;
      std::vector<int> *GenMotherId;

      float cosTheta1, cosTheta2, cosThetaStar, Phi, Phi1;
      double massZ1, massZ2, mass4l;
      int finalState;
      int passedFiducial;

      float LepPt[4], LepEta[4], LepPhi[4], LepMass[4];
      int LepId[4];
      float GenLepPt[4], GenLepEta[4], GenLepPhi[4], GenLepMass[4];
      int GenLepId[4];

      double MaxDijetM, Dijet01M;

      const constexpr static double electronetacut = 2.5;
      const constexpr static double electronpTcut = 7.;
      const constexpr static double muonetacut = 2.4;
      const constexpr static double muonpTcut = 5.;

      int isSelected;

      enum EnumVBForVH {VBF, VH} VBForVH;
};

PlotHiggsMass::PlotHiggsMass(const edm::ParameterSet& iConfig) :
  iC(consumesCollector())
{
    if (iConfig.getParameter<std::string>("VBForVH") == "VBF")
      VBForVH = VBF;
    else if (iConfig.getParameter<std::string>("VBForVH") == "VH")
      VBForVH = VBF;
    else
      throw cms::Exception("BadInput") << "VBForVH should be VBF or VH";

    edm::Service<TFileService> fs;
    h_nleptons    = fs->make<TH1F>( "h_nleptons", "", 1,  0., 1. );

    GenEventsTree = new TTree("GenEvents","GenEvents");

    //now do what ever initialization is needed

    particleCollectionToken_ = iC.consumes<reco::GenParticleCollection> (edm::InputTag("genParticles"));
    jetCollectionToken_ = iC.consumes<reco::GenJetCollection> (edm::InputTag("ak5GenJets"));

    AssociatedParticlePt = new std::vector<double>;
    AssociatedParticleEta = new std::vector<double>;
    AssociatedParticlePhi = new std::vector<double>;
    AssociatedParticleMass = new std::vector<double>;
    AssociatedParticleId = new std::vector<int>;

    GenAssociatedParticlePt = new std::vector<double>;
    GenAssociatedParticleEta = new std::vector<double>;
    GenAssociatedParticlePhi = new std::vector<double>;
    GenAssociatedParticleMass = new std::vector<double>;
    GenAssociatedParticleId = new std::vector<int>;

    GenMotherPz = new std::vector<double>;
    GenMotherId = new std::vector<int>;

}


PlotHiggsMass::~PlotHiggsMass()
{
    delete AssociatedParticlePt;
    delete AssociatedParticleEta;
    delete AssociatedParticlePhi;
    delete AssociatedParticleMass;
    delete AssociatedParticleId;

    delete GenAssociatedParticlePt;
    delete GenAssociatedParticleEta;
    delete GenAssociatedParticlePhi;
    delete GenAssociatedParticleMass;
    delete GenAssociatedParticleId;

    delete GenMotherPz;
    delete GenMotherId;
}


//
// member functions
//

// ------------ method called for each event  ------------
void
PlotHiggsMass::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    Run = iEvent.id().run();
    Event = iEvent.id().event();
    LumiSect = iEvent.id().luminosityBlock();

    // GEN collection
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByToken(particleCollectionToken_, genParticles);

    edm::Handle<reco::GenJetCollection> genJets;
    iEvent.getByToken(jetCollectionToken_, genJets);

    std::vector<reco::GenParticle> leptons;
    std::vector<reco::GenParticle> leptonsS1;
    std::vector<reco::GenParticle> genleptons;
    std::vector<reco::GenParticle> Zs;
    std::vector<reco::GenJet> jets;

    nleptons=0;
    njets=0;
    pT_H = -1.0; mH=1.0;

    p4_lep->Clear(); p4_lepS1->Clear(); p4_jet->Clear(); p4_4l->Clear(); p4_4lS1->Clear(); p4_Z->Clear(); p4_ZZ->Clear();
    id_lep.clear(); status_lep.clear(); motherid_lep.clear();
    AssociatedParticlePt->clear(); AssociatedParticleEta->clear(); AssociatedParticlePhi->clear(); AssociatedParticleMass->clear(); AssociatedParticleId->clear();
    GenAssociatedParticlePt->clear(); GenAssociatedParticleEta->clear(); GenAssociatedParticlePhi->clear(); GenAssociatedParticleMass->clear(); GenAssociatedParticleId->clear();
    GenMotherPz->clear(); GenMotherId->clear();

    LepPt[0] = LepEta[0] = LepPhi[0] = LepMass[0] = LepId[0] = 0;
    LepPt[1] = LepEta[1] = LepPhi[1] = LepMass[1] = LepId[1] = 0;
    LepPt[2] = LepEta[2] = LepPhi[2] = LepMass[2] = LepId[2] = 0;
    LepPt[3] = LepEta[3] = LepPhi[3] = LepMass[3] = LepId[3] = 0;

    massZ1=-1.0; massZ2=-1.0; mass4l=-1.0;
    cosTheta1=9999.0; cosTheta2=9999.0; cosThetaStar=9999.0; Phi=9999.0; Phi1=9999.0;
    finalState=-1;
    passedFiducial=0;

    MaxDijetM=-1.0;
    Dijet01M=-1.0;

    reco::GenParticle higgs;
    vector<const reco::Candidate*> higgsmothers;

    for(auto genPart : *genParticles) {

        if (genPart.pdgId()==25) {
            higgs = genPart;
            pT_H = genPart.pt();
            mH = genPart.mass();
            const reco::Candidate *higgsptr = &genPart;
            if (higgsptr->mother(0)->pdgId()!=25) {
              for (unsigned int i = 0; i < genPart.numberOfMothers(); i++) {
                edm::LogInfo("HiggsMother") << i << " " << higgsptr->mother(i)->pt() << " " << higgsptr->mother(i)->pdgId();
                higgsmothers.push_back(higgsptr->mother(i));
                GenMotherId->push_back(higgsptr->mother(i)->pdgId());
                GenMotherPz->push_back(higgsptr->mother(i)->pz());
              }
              std::sort(higgsmothers.begin(), higgsmothers.end(), sortPointerByPz);
              break;
            }
        }
    }
    for(auto genPart : *genParticles) {

        if (genPart.pdgId()==23) {
            Zs.push_back(genPart);
        }

        if (abs(genPart.pdgId()) < 6 || (abs(genPart.pdgId()) >= 11 && abs(genPart.pdgId()) <= 16)) {
            if( genPart.mother()->numberOfMothers() >= 1 && genPart.mother()->mother()->pdgId()==25 ) {
                edm::LogInfo("Leptons") << "gen  " << genPart.pdgId() << " " << genPart.pt() << " " << genPart.eta();
                genleptons.push_back(genPart);
            } else if (genPart.numberOfMothers() == higgsmothers.size() && VBForVH == VBF){
                vector<const reco::Candidate*> mothers;
                for (unsigned int i = 0; i < genPart.numberOfMothers(); i++) {
                  edm::LogInfo("JetMother") << i << " " << genPart.mother(i)->pt() << " " << genPart.mother(i)->pdgId();
                  mothers.push_back(genPart.mother(i));
                }
                std::sort(mothers.begin(), mothers.end(), sortPointerByPz);

                bool matches = true;
                for (unsigned int i = 0; i < mothers.size(); i++) {
                    if (! (mothers[i]->p4() == higgsmothers[i]->p4() && mothers[i]->pdgId() == higgsmothers[i]->pdgId()))
                        matches = false;
                }
                if (matches) {
                    GenAssociatedParticlePt->push_back(genPart.pt());
                    GenAssociatedParticleEta->push_back(genPart.eta());
                    GenAssociatedParticlePhi->push_back(genPart.phi());
                    GenAssociatedParticleMass->push_back(genPart.mass());
                    GenAssociatedParticleId->push_back(genPart.pdgId());
                }
            } else if (genPart.numberOfMothers() == 1 && (genPart.mother(0)->pdgId()==23 || abs(genPart.mother(0)->pdgId())==24)) {
                GenAssociatedParticlePt->push_back(genPart.pt());
                GenAssociatedParticleEta->push_back(genPart.eta());
                GenAssociatedParticlePhi->push_back(genPart.phi());
                GenAssociatedParticleMass->push_back(genPart.mass());
                GenAssociatedParticleId->push_back(genPart.pdgId());
            }
        }

        if( abs(genPart.pdgId())==15 || abs(genPart.pdgId()) == 13 || abs(genPart.pdgId()) == 11 ) {

            bool passcut = true;
            if (abs(genPart.pdgId())==11 && (genPart.pt()<electronpTcut || abs(genPart.eta()) > electronetacut)) passcut = false;
            if (abs(genPart.pdgId())==13 && (genPart.pt()<    muonpTcut || abs(genPart.eta()) >     muonetacut)) passcut = false;
            if (abs(genPart.pdgId())==15) passcut = false;

            if(genPart.status()==1) {
                edm::LogInfo("Leptons") << "reco " << passcut << " " << genPart.pdgId() << " " << genPart.pt() << " " << genPart.eta();
                if (passcut) {
                    nleptons++;
                    leptons.push_back(genPart);
                    leptonsS1.push_back(genPart);
                }
            }
        }
    }

    if (genleptons.size() != 4)
      throw cms::Exception("GenLeptons") << genleptons.size() << " gen leptons in the event";

    for (unsigned int i = 0; i < 4; i++) {
        GenLepPt[i] = genleptons[i].pt();
        GenLepEta[i] = genleptons[i].eta();
        GenLepPhi[i] = genleptons[i].phi();
        GenLepMass[i] = genleptons[i].mass();
        GenLepId[i] = genleptons[i].pdgId();
    }

    int ep=0, em=0, mup=0, mum=0, taup=0, taum=0;
    for (auto l : leptons) {
      if (l.pdgId() == +11) ep++;
      if (l.pdgId() == -11) em++;
      if (l.pdgId() == +13) mup++;
      if (l.pdgId() == -13) mum++;
      if (l.pdgId() == +15) taup++;
      if (l.pdgId() == -15) taum++;
    }
    isSelected = (min(ep, em) + min(mup, mum) + min(taup, taum) >= 2);

    std::sort(Zs.begin(), Zs.end(), sortByM);
    for (unsigned int i=0; i<Zs.size(); ++i) {
        new ( (*p4_Z)[i] ) TLorentzVector(Zs[i].px(),Zs[i].py(),Zs[i].pz(),Zs[i].energy());
    }
    if (Zs.size()>=2) {
        TLorentzVector *Z_1 = (TLorentzVector*) p4_Z->At(0);
        TLorentzVector *Z_2 = (TLorentzVector*) p4_Z->At(1);
        TLorentzVector p4_Z1Z2 = (*Z_1)+(*Z_2);
        new ( (*p4_ZZ)[0] ) TLorentzVector(p4_Z1Z2.Px(), p4_Z1Z2.Py(), p4_Z1Z2.Pz(), p4_Z1Z2.E());
    }

    std::sort(leptons.begin(), leptons.end(), sortByPt);
    for (unsigned int i=0; i<leptons.size(); ++i) {
        new ( (*p4_lep)[i] ) TLorentzVector(leptons[i].px(),leptons[i].py(),leptons[i].pz(),leptons[i].energy());
        id_lep.push_back(leptons[i].pdgId());
        status_lep.push_back(leptons[i].status());
        motherid_lep.push_back(leptons[i].mother()->pdgId());
    }

    unsigned int Nlep = p4_lep->GetLast()+1;
    if (isSelected) {

        TLorentzVector *L_1 = (TLorentzVector*) p4_lep->At(0);
        TLorentzVector *L_2 = (TLorentzVector*) p4_lep->At(1);
        TLorentzVector *L_3 = (TLorentzVector*) p4_lep->At(2);
        TLorentzVector *L_4 = (TLorentzVector*) p4_lep->At(3);

        TLorentzVector p4_4lep = (*L_1)+(*L_2)+(*L_3)+(*L_4);
        new ( (*p4_4l)[0] ) TLorentzVector(p4_4lep.Px(), p4_4lep.Py(), p4_4lep.Pz(), p4_4lep.E());
        mass4l=p4_4lep.M();

        double offshell=99999.0;
        int L1=10, L2=10, L3=10, L4=10;

        for(unsigned int i=0; i<Nlep; i++){
            for(unsigned int j=i+1; j<Nlep; j++){

                if((id_lep[i]+id_lep[j])!=0) continue;

                TLorentzVector *li, *lj;
                li = (TLorentzVector*) p4_lep->At(i);
                lj = (TLorentzVector*) p4_lep->At(j);

                TLorentzVector mll = (*li)+(*lj);
                if(abs(mll.M()-91.1876)<offshell){
                    double mZ1 = mll.M();
                    L1 = i; L2 = j; offshell = abs(mZ1-91.1876);
                }
            }
        }

        TLorentzVector *Z11, *Z12;
        Z11 = (TLorentzVector*) p4_lep->At(L1);
        Z12 = (TLorentzVector*) p4_lep->At(L2);
        TLorentzVector Z1 = (*Z11)+(*Z12);

        bool passZ1=false;
        if(Z1.M()>40 && Z1.M()<120) passZ1 = true;
        massZ1=Z1.M();

        double pTL34 = 0.0; bool passZ2 = false;

        for(unsigned int i=0; i<Nlep; i++){
            if((int)i==L1 || (int)i==L2) continue; // can not be the lep from Z1
            for(unsigned int j=i+1; j<Nlep; j++){
                if((int)j==L1 || (int)j==L2) continue; // can not be the lep from Z1
                if((id_lep[i]+id_lep[j])!=0) continue;

                TLorentzVector *li, *lj;
                li = (TLorentzVector*) p4_lep->At(i);
                lj = (TLorentzVector*) p4_lep->At(j);
                TLorentzVector Z2 = (*li)+(*lj);

                if ( ( (*li).Pt()+(*lj).Pt() ) >=pTL34 ) { // choose high sum pT pair satisfy the following selection
                    double mZ2 = Z2.M();
                    if(mZ2>12.0 && mZ2<120.0){
                        L3 = i; L4 = j; passZ2 = true; pTL34 = (*li).Pt()+(*lj).Pt();
                    } else {
                        // still assign L3 and L4 to this pair if we don't have a passing Z2 yet
                        if (passZ2 == false) {L3 = i; L4 = j;}
                    }
                }

            }
        }

        LepPt[0] = ((TLorentzVector*)p4_lep->At(L1))->Pt();
        LepEta[0] = ((TLorentzVector*)p4_lep->At(L1))->Eta();
        LepPhi[0] = ((TLorentzVector*)p4_lep->At(L1))->Phi();
        LepMass[0] = ((TLorentzVector*)p4_lep->At(L1))->M();

        LepPt[1] = ((TLorentzVector*)p4_lep->At(L2))->Pt();
        LepEta[1] = ((TLorentzVector*)p4_lep->At(L2))->Eta();
        LepPhi[1] = ((TLorentzVector*)p4_lep->At(L2))->Phi();
        LepMass[1] = ((TLorentzVector*)p4_lep->At(L2))->M();

        LepPt[2] = ((TLorentzVector*)p4_lep->At(L3))->Pt();
        LepEta[2] = ((TLorentzVector*)p4_lep->At(L3))->Eta();
        LepPhi[2] = ((TLorentzVector*)p4_lep->At(L3))->Phi();
        LepMass[2] = ((TLorentzVector*)p4_lep->At(L3))->M();

        LepPt[3] = ((TLorentzVector*)p4_lep->At(L4))->Pt();
        LepEta[3] = ((TLorentzVector*)p4_lep->At(L4))->Eta();
        LepPhi[3] = ((TLorentzVector*)p4_lep->At(L4))->Phi();
        LepMass[3] = ((TLorentzVector*)p4_lep->At(L4))->M();

        LepId[0] = leptons[L1].pdgId();
        LepId[1] = leptons[L2].pdgId();
        LepId[2] = leptons[L3].pdgId();
        LepId[3] = leptons[L4].pdgId();

        if (passZ1 && passZ2) {

            passedFiducial=1;

            TLorentzVector *Z21, *Z22;
            Z21 = (TLorentzVector*) p4_lep->At(L3);
            Z22 = (TLorentzVector*) p4_lep->At(L4);
            TLorentzVector Z2 = (*Z21)+(*Z22);
            massZ2=Z2.M();

            // prepare TLorentzVector(s)
            TLorentzVector L11P4,L12P4,L21P4,L22P4;
            // convention for leptons - example for 2mu2e event:    pp -> H -> ZZ -> mu-(p3)+mu+(p4)+e-(p5)+e+(p6)
            L11P4.SetPxPyPzE(Z11->Px(), Z11->Py(), Z11->Pz(), Z11->E());
            L12P4.SetPxPyPzE(Z12->Px(), Z12->Py(), Z12->Pz(), Z12->E());
            L21P4.SetPxPyPzE(Z21->Px(), Z21->Py(), Z21->Pz(), Z21->E());
            L22P4.SetPxPyPzE(Z22->Px(), Z22->Py(), Z22->Pz(), Z22->E());

            // prepare vectors of TLorentzVector(s)
            std::vector<TLorentzVector> partP;
            std::vector<int> partId;
            partP.push_back(L11P4); partP.push_back(L12P4); partP.push_back(L21P4); partP.push_back(L22P4);
            partId.push_back(id_lep[L1]); partId.push_back(id_lep[L2]); partId.push_back(id_lep[L3]); partId.push_back(id_lep[L4]);

        }
    }

    std::sort(leptonsS1.begin(), leptonsS1.end(), sortByPt);
    for (unsigned int i=0; i<leptonsS1.size(); ++i) {
        new ( (*p4_lepS1)[i] ) TLorentzVector(leptonsS1[i].px(),leptonsS1[i].py(),leptonsS1[i].pz(),leptonsS1[i].energy());
    }

    unsigned int NlepS1 = p4_lepS1->GetLast()+1;
    if (NlepS1>=4) {

        TLorentzVector *LS1_1 = (TLorentzVector*) p4_lepS1->At(0);
        TLorentzVector *LS1_2 = (TLorentzVector*) p4_lepS1->At(1);
        TLorentzVector *LS1_3 = (TLorentzVector*) p4_lepS1->At(2);
        TLorentzVector *LS1_4 = (TLorentzVector*) p4_lepS1->At(3);

        TLorentzVector p4_4lepS1 = (*LS1_1)+(*LS1_2)+(*LS1_3)+(*LS1_4);
        new ( (*p4_4lS1)[0] ) TLorentzVector(p4_4lepS1.Px(), p4_4lepS1.Py(), p4_4lepS1.Pz(), p4_4lepS1.E());
    }

    reco::GenJetCollection::const_iterator  genJet;
    for(genJet = genJets->begin(); genJet != genJets->end(); genJet++) {
        double min_dR = 9999.0;
        for (unsigned int i=0; i<leptonsS1.size(); ++i) {
            double dR = deltaR(leptonsS1[i].eta(), leptonsS1[i].phi(), genJet->eta(), genJet->phi());
            if (dR<min_dR) min_dR = dR;
        }
        if (min_dR>0.5 && genJet->pt()>30.0) {
            njets++;
            jets.push_back(*genJet);
        }
    }
    std::sort(jets.begin(), jets.end(), sortJetsByPt);
    for (unsigned int i=0; i<jets.size(); ++i) {
        new ( (*p4_jet)[i] ) TLorentzVector(jets[i].px(),jets[i].py(),jets[i].pz(),jets[i].energy());
    }

    for (unsigned int i=0; i<jets.size(); ++i) {
        AssociatedParticlePt->push_back(jets[i].pt());
        AssociatedParticleEta->push_back(jets[i].eta());
        AssociatedParticlePhi->push_back(jets[i].phi());
        AssociatedParticleMass->push_back(jets[i].mass());
        AssociatedParticleId->push_back(jets[i].pdgId());
        for (unsigned int j=0; j<jets.size(); ++j) {
            if (i==j) continue;
            const TLorentzVector jet1(jets[i].px(), jets[i].py(), jets[i].pz(), jets[i].energy());
            const TLorentzVector jet2(jets[j].px(), jets[j].py(), jets[j].pz(), jets[j].energy());
            TLorentzVector dijet = jet1+jet2;
            if (dijet.M()>MaxDijetM) MaxDijetM=dijet.M();
            if (i==0 and j==1) Dijet01M=dijet.M();
        }
    }

    h_nleptons->Fill( nleptons );

    GenEventsTree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void
PlotHiggsMass::beginJob()
{

    bookPassedEventTree("GenEvents", GenEventsTree);

}

void PlotHiggsMass::bookPassedEventTree(TString treeName, TTree *tree)
{


    tree->Branch("Run",&Run,"Run/l");
    tree->Branch("Event",&Event,"Event/l");
    tree->Branch("LumiSect",&LumiSect,"LumiSect/l");

    tree->Branch("isSelected", &isSelected, "isSelected/I");

    tree->Branch("nleptons",&nleptons,"nleptons/I");
    tree->Branch("njets",&njets,"njets/I");
    tree->Branch("pT_H",&pT_H,"pT_H/D");
    tree->Branch("mH",&mH,"mH/D");

    p4_lep = new TClonesArray("TLorentzVector", 10);
    tree->Branch("p4_lep","TClonesArray", &p4_lep,128000, 0);
    tree->Branch("id_lep",&id_lep);
    tree->Branch("status_lep",&status_lep);
    tree->Branch("motherid_lep",&motherid_lep);

    p4_4l = new TClonesArray("TLorentzVector", 10);
    tree->Branch("p4_4l","TClonesArray", &p4_4l,128000, 0);

    p4_lepS1 = new TClonesArray("TLorentzVector", 10);
    tree->Branch("p4_lepS1","TClonesArray", &p4_lepS1,128000, 0);

    p4_4lS1 = new TClonesArray("TLorentzVector", 10);
    tree->Branch("p4_4lS1","TClonesArray", &p4_4lS1,128000, 0);

    p4_Z = new TClonesArray("TLorentzVector", 10);
    tree->Branch("p4_Z","TClonesArray", &p4_Z,128000, 0);

    p4_ZZ = new TClonesArray("TLorentzVector", 10);
    tree->Branch("p4_ZZ","TClonesArray", &p4_ZZ,128000, 0);

    p4_jet = new TClonesArray("TLorentzVector", 10);
    tree->Branch("p4_jet","TClonesArray", &p4_jet,128000, 0);

    tree->Branch("AssociatedParticlePt", &AssociatedParticlePt);
    tree->Branch("AssociatedParticleEta", &AssociatedParticleEta);
    tree->Branch("AssociatedParticlePhi", &AssociatedParticlePhi);
    tree->Branch("AssociatedParticleMass", &AssociatedParticleMass);
    tree->Branch("AssociatedParticleId", &AssociatedParticleId);

    tree->Branch("GenAssociatedParticlePt", &GenAssociatedParticlePt);
    tree->Branch("GenAssociatedParticleEta", &GenAssociatedParticleEta);
    tree->Branch("GenAssociatedParticlePhi", &GenAssociatedParticlePhi);
    tree->Branch("GenAssociatedParticleMass", &GenAssociatedParticleMass);
    tree->Branch("GenAssociatedParticleId", &GenAssociatedParticleId);

    tree->Branch("GenMotherPz", &GenMotherPz);
    tree->Branch("GenMotherId", &GenMotherId);

    tree->Branch("cosTheta1",&cosTheta1,"cosTheta1/F");
    tree->Branch("cosTheta2",&cosTheta2,"cosTheta2/F");
    tree->Branch("cosThetaStar",&cosThetaStar,"cosThetaStar/F");
    tree->Branch("Phi",&Phi,"Phi/F");
    tree->Branch("Phi1",&Phi1,"Phi1/F");
    tree->Branch("massZ1",&massZ1,"massZ1/D");
    tree->Branch("massZ2",&massZ2,"massZ2/D");
    tree->Branch("mass4l",&mass4l,"mass4l/D");
    tree->Branch("finalState",&finalState,"finalState/I");
    tree->Branch("passedFiducial",&passedFiducial,"passedFiducial/I");

    tree->Branch("MaxDijetM",&MaxDijetM,"MaxDijetM/D");
    tree->Branch("Dijet01M",&Dijet01M,"Dijet01M/D");

    tree->Branch("Lep1Pt", &LepPt[0], "Lep1Pt/F");
    tree->Branch("Lep1Eta", &LepEta[0], "Lep1Eta/F");
    tree->Branch("Lep1Phi", &LepPhi[0], "Lep1Phi/F");
    tree->Branch("Lep1Mass", &LepMass[0], "Lep1Mass/F");
    tree->Branch("Lep1Id", &LepId[0], "Lep1Id/I");

    tree->Branch("Lep2Pt", &LepPt[1], "Lep2Pt/F");
    tree->Branch("Lep2Eta", &LepEta[1], "Lep2Eta/F");
    tree->Branch("Lep2Phi", &LepPhi[1], "Lep2Phi/F");
    tree->Branch("Lep2Mass", &LepMass[1], "Lep2Mass/F");
    tree->Branch("Lep2Id", &LepId[1], "Lep2Id/I");

    tree->Branch("Lep3Pt", &LepPt[2], "Lep3Pt/F");
    tree->Branch("Lep3Eta", &LepEta[2], "Lep3Eta/F");
    tree->Branch("Lep3Phi", &LepPhi[2], "Lep3Phi/F");
    tree->Branch("Lep3Mass", &LepMass[2], "Lep3Mass/F");
    tree->Branch("Lep3Id", &LepId[2], "Lep3Id/I");

    tree->Branch("Lep4Pt", &LepPt[3], "Lep4Pt/F");
    tree->Branch("Lep4Eta", &LepEta[3], "Lep4Eta/F");
    tree->Branch("Lep4Phi", &LepPhi[3], "Lep4Phi/F");
    tree->Branch("Lep4Mass", &LepMass[3], "Lep4Mass/F");
    tree->Branch("Lep4Id", &LepId[3], "Lep4Id/I");

    tree->Branch("GenLep1Pt", &GenLepPt[0], "GenLep1Pt/F");
    tree->Branch("GenLep1Eta", &GenLepEta[0], "GenLep1Eta/F");
    tree->Branch("GenLep1Phi", &GenLepPhi[0], "GenLep1Phi/F");
    tree->Branch("GenLep1Mass", &GenLepMass[0], "GenLep1Mass/F");
    tree->Branch("GenLep1Id", &GenLepId[0], "GenLep1Id/I");

    tree->Branch("GenLep2Pt", &GenLepPt[1], "GenLep2Pt/F");
    tree->Branch("GenLep2Eta", &GenLepEta[1], "GenLep2Eta/F");
    tree->Branch("GenLep2Phi", &GenLepPhi[1], "GenLep2Phi/F");
    tree->Branch("GenLep2Mass", &GenLepMass[1], "GenLep2Mass/F");
    tree->Branch("GenLep2Id", &GenLepId[1], "GenLep2Id/I");

    tree->Branch("GenLep3Pt", &GenLepPt[2], "GenLep3Pt/F");
    tree->Branch("GenLep3Eta", &GenLepEta[2], "GenLep3Eta/F");
    tree->Branch("GenLep3Phi", &GenLepPhi[2], "GenLep3Phi/F");
    tree->Branch("GenLep3Mass", &GenLepMass[2], "GenLep3Mass/F");
    tree->Branch("GenLep3Id", &GenLepId[2], "GenLep3Id/I");

    tree->Branch("GenLep4Pt", &GenLepPt[3], "GenLep4Pt/F");
    tree->Branch("GenLep4Eta", &GenLepEta[3], "GenLep4Eta/F");
    tree->Branch("GenLep4Phi", &GenLepPhi[3], "GenLep4Phi/F");
    tree->Branch("GenLep4Mass", &GenLepMass[3], "GenLep4Mass/F");
    tree->Branch("GenLep4Id", &GenLepId[3], "GenLep4Id/I");

}

// ------------ method called once each job just after ending the event loop  ------------
void
PlotHiggsMass::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
/*
void
PlotHiggsMass::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
PlotHiggsMass::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
PlotHiggsMass::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
PlotHiggsMass::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PlotHiggsMass::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PlotHiggsMass);
