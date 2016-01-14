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

// user include files
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

#include "computeAngles.cc"

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

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      void bookPassedEventTree(TString treeName, TTree *tree);

      TH1F* h_nleptons;

      TTree *GenEventsTree;

      ULong64_t Run, Event, LumiSect;    

      int nleptons, njets;
      double pT_H, mH, pT_jet0; 
    
      TClonesArray *p4_lep; TClonesArray *p4_lepS1; TClonesArray *p4_jet; 
      TClonesArray *p4_4l; TClonesArray *p4_4lS1; TClonesArray *p4_Z; TClonesArray *p4_ZZ;

      std::vector<int> id_lep; std::vector<int> status_lep; std::vector<int> motherid_lep;

      float cosTheta1, cosTheta2, cosThetaStar, Phi, Phi1;
      double massZ1, massZ2, mass4l;
      int finalState;
      int passedFiducial;

      double MaxDijetM, Dijet01M;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PlotHiggsMass::PlotHiggsMass(const edm::ParameterSet& iConfig)

{
    edm::Service<TFileService> fs;
    h_nleptons    = fs->make<TH1F>( "h_nleptons", "", 1,  0., 1. );

    GenEventsTree = new TTree("GenEvents","GenEvents");
    
    //now do what ever initialization is needed

}


PlotHiggsMass::~PlotHiggsMass()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

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

    //std::cout<<"Run: "<<Run<<" Event: "<<Event<<" Lumi: "<<LumiSect<<std::endl;
    // GEN collection
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel("genParticles", genParticles);
   
    edm::Handle<reco::GenJetCollection> genJets;
    iEvent.getByLabel("ak5GenJets", genJets);

    std::vector<reco::GenParticle> leptons;
    std::vector<reco::GenParticle> leptonsS1;
    std::vector<reco::GenParticle> Zs;
    std::vector<reco::GenJet> jets;

    nleptons=0;
    njets=0;
    pT_H = -1.0; mH=1.0;
    p4_lep->Clear(); p4_lepS1->Clear(); p4_jet->Clear(); p4_4l->Clear(); p4_4lS1->Clear(); p4_Z->Clear(); p4_ZZ->Clear();
    id_lep.clear(); status_lep.clear(); motherid_lep.clear();

    massZ1=-1.0; massZ2=-1.0; mass4l=-1.0; 
    cosTheta1=9999.0; cosTheta2=9999.0; cosThetaStar=9999.0; Phi=9999.0; Phi1=9999.0;
    finalState=-1;
    passedFiducial=0;

    MaxDijetM=-1.0;
    Dijet01M=-1.0;

    reco::GenParticleCollection::const_iterator  genPart;
    for(genPart = genParticles->begin(); genPart != genParticles->end(); genPart++) {
        
        if (genPart->pdgId()==25) {
            pT_H = genPart->pt();
            mH = genPart->mass();
            //std::cout<<"H status: "<<genPart->status()<<" pt: "<<genPart->pt()<<std::endl;
        }
        if (genPart->pdgId()==23) {
            Zs.push_back(*genPart);
        }

        if( abs(genPart->pdgId())==15 || abs(genPart->pdgId()) == 13 || abs(genPart->pdgId()) == 11 ) {

            if( genPart->mother()->pdgId()==23 or abs(genPart->mother()->pdgId())==24 ) {
                //std::cout<<"leptons from W or Z:"<<std::endl;
                //std::cout<<"lepton id: "<<genPart->pdgId()<<" pt: "<<genPart->pt()<<" eta: "<<genPart->eta()<<" status: "<<genPart->status()<<" mother: "<<genPart->mother()->pdgId()<<std::endl;
                leptons.push_back(*genPart);
                nleptons++;
            }        

            //if( genPart->status()==1 && (genPart->mother()->pdgId()==23 or abs(genPart->mother()->pdgId())==24 or genPart->mother()->pdgId()==genPart->pdgId())) {
            if(genPart->status()==1) {
                //std::cout<<"status 1 leptons:"<<std::endl;
                //std::cout<<"lepton id: "<<genPart->pdgId()<<" pt: "<<genPart->pt()<<" eta: "<<genPart->eta()<<" status: "<<genPart->status()<<" mother: "<<genPart->mother()->pdgId()<<std::endl;
                leptonsS1.push_back(*genPart);
                
            }        

        }
    }

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
        //std::cout<<"lepton id: "<<leptons[i].pdgId()<<" pt: "<<leptons[i].pt()<<" eta: "<<leptons[i].eta()<<" status: "<<leptons[i].status()<<" mother: "<<leptons[i].mother()->pdgId()<<std::endl;
        new ( (*p4_lep)[i] ) TLorentzVector(leptons[i].px(),leptons[i].py(),leptons[i].pz(),leptons[i].energy());
        id_lep.push_back(leptons[i].pdgId());
        status_lep.push_back(leptons[i].status());       
        motherid_lep.push_back(leptons[i].mother()->pdgId());
    }

    unsigned int Nlep = p4_lep->GetLast()+1;
    if (Nlep>=4) {

        TLorentzVector *L_1 = (TLorentzVector*) p4_lep->At(0);
        TLorentzVector *L_2 = (TLorentzVector*) p4_lep->At(1);
        TLorentzVector *L_3 = (TLorentzVector*) p4_lep->At(2);
        TLorentzVector *L_4 = (TLorentzVector*) p4_lep->At(3);

        TLorentzVector p4_4lep = (*L_1)+(*L_2)+(*L_3)+(*L_4);
        new ( (*p4_4l)[0] ) TLorentzVector(p4_4lep.Px(), p4_4lep.Py(), p4_4lep.Pz(), p4_4lep.E());
        //std::cout<< "p4_4l[0].Pt()" << p4_4l[0].Pt() << std::endl;
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

        if (abs(id_lep[L1])==11 && abs(id_lep[L3])==11) finalState=1; // 4e
        else if (abs(id_lep[L1])==13 && abs(id_lep[L3])==13) finalState=2; // 4mu
        else if (abs(id_lep[L1])==15 && abs(id_lep[L3])==15) finalState=3; // 4tau 
        else if ( (abs(id_lep[L1])+abs(id_lep[L3]))==24 ) finalState=4; // 2e2mu
        else if ( (abs(id_lep[L1])+abs(id_lep[L3]))==26 ) finalState=5; // 2e2tau
        else if ( (abs(id_lep[L1])+abs(id_lep[L3]))==28 ) finalState=6; // 2mu2tau

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
            vector<TLorentzVector> partP;
            vector<int> partId;
            partP.push_back(L11P4); partP.push_back(L12P4); partP.push_back(L21P4); partP.push_back(L22P4);
            partId.push_back(id_lep[L1]); partId.push_back(id_lep[L2]); partId.push_back(id_lep[L3]); partId.push_back(id_lep[L4]);

            mela::computeAngles(partP[0],partId[0],partP[1],partId[1],partP[2],partId[2],partP[3],partId[3],cosThetaStar,cosTheta1,cosTheta2,Phi,Phi1);
            
        }
    }
        
    std::sort(leptonsS1.begin(), leptonsS1.end(), sortByPt);
    for (unsigned int i=0; i<leptonsS1.size(); ++i) {
        //std::cout<<"lepton id: "<<leptons[i].pdgId()<<" pt: "<<leptons[i].pt()<<" eta: "<<leptons[i].eta()<<" status: "<<leptons[i].status()<<" mother: "<<leptons[i].mother()->pdgId()<<std::endl;
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
        //std::cout<< "p4_4l[0].Pt()" << p4_4l[0].Pt() << std::endl;
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
            //std::cout<<"gen jet: "<<genJet->pdgId()<<" pt: "<<genJet->pt()<<" eta: "<<genJet->eta()<<std::endl;
        }
    }
    std::sort(jets.begin(), jets.end(), sortJetsByPt);
    for (unsigned int i=0; i<jets.size(); ++i) {
        new ( (*p4_jet)[i] ) TLorentzVector(jets[i].px(),jets[i].py(),jets[i].pz(),jets[i].energy());
    }

    for (unsigned int i=0; i<jets.size(); ++i) {
        for (unsigned int j=0; j<jets.size(); ++j) {
            if (i==j) continue;
            TLorentzVector *jet1 = new TLorentzVector(jets[i].px(),jets[i].py(),jets[i].pz(),jets[i].energy());
            TLorentzVector *jet2 = new TLorentzVector(jets[j].px(),jets[j].py(),jets[j].pz(),jets[j].energy());
            TLorentzVector dijet = (*jet1)+(*jet2);
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
