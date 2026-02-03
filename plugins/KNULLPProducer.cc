// -*- C++ -*-
//
// Package:    PhysicsTools/EXOnanoAOD
// Class:      KNULLPProducer
//
/**\class KNULLPProducer

 Description: Simple template producer as a guideline for createting custom nanoAOD producers.
 Example included for BeamSpot object.

*/
//
// Original Author:  Lovisa Rygaard
//         Created:  Fri, 14 Feb 2025 08:13:36 GMT
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

// nanoAOD include files
#include "DataFormats/NanoAOD/interface/FlatTable.h"

// object specific include files
#include "DataFormats/BeamSpot/interface/BeamSpot.h"  // replace with inclusions 

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "TLorentzVector.h"

#define OBJECTARRAYSIZE 1000

class KNULLPProducer : public edm::stream::EDProducer<> {
public:
  explicit KNULLPProducer(edm::ParameterSet const& cfg)
      : src_(cfg.getParameter<edm::InputTag>("src")),
        token_(consumes<edm::View<pat::Photon>>(src_)),
        v_photonsInputTag_(cfg.getParameter<std::vector<edm::InputTag>>("photons")) 
  {
    // instance label "Photon" -> branches go under Photon_*
    produces<nanoaod::FlatTable>("Photon");

    produces<nanoaod::FlatTable>("pho");

    for (auto const& tag : v_photonsInputTag_) {
        v_photonsToken_.push_back(consumes<pat::PhotonCollection>(tag));
    }
  }

  void produce(edm::Event& iEvent, edm::EventSetup const&) override {
    edm::Handle<edm::View<pat::Photon>> photons;
    iEvent.getByToken(token_, photons);

    const size_t n = photons->size();

    // extension=true so it appends to existing Photon table
    auto tab = std::make_unique<nanoaod::FlatTable>(n, "Photon", /*singleton=*/false, /*extension=*/true);

    std::vector<float> pt(n, -999.f);
    std::vector<float> sMajor(n, -999.f);
    std::vector<float> sMinor(n, -999.f);

    for (size_t i = 0; i < n; ++i) {
      const auto& pho = photons->at(i);

      pt[i] = pho.pt();

      const auto& ss = pho.full5x5_showerShapeVariables();
      sMajor[i] = ss.smMajor;
      sMinor[i] = ss.smMinor;
    }

    tab->addColumn<float>("pt_TEST", pt, "pT/E (custom)", 10);  // just to check if the same photon object is used comparing original EXOnanoAOD
    tab->addColumn<float>("sMajor", sMajor, "ECAL shower shape sMajor (full5x5)", 10);
    tab->addColumn<float>("sMinor", sMinor, "ECAL shower shape sMinor (full5x5)", 10);


    // slimmedPhoton and slimmedOOTPhoton
   
    //float metType1Px = -99.0;
    //float metType1Py = -99.0;
    //float metType1Pt = -99.0;
    //float metType1Phi = -99.0;

    //int nPhotons_overlap = 0;

    edm::Handle<pat::PhotonCollection> itPhotons;
    edm::Handle<pat::PhotonCollection> ootPhotons;

    iEvent.getByToken(v_photonsToken_[0], itPhotons);
    iEvent.getByToken(v_photonsToken_[1], ootPhotons);

    std::vector<char> keepIT(itPhotons->size(), 1);
    std::vector<char> keepOOT(ootPhotons->size(), 1);

    for (size_t i = 0; i < itPhotons->size(); ++i) {
      if (!keepIT[i]) continue;
      const auto& it = (*itPhotons)[i];

      for (size_t j = 0; j < ootPhotons->size(); ++j) {
	if (!keepOOT[j]) continue;
	const auto& oot = (*ootPhotons)[j];

	if (deltaR(it.eta(), it.phi(), oot.eta(), oot.phi()) > 0.3) continue;

	if (it.pt() < oot.pt()) {
	  keepIT[i] = 0;
          // TODO correct MET
	  break;              // IT removed, stop checking this IT photon
	} else {
	  keepOOT[j] = 0;     // OOT removed, continue checking IT vs others
          // TODO correct MET
	}
      }
    } 

    // photons(in-time or out-of-time) without overlap
    std::vector<const pat::Photon*> final_photons;
    std::vector<uint8_t> isOOT;

    for (size_t i = 0; i < itPhotons->size(); ++i) {
      if (!keepIT[i]) continue;
      final_photons.push_back(&(*itPhotons)[i]);
      isOOT.push_back(0);
    }

    for (size_t j = 0; j < ootPhotons->size(); ++j) {
      if (!keepOOT[j]) continue;
      final_photons.push_back(&(*ootPhotons)[j]);
      isOOT.push_back(1);
    }
    
    const size_t n_final_photons = final_photons.size();
    auto final_photons_tab = std::make_unique<nanoaod::FlatTable>(n_final_photons, "pho", false, false);
    std::vector<float> final_photons_pt(n_final_photons, -999.f);
    std::vector<float> final_photons_sMajor(n_final_photons, -999.f);
    std::vector<float> final_photons_sMinor(n_final_photons, -999.f);
  
    for (size_t r = 0; r < n_final_photons; ++r) {
      const auto& p = *final_photons[r];
      final_photons_pt[r] = p.pt();

      const auto& ss = p.full5x5_showerShapeVariables();
      final_photons_sMajor[r] = ss.smMajor;
      final_photons_sMinor[r] = ss.smMinor;
    }

    final_photons_tab->addColumn<float>("pt", final_photons_pt, "pT", 10);
    final_photons_tab->addColumn<float>("sMajor", final_photons_sMajor, "ECAL shower shape sMajor (full5x5)", 10);
    final_photons_tab->addColumn<float>("sMinor", final_photons_sMinor, "ECAL shower shape sMinor (full5x5)", 10);
    final_photons_tab->addColumn<uint8_t>("isOOT", isOOT, "1 if OOT", -1);
    
    iEvent.put(std::move(tab), "Photon");  // this is just to check
    iEvent.put(std::move(final_photons_tab), "pho");
  }

private:
  edm::InputTag src_;
  edm::EDGetTokenT<edm::View<pat::Photon>> token_;

  std::vector<edm::InputTag> v_photonsInputTag_;
  std::vector<edm::EDGetTokenT<pat::PhotonCollection>> v_photonsToken_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(KNULLPProducer);
