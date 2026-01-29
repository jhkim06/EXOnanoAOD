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

class KNULLPProducer : public edm::stream::EDProducer<> {
public:
  explicit KNULLPProducer(edm::ParameterSet const& cfg)
      : src_(cfg.getParameter<edm::InputTag>("src")),
        token_(consumes<edm::View<pat::Photon>>(src_)) {
    // instance label "Photon" -> branches go under Photon_*
    produces<nanoaod::FlatTable>("Photon");
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

    iEvent.put(std::move(tab), "Photon");
  }

private:
  edm::InputTag src_;
  edm::EDGetTokenT<edm::View<pat::Photon>> token_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(KNULLPProducer);
