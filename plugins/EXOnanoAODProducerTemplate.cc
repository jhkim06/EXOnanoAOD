// -*- C++ -*-
//
// Package:    PhysicsTools/EXOnanoAOD
// Class:      EXOnanoAODProducerTemplate
//
/**\class EXOnanoAODProducerTemplate

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

class EXOnanoAODProducerTemplate : public edm::stream::EDProducer<> {
protected:
  const edm::InputTag inputLabel;
  const edm::EDGetTokenT<reco::BeamSpot> inputTag_;  // add your object class
  
public:
  EXOnanoAODProducerTemplate(edm::ParameterSet const& params)
    :
    inputLabel(params.getParameter<edm::InputTag>("inputExample")), // add your input name
    inputTag_(consumes<reco::BeamSpot>(inputLabel))  // add your object class
    {
    produces<nanoaod::FlatTable>("tableName");  // add your output name
  }

  ~EXOnanoAODProducerTemplate() override {}

  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override {
    
    edm::Handle<reco::BeamSpot> inputHandle;  // add your object class
    iEvent.getByToken(inputTag_, inputHandle);
    auto outputTab = std::make_unique<nanoaod::FlatTable>(1, "tableName", false, false);  // add your output name

    const auto& input = *inputHandle;  // example of using getting beamspot object as input

    // add variables as vectors to store in table
    std::vector<float> x, y, z;
    x.push_back(input.position().x());
    y.push_back(input.position().y());
    z.push_back(input.position().z());

    // add columns to table
    outputTab->addColumn<float>("x", x, "");
    outputTab->addColumn<float>("y", y, "");
    outputTab->addColumn<float>("z", z, "");

    iEvent.put(std::move(outputTab), "tableName");  // add your output name
  }

};
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(EXOnanoAODProducerTemplate);