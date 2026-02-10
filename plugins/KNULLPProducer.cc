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
#include "PhysicsTools/EXOnanoAOD/interface/PhoVars.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#define OBJECTARRAYSIZE 1000
#define MAX_NPV 200

class KNULLPProducer : public edm::stream::EDProducer<> {
	public:
	explicit KNULLPProducer(edm::ParameterSet const& cfg):
	src_(cfg.getParameter<edm::InputTag>("src")),
	token_(consumes<edm::View<pat::Photon>>(src_)),
	v_photonsInputTag_(cfg.getParameter<std::vector<edm::InputTag>>("photons")),
	verticesToken_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vertices"))),
	packedPFCandsToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("packedPfCands"))),
        tracksToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("tracks")))
	{
		produces<nanoaod::FlatTable>("Photon");
		produces<nanoaod::FlatTable>("pho");
		produces<nanoaod::FlatTable>("Razor");
		for (auto const& tag : v_photonsInputTag_) {
			v_photonsToken_.push_back(consumes<pat::PhotonCollection>(tag));
		}
    }

	void fillPFIsoVar(const pat::Photon& pho,
		const int photon_index,
		const pat::PackedCandidateCollection& packedPFCands,
		const reco::VertexCollection& vertices) {

		// https://github.com/cms-lpc-llp/DelayedPhotonTuplizer/blob/master/plugins/RazorTuplizer.cc#L2022
	    // Compute PF isolation

	    const float coneSizeDR = 0.3;
	    const float dxyMax = 0.1;
	    const float dzMax = 0.2;
	    float chargedIsoSumAllVertices[MAX_NPV];
	    for (int q=0;q<MAX_NPV;++q) chargedIsoSumAllVertices[q] = 0.0;
	    float chargedIsoSum = 0;
	    float neutralHadronIsoSum = 0;
	    float photonIsoSum = 0;
	    // First, find photon direction with respect to the good PV
	    math::XYZVector photon_directionWrtVtx(
	    	pho.superCluster()->x() - myPV->x(),
	    	pho.superCluster()->y() - myPV->y(),
	    	pho.superCluster()->z() - myPV->z()
	    	);
		// Loop over all PF candidates
	    for (const pat::PackedCandidate &candidate : packedPFCands) {
	    	// Check if this candidate is within the isolation cone
	    	float dR=deltaR(
	    		photon_directionWrtVtx.Eta(),photon_directionWrtVtx.Phi(),
	    		candidate.eta(), candidate.phi()
	    		);
	    	if ( dR > coneSizeDR ) continue;
	    	// Check if this candidate is not in the footprint

	    	bool inFootprint = false;
	    	for (auto itr : pho.associatedPackedPFCandidates()) {
	    		if ( &(*itr) == &candidate) {
	    			inFootprint = true;
	    		}
	    	}
	    	if ( inFootprint ) continue;

	    	// Find candidate type
	    	reco::PFCandidate::ParticleType thisCandidateType = reco::PFCandidate::X;

			// the neutral hadrons and charged hadrons can be of pdgId types
			// only 130 (K0L) and +-211 (pi+-) in packed candidates

	    	const int pdgId = candidate.pdgId();
	    	if( pdgId == 22 )
	    		thisCandidateType = reco::PFCandidate::gamma;
	    	else if( abs(pdgId) == 130) // PDG ID for K0L
	    		thisCandidateType = reco::PFCandidate::h0;
	    	else if( abs(pdgId) == 211) // PDG ID for pi+-
	    		thisCandidateType = reco::PFCandidate::h;

	    	// Increment the appropriate isolation sum
			if ( thisCandidateType == reco::PFCandidate::h && candidate.hasTrackDetails() ) {
				// for charged hadrons, additionally check consistency
				// with the PV
				float dxy = -999, dz = -999;

				//For the primary vertex
				dz = candidate.pseudoTrack().dz(myPV->position());
				dxy =candidate.pseudoTrack().dxy(myPV->position());
				if (fabs(dz) <= dzMax && fabs(dxy) <= dxyMax) {
					chargedIsoSum += candidate.pt();
				}

				//loop over all vertices
				for(int q = 0; q < nPVAll; q++) {
					if(!(vertices.at(q).isValid() && !vertices.at(q).isFake())) continue;
					dz = candidate.pseudoTrack().dz(vertices.at(q).position());
					dxy =candidate.pseudoTrack().dxy(vertices.at(q).position());
					if (fabs(dz) > dzMax) continue;
					if(fabs(dxy) > dxyMax) continue;
					// The candidate is eligible, increment the isolation
					chargedIsoSumAllVertices[q] += candidate.pt();
				}
			}
	    	if( thisCandidateType == reco::PFCandidate::h0 )
	    		neutralHadronIsoSum += candidate.pt();
	    	if( thisCandidateType == reco::PFCandidate::gamma )
	    		photonIsoSum += candidate.pt();
	    }
		//fill the proper variables

		for(int q = 0; q < nPVAll; q++) {
			// FIXME
			//vars_.sumChargedHadronPtAllVertices[i][q] = chargedIsoSumAllVertices[q];
			//vars_.sumChargedHadronPtAllVertices[i*nPVAll + q] = chargedIsoSumAllVertices[q]; // flatten?
		}
		vars_.sumChargedHadronPt[photon_index] = chargedIsoSum;
		vars_.sumNeutralHadronEt[photon_index] = neutralHadronIsoSum;
		vars_.sumPhotonEt[photon_index] = photonIsoSum;
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

	    // TODO think if it is better to make a separate methode
	    // slimmedPhoton and slimmedOOTPhoton
	    edm::Handle<pat::PhotonCollection> itPhotons;
	    edm::Handle<pat::PhotonCollection> ootPhotons;
	    edm::Handle<reco::VertexCollection> vertices;
	    edm::Handle<pat::PackedCandidateCollection> packedPFCands;
	    edm::Handle<reco::TrackCollection> tracks;

	    iEvent.getByToken(v_photonsToken_[0], itPhotons);
	    iEvent.getByToken(v_photonsToken_[1], ootPhotons);

	    iEvent.getByToken(verticesToken_, vertices);
	    iEvent.getByToken(packedPFCandsToken_, packedPFCands);
            iEvent.getByToken(tracksToken_, tracks);

	    //select the primary vertex, if any
	    nPV = 0; 
	    myPV = &(vertices->front());
	    nPVAll = std::min(int(vertices->size()),int(MAX_NPV));
	    bool foundPV = false;
	    for(unsigned int i = 0; i < vertices->size(); i++) {
	    	if(vertices->at(i).isValid() && !vertices->at(i).isFake()) {
	    	        if (!foundPV) {
                            myPV = &(vertices->at(i));
	    		    foundPV = true;
	    		}
	    		nPV++;
	    	}
	    }

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
	    		}
	    		else {
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

	    // fill table
	    const size_t n_final_photons = final_photons.size();
	    auto final_photons_tab = std::make_unique<nanoaod::FlatTable>(n_final_photons, "pho", false, false);
	    vars_.resize(n_final_photons);

	    for (size_t i = 0; i < n_final_photons; ++i) {
	    	const auto& pho = *final_photons[i];
	    	vars_.fillFromPho(pho, i);
	    	vars_.isOOT[i] = isOOT[i];

		//**************************************************************************
		//Veto photons that overlap with tracks
		//**************************************************************************
		bool closeToTrack = false;
		for (const auto & track : *tracks)
		{
		    if (track.pt() < 5) continue;
		    if (reco::deltaR(pho, track) < 0.2) closeToTrack = true;
		}
		vars_.trackMatching[i] = closeToTrack;

	    	fillPFIsoVar(pho, i, *packedPFCands, *vertices);
	    } // n_final_photons

	    addPhoColumns(*final_photons_tab, vars_);

	    // test
	    auto razor_event = std::make_unique<nanoaod::FlatTable>(1, "Razor", /*singleton=*/true, /*extension=*/false);
	    razor_event->addColumnValue<float>("pv_x", myPV->x(), "pv_x", 10);
	    razor_event->addColumnValue<float>("pv_y", myPV->y(), "pv_y", 10);
	    razor_event->addColumnValue<float>("pv_z", myPV->z(), "pv_z", 10);

	    iEvent.put(std::move(tab), "Photon");  // this is just to check
	    iEvent.put(std::move(final_photons_tab), "pho");
	    iEvent.put(std::move(razor_event), "Razor");
	}

    private:
	edm::InputTag src_;
	edm::EDGetTokenT<edm::View<pat::Photon>> token_;

	std::vector<edm::InputTag> v_photonsInputTag_;
	std::vector<edm::EDGetTokenT<pat::PhotonCollection>> v_photonsToken_;

	edm::EDGetTokenT<reco::VertexCollection> verticesToken_;
	edm::EDGetTokenT<pat::PackedCandidateCollection> packedPFCandsToken_;
        edm::EDGetTokenT<reco::TrackCollection> tracksToken_;

	int nPV;
	int nPVAll;
	const reco::Vertex *myPV;
	PhoVars vars_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(KNULLPProducer);
