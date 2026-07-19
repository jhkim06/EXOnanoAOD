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
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <vector>

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
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"

// ECAL Record info (Pedestals)
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"

//ECAL conditions
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"

#define OBJECTARRAYSIZE 1000
#define MAX_NPV 200

class KNULLPProducer : public edm::stream::EDProducer<> {
	public:
	explicit KNULLPProducer(edm::ParameterSet const& cfg):
	src_(cfg.getParameter<edm::InputTag>("src")),
	token_(consumes<edm::View<pat::Photon>>(src_)),

	v_photonsInputTag_(cfg.getParameter<std::vector<edm::InputTag>>("photons")),
	metToken_(consumes<pat::METCollection>(cfg.getParameter<edm::InputTag>("mets"))),
    electronsToken_(consumes<edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("electrons"))),
	jetsToken_(consumes<pat::JetCollection>(cfg.getParameter<edm::InputTag>("jets"))),
	verticesToken_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vertices"))),
	packedPFCandsToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("packedPfCands"))),
	tracksToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
	conversionsToken_(consumes<std::vector<reco::Conversion>>(cfg.getParameter<edm::InputTag>("conversions"))),
	singleLegConversionsToken_(consumes<std::vector<reco::Conversion> >(cfg.getParameter<edm::InputTag>("singleLegConversions"))),
	beamSpotToken_(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamSpot"))),
	triggerObjectsToken_(consumes<pat::TriggerObjectStandAloneCollection>(cfg.getParameter<edm::InputTag>("triggerObjects"))),
	ebRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(cfg.getParameter<edm::InputTag>("ebRecHits"))),
	eeRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(cfg.getParameter<edm::InputTag>("eeRecHits"))),
	caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
	laserToken_(esConsumes<EcalLaserDbService, EcalLaserDbRecord>()),
	pedestalToken_(esConsumes<EcalPedestals, EcalPedestalsRcd>())
	{
		produces<nanoaod::FlatTable>("Photon");
		produces<nanoaod::FlatTable>("pho");
		produces<nanoaod::FlatTable>("phoVtx");
		produces<nanoaod::FlatTable>("phoECALIdx"); 
		produces<nanoaod::FlatTable>("phoSeedIdx"); 
		produces<nanoaod::FlatTable>("ecalRechit");
		produces<nanoaod::FlatTable>("Razor");
		for (auto const& tag : v_photonsInputTag_) {
			v_photonsToken_.push_back(consumes<pat::PhotonCollection>(tag));
		}
	}

    void fillElctrons(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
        for (size_t i = 0; i < electrons->size(); ++i){
			const auto ele = electrons->ptrAt(i);

			ele_SeedRechitID.push_back(ele->superCluster()->seed()->seed().rawId());

			//*************************************************
			//Find relevant rechits
			//*************************************************
			std::vector<uint> rechits; rechits.clear();
			const std::vector< std::pair<DetId, float>>& v_id =ele->superCluster()->seed()->hitsAndFractions();
			for ( size_t i = 0; i < v_id.size(); ++i ) {
				ecalRechitID_ToBeSaved.push_back(v_id[i].first.rawId());
				rechits.push_back(v_id[i].first.rawId());
			}
			ecalRechitEtaPhi_ToBeSaved.push_back( std::pair<double,double>( ele->superCluster()->eta(), ele->superCluster()->phi() ));
			ele_EcalRechitID.push_back(rechits);

        }
    }

	void fillPFIsoVar(const pat::Photon& pho, const int photon_index) {
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
		for (const pat::PackedCandidate &candidate : *packedPFCands) {
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
					if(!(vertices->at(q).isValid() && !vertices->at(q).isFake())) continue;
					dz = candidate.pseudoTrack().dz(vertices->at(q).position());
					dxy =candidate.pseudoTrack().dxy(vertices->at(q).position());
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
		} // PF candidates loop
		//fill the proper variables
		for(int q = 0; q < nPVAll; q++) {
			sumChargedHadronPtAllVertices.push_back(chargedIsoSumAllVertices[q]);
		}
		vars_.sumChargedHadronPt[photon_index] = chargedIsoSum;
		vars_.sumNeutralHadronEt[photon_index] = neutralHadronIsoSum;
		vars_.sumPhotonEt[photon_index] = photonIsoSum;

		//*****************************************************************
		//Compute Worst Isolation Looping over all vertices
		//*****************************************************************
		const double ptMin = 0.0;
		const float dRvetoBarrel = 0.0;
		const float dRvetoEndcap = 0.0;
		float dRveto = 0;
		if (pho.isEB())
			dRveto = dRvetoBarrel;
		else
			dRveto = dRvetoEndcap;

		float worstIsolation = 999;
		std::vector<float> allIsolations;
		for(unsigned int ivtx=0; ivtx<vertices->size(); ++ivtx) {

			// Shift the photon according to the vertex
			reco::VertexRef vtx(vertices, ivtx);
			math::XYZVector photon_directionWrtVtx(pho.superCluster()->x() - vtx->x(),
			pho.superCluster()->y() - vtx->y(),
			pho.superCluster()->z() - vtx->z());

			float sum = 0;
			// Loop over all PF candidates
			for (const pat::PackedCandidate &candidate : *packedPFCands) {
				//require that PFCandidate is a charged hadron
                const int pdgId = candidate.pdgId();
                if( abs(pdgId) != 211) continue;

                if (candidate.pt() < ptMin) continue;

                float dxy = -999, dz = -999;
                dz = candidate.dz(myPV->position());
                dxy =candidate.dxy(myPV->position());
                if ( fabs(dxy) > dxyMax) continue;
                if ( fabs(dz) > dzMax) continue;

                float dR = deltaR(
                	photon_directionWrtVtx.Eta(), photon_directionWrtVtx.Phi(),
                	candidate.eta(), candidate.phi()
                	);
                if (dR > coneSizeDR || dR < dRveto) continue;
				sum += candidate.pt();
			}
			allIsolations.push_back(sum);
		}

		if( allIsolations.size()>0 )
			worstIsolation = * std::max_element( allIsolations.begin(), allIsolations.end() );

		vars_.sumWorstVertexChargedHadronPt[photon_index] = worstIsolation;
	}

	bool fillTriggerVars(const std::vector<const pat::Photon*>& final_photons) {
		bool hasHLTCscClusterVeryLoose = false;

		for (const auto& trigObj : *triggerObjects) {
			if (trigObj.hasFilterLabel("hltCscClusterVeryLoose")) {
				hasHLTCscClusterVeryLoose = true;
			}

			const bool passEG20HE = trigObj.hasFilterLabel("hltEG20HEFilterUnseeded");
			const bool passEG30HE = trigObj.hasFilterLabel("hltEG30HEFilterUnseeded");
			const bool passEG20Et = trigObj.hasFilterLabel("hltEG20EtFilterUnseeded");
			const bool passEG30Et = trigObj.hasFilterLabel("hltEG30EtFilterUnseeded");
			if (!(passEG20HE || passEG30HE || passEG20Et || passEG30Et)) continue;

			for (size_t ipho = 0; ipho < final_photons.size(); ++ipho) {
				const auto& pho = *final_photons[ipho];
				const float dr = reco::deltaR(pho.eta(), pho.phi(), trigObj.eta(), trigObj.phi());
				if (passEG20HE && dr < vars_.drHLTEG20HEUnseeded[ipho]) vars_.drHLTEG20HEUnseeded[ipho] = dr;
				if (passEG30HE && dr < vars_.drHLTEG30HEUnseeded[ipho]) vars_.drHLTEG30HEUnseeded[ipho] = dr;
				if (passEG20Et && dr < vars_.drHLTEG20EtUnseeded[ipho]) vars_.drHLTEG20EtUnseeded[ipho] = dr;
				if (passEG30Et && dr < vars_.drHLTEG30EtUnseeded[ipho]) vars_.drHLTEG30EtUnseeded[ipho] = dr;
			}
		}

		const float matchCone = 0.3;
		for (size_t ipho = 0; ipho < final_photons.size(); ++ipho) {
			vars_.passHLTEG20HEUnseeded[ipho] = vars_.drHLTEG20HEUnseeded[ipho] < matchCone;
			vars_.passHLTEG30HEUnseeded[ipho] = vars_.drHLTEG30HEUnseeded[ipho] < matchCone;
			vars_.passHLTEG20EtUnseeded[ipho] = vars_.drHLTEG20EtUnseeded[ipho] < matchCone;
			vars_.passHLTEG30EtUnseeded[ipho] = vars_.drHLTEG30EtUnseeded[ipho] < matchCone;
		}

		return hasHLTCscClusterVeryLoose;
	}

	void fillPhoVtxVars(const std::vector<const pat::Photon*>& final_photons) {
		const size_t nRows = final_photons.size() * nPVAll;
		phoVtxPhotonIdx.assign(nRows, -1);
		phoVtxVertexIdx.assign(nRows, -1);
		phoVtxSumPx.assign(nRows, 0.f);
		phoVtxSumPy.assign(nRows, 0.f);

		for (size_t ipho = 0; ipho < final_photons.size(); ++ipho) {
			for (int ipv = 0; ipv < nPVAll; ++ipv) {
				const size_t idx = ipho * nPVAll + ipv;
				phoVtxPhotonIdx[idx] = static_cast<int32_t>(ipho);
				phoVtxVertexIdx[idx] = static_cast<int32_t>(ipv);
			}
		}

		for (const auto& pfcand : *packedPFCands) {
			if (pfcand.charge() == 0) continue;

			double mindz = std::numeric_limits<double>::max();
			int ipvmin = -1;
			for (int ipv = 0; ipv < nPVAll; ++ipv) {
				const reco::Vertex& vtx = vertices->at(ipv);
				const double dz = std::abs(pfcand.dz(vtx.position()));
				if (dz < mindz) {
					mindz = dz;
					ipvmin = ipv;
				}
			}

			if (mindz >= 0.2 || ipvmin < 0) continue;
			const reco::Vertex& vtx = vertices->at(ipvmin);
			for (size_t ipho = 0; ipho < final_photons.size(); ++ipho) {
				const auto& pho = *final_photons[ipho];
				math::XYZVector phodir(pho.superCluster()->x() - vtx.x(),
				                       pho.superCluster()->y() - vtx.y(),
				                       pho.superCluster()->z() - vtx.z());
				if (reco::deltaR(phodir.Eta(), phodir.Phi(), pfcand.eta(), pfcand.phi()) >= 0.05) continue;

				const size_t idx = ipho * nPVAll + ipvmin;
				phoVtxSumPx[idx] += pfcand.px();
				phoVtxSumPy[idx] += pfcand.py();
			}
		}
	}

	bool fillEcalRechits(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
		edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> ebRecHits;
		edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> eeRecHits;
		iEvent.getByToken(ebRecHitsToken_,ebRecHits);
		iEvent.getByToken(eeRecHitsToken_,eeRecHits);

		// geometry (from ECAL ELF)
        //edm::ESHandle<CaloGeometry> geoHandle;
        //iSetup.get<CaloGeometryRecord>().get(geoHandle);

        //const CaloSubdetectorGeometry *barrelGeometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
        //const CaloSubdetectorGeometry *endcapGeometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);

        const CaloGeometry& geo = iSetup.getData(caloGeomToken_);
        const CaloSubdetectorGeometry* barrelGeometry = geo.getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
        const CaloSubdetectorGeometry* endcapGeometry = geo.getSubdetectorGeometry(DetId::Ecal, EcalEndcap);

		std::map<uint, uint> mapRecHitIdToIndex; mapRecHitIdToIndex.clear();
		uint rechitIndex = 0;

        // ECAL conditions
        const auto& laser = iSetup.getData(laserToken_);

        // ECAL Pedestal
        const auto& pedestals = iSetup.getData(pedestalToken_);

        //ECAL conditions
        //edm::ESHandle<EcalLaserDbService> laser_;
        //iSetup.get<EcalLaserDbRecord>().get(laser_);

        //ECAL Pedestal
        //edm::ESHandle<EcalPedestals> pedestalsH;
        //iSetup.get<EcalPedestalsRcd>().get(pedestalsH);

        auto saveRecHit = [&](const EcalRecHitCollection::const_iterator& recHit,
                              const DetId& recHitId,
                              const GlobalPoint& recHitPos) {
        	mapRecHitIdToIndex[recHitId.rawId()] = rechitIndex;
        	ecalRechitID.push_back(static_cast<int32_t>(recHitId.rawId()));
        	ecalRechitEta.push_back(recHitPos.eta());
        	ecalRechitPhi.push_back(recHitPos.phi());
        	ecalRechitX.push_back(recHitPos.x());
        	ecalRechitY.push_back(recHitPos.y());
        	ecalRechitZ.push_back(recHitPos.z());
        	ecalRechitEnergy.push_back(recHit->energy());
        	ecalRechitTime.push_back(recHit->time());
        	ecalRechitFlagOOT.push_back(recHit->checkFlag(EcalRecHit::kOutOfTime));
        	ecalRechitGainSwitch1.push_back(recHit->checkFlag(EcalRecHit::kHasSwitchToGain1));
        	ecalRechitGainSwitch6.push_back(recHit->checkFlag(EcalRecHit::kHasSwitchToGain6));
        	ecalRechitTranspCorr.push_back(laser.getLaserCorrection(recHitId, iEvent.eventAuxiliary().time()));

        	float pedrms12 = -999.f;
        	float pedrms6 = -999.f;
        	float pedrms1 = -999.f;
        	const auto& pediter = pedestals.find(recHitId);
        	if (pediter != pedestals.end()) {
        		const auto& ped = (*pediter);
        		pedrms12 = ped.rms(1);
        		pedrms6 = ped.rms(2);
        		pedrms1 = ped.rms(3);
        	}
        	ecalRechitPedrms12.push_back(pedrms12);
        	ecalRechitPedrms6.push_back(pedrms6);
        	ecalRechitPedrms1.push_back(pedrms1);
        	rechitIndex++;
        };

        //Barrel Rechits
        for (EcalRecHitCollection::const_iterator recHit = ebRecHits->begin(); recHit != ebRecHits->end(); ++recHit) {
        	// first get detector id
        	const DetId recHitId = recHit->detid();
			//const uint32_t rhId  = recHitId.rawId();

			//Find rechits by ID that are explicitly marked to be saved.
            bool matchedRechit = false;
            std::vector<uint>::iterator it;
            it = find (ecalRechitID_ToBeSaved.begin(), ecalRechitID_ToBeSaved.end(), recHitId.rawId());
            if (it != ecalRechitID_ToBeSaved.end()) {
            	matchedRechit = true;
            }

            const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();

            // Find rechits by deltaR proximity
            bool dRProximity = false;
            for (int k=0; k<int(ecalRechitEtaPhi_ToBeSaved.size()); ++k) {
                if ( deltaR(ecalRechitEtaPhi_ToBeSaved[k].first,
                	ecalRechitEtaPhi_ToBeSaved[k].second,
                	recHitPos.eta(), recHitPos.phi()) < 0.5 ) {
                	dRProximity = true;
                	}
            }
            bool dRJetProximity = false;
            for (int k=0; k<int(ecalRechitJetEtaPhi_ToBeSaved.size()); ++k) {
                if ( deltaR(ecalRechitJetEtaPhi_ToBeSaved[k].first, ecalRechitJetEtaPhi_ToBeSaved[k].second,
                            recHitPos.eta(), recHitPos.phi()) < 0.7) {
                    dRJetProximity = true;
                }
            }

            //skip rechits that are not relevant
            if (!(matchedRechit || dRProximity || dRJetProximity)) {
                continue;
            }

        	saveRecHit(recHit, recHitId, recHitPos);
        }

		//Endcap Rechits
		for (EcalRecHitCollection::const_iterator recHit = eeRecHits->begin(); recHit != eeRecHits->end(); ++recHit) {
			// first get detector id

			const DetId recHitId = recHit->detid();
			// const uint32_t rhId  = recHitId.rawId();

            bool matchedRechit = false;
            std::vector<uint>::iterator it;
            it = find (ecalRechitID_ToBeSaved.begin(), ecalRechitID_ToBeSaved.end(), recHitId.rawId());
            if (it != ecalRechitID_ToBeSaved.end()) {
                matchedRechit = true;
            }

            const auto recHitPos = endcapGeometry->getGeometry(recHitId)->getPosition();

            //Find rechits by deltaR proximity
            bool dRProximity = false;
            for (int k=0; k<int(ecalRechitEtaPhi_ToBeSaved.size()); ++k) {
                if ( deltaR(ecalRechitEtaPhi_ToBeSaved[k].first, ecalRechitEtaPhi_ToBeSaved[k].second,
                            recHitPos.eta(), recHitPos.phi()) < 0.5) {
                    dRProximity = true;
                }
            }

            bool dRJetProximity = false;
            for (int k=0; k<int(ecalRechitJetEtaPhi_ToBeSaved.size()); ++k) {
                if ( deltaR(ecalRechitJetEtaPhi_ToBeSaved[k].first, ecalRechitJetEtaPhi_ToBeSaved[k].second,
                            recHitPos.eta(), recHitPos.phi()) < 0.7) {
                    dRJetProximity = true;
                }
            }

            //skip rechits that are not relevant
            if (!(matchedRechit || dRProximity || dRJetProximity)) {
                continue;
            }

            saveRecHit(recHit, recHitId, recHitPos);
		}

        for (uint k=0; k<ele_EcalRechitID.size(); k++) {
        	for (uint l=0; l<ele_EcalRechitID[k].size(); l++) {
        		const auto mapIt = mapRecHitIdToIndex.find(ele_EcalRechitID[k][l]);
        		ele_EcalRechitIndex.push_back(mapIt != mapRecHitIdToIndex.end() ? static_cast<int32_t>(mapIt->second) : -1); // flatten
        	}
        	const auto seedMapIt = mapRecHitIdToIndex.find(ele_SeedRechitID[k]);
        	ele_SeedRechitIndex.push_back(seedMapIt != mapRecHitIdToIndex.end() ? static_cast<int32_t>(seedMapIt->second) : -1);
        }

        for (uint k=0; k<EcalRechitID.size(); k++) {
        	vars_.ecalRechitIdxStart[k] = static_cast<int32_t>(EcalRechitIndex.size());
        	vars_.nEcalRechits[k] = static_cast<int32_t>(EcalRechitID[k].size());
        	for (uint l=0; l<EcalRechitID[k].size(); l++) {
        		const auto mapIt = mapRecHitIdToIndex.find(EcalRechitID[k][l]);
        		EcalRechitIndex.push_back(mapIt != mapRecHitIdToIndex.end() ? static_cast<int32_t>(mapIt->second) : -1); // flatten
        	}
        	const auto seedMapIt = mapRecHitIdToIndex.find(SeedRechitID[k]);
        	const int32_t seedIndex = seedMapIt != mapRecHitIdToIndex.end() ? static_cast<int32_t>(seedMapIt->second) : -1;
        	SeedRechitIndex.push_back(seedIndex);
        	vars_.seedRechitIdx[k] = seedIndex;
        	vars_.seedRechitID[k] = static_cast<int32_t>(SeedRechitID[k]);
        	if (seedIndex >= 0) {
        		const size_t idx = seedIndex;
        		vars_.seedRechitEnergy[k] = ecalRechitEnergy[idx];
        		vars_.seedRechitTime[k] = ecalRechitTime[idx];
        		vars_.seedRechitEta[k] = ecalRechitEta[idx];
        		vars_.seedRechitPhi[k] = ecalRechitPhi[idx];
        		vars_.seedRechitX[k] = ecalRechitX[idx];
        		vars_.seedRechitY[k] = ecalRechitY[idx];
        		vars_.seedRechitZ[k] = ecalRechitZ[idx];
        		vars_.seedRechitFlagOOT[k] = ecalRechitFlagOOT[idx];
        		vars_.seedRechitGainSwitch1[k] = ecalRechitGainSwitch1[idx];
        		vars_.seedRechitGainSwitch6[k] = ecalRechitGainSwitch6[idx];
        	}
        }
		return true;
	}
    void fillPhoCovVar(const pat::Photon& pho, const int photon_index) {
        //conversion matching for beamspot pointing
		const reco::Conversion *convmatch = 0;
		double drmin = std::numeric_limits<double>::max();
        //double leg conversions
        for (const reco::Conversion &conv : *conversions) {
            if (conv.refittedPairMomentum().rho()<10.) continue;
            if (!conv.conversionVertex().isValid()) continue;
            if (TMath::Prob(conv.conversionVertex().chi2(),  conv.conversionVertex().ndof())<1e-6) continue;
            math::XYZVector mom(conv.refittedPairMomentum());
            math::XYZPoint scpos(pho.superCluster()->position());
            math::XYZPoint cvtx(conv.conversionVertex().position());
            math::XYZVector cscvector = scpos - cvtx;

            double dr = reco::deltaR(mom,cscvector);

            if (dr<drmin && dr<0.1) {
            	drmin = dr;
            	convmatch = &conv;
            }

            if (!convmatch) {
                drmin = std::numeric_limits<double>::max();
                //single leg conversions
                for (const reco::Conversion &conv : *singleLegConversions) {
                    math::XYZVector mom(conv.tracksPin()[0]);
                    math::XYZPoint scpos(pho.superCluster()->position());
                    math::XYZPoint cvtx(conv.conversionVertex().position());
                    math::XYZVector cscvector = scpos - cvtx;

                    double dr = reco::deltaR(mom,cscvector);

                    if (dr<drmin && dr<0.1) {
                        drmin = dr;
                        convmatch = &conv;
                    }
                }
            }

            //std::cout << "[DEBUG photon conversion] "
            //<< "convmatch=" << convmatch
            //<< " hasConversionTracks=" << pho.hasConversionTracks()
            //<< " nConversions=" << conversions->size()
            //<< std::endl;

            //matched conversion, compute conversion type
            //and extrapolation to beamline
            //*FIXME* Both of these additional two requirements are inconsistent and make the conversion
            //selection depend on poorly defined criteria, but we keep them for sync purposes
            if (convmatch && pho.hasConversionTracks() && conversions->size()>0) {
                int ntracks = convmatch->nTracks();

                math::XYZVector mom(ntracks==2 ? convmatch->refittedPairMomentum() : convmatch->tracksPin()[0]);
                math::XYZPoint scpos(pho.superCluster()->position());
                math::XYZPoint cvtx(convmatch->conversionVertex().position());
                math::XYZVector cscvector = scpos - cvtx;

                double z = cvtx.z();
                double rho = cvtx.rho();

                int legtype = ntracks==2 ? 0 : 1;
                int dettype = pho.isEB() ? 0 : 1;
                int postype =0;

                if (pho.isEB()) {
                    if (rho<15.) {
                        postype = 0;
                    }
                    else if (rho>=15. && rho<60.) {
                        postype = 1;
                    }
                    else {
                        postype = 2;
                    }
                }
                else {
                    if (std::abs(z) < 50.) {
                        postype = 0;
                    }
                    else if (std::abs(z) >= 50. && std::abs(z) < 100.) {
                        postype = 1;
                    }
                    else {
                        postype = 2;
                    } 
                }
                vars_.convType[photon_index] = legtype + 2*dettype + 4*postype;
                vars_.convTrkZ[photon_index] = cvtx.z() - ((cvtx.x()-beamSpot->x0())*mom.x()+(cvtx.y()-beamSpot->y0())*mom.y())/mom.rho() * mom.z()/mom.rho();
                vars_.convTrkClusZ[photon_index] = cvtx.z() - ((cvtx.x()-beamSpot->x0())*cscvector.x()+(cvtx.y()-beamSpot->y0())*cscvector.y())/cscvector.rho() * cscvector.z()/cscvector.rho();
            }
        }
    }

	void produce(edm::Event& iEvent, edm::EventSetup const& iSetup) override {
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
		iEvent.getByToken(metToken_, mets);
        iEvent.getByToken(electronsToken_, electrons);
        iEvent.getByToken(v_photonsToken_[0], itPhotons);
        iEvent.getByToken(v_photonsToken_[1], ootPhotons);
		iEvent.getByToken(jetsToken_, jets);
        iEvent.getByToken(verticesToken_, vertices);
        iEvent.getByToken(packedPFCandsToken_, packedPFCands);
        iEvent.getByToken(tracksToken_, tracks);
        iEvent.getByToken(conversionsToken_,conversions);
        iEvent.getByToken(singleLegConversionsToken_,singleLegConversions);
        iEvent.getByToken(beamSpotToken_,beamSpot);
        iEvent.getByToken(triggerObjectsToken_, triggerObjects);

        sumChargedHadronPtAllVertices.clear();
        phoVtxPhotonIdx.clear();
        phoVtxVertexIdx.clear();
        phoVtxSumPx.clear();
        phoVtxSumPy.clear();
        ele_SeedRechitID.clear();
        SeedRechitID.clear();
		ele_EcalRechitID.clear();
	    EcalRechitID.clear();
	    ecalRechitID_ToBeSaved.clear();
	    ecalRechitEtaPhi_ToBeSaved.clear();;
	    ecalRechitID.clear();
	    ecalRechitEta.clear();
	    ecalRechitPhi.clear();
	    ecalRechitX.clear();
	    ecalRechitY.clear();
	    ecalRechitZ.clear();
	    ecalRechitEnergy.clear();
	    ecalRechitTime.clear();
	    ecalRechitFlagOOT.clear();
	    ecalRechitGainSwitch1.clear();
	    ecalRechitGainSwitch6.clear();
	    ecalRechitTranspCorr.clear();
	    ecalRechitPedrms12.clear();
	    ecalRechitPedrms6.clear();
	    ecalRechitPedrms1.clear();

		const pat::MET &Met = mets->front();

        // loop over electrons for 
        fillElctrons(iEvent, iSetup);

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

		ecalRechitJetEtaPhi_ToBeSaved.clear();
		for (const pat::Jet &j : *jets) {
			if (j.pt() < 10) continue;
			ecalRechitJetEtaPhi_ToBeSaved.push_back( std::pair<double,double>( j.eta(), j.phi() ));
		}

		float metType1Px = Met.px();
		float metType1Py = Met.py();
		float metType1Pt = Met.pt();
		float metType1Phi = Met.phi(); 
		
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
					// correct MET
					metType1Px = metType1Px + it.px();
					metType1Py = metType1Py + it.py();
					TVector3 vec_met_temp(metType1Px, metType1Py, 0);			
					metType1Pt = sqrt(metType1Px*metType1Px+metType1Py*metType1Py);
					metType1Phi = vec_met_temp.Phi();
					break;              // IT removed, stop checking this IT photon
				}
				else {
					keepOOT[j] = 0;     // OOT removed, continue checking IT vs others
					// correct MET
					metType1Px = metType1Px - it.px();
					metType1Py = metType1Py - it.py();
					TVector3 vec_met_temp(metType1Px, metType1Py, 0);
					metType1Pt = sqrt(metType1Px*metType1Px+metType1Py*metType1Py);
					metType1Phi = vec_met_temp.Phi();
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

	    	SeedRechitID.push_back(pho.superCluster()->seed()->seed().rawId());
	    	const std::vector< std::pair<DetId, float>>& v_id =pho.seed()->hitsAndFractions();
	    	std::vector<uint> rechits;
	    	rechits.clear();
	    	for ( size_t i = 0; i < v_id.size(); ++i ) {
	    		//EcalRecHitCollection::const_iterator it = ebRecHits->find( v_id[i].first );

	    		ecalRechitID_ToBeSaved.push_back(v_id[i].first.rawId());
	    		rechits.push_back(v_id[i].first.rawId());

		    /*
                if (it != ebRecHits->end()) {
					float energy = it->energy() * v_id[i].second;
                    if (it->checkFlag(EcalRecHit::kHasSwitchToGain6)) anySwitchToGain6 = true;
                    if (it->checkFlag(EcalRecHit::kHasSwitchToGain1)) anySwitchToGain1 = true;
                    if ( energy > max ) {
                        max = energy;
                        maxSwitchToGain6 = it->checkFlag(EcalRecHit::kHasSwitchToGain6);
                        maxSwitchToGain1 = it->checkFlag(EcalRecHit::kHasSwitchToGain1);
                    }
                }
                else {
                    //cout << "rechit not found\n";
                }
	            */
	    	}
	    	EcalRechitID.push_back(rechits);
	    	ecalRechitEtaPhi_ToBeSaved.push_back(
	    		std::pair<double,double>(pho.superCluster()->eta(),
	    		pho.superCluster()->phi())
	    		);

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

            fillPhoCovVar(pho, i);
	    	fillPFIsoVar(pho, i);
	    } // n_final_photons

	    fillPhoVtxVars(final_photons);
	    const bool hasHLTCscClusterVeryLoose = fillTriggerVars(final_photons);

		ele_EcalRechitIndex.clear();
		EcalRechitIndex.clear();
	    ele_SeedRechitIndex.clear(); 
	    SeedRechitIndex.clear(); 
	    fillEcalRechits(iEvent, iSetup);

	    addPhoColumns(*final_photons_tab, vars_);

        // TODO flattened variables needs reference index
	    auto final_photons_sumChargedHadronPtAllVertices_tab = std::make_unique<nanoaod::FlatTable>(sumChargedHadronPtAllVertices.size(), "phoVtx", false, false);
		final_photons_sumChargedHadronPtAllVertices_tab->addColumn<int32_t>("phoIdx", phoVtxPhotonIdx, "photon index in pho table", -1);
		final_photons_sumChargedHadronPtAllVertices_tab->addColumn<int32_t>("vtxIdx", phoVtxVertexIdx, "vertex index in offlineSlimmedPrimaryVertices", -1);
		final_photons_sumChargedHadronPtAllVertices_tab->addColumn<float>("sumChargedHadronPtAllVertices", sumChargedHadronPtAllVertices, "flattened sumChargedHadronPtAllVertices", 10);
		final_photons_sumChargedHadronPtAllVertices_tab->addColumn<float>("vtxSumPx", phoVtxSumPx, "Run2-style charged PF candidate px sum near photon for this vertex", 10);
		final_photons_sumChargedHadronPtAllVertices_tab->addColumn<float>("vtxSumPy", phoVtxSumPy, "Run2-style charged PF candidate py sum near photon for this vertex", 10);

	    auto final_photons_EcalRechitIndex_tab = std::make_unique<nanoaod::FlatTable>(EcalRechitIndex.size(), "phoECALIdx", false, false);
		final_photons_EcalRechitIndex_tab->addColumn<int32_t>("EcalRechitIndex", EcalRechitIndex, "flattend EcalRechitIndex", -1);

	    auto final_photons_SeedRechitIndex_tab = std::make_unique<nanoaod::FlatTable>(SeedRechitIndex.size(), "phoSeedIdx", false, false);
		final_photons_SeedRechitIndex_tab->addColumn<int32_t>("SeedRechitIndex", SeedRechitIndex, "flattend SeedRechitIndex", -1);

	    auto ecal_rechit_tab = std::make_unique<nanoaod::FlatTable>(ecalRechitID.size(), "ecalRechit", false, false);
	    ecal_rechit_tab->addColumn<int32_t>("ID", ecalRechitID, "raw ECAL rechit DetId", -1);
	    ecal_rechit_tab->addColumn<float>("eta", ecalRechitEta, "ECAL rechit eta", 10);
	    ecal_rechit_tab->addColumn<float>("phi", ecalRechitPhi, "ECAL rechit phi", 10);
	    ecal_rechit_tab->addColumn<float>("x", ecalRechitX, "ECAL rechit x position", 10);
	    ecal_rechit_tab->addColumn<float>("y", ecalRechitY, "ECAL rechit y position", 10);
	    ecal_rechit_tab->addColumn<float>("z", ecalRechitZ, "ECAL rechit z position", 10);
	    ecal_rechit_tab->addColumn<float>("energy", ecalRechitEnergy, "ECAL rechit energy", 10);
	    ecal_rechit_tab->addColumn<float>("time", ecalRechitTime, "ECAL rechit time", 10);
	    ecal_rechit_tab->addColumn<uint8_t>("flagOOT", ecalRechitFlagOOT, "ECAL rechit kOutOfTime flag", -1);
	    ecal_rechit_tab->addColumn<uint8_t>("gainSwitch1", ecalRechitGainSwitch1, "ECAL rechit kHasSwitchToGain1 flag", -1);
	    ecal_rechit_tab->addColumn<uint8_t>("gainSwitch6", ecalRechitGainSwitch6, "ECAL rechit kHasSwitchToGain6 flag", -1);
	    ecal_rechit_tab->addColumn<float>("transpCorr", ecalRechitTranspCorr, "ECAL laser transparency correction", 10);
	    ecal_rechit_tab->addColumn<float>("pedrms12", ecalRechitPedrms12, "ECAL pedestal RMS gain 12", 10);
	    ecal_rechit_tab->addColumn<float>("pedrms6", ecalRechitPedrms6, "ECAL pedestal RMS gain 6", 10);
	    ecal_rechit_tab->addColumn<float>("pedrms1", ecalRechitPedrms1, "ECAL pedestal RMS gain 1", 10);

	    // Lets use this table for corrected MET
	    auto razor_event = std::make_unique<nanoaod::FlatTable>(1, "Razor", /*singleton=*/true, /*extension=*/false);
	    razor_event->addColumnValue<float>("pv_x", myPV->x(), "pv_x", 10);
	    razor_event->addColumnValue<float>("pv_y", myPV->y(), "pv_y", 10);
	    razor_event->addColumnValue<float>("pv_z", myPV->z(), "pv_z", 10);

	    razor_event->addColumnValue<float>("metType1Px", metType1Px, "metType1Px", 10);
	    razor_event->addColumnValue<float>("metType1Py", metType1Py, "metType1Py", 10);
	    razor_event->addColumnValue<float>("metType1Pt", metType1Pt, "metType1Pt", 10);
	    razor_event->addColumnValue<float>("metType1Phi", metType1Phi, "metType1Phi", 10);
	    razor_event->addColumnValue<uint8_t>("hasHLTCscClusterVeryLoose", static_cast<uint8_t>(hasHLTCscClusterVeryLoose), "1 if any trigger object has hltCscClusterVeryLoose filter label", -1);

	    iEvent.put(std::move(tab), "Photon");  // this is just to check
	    iEvent.put(std::move(final_photons_tab), "pho");
	    iEvent.put(std::move(final_photons_sumChargedHadronPtAllVertices_tab), "phoVtx");
	    iEvent.put(std::move(final_photons_EcalRechitIndex_tab), "phoECALIdx");
	    iEvent.put(std::move(final_photons_SeedRechitIndex_tab), "phoSeedIdx");
	    iEvent.put(std::move(ecal_rechit_tab), "ecalRechit");
	    iEvent.put(std::move(razor_event), "Razor");
	}

    private:
    edm::InputTag src_;
	edm::EDGetTokenT<edm::View<pat::Photon>> token_;

	std::vector<edm::InputTag> v_photonsInputTag_;
	std::vector<edm::EDGetTokenT<pat::PhotonCollection>> v_photonsToken_;

	edm::EDGetTokenT<pat::METCollection> metToken_;
    edm::EDGetTokenT<edm::View<reco::GsfElectron> > electronsToken_;
	edm::EDGetTokenT<pat::JetCollection> jetsToken_;
	edm::EDGetTokenT<reco::VertexCollection> verticesToken_;
	edm::EDGetTokenT<pat::PackedCandidateCollection> packedPFCandsToken_;
	edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
	edm::EDGetTokenT<std::vector<reco::Conversion>> conversionsToken_;
	edm::EDGetTokenT<std::vector<reco::Conversion> > singleLegConversionsToken_;
	edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
	edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsToken_;
	edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> ebRecHitsToken_;
	edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> eeRecHitsToken_;
	edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
	edm::ESGetToken<EcalLaserDbService, EcalLaserDbRecord> laserToken_;
	edm::ESGetToken<EcalPedestals, EcalPedestalsRcd> pedestalToken_;

	edm::Handle<pat::METCollection> mets;
    edm::Handle<pat::PhotonCollection> itPhotons;
    edm::Handle<pat::PhotonCollection> ootPhotons;
    edm::Handle<edm::View<reco::GsfElectron> > electrons;
	edm::Handle<pat::JetCollection> jets;
    edm::Handle<reco::VertexCollection> vertices;
    edm::Handle<pat::PackedCandidateCollection> packedPFCands;
    edm::Handle<reco::TrackCollection> tracks;
    edm::Handle<std::vector<reco::Conversion> > conversions;
    edm::Handle<std::vector<reco::Conversion>> singleLegConversions;
    edm::Handle<reco::BeamSpot> beamSpot;
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;

	int nPV;
	int nPVAll;
	const reco::Vertex *myPV;
	std::vector<float> sumChargedHadronPtAllVertices;
	std::vector<int32_t> phoVtxPhotonIdx;
	std::vector<int32_t> phoVtxVertexIdx;
	std::vector<float> phoVtxSumPx;
	std::vector<float> phoVtxSumPy;
	std::vector<std::vector<uint>> ele_EcalRechitID;
	std::vector<std::vector<uint>> EcalRechitID;
	std::vector<uint> ele_SeedRechitID;
	std::vector<uint> SeedRechitID;
	std::vector<int32_t> ele_EcalRechitIndex;  // flatten
	std::vector<int32_t> EcalRechitIndex;  // flatten
	std::vector<int32_t> ele_SeedRechitIndex;
	std::vector<int32_t> SeedRechitIndex;
	std::vector<uint> ecalRechitID_ToBeSaved;
	std::vector<std::pair<double,double>> ecalRechitEtaPhi_ToBeSaved;  // FIXME https://github.com/cms-lpc-llp/DelayedPhotonTuplizer/blob/master/plugins/RazorTuplizer.cc#L1687
	std::vector<std::pair<double,double>> ecalRechitJetEtaPhi_ToBeSaved;
	std::vector<int32_t> ecalRechitID;
	std::vector<float> ecalRechitEta;
	std::vector<float> ecalRechitPhi;
	std::vector<float> ecalRechitX;
	std::vector<float> ecalRechitY;
	std::vector<float> ecalRechitZ;
	std::vector<float> ecalRechitEnergy;
	std::vector<float> ecalRechitTime;
	std::vector<uint8_t> ecalRechitFlagOOT;
	std::vector<uint8_t> ecalRechitGainSwitch1;
	std::vector<uint8_t> ecalRechitGainSwitch6;
	std::vector<float> ecalRechitTranspCorr;
	std::vector<float> ecalRechitPedrms12;
	std::vector<float> ecalRechitPedrms6;
	std::vector<float> ecalRechitPedrms1;
	PhoVars vars_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(KNULLPProducer);
