#ifndef PHOVARS_H
#define PHOVARS_H

#include <vector>
#include <cstdint>

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

struct PhoVars {
  // --- counters/bookkeeping ---
  std::vector<uint8_t> isStandardPhoton;  // you set false for OOT group in old code (if you want)
  std::vector<uint8_t> passEleVeto;
  std::vector<uint8_t> isConversion;
  std::vector<uint8_t> hasPixelSeed;

  // --- kinematics ---
  std::vector<float> E;
  std::vector<float> pt;
  std::vector<float> eta;
  std::vector<float> phi;
  std::vector<float> px;
  std::vector<float> py;
  std::vector<float> pz;

  // --- shower shapes / ID-like ---
  std::vector<float> sigmaIetaIeta;          // pho.see()
  std::vector<float> full5x5SigmaIetaIeta;   // pho.full5x5_sigmaIetaIeta()
  std::vector<float> r9;                     // pho.full5x5_r9()
  std::vector<float> hOverE;                 // pho.hadTowOverEm()

  std::vector<uint8_t> isOOT;

  // --- miniAOD isolation quantities ---
  std::vector<float> pfIsoChargedHadronIso;          // pho.chargedHadronIso()
  std::vector<float> pfIsoChargedHadronIsoWrongVtx;  // pho.chargedHadronIsoWrongVtx()
  std::vector<float> pfIsoNeutralHadronIso;          // pho.neutralHadronIso()
  std::vector<float> pfIsoPhotonIso;                 // pho.photonIso()
  std::vector<float> pfIsoModFrixione;               // pho.getPflowIsolationVariables().modFrixione
  std::vector<float> pfIsoSumPUPt;                   // pho.sumPUPt()

  // --- PF cluster / track iso from pat::Photon ---
  std::vector<float> ecalPFClusterIso;               // pho.ecalPFClusterIso()
  std::vector<float> hcalPFClusterIso;               // pho.hcalPFClusterIso()
  std::vector<float> trkSumPtHollowConeDR03;         // pho.trkSumPtHollowConeDR03()

  // --- regression energy ---
  std::vector<float> regressionE;                    // pho.getCorrectedEnergy( pho.getCandidateP4type() )
  std::vector<float> regressionEUnc;                 // pho.getCorrectedEnergyError( pho.getCandidateP4type() )

  // --- supercluster (through pho.superCluster()) ---
  std::vector<float> scEnergy;
  std::vector<float> scRawEnergy;
  std::vector<float> scEta;
  std::vector<float> scPhi;
  std::vector<float> scX;
  std::vector<float> scY;
  std::vector<float> scZ;
  std::vector<uint32_t> scSeedRawId;                 // pho.superCluster()->seed()->seed().rawId()

  // --- optional: cutbased ID + MVA (stored on pat::Photon as IDs/user data) ---
  std::vector<uint8_t> cutBasedID_loose;
  std::vector<uint8_t> cutBasedID_medium;
  std::vector<uint8_t> cutBasedID_tight;
  std::vector<float> mvaValue;                       // pho.userFloat(...)
  std::vector<int32_t> mvaCategory;                  // pho.userInt(...)

  explicit PhoVars(size_t n)
      : isStandardPhoton(n, 1),
        passEleVeto(n, 0),
        isConversion(n, 0),
        hasPixelSeed(n, 0),

        E(n, -999.f),
        pt(n, -999.f),
        eta(n, -999.f),
        phi(n, -999.f),
        px(n, -999.f),
        py(n, -999.f),
        pz(n, -999.f),

        sigmaIetaIeta(n, -999.f),
        full5x5SigmaIetaIeta(n, -999.f),
        r9(n, -999.f),
        hOverE(n, -999.f),

	isOOT(n, 0),

        pfIsoChargedHadronIso(n, -999.f),
        pfIsoChargedHadronIsoWrongVtx(n, -999.f),
        pfIsoNeutralHadronIso(n, -999.f),
        pfIsoPhotonIso(n, -999.f),
        pfIsoModFrixione(n, -999.f),
        pfIsoSumPUPt(n, -999.f),

        ecalPFClusterIso(n, -999.f),
        hcalPFClusterIso(n, -999.f),
        trkSumPtHollowConeDR03(n, -999.f),

        regressionE(n, -999.f),
        regressionEUnc(n, -999.f),

        scEnergy(n, -999.f),
        scRawEnergy(n, -999.f),
        scEta(n, -999.f),
        scPhi(n, -999.f),
        scX(n, -999.f),
        scY(n, -999.f),
        scZ(n, -999.f),
        scSeedRawId(n, 0u),

        cutBasedID_loose(n, 0),
        cutBasedID_medium(n, 0),
        cutBasedID_tight(n, 0),
        mvaValue(n, -999.f),
        mvaCategory(n, -999) {}

  inline void fillFromPho(const pat::Photon& pho, size_t i) {
    // kinematics
    E[i]   = pho.energy();
    pt[i]  = pho.pt();
    eta[i] = pho.eta();
    phi[i] = pho.phi();
    px[i]  = pho.px();
    py[i]  = pho.py();
    pz[i]  = pho.pz();

    // shapes / flags
    sigmaIetaIeta[i]        = pho.see();
    full5x5SigmaIetaIeta[i] = pho.full5x5_sigmaIetaIeta();
    r9[i]                   = pho.full5x5_r9();
    hOverE[i]               = pho.hadTowOverEm();

    isConversion[i] = pho.hasConversionTracks();
    passEleVeto[i]  = pho.passElectronVeto();
    hasPixelSeed[i] = pho.hasPixelSeed();

    // miniAOD iso quantities
    pfIsoChargedHadronIso[i]         = pho.chargedHadronIso();
    //pfIsoChargedHadronIsoWrongVtx[i] = pho.chargedHadronIsoWrongVtx();  // FIXME
    pfIsoNeutralHadronIso[i]         = pho.neutralHadronIso();
    pfIsoPhotonIso[i]                = pho.photonIso();
    //pfIsoModFrixione[i]              = pho.getPflowIsolationVariables().modFrixione;  // FIXME
    //pfIsoSumPUPt[i]                  = pho.sumPUPt();  // FIXME

    // pat::Photon PF cluster / trk iso (if available in your release)
    ecalPFClusterIso[i]       = pho.ecalPFClusterIso();
    hcalPFClusterIso[i]       = pho.hcalPFClusterIso();
    trkSumPtHollowConeDR03[i] = pho.trkSumPtHollowConeDR03();

    // regression energy
    regressionE[i]    = pho.getCorrectedEnergy( pho.getCandidateP4type() );
    regressionEUnc[i] = pho.getCorrectedEnergyError( pho.getCandidateP4type() );

    // supercluster
    if (pho.superCluster().isNonnull()) {
      const auto& sc = *pho.superCluster();
      scEnergy[i]    = sc.energy();
      scRawEnergy[i] = sc.rawEnergy();
      scEta[i]       = sc.eta();
      scPhi[i]       = sc.phi();
      scX[i]         = sc.x();
      scY[i]         = sc.y();
      scZ[i]         = sc.z();

      // seed rawId (guard all pointers)
      if (sc.seed().isNonnull()) {
        scSeedRawId[i] = sc.seed()->seed().rawId();
      }
    }
  }
};

// ---- NanoAOD column booking ----

inline void addPhoColumns(nanoaod::FlatTable& tab, const PhoVars& v) {
  // --- flags / bookkeeping ---
  tab.addColumn<uint8_t>("isStandardPhoton", v.isStandardPhoton, "1 if in-time (standard) photon", -1);
  tab.addColumn<uint8_t>("passEleVeto",      v.passEleVeto,      "1 if passElectronVeto()", -1);
  tab.addColumn<uint8_t>("isConversion",     v.isConversion,     "1 if hasConversionTracks()", -1);
  tab.addColumn<uint8_t>("hasPixelSeed",     v.hasPixelSeed,     "1 if hasPixelSeed()", -1);

  // --- kinematics ---
  tab.addColumn<float>("E",   v.E,   "energy()", 10);
  tab.addColumn<float>("pt",  v.pt,  "pT", 10);
  tab.addColumn<float>("eta", v.eta, "eta", 10);
  tab.addColumn<float>("phi", v.phi, "phi", 10);
  tab.addColumn<float>("px",  v.px,  "px()", 10);
  tab.addColumn<float>("py",  v.py,  "py()", 10);
  tab.addColumn<float>("pz",  v.pz,  "pz()", 10);

  // --- shower shapes / ID-like ---
  tab.addColumn<float>("sigmaIetaIeta",        v.sigmaIetaIeta,        "see()", 10);
  tab.addColumn<float>("full5x5SigmaIetaIeta", v.full5x5SigmaIetaIeta, "full5x5_sigmaIetaIeta()", 10);
  tab.addColumn<float>("r9",                   v.r9,                   "full5x5_r9()", 10);
  tab.addColumn<float>("hOverE",               v.hOverE,               "hadTowOverEm()", 10);

  tab.addColumn<uint8_t>("isOOT", v.isOOT, "1 if OOT", -1);

  // --- miniAOD isolation quantities ---
  tab.addColumn<float>("pfIsoChargedHadronIso",         v.pfIsoChargedHadronIso,         "chargedHadronIso()", 10);
  tab.addColumn<float>("pfIsoChargedHadronIsoWrongVtx", v.pfIsoChargedHadronIsoWrongVtx, "chargedHadronIsoWrongVtx()", 10);
  tab.addColumn<float>("pfIsoNeutralHadronIso",         v.pfIsoNeutralHadronIso,         "neutralHadronIso()", 10);
  tab.addColumn<float>("pfIsoPhotonIso",                v.pfIsoPhotonIso,                "photonIso()", 10);
  tab.addColumn<float>("pfIsoModFrixione",              v.pfIsoModFrixione,              "getPflowIsolationVariables().modFrixione", 10);
  tab.addColumn<float>("pfIsoSumPUPt",                  v.pfIsoSumPUPt,                  "sumPUPt()", 10);

  // --- pat::Photon PF cluster / track iso (as used in your fillPhotons) ---
  tab.addColumn<float>("ecalPFClusterIso",       v.ecalPFClusterIso,       "ecalPFClusterIso()", 10);
  tab.addColumn<float>("hcalPFClusterIso",       v.hcalPFClusterIso,       "hcalPFClusterIso()", 10);
  tab.addColumn<float>("trkSumPtHollowConeDR03", v.trkSumPtHollowConeDR03, "trkSumPtHollowConeDR03()", 10);

  // --- regression energy ---
  tab.addColumn<float>("regressionE",    v.regressionE,    "getCorrectedEnergy(getCandidateP4type())", 10);
  tab.addColumn<float>("regressionEUnc", v.regressionEUnc, "getCorrectedEnergyError(getCandidateP4type())", 10);

  // --- supercluster ---
  tab.addColumn<float>("superClusterEnergy",    v.scEnergy,    "superCluster()->energy()", 10);
  tab.addColumn<float>("superClusterRawEnergy", v.scRawEnergy, "superCluster()->rawEnergy()", 10);
  tab.addColumn<float>("superClusterEta",       v.scEta,       "superCluster()->eta()", 10);
  tab.addColumn<float>("superClusterPhi",       v.scPhi,       "superCluster()->phi()", 10);
  tab.addColumn<float>("superClusterX",         v.scX,         "superCluster()->x()", 10);
  tab.addColumn<float>("superClusterY",         v.scY,         "superCluster()->y()", 10);
  tab.addColumn<float>("superClusterZ",         v.scZ,         "superCluster()->z()", 10);
  tab.addColumn<uint32_t>("superClusterSeedRawId", v.scSeedRawId, "superCluster()->seed()->seed().rawId()", -1);

  // --- optional IDs / MVA (kept even if defaults when missing) ---
  tab.addColumn<uint8_t>("cutBasedID_loose",  v.cutBasedID_loose,  "cutBased loose ID (if present)", -1);
  tab.addColumn<uint8_t>("cutBasedID_medium", v.cutBasedID_medium, "cutBased medium ID (if present)", -1);
  tab.addColumn<uint8_t>("cutBasedID_tight",  v.cutBasedID_tight,  "cutBased tight ID (if present)", -1);
  tab.addColumn<float>("mvaValue",            v.mvaValue,          "userFloat PhotonMVAEstimator... (if present)", 10);
  tab.addColumn<int32_t>("mvaCategory",       v.mvaCategory,       "userInt PhotonMVAEstimator... (if present)", -1);
}

#endif
