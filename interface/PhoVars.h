#ifndef PHOVARS_H
#define PHOVARS_H

#include <vector>
#include <cstdint>

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

// FIXME may be RazorVars better
// photon table
// photon related table (different size)
struct PhoVars {
  // --- counters/bookkeeping ---
  std::vector<uint8_t> isStandardPhoton;
  std::vector<uint8_t> passEleVeto;
  std::vector<uint8_t> isConversion;
  std::vector<uint8_t> hasPixelSeed;

  std::vector<uint8_t> trackMatching;

  // --- kinematics ---
  std::vector<float> E;
  std::vector<float> pt;
  std::vector<float> eta;
  std::vector<float> phi;

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
  //std::vector<float> sumChargedHadronPtAllVertices;
  std::vector<float> sumChargedHadronPt;               // 
  std::vector<float> sumNeutralHadronEt;               // 
  std::vector<float> sumPhotonEt;               // 
  std::vector<float> ecalPFClusterIso;               // pho.ecalPFClusterIso()
  std::vector<float> hcalPFClusterIso;               // pho.hcalPFClusterIso()
  std::vector<float> trkSumPtHollowConeDR03;         // pho.trkSumPtHollowConeDR03()
  std::vector<float> sumWorstVertexChargedHadronPt;               // 

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

  // --- optional: cutbased ID + MVA (stored on pat::Photon as IDs/user data) ---
  std::vector<uint8_t> cutBasedID_loose;
  std::vector<uint8_t> cutBasedID_medium;
  std::vector<uint8_t> cutBasedID_tight;
  std::vector<float> mvaValue;                       // pho.userFloat(...)
  std::vector<int32_t> mvaCategory;                  // pho.userInt(...)

  std::vector<float> energy_scale;
  std::vector<float> energy_scale_up;
  std::vector<float> energy_scale_down;
  std::vector<float> energy_smear;

  std::vector<int32_t> convType;                  // 
  std::vector<float> convTrkZ;
  std::vector<float> convTrkClusZ;

  PhoVars() = default;

  explicit PhoVars(size_t n) {
    resize(n);
  }

  void resize(size_t n) {
    // bool / uint8-like flags
    isStandardPhoton.assign(n, 1);
    passEleVeto.assign(n, 0);
    isConversion.assign(n, 0);
    hasPixelSeed.assign(n, 0);

    trackMatching.assign(n, 0);

    // kinematics
    E.assign(n, -999.f);
    pt.assign(n, -999.f);
    eta.assign(n, -999.f);
    phi.assign(n, -999.f);

    // shapes / ID-like floats
    sigmaIetaIeta.assign(n, -999.f);
    full5x5SigmaIetaIeta.assign(n, -999.f);
    r9.assign(n, -999.f);
    hOverE.assign(n, -999.f);

    // oot flag
    isOOT.assign(n, 0);

    // miniAOD PF isolations
    pfIsoChargedHadronIso.assign(n, -999.f);
    pfIsoChargedHadronIsoWrongVtx.assign(n, -999.f);
    pfIsoNeutralHadronIso.assign(n, -999.f);
    pfIsoPhotonIso.assign(n, -999.f);
    pfIsoModFrixione.assign(n, -999.f);
    pfIsoSumPUPt.assign(n, -999.f);

    // 2D (flattened) optional
    // sumChargedHadronPtAllVertices.assign(n * npv, -999.f);

    // other isolations
    sumChargedHadronPt.assign(n, -999.f);
    sumNeutralHadronEt.assign(n, -999.f);
    sumPhotonEt.assign(n, -999.f);
    ecalPFClusterIso.assign(n, -999.f);
    hcalPFClusterIso.assign(n, -999.f);
    trkSumPtHollowConeDR03.assign(n, -999.f);
    sumWorstVertexChargedHadronPt.assign(n, -999.f);

    // regression
    regressionE.assign(n, -999.f);
    regressionEUnc.assign(n, -999.f);

    // supercluster
    scEnergy.assign(n, -999.f);
    scRawEnergy.assign(n, -999.f);
    scEta.assign(n, -999.f);
    scPhi.assign(n, -999.f);
    scX.assign(n, -999.f);
    scY.assign(n, -999.f);
    scZ.assign(n, -999.f);

    // VID / MVA
    cutBasedID_loose.assign(n, 0);
    cutBasedID_medium.assign(n, 0);
    cutBasedID_tight.assign(n, 0);
    mvaValue.assign(n, -999.f);
    mvaCategory.assign(n, -999);

    //  scales/smears
    energy_scale.assign(n, -999.f);
    energy_scale_up.assign(n, -999.f);
    energy_scale_down.assign(n, -999.f);
    energy_smear.assign(n, -999.f);

    convType.assign(n, -99);
    convTrkZ.assign(n, -99.f);
    convTrkClusZ.assign(n, -99.f);
  }

  void clear() {

    // bool / uint8-like flags
    isStandardPhoton.clear();
    passEleVeto.clear();
    isConversion.clear();
    hasPixelSeed.clear();

    trackMatching.clear();

    // kinematics
    E.clear();
    pt.clear();
    eta.clear();
    phi.clear();

    // shapes / ID-like floats
    sigmaIetaIeta.clear();
    full5x5SigmaIetaIeta.clear();
    r9.clear();
    hOverE.clear();

    // oot flag
    isOOT.clear();

    // miniAOD PF isolations
    pfIsoChargedHadronIso.clear();
    pfIsoChargedHadronIsoWrongVtx.clear();
    pfIsoNeutralHadronIso.clear();
    pfIsoPhotonIso.clear();
    pfIsoModFrixione.clear();
    pfIsoSumPUPt.clear();

    // optional flattened 2D
    // sumChargedHadronPtAllVertices.clear();

    // other isolations
    sumChargedHadronPt.clear();
    sumNeutralHadronEt.clear();
    sumPhotonEt.clear();
    ecalPFClusterIso.clear();
    hcalPFClusterIso.clear();
    trkSumPtHollowConeDR03.clear();
    sumWorstVertexChargedHadronPt.clear();

    // regression
    regressionE.clear();
    regressionEUnc.clear();

    // supercluster
    scEnergy.clear();
    scRawEnergy.clear();
    scEta.clear();
    scPhi.clear();
    scX.clear();
    scY.clear();
    scZ.clear();

    // VID / MVA
    cutBasedID_loose.clear();
    cutBasedID_medium.clear();
    cutBasedID_tight.clear();
    mvaValue.clear();
    mvaCategory.clear();

    // scales/smears
    energy_scale.clear();
    energy_scale_up.clear();
    energy_scale_down.clear();
    energy_smear.clear();

    convType.clear();
    convTrkZ.clear();
    convTrkClusZ.clear();
  }

  inline void fillFromPho(const pat::Photon& pho, size_t i) {
    // kinematics
    E[i]   = pho.energy();
    pt[i]  = pho.pt();
    eta[i] = pho.eta();
    phi[i] = pho.phi();

    // shapes / flags
    sigmaIetaIeta[i]        = pho.see();  // FIXME this is always 0: slimmedPhotons_cfi.py:    dropRegressionData = cms.string("1"),
    full5x5SigmaIetaIeta[i] = pho.full5x5_sigmaIetaIeta();
    r9[i]                   = pho.full5x5_r9();
    hOverE[i]               = pho.hadTowOverEm();

    isConversion[i] = pho.hasConversionTracks();
    passEleVeto[i]  = pho.passElectronVeto();
    hasPixelSeed[i] = pho.hasPixelSeed();

    // miniAOD iso quantities
    pfIsoChargedHadronIso[i]         = pho.chargedHadronIso();
    //pfIsoChargedHadronIsoWrongVtx[i] = pho.chargedHadronIsoWrongVtx();  // FIXME deprecated
    pfIsoNeutralHadronIso[i]         = pho.neutralHadronIso();
    pfIsoPhotonIso[i]                = pho.photonIso();
    //pfIsoModFrixione[i]              = pho.getPflowIsolationVariables().modFrixione;  // FIXME deprecated
    //pfIsoSumPUPt[i]                  = pho.sumPUPt();  // FIXME deprecated

    // pat::Photon PF cluster / trk iso (if available in your release)
    ecalPFClusterIso[i]       = pho.ecalPFClusterIso();
    hcalPFClusterIso[i]       = pho.hcalPFClusterIso();
    trkSumPtHollowConeDR03[i] = pho.trkSumPtHollowConeDR03();

    // regression energy
    regressionE[i]    = pho.getCorrectedEnergy( pho.getCandidateP4type() );
    regressionEUnc[i] = pho.getCorrectedEnergyError( pho.getCandidateP4type() );

    //std::cout << "photon ID: " << std::endl;
    //for (auto const& n : pho.userFloatNames()) edm::LogPrint("KNULLP") << "  " << n;
    //for (auto const& n : pho.userIntNames()) edm::LogPrint("KNULLP") << "  " << n;

    try {// as saved in slimmedPhoton
      cutBasedID_loose[i] = pho.photonID("cutBasedPhotonID-RunIIIWinter22-122X-V1-loose");  //
      cutBasedID_medium[i] = pho.photonID("cutBasedPhotonID-RunIIIWinter22-122X-V1-medium");
      cutBasedID_tight[i] = pho.photonID("cutBasedPhotonID-RunIIIWinter22-122X-V1-tight");

      mvaValue[i] = pho.userFloat("PhotonMVAEstimatorRunIIIWinter22v1Values");
      mvaCategory[i] = pho.userInt("PhotonMVAEstimatorRunIIIWinter22v1Categories");
    } catch (...) {
      std::cout << "No Photon ID / MVA found." << std::endl;
    }

    try {
      // FIXME https://egammapog.docs.cern.ch/Run3/SaS/
      // or Photon_energyErr?
      energy_scale[i] = pho.userFloat("energyScaleValue");
      energy_scale_up[i] = pho.userFloat("energyScaleUp");
      energy_scale_down[i] = pho.userFloat("energyScaleDown");
      energy_smear[i] = pho.userFloat("energySigmaValue");
    } catch (...) {
      std::cout << "No Photon scale found. Set it to 1." << std::endl;
      energy_scale[i] = 1.0;
      energy_scale_up[i] = 1.0;
      energy_scale_down[i] = 1.0;
      energy_smear[i] = 1.0;
    }

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
    }
    /*
    TODO
    float pho_vtxSumPx[OBJECTARRAYSIZE][MAX_NPV]; CALCULATE as in Razor, Note 2D vector not supported in FlatTable
    float pho_vtxSumPy[OBJECTARRAYSIZE][MAX_NPV]; CALCULATE as in Razor, Note 2D vector not supported in FlatTable

    bool pho_passHLTFilter[OBJECTARRAYSIZE][MAX_PhotonHLTFilters];  Note 2D vector not supported in FlatTable
    */
  }
};

// ---- NanoAOD column booking ----
inline void addPhoColumns(nanoaod::FlatTable& tab, const PhoVars& v) {

  // --- flags / bookkeeping ---
  tab.addColumn<uint8_t>("isStandardPhoton", v.isStandardPhoton, "1 if in-time (standard) photon", -1);
  tab.addColumn<uint8_t>("passEleVeto",      v.passEleVeto,      "1 if passElectronVeto()", -1);
  tab.addColumn<uint8_t>("isConversion",     v.isConversion,     "1 if hasConversionTracks()", -1);
  tab.addColumn<uint8_t>("hasPixelSeed",     v.hasPixelSeed,     "1 if hasPixelSeed()", -1);

  tab.addColumn<uint8_t>("trackMatching", v.trackMatching, "1 there is close track", -1);

  // --- kinematics ---
  tab.addColumn<float>("E",   v.E,   "energy()", 10);
  tab.addColumn<float>("pt",  v.pt,  "pT", 10);
  tab.addColumn<float>("eta", v.eta, "eta", 10);
  tab.addColumn<float>("phi", v.phi, "phi", 10);

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
  //tab.addColumn<float>("sumChargedHadronPtAllVertices",       v.sumChargedHadronPtAllVertices,       "", 10);  // 2d array not supported  in FlatTable
  tab.addColumn<float>("sumChargedHadronPt",       v.sumChargedHadronPt,       "", 10);
  tab.addColumn<float>("sumNeutralHadronEt",       v.sumNeutralHadronEt,       "", 10);
  tab.addColumn<float>("sumPhotonEt",       v.sumPhotonEt,       "", 10);
  tab.addColumn<float>("ecalPFClusterIso",       v.ecalPFClusterIso,       "ecalPFClusterIso()", 10);
  tab.addColumn<float>("hcalPFClusterIso",       v.hcalPFClusterIso,       "hcalPFClusterIso()", 10);
  tab.addColumn<float>("trkSumPtHollowConeDR03", v.trkSumPtHollowConeDR03, "trkSumPtHollowConeDR03()", 10);
  tab.addColumn<float>("sumWorstVertexChargedHadronPt",       v.sumWorstVertexChargedHadronPt,       "", 10);

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

  // --- optional IDs / MVA (kept even if defaults when missing) ---
  tab.addColumn<uint8_t>("cutBasedID_loose",  v.cutBasedID_loose,  "cutBased loose ID (if present)", -1);
  tab.addColumn<uint8_t>("cutBasedID_medium", v.cutBasedID_medium, "cutBased medium ID (if present)", -1);
  tab.addColumn<uint8_t>("cutBasedID_tight",  v.cutBasedID_tight,  "cutBased tight ID (if present)", -1);
  tab.addColumn<float>("mvaValue",            v.mvaValue,          "userFloat PhotonMVAEstimator... (if present)", 10);
  tab.addColumn<int32_t>("mvaCategory",       v.mvaCategory,       "userInt PhotonMVAEstimator... (if present)", -1);

  tab.addColumn<float>("energy_scale",       v.energy_scale,       "userFloat energy_scale... (if present)", 10);
  tab.addColumn<float>("energy_scale_up",       v.energy_scale_up,       "userFloat energy_scale_up... (if present)", 10);
  tab.addColumn<float>("energy_scale_down",       v.energy_scale_down,       "userFloat energy_scale_down... (if present)", 10);
  tab.addColumn<float>("energy_smear",       v.energy_smear,       "userFloat energy_smear... (if present)", 10);

  tab.addColumn<int32_t>("convType",  v.convType,       "", -1);
  tab.addColumn<float>("convTrkZ",    v.convTrkZ,    "", 10);
  tab.addColumn<float>("convTrkClusZ", v.convTrkClusZ,    "", 10);
}

#endif
