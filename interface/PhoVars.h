#ifndef PHOVARS_H
#define PHOVARS_H

#include <vector>
#include <cstdint>

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

struct PhoVars {
  // --- kinematics ---
  std::vector<float> pt;
  std::vector<float> eta;
  std::vector<float> phi;
  std::vector<float> energy;

  // --- supercluster ---
  std::vector<float> scEnergy;
  std::vector<float> scEta;
  std::vector<float> scPhi;
  std::vector<float> scRawEnergy;
  std::vector<float> scEtaWidth;
  std::vector<float> scPhiWidth;

  // --- shower shapes (full5x5) ---
  std::vector<float> r9;
  std::vector<float> sigmaIetaIeta;
  std::vector<float> sigmaIetaIphi;
  std::vector<float> sigmaIphiIphi;
  std::vector<float> e1x5;
  std::vector<float> e2x5;
  std::vector<float> e3x3;
  std::vector<float> e5x5;
  std::vector<float> sMajor;
  std::vector<float> sMinor;

  // --- H/E ---
  std::vector<float> hOverE;

  // --- isolation (PAT-level) ---
  std::vector<float> ecalIso;
  std::vector<float> hcalIso;
  std::vector<float> trkIsoHollow;
  std::vector<float> trkIsoSolid;

  // --- flags ---
  std::vector<uint8_t> isEB;
  std::vector<uint8_t> isEE;
  std::vector<uint8_t> passEleVeto;
  std::vector<uint8_t> hasConversion;
  std::vector<uint8_t> isOOT;

  explicit PhoVars(size_t n)
      : pt(n, -999.f),
        eta(n, -999.f),
        phi(n, -999.f),
        energy(n, -999.f),
        scEnergy(n, -999.f),
        scEta(n, -999.f),
        scPhi(n, -999.f),
        scRawEnergy(n, -999.f),
        scEtaWidth(n, -999.f),
        scPhiWidth(n, -999.f),
        r9(n, -999.f),
        sigmaIetaIeta(n, -999.f),
        sigmaIetaIphi(n, -999.f),
        sigmaIphiIphi(n, -999.f),
        e1x5(n, -999.f),
        e2x5(n, -999.f),
        e3x3(n, -999.f),
        e5x5(n, -999.f),
        sMajor(n, -999.f),
        sMinor(n, -999.f),
        hOverE(n, -999.f),
        ecalIso(n, -999.f),
        hcalIso(n, -999.f),
        trkIsoHollow(n, -999.f),
        trkIsoSolid(n, -999.f),
        isEB(n, 0),
        isEE(n, 0),
        passEleVeto(n, 0),
        hasConversion(n, 0),
        isOOT(n, 0) {}
};

// -------------------- helpers --------------------

inline void fillPhoVars(const pat::Photon& p,
                        PhoVars& v,
                        size_t i) {
  // kinematics
  v.pt[i]     = p.pt();
  v.eta[i]    = p.eta();
  v.phi[i]    = p.phi();
  v.energy[i] = p.energy();

  // supercluster
  if (p.superCluster().isNonnull()) {
    const auto& sc = *p.superCluster();
    v.scEnergy[i]    = sc.energy();
    v.scEta[i]       = sc.eta();
    v.scPhi[i]       = sc.phi();
    v.scRawEnergy[i] = sc.rawEnergy();
    v.scEtaWidth[i]  = sc.etaWidth();
    v.scPhiWidth[i]  = sc.phiWidth();
  }

  // shower shapes
  const auto& ss = p.full5x5_showerShapeVariables();
//  v.r9[i]               = ss.r9;
//  v.sigmaIetaIeta[i]    = ss.sigmaIetaIeta;
//  v.sigmaIetaIphi[i]    = ss.sigmaIetaIphi;
//  v.sigmaIphiIphi[i]    = ss.sigmaIphiIphi;
//  v.e1x5[i]             = ss.e1x5;
//  v.e2x5[i]             = ss.e2x5;
//  v.e3x3[i]             = ss.e3x3;
//  v.e5x5[i]             = ss.e5x5;
  v.sMajor[i]           = ss.smMajor;
  v.sMinor[i]           = ss.smMinor;

  // H/E
  v.hOverE[i] = p.hadronicOverEm();

  // isolation
  v.ecalIso[i]      = p.ecalRecHitSumEtConeDR04();
  v.hcalIso[i]      = p.hcalTowerSumEtConeDR04();
  v.trkIsoHollow[i] = p.trkSumPtHollowConeDR04();
  v.trkIsoSolid[i]  = p.trkSumPtSolidConeDR04();

  // flags
  v.isEB[i]          = p.isEB();
  v.isEE[i]          = p.isEE();
  v.passEleVeto[i]   = p.passElectronVeto();
  v.hasConversion[i] = p.hasConversionTracks();
}

// Book columns into a FlatTable (kept outside the struct = cleaner separation)
inline void addPhoColumns(nanoaod::FlatTable& tab, const PhoVars& v) {
  tab.addColumn<float>("pt", v.pt, "pT", 10);
  tab.addColumn<float>("eta", v.eta, "eta", 10);
  tab.addColumn<float>("phi", v.phi, "phi", 10);
  tab.addColumn<float>("energy", v.energy, "energy", 10);

  tab.addColumn<float>("superclusterEnergy", v.scEnergy, "supercluster energy", 10);
  tab.addColumn<float>("superclusterEta", v.scEta, "supercluster eta", 10);
  tab.addColumn<float>("superclusterPhi", v.scPhi, "supercluster phi", 10);
  tab.addColumn<float>("superclusterRawEnergy", v.scRawEnergy, "supercluster raw energy", 10);
  tab.addColumn<float>("superclusterEtaWidth", v.scEtaWidth, "supercluster eta width", 10);
  tab.addColumn<float>("superclusterPhiWidth", v.scPhiWidth, "supercluster phi width", 10);

  //tab.addColumn<float>("r9", v.r9, "R9 (full5x5)", 10);
  //tab.addColumn<float>("sieie", v.sigmaIetaIeta, "sigmaIetaIeta (full5x5)", 10);
  //tab.addColumn<float>("sieip", v.sigmaIetaIphi, "sigmaIetaIphi (full5x5)", 10);
  //tab.addColumn<float>("sipip", v.sigmaIphiIphi, "sigmaIphiIphi (full5x5)", 10);
  //tab.addColumn<float>("e1x5", v.e1x5, "E1x5 (full5x5)", 10);
  //tab.addColumn<float>("e2x5", v.e2x5, "E2x5 (full5x5)", 10);
  //tab.addColumn<float>("e3x3", v.e3x3, "E3x3 (full5x5)", 10);
  //tab.addColumn<float>("e5x5", v.e5x5, "E5x5 (full5x5)", 10);
  tab.addColumn<float>("sMajor", v.sMajor, "shower shape sMajor (full5x5)", 10);
  tab.addColumn<float>("sMinor", v.sMinor, "shower shape sMinor (full5x5)", 10);

  tab.addColumn<float>("hOverE", v.hOverE, "hadronic over EM", 10);

  tab.addColumn<float>("ecalIso", v.ecalIso, "ecalRecHitSumEtConeDR04", 10);
  tab.addColumn<float>("hcalIso", v.hcalIso, "hcalTowerSumEtConeDR04", 10);
  tab.addColumn<float>("trkIsoHollow", v.trkIsoHollow, "trkSumPtHollowConeDR04", 10);
  tab.addColumn<float>("trkIsoSolid", v.trkIsoSolid, "trkSumPtSolidConeDR04", 10);

  tab.addColumn<uint8_t>("isEB", v.isEB, "1 if EB", -1);
  tab.addColumn<uint8_t>("isEE", v.isEE, "1 if EE", -1);
  tab.addColumn<uint8_t>("passEleVeto", v.passEleVeto, "1 if pass electron veto", -1);
  tab.addColumn<uint8_t>("hasConversion", v.hasConversion, "1 if has conversion tracks", -1);

  tab.addColumn<uint8_t>("isOOT", v.isOOT, "1 if OOT", -1);
}

#endif
