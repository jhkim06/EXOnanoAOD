import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *


myPhotonTable = cms.EDProducer(
    "KNULLPProducer",
    src = cms.InputTag("linkedObjects", "photons"),

    photons = cms.VInputTag(
        cms.InputTag("slimmedPhotons"),
        cms.InputTag("slimmedOOTPhotons"),
    ),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    packedPfCands = cms.InputTag("packedPFCandidates"),
    tracks = cms.InputTag("unpackedTracksAndVertices"),
    conversions = cms.InputTag("reducedEgamma", "reducedConversions"),  # "PAT" used in Razor
    singleLegConversions = cms.InputTag("reducedEgamma", "reducedSingleLegConversions"),  # "PAT" used in Razor
    beamSpot = cms.InputTag("offlineBeamSpot"),
)

def add_customTables_template(process):

    process.outputTable = myPhotonTable
    process.exonanoaodTask = cms.Task(process.outputTable)
    process.nanoTableTaskCommon.add(process.exonanoaodTask)

    return process
