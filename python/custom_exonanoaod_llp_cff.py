import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

myPhotonTable = cms.EDProducer(
    "KNULLPProducer",
    src = cms.InputTag("linkedObjects", "photons")
)

def add_customTables_template(process):

    process.outputTable = myPhotonTable
    process.exonanoaodTask = cms.Task(process.outputTable)
    process.nanoTableTaskCommon.add(process.exonanoaodTask)

    return process
