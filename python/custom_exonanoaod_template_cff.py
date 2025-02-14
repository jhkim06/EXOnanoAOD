import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

outputTable = cms.EDProducer("EXOnanoAODProducerTemplate",
    inputExample = cms.InputTag("offlineBeamSpot")
)

def add_customTables_template(process):

    process.outputTable = outputTable
    process.exonanoaodTask = cms.Task(process.outputTable)
    process.nanoTableTaskCommon.add(process.exonanoaodTask)

    return process
