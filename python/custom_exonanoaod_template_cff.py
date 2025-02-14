import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *


outputTable = cms.EDProducer("EXOnanoAODProducerTemplate",
    inputExample = cms.InputTag("offlineBeamSpot")
)

def add_customTables_template(process):

    process.outputTable = outputTable
    process.exonanoaodTask = cms.Task(process.outputTable)
    process.nanoTableTaskCommon.add(process.exonanoaodTask)

## To be filled in once we converge on final EXO nanoAOD content

exonanoTable = cms.EDProducer("EXOnanoAODProducer",
    inputExample = cms.InputTag("input")
)

def add_exonanoTables(process):

    process.exonanoTable = exonanoTable
    process.exonanoTask = cms.Task(process.exonanoTable)
    process.nanoTableTaskCommon.add(process.exonanoTask)


    return process
