import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *


myPhotonTable = cms.EDProducer(
    "KNULLPProducer",
    src = cms.InputTag("linkedObjects", "photons"),
    photons = cms.VInputTag(
        cms.InputTag("slimmedPhotons"),
        cms.InputTag("slimmedOOTPhotons"),
    )
)

def add_customTables_template(process):

    ##------ Electron MVA Setup ------#
    ## define which IDs we want to produce
    #dataFormat = DataFormat.MiniAOD
    #switchOnVIDPhotonIdProducer(process, dataFormat)
    #my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_RunIIIWinter22_122X_V1_cff', 'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Winter22_122X_V1_cff']

    ##add them to the VID producer
    #for idmod in my_id_modules:
    #    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

    #process.nanoAOD_step.associate(process.egmPhotonIDTask)

    process.outputTable = myPhotonTable
    process.exonanoaodTask = cms.Task(process.outputTable)
    process.nanoTableTaskCommon.add(process.exonanoaodTask)

    return process
