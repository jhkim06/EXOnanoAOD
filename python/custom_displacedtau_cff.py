import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.nano_eras_cff import *
from PhysicsTools.NanoAOD.genparticles_cff import *
from PhysicsTools.NanoAOD.taus_cff import *
from PhysicsTools.NanoAOD.muons_cff import *
from PhysicsTools.NanoAOD.jetsAK4_CHS_cff import *
from PhysicsTools.NanoAOD.jetsAK4_Puppi_cff import *
from PhysicsTools.NanoAOD.custom_btv_cff import *
from Configuration.Eras.Modifier_run3_common_cff import run3_common
import os

def add_displacedtauTables(process, isMC, useCHSJets = True):

    if isMC:
        # GenParticles                                                                                                                            
        genParticleTable.variables.vertexX        = Var("vertex.X"      , float)
        genParticleTable.variables.vertexY        = Var("vertex.Y"      , float)
        genParticleTable.variables.vertexZ        = Var("vertex.Z"      , float)
        genParticleTable.variables.vertexRho      = Var("vertex.Rho"    , float)
        genParticleTable.variables.vertexR        = Var("vertex.R"      , float)


    file = "NanoProd/data/particlenet_v1_a27159734e304ea4b7f9e0042baa9e22.pb"
    if os.path.exists(file):
        file_string = file
    elif os.path.exists( os.path.basename(file) ):
        file_string = os.path.basename(file)
    else:
        #file_string = "data/particlenet_v1_a27159734e304ea4b7f9e0042baa9e22.pb"
        file_string = "/afs/cern.ch/work/p/ppalit2/public/tau_pog_reco/exonanoaod_v2/CMSSW_15_0_0_pre3/src/PhysicsTools/EXOnanoAOD/data/particlenet_v1_a27159734e304ea4b7f9e0042baa9e22.pb"

    process.disTauTag = cms.EDProducer(
          "DisTauTag",
        graphPath = cms.string(file_string),
        jets = process.jetTable.src,
        pfCandidates = cms.InputTag('packedPFCandidates'),
        save_inputs  = cms.bool(False)
    )
    
    d_disTauTagVars = {
          "disTauTag_score0":     ExtVar("disTauTag:score0"       , float, doc = "Score 0"),
          "disTauTag_score1":     ExtVar("disTauTag:score1"       , float, doc = "Score 1"),
    }

    # Create the task                                                                                                                             
    print ('adding disTau edproducer')
    if useCHSJets:
      process.jetTable.externalVariables = process.jetTable.externalVariables.clone(**d_disTauTagVars)
    ## for puppi jets, use this!                                                                                                                  
    else:
      process.jetPuppiTable.externalVariables = process.jetPuppiTable.externalVariables.clone(**d_disTauTagVars)

    ##//process.custom_nanoaod_task = cms.Task(  
    process.exonanoaodTask = cms.Task(
        process.disTauTag,
    )

    ##//process.schedule.associate(process.custom_nanoaod_task)
    process.nanoTableTaskCommon.add(process.exonanoaodTask)
    
    return process


def BTVCustomNanoAODStaus(process, useCHSJets = True):
    addPFCands(process,btvNano_switch.btvNano_addallPF_switch,btvNano_switch.btvNano_addAK4_switch,btvNano_switch.btvNano_addAK8_switch)
    
    ### for MC   
    process.load("PhysicsTools.NanoAOD.btvMC_cff")
    process.nanoSequenceMC+=ak4onlyPFCandsMCSequence

    if useCHSJets:
        process.finalJetsAK4Constituents.src = src = cms.InputTag("finalJets")
        process.customAK4ConstituentsTable.jets = cms.InputTag("finalJets")
    
    return process


def run3modifier_chs(process, useCHSJets = True):
    if useCHSJets:
        run3_common.toModify(
            process.linkedObjects, jets="finalJets" ## run 2                                                                                          
        )
        _nanoTableTaskCommonRun3 = process.nanoTableTaskCommon.copy()
        _nanoTableTaskCommonRun3.add(process.jetTask)
        _nanoTableTaskCommonRun3.add(process.jetForMETTask)
        ## remove puppi table otherwise it tries to save the jet table twice                                                                          
        process.jetPuppiForMETTask.remove(process.corrT1METJetPuppiTable)
        
        process.jetTablesTask.remove(process.bjetNN)
        process.jetTablesTask.remove(process.cjetNN)
        _nanoTableTaskCommonRun3.replace(process.jetPuppiTablesTask, process.jetTablesTask)
        
        process.jetTable.externalVariables = cms.PSet()
        run3_common.toReplaceWith(
            process.nanoTableTaskCommon, _nanoTableTaskCommonRun3
        )
        run3_common.toModify(
            process.ptRatioRelForEle, srcJet="updatedJets"
        )
        run3_common.toModify(
            process.ptRatioRelForMu, srcJet="updatedJets"
        )


    return process

outputTable = cms.EDProducer("EXOnanoAODProducerTemplate",
    inputExample = cms.InputTag("offlineBeamSpot")
)

def add_customTables_template(process):

    process.outputTable = outputTable
    process.exonanoaodTask = cms.Task(process.outputTable)
    process.nanoTableTaskCommon.add(process.exonanoaodTask)

    return process
