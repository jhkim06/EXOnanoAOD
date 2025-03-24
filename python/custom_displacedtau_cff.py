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

'''
def addPFCandsDxyz(process, allPF = False, addAK4=False, addAK8=False):
    process.customizedPFCandsTask = cms.Task()
    process.schedule.associate(process.customizedPFCandsTask)
    
    process.finalJetsAK8Constituents = cms.EDProducer("PatJetConstituentPtrSelector",
                                                      src = cms.InputTag("finalJetsAK8"),
                                                      cut = cms.string("")
                                                      )
    process.finalJetsAK4Constituents = cms.EDProducer("PatJetConstituentPtrSelector",
                                                      src = cms.InputTag("finalJetsPuppi"),
                                                      cut = cms.string("")
                                                      )
    if allPF:
        candInput = cms.InputTag("packedPFCandidates")
    elif not addAK8:
        candList = cms.VInputTag(cms.InputTag("finalJetsAK4Constituents", "constituents"))
        process.customizedPFCandsTask.add(process.finalJetsAK4Constituents)
        process.finalJetsConstituentsTable = cms.EDProducer("PackedCandidatePtrMerger", src = candList, skipNulls = cms.bool(True), warnOnSkip = cms.bool(True))
        candInput = cms.InputTag("finalJetsConstituentsTable")
    elif not addAK4:
        candList = cms.VInputTag(cms.InputTag("finalJetsAK8Constituents", "constituents"))
        process.customizedPFCandsTask.add(process.finalJetsAK8Constituents)
        process.finalJetsConstituentsTable = cms.EDProducer("PackedCandidatePtrMerger", src = candList, skipNulls = cms.bool(True), warnOnSkip = cms.bool(True))
        candInput = cms.InputTag("finalJetsConstituentsTable")
    else:
        candList = cms.VInputTag(cms.InputTag("finalJetsAK4Constituents", "constituents"), cms.InputTag("finalJetsAK8Constituents", "constituents"))
        process.customizedPFCandsTask.add(process.finalJetsAK4Constituents)
        process.customizedPFCandsTask.add(process.finalJetsAK8Constituents)
        process.finalJetsConstituentsTable = cms.EDProducer("PackedCandidatePtrMerger", src = candList, skipNulls = cms.bool(True), warnOnSkip = cms.bool(True))
        candInput = cms.InputTag("finalJetsConstituentsTable")
        
    process.customConstituentsExtTable = cms.EDProducer("SimplePATCandidateFlatTableProducer",
                                                        src = candInput,
                                                        cut = cms.string(""), #we should not filter after pruning
                                                        name = cms.string("PFCands"),
                                                        doc = cms.string("interesting particles from AK4 and AK8 jets"),
                                                        singleton = cms.bool(False), # the number of entries is variable
                                                        extension = cms.bool(False), # this is the extension table for the AK8 constituents
                                                        variables = cms.PSet(##//CandVars,
                                                            dz = Var("?hasTrackDetails()?dz():-1", float, doc="pf dz", precision=10),
                                                            dzErr = Var("?hasTrackDetails()?dzError():-1", float, doc="pf dz err", precision=10),
                                                            d0 = Var("?hasTrackDetails()?dxy():-1", float, doc="pf d0", precision=10),
                                                            d0Err = Var("?hasTrackDetails()?dxyError():-1", float, doc="pf d0 err", precision=10),
                                                            trkP = Var("?hasTrackDetails()?pseudoTrack().p():-1", float, doc="track momemtum", precision=-1),
                                                            trkPt = Var("?hasTrackDetails()?pseudoTrack().pt():-1", float, doc="track pt", precision=-1),
                                                            trkEta = Var("?hasTrackDetails()?pseudoTrack().eta():-1", float, doc="track pt", precision=12),
                                                            trkPhi = Var("?hasTrackDetails()?pseudoTrack().phi():-1", float, doc="track phi", precision=12),
                                                            ),
                                                        externalVariables = cms.PSet(
                                                            rankByPt = cms.VPSet(
                                                                cms.PSet(
                                                                    tag = cms.string("rank"),
                                                                    quantity = cms.string("-pt()")  # Sorting by descending pt
                                                                )
                                                            )
                                                        ),
                                                        
                                                        maxLen = cms.uint32(50) 
                                                        

                                                        )
    kwargs = { }
    import os
    sv_sort = os.getenv('CMSSW_NANOAOD_SV_SORT')
    if sv_sort is not None: kwargs['sv_sort'] = cms.untracked.string(sv_sort)
    pf_sort = os.getenv('CMSSW_NANOAOD_PF_SORT')
    if pf_sort is not None: kwargs['pf_sort'] = cms.untracked.string(pf_sort)
    process.customAK8ConstituentsTable = cms.EDProducer("PatJetConstituentTableProducer",
                                                        candidates = candInput,
                                                        jets = cms.InputTag("finalJetsAK8"),
                                                        jet_radius = cms.double(0.8),
                                                        name = cms.string("FatJetPFCands"),
                                                        idx_name = cms.string("pFCandsIdx"),
                                                        nameSV = cms.string("FatJetSVs"),
                                                        idx_nameSV = cms.string("sVIdx"),
                                                        **kwargs,
                                                        )
    process.customAK4ConstituentsTable = cms.EDProducer("PatJetConstituentTableProducer",
                                                        candidates = candInput,
                                                        jets = cms.InputTag("finalJetsPuppi"), # was finalJets before
                                                        jet_radius = cms.double(0.4),
                                                        name = cms.string("JetPFCands"),
                                                        idx_name = cms.string("pFCandsIdx"),
                                                        nameSV = cms.string("JetSVs"),
                                                        idx_nameSV = cms.string("sVIdx"),
                                                        **kwargs,
                                                        )
    process.customizedPFCandsTask.add(process.customConstituentsExtTable)

    if not allPF:
        process.customizedPFCandsTask.add(process.finalJetsConstituentsTable)
    # linkedObjects are WIP for Run3
    if addAK8:
        process.customizedPFCandsTask.add(process.customAK8ConstituentsTable)
    if addAK4: 
        process.customizedPFCandsTask.add(process.customAK4ConstituentsTable)
    
        
    return process
'''


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
    
    process.jetImpactParameters = cms.EDProducer(
        "JetImpactParameters",
        jets = process.jetTable.src,
        pfCandidates = cms.InputTag('packedPFCandidates'),
        deltaRMax = cms.double(0.4)
    )

    
    
    d_disTauTagVars = {
        "disTauTag_score0":     ExtVar("disTauTag:score0"       , float, doc = "Score 0"),
        "disTauTag_score1":     ExtVar("disTauTag:score1"       , float, doc = "Score 1"),
        "dxy": ExtVar("jetImpactParameters:jetDxy", float, doc = "leadingPtPFCand_dxy which is within dR=0.4 and charged/hasTrackDetails"),
        "dz": ExtVar("jetImpactParameters:jetDz", float, doc = "leadingPtPFCand_dz which is within dR=0.4 and charged/hasTrackDetails"),
        "dxyerror": ExtVar("jetImpactParameters:jetDxyError", float, doc = "leadingPtPFCand_dxyerror which is within dR=0.4 and charged/hasTrackDetails"),
        "dzerror": ExtVar("jetImpactParameters:jetDzError", float, doc = "leadingPtPFCand_dzerror which is within dR=0.4 and charged/hasTrackDetails"), 
    }

    print ('adding disTau edproducer')
    if useCHSJets:
        process.jetTable.externalVariables = process.jetTable.externalVariables.clone(**d_disTauTagVars)
        ## for puppi jets, use this!
        
    else:
        process.jetPuppiTable.externalVariables = process.jetPuppiTable.externalVariables.clone(**d_disTauTagVars)
      
    ##//process.custom_nanoaod_task = cms.Task(  
    process.disTauTagTask = cms.Task(process.disTauTag)
    process.jetImpactParametersTask = cms.Task(process.jetImpactParameters)
    ##process.schedule.associate(process.custom_nanoaod_task)
    ##//process.schedule.associate(process.exonanoaodTask1, process.exonanoaodTask2)
    ##//process.nanoTableTaskCommon.add(process.exonanoaodTask)
    ##//process.nanoTableTaskCommon.add(process.exonanoaodTask1)
    process.nanoTableTaskCommon.add(process.disTauTagTask, process.jetImpactParametersTask)
    
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


def BTVCustomNanoAODStausDxyzInfo(process, useCHSJets = True):
    addPFCandsDxyz(process,btvNano_switch.btvNano_addallPF_switch,btvNano_switch.btvNano_addAK4_switch,btvNano_switch.btvNano_addAK8_switch)
    
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
