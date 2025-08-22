import FWCore.ParameterSet.Config as cms

from Configuration.ProcessModifiers.Era_Run3_CTPPS_directSim_cff import *
process = cms.Process('CTPPSTest', Run3_CTPPS_directSim)

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Validation.CTPPS.ctppsLHCInfoPlotter_cfi')
process.load('Configuration.Generator.randomXiThetaGunProducer_cfi')
process.load("CondCore.CondDB.CondDB_cfi")
process.load('SimPPS.DirectSimProducer.ctppsGregDucer_cfi')


# minimal logger settings
process.MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cout'),
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('INFO')
    )
)

# particle generator
process.generator.xi_max = 0.25
process.generator.theta_x_sigma = 60.e-6
process.generator.theta_y_sigma = 60.e-6

# default source
process.source = cms.Source("EmptySource",
    firstRun = cms.untracked.uint32(1),
)

process.CondDB.connect = 'frontier://FrontierProd/CMS_CONDITIONS'
process.PoolDBESSource = cms.ESSource("PoolDBESSource",
    process.CondDB,
    toGet = cms.VPSet(cms.PSet(
        record = cms.string('CTPPSPixelAnalysisMaskRcd'),
        tag = cms.string("CTPPSPixelAnalysisMask_Run3_v1_hlt"))
        ))

# random seeds
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    sourceSeed = cms.PSet(initialSeed = cms.untracked.uint32(98765)),
    generator = cms.PSet(initialSeed = cms.untracked.uint32(98766)),
    beamDivergenceVtxGenerator = cms.PSet(initialSeed = cms.untracked.uint32(3849)),
    ppsDirectProtonSimulation = cms.PSet(initialSeed = cms.untracked.uint32(4981))
)

# number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(int(100000))
)

# LHCInfo plotter
process.ctppsLHCInfoPlotter.outputFile = "simu_2025_lhcInfo.root"

# track distribution plotter
process.ctppsTrackDistributionPlotter = cms.EDAnalyzer("CTPPSTrackDistributionPlotter",
    tagTracks = cms.InputTag("ctppsLocalTrackLiteProducer"),
    outputFile = cms.string("simu_2025_tracks.root")
)

# reconstruction plotter
process.ctppsProtonReconstructionPlotter = cms.EDAnalyzer("CTPPSProtonReconstructionPlotter",
    tagTracks = cms.InputTag("ctppsLocalTrackLiteProducer"),
    tagRecoProtonsSingleRP = cms.InputTag("ctppsProtons", "singleRP"),
    tagRecoProtonsMultiRP = cms.InputTag("ctppsProtons", "multiRP"),
    outputFile = cms.string("simu_2025_protons.root")
)


# Greg plotter 1 - Unfiltered
process.ctppsGregPlotter = cms.EDAnalyzer("CTPPSGregPlotter",
    tagTracks = cms.InputTag("GenParticles"),
    outputFile = cms.string("simu_2018_Greg.root")#,
)

# Greg producer 1 - Calibration
process.ctppsGregDucer1 = cms.EDProducer("CTPPSGregDucer",
    tagTracks = cms.InputTag("GenParticlesNew"),
    hepMCTag = cms.InputTag("generator", "unsmeared"),
    # filename = cms.string("/afs/cern.ch/user/g/gjedrzej/private/mainTask/CMSSW_15_0_11/src/SimPPS/DirectSimProducer/cutFiles/thetaphilimits_-160urad_18cm_60cm_calib-nodet-xrphd-64bins.out")
    filename = cms.string("SimPPS/DirectSimProducer/cutFiles/thetaphilimits_-160urad_18cm_60cm_calib-nodet-xrphd-64bins.out")
)

# Greg plotter 2 - Calibration
process.ctppsGregPlotter2 = cms.EDAnalyzer("CTPPSGregPlotter",
    tagTracks = cms.InputTag("ctppsGregDucer1", "selectedProtons"),
    hepMCTag = cms.InputTag("ctppsGregDucer1", "selectedProtons"),
    outputFile = cms.string("simu_2018_64binsCalib.root")

)

# Greg producer 2 - Physics
process.ctppsGregDucer2 = cms.EDProducer("CTPPSGregDucer",
    tagTracks = cms.InputTag("GenParticlesNew"),
    hepMCTag = cms.InputTag("generator", "unsmeared"),
    # filename = cms.string("/afs/cern.ch/user/g/gjedrzej/private/mainTask/CMSSW_15_0_11/src/SimPPS/DirectSimProducer/cutFiles/thetaphilimits_-160urad_18cm_60cm_phys-nodet-xrphd-64bins.out")
    filename = cms.string("SimPPS/DirectSimProducer/cutFiles/thetaphilimits_-160urad_18cm_60cm_phys-nodet-xrphd-64bins.out")

)


# Greg plotter 3 - Physics
process.ctppsGregPlotter3 = cms.EDAnalyzer("CTPPSGregPlotter",
    tagTracks = cms.InputTag("ctppsGregDucer2", "selectedProtons"),
    hepMCTag = cms.InputTag("ctppsGregDucer2", "selectedProtons"),
    outputFile = cms.string("simu_2018_64binsPhys.root")
)



process.generation = cms.Path(process.generator)

process.validation = cms.Path(
    process.ctppsLHCInfoPlotter
    * process.ctppsTrackDistributionPlotter
    * process.ctppsProtonReconstructionPlotter
    * process.ctppsGregPlotter 
)

# Calibration
process.cutAndValidate = cms.Path(
    process.ctppsGregDucer1 * process.ctppsGregPlotter2
)

# Physics
process.cutAndValidate2 = cms.Path(
    process.ctppsGregDucer2 * process.ctppsGregPlotter3)


# processing path
process.schedule = cms.Schedule(
    process.generation,
    process.validation,
    process.cutAndValidate,
    process.cutAndValidate2,
)


# Built manually below
# from SimPPS.Configuration.Utils import setupPPSDirectSim
# setupPPSDirectSim(process)


# Import the direct sim manually
process.load('SimPPS.DirectSimProducer.ppsDirectProtonSimulation_cff')
process.ppsDirectProtonSimulation.verbosity = 0
process.ppsDirectProtonSimulation.useEmpiricalApertures = cms.bool(False) # TODO: disable this for now

process.directSimPPSTask = cms.Task(
    process.beamDivergenceVtxGenerator,
    process.ppsDirectProtonSimulation
)

process.directSimPPS = cms.Sequence(process.directSimPPSTask)

# Manual conditions for 2025
from CalibPPS.ESProducers.ctppsCompositeESSource_cfi import ctppsCompositeESSource as _esComp
from CalibPPS.ESProducers.ppsAssociationCuts_non_DB_cff import use_single_infinite_iov_entry, p2022
from CalibPPS.ESProducers.ppsAssociationCuts_non_DB_cff import ppsAssociationCutsESSource as _esAssCuts
from Geometry.VeryForwardGeometry.commons_cff import cloneGeometry

# TODO: probably skipping some RP IDs here (like the new timing RP)
from SimPPS.DirectSimProducer.simPPS2017_cfi import rpIds 
process.rpIds = rpIds.clone()

from SimPPS.DirectSimProducer.profile_base_cff import profile_base as _base


# take optics from /eos/cms/store/group/phys_pps/reconstruction/optical_functions/2025/version_0/
optics_2025 = cms.PSet(
  validityRange = cms.EventRange("000001:min - 999999:max"),

  opticalFunctions = cms.VPSet(
    cms.PSet( xangle = cms.double(80.0), fileName = cms.FileInPath("Validation/CTPPS/tmp_optical_functions/80urad.root") ),
    cms.PSet( xangle = cms.double(120.0), fileName = cms.FileInPath("Validation/CTPPS/tmp_optical_functions/120urad.root") ),
    cms.PSet( xangle = cms.double(140.0), fileName = cms.FileInPath("Validation/CTPPS/tmp_optical_functions/140urad.root") ),
    cms.PSet( xangle = cms.double(160.0), fileName = cms.FileInPath("Validation/CTPPS/tmp_optical_functions/160urad.root") )
  ),

  scoringPlanes = cms.VPSet(
    # z in cm
    cms.PSet( rpId = cms.uint32(2014838784), dirName = cms.string("XRPH_D6L5_B2"), z = cms.double(-21255.0) ),  # RP 003, pixel
    cms.PSet( rpId = cms.uint32(2056257536), dirName = cms.string("XRPH_A6L5_B2"), z = cms.double(-21507.8) ),  # RP 022, diamond
    cms.PSet( rpId = cms.uint32(2054160384), dirName = cms.string("XRPH_E6L5_B2"), z = cms.double(-21570.0) ),  # RP 016, diamond
    cms.PSet( rpId = cms.uint32(2023227392), dirName = cms.string("XRPH_B6L5_B2"), z = cms.double(-21955.0) ),  # RP 023, pixel

    cms.PSet( rpId = cms.uint32(2031616000), dirName = cms.string("XRPH_D6R5_B1"), z = cms.double(+21255.0) ),  # RP 103, pixel
    cms.PSet( rpId = cms.uint32(2073034752), dirName = cms.string("XRPH_A6R5_B1"), z = cms.double(+21507.8) ),  # RP 122, diamond
    cms.PSet( rpId = cms.uint32(2070937600), dirName = cms.string("XRPH_E6R5_B1"), z = cms.double(+21570.0) ),  # RP 116, diamond
    cms.PSet( rpId = cms.uint32(2040004608), dirName = cms.string("XRPH_B6R5_B1"), z = cms.double(+21955.0) ),  # RP 123, pixel
  )
)

profile_2025 = _base.clone(
    L_int = 1.,
    ctppsOpticalFunctions = _base.ctppsOpticalFunctions.clone(
        opticalFunctions = optics_2025.opticalFunctions,
        scoringPlanes = optics_2025.scoringPlanes,
    ),
    ctppsLHCInfo = _base.ctppsLHCInfo.clone(
        # NB: until a dedicated 2022 distributions are issued, it is OK to use 2021 ones here
        xangle = cms.double(160.),# TODO: figure out is this should be positive or negative (I think positive)
        betaStar = cms.double(0.68), # TODO: figure out if this is in the right units and if it should be beta_x or beta_y
        beamEnergy = cms.double(6.8e3),
        # xangle-beta* histogram is ignored since xangle is > 0
        xangleBetaStarHistogramFile = cms.string("CalibPPS/ESProducers/data/xangle_beta_distributions/version1.root"),
        xangleBetaStarHistogramObject = cms.string("")
    ),
    ctppsRPAlignmentCorrectionsDataXML = _base.ctppsRPAlignmentCorrectionsDataXML.clone(
        MisalignedFiles = ["Validation/CTPPS/alignment/null.xml"], 
        RealFiles = ["Validation/CTPPS/alignment/null.xml"] 
    ),
    ctppsDirectSimuData = _base.ctppsDirectSimuData.clone( 
        timeResolutionDiamonds45 = "0.200", # TODO: well, we got bigger problems than this at the moment
        timeResolutionDiamonds56 = "0.200", # TODO: well, we got bigger problems than this at the moment
        empiricalAperture45 = cms.string("999"), # TODO: huge apertures, equivalent to no cuts
        empiricalAperture56 = cms.string("999"), # TODO: huge apertures, equivalent to no cuts
        efficienciesPerRP = cms.VPSet(),
        efficienciesPerPlane = cms.VPSet()
    )
)


process.ppsAssociationCutsESSource = _esAssCuts.clone()
use_single_infinite_iov_entry(process.ppsAssociationCutsESSource, p2022)
process.XMLIdealGeometryESSource_CTPPS, _ctppsGeometryESModule = cloneGeometry('Geometry.VeryForwardGeometry.geometryRPFromDD_2025_cfi')

process.ctppsCompositeESSource = _esComp.clone(
    generateEveryNEvents = 100,
    periods = [profile_2025],
    compactViewTag = _ctppsGeometryESModule.compactViewTag,
    isRun2 = _ctppsGeometryESModule.isRun2
)


# handle clashes between simulation and GT conditions (if GT is there)
process.es_prefer_composrc = cms.ESPrefer('CTPPSCompositeESSource', 'ctppsCompositeESSource')
process.es_prefer_pixtopo = cms.ESPrefer('PPSPixelTopologyESSource', 'ppsPixelTopologyESSource')
process.es_prefer_lhcinfo = cms.ESPrefer('CTPPSBeamParametersFromLHCInfoESSource', 'ctppsBeamParametersFromLHCInfoESSource')
process.es_prefer_assocuts = cms.ESPrefer('PPSAssociationCutsESSource', 'ppsAssociationCutsESSource')


process.load('RecoPPS.Configuration.recoCTPPS_cff')

process.totemRPUVPatternFinder.tagRecHit = cms.InputTag('ppsDirectProtonSimulation')
process.ctppsPixelLocalTracks.tag = cms.InputTag('ppsDirectProtonSimulation')
process.ctppsDiamondLocalTracks.recHitsTag = cms.InputTag('ppsDirectProtonSimulation')


process.ppsDirectSim = cms.Path(process.directSimPPS * process.recoDirectSimPPS)
process.schedule.append(process.ppsDirectSim)




# Set this in >2022 because we don't know the exact values yet
# TODO: update this to the correct values 
process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetX45 = 0.
process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetY45 = 0.
process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetZ45 = 0.


# Needed for proper generation of settings
process.source.numberEventsInLuminosityBlock = process.ctppsCompositeESSource.generateEveryNEvents

# Set RP IDs for track plotter
process.ctppsTrackDistributionPlotter.rpId_45_F = process.rpIds.rp_45_F
process.ctppsTrackDistributionPlotter.rpId_45_N = process.rpIds.rp_45_N
process.ctppsTrackDistributionPlotter.rpId_56_N = process.rpIds.rp_56_N
process.ctppsTrackDistributionPlotter.rpId_56_F = process.rpIds.rp_56_F

# Set RP IDs for proton plotter
process.ctppsProtonReconstructionPlotter.rpId_45_F = process.rpIds.rp_45_F
process.ctppsProtonReconstructionPlotter.rpId_45_N = process.rpIds.rp_45_N
process.ctppsProtonReconstructionPlotter.rpId_56_N = process.rpIds.rp_56_N
process.ctppsProtonReconstructionPlotter.rpId_56_F = process.rpIds.rp_56_F

# import sys
# print(process.es_sources_())
# print(process.es_producers_())

# import sys
# sys.exit(1)