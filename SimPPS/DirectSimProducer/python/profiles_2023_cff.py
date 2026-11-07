import FWCore.ParameterSet.Config as cms
from SimPPS.DirectSimProducer.profile_base_cff import profile_base as _base
from CalibPPS.ESProducers.ctppsOpticalFunctions_non_DB_cff import optics_2022 as _optics

# base profile settings for 2023
_base_2023 = _base.clone(
    ctppsOpticalFunctions = _base.ctppsOpticalFunctions.clone(
        opticalFunctions = _optics.opticalFunctions,
        scoringPlanes = _optics.scoringPlanes,
    ),
    ctppsDirectSimuData = _base.ctppsDirectSimuData.clone(
        empiricalAperture45 = "(-7.5646421315777661e-08*[xangle]+3.3983287691553965e-07)+(([xi]<(0.00035395137081931403*[xangle]+0.10361157845282201))*(1.0455353725918569e-05*[xangle]+-0.0060158593045788787)+([xi]>=(0.00035395137081931403*[xangle]+0.10361157845282201))*(1.5555152274688385e-05*[xangle]+-0.0093198916512878573))*([xi]-(0.00035395137081931403*[xangle]+0.10361157845282201))",
        empiricalAperture56 = "(-7.5792161116942461e-07*[xangle]+0.00010386977124192235)+(([xi]<(0.00061722482073612542*[xangle]+0.056131117510882927))*(-3.3945438745896078e-06*[xangle]+-0.0042577282646826609)+([xi]>=(0.00061722482073612542*[xangle]+0.056131117510882927))*(7.7067899693672399e-05*[xangle]+-0.018977291288112909))*([xi]-(0.00061722482073612542*[xangle]+0.056131117510882927))"
    )
)

profile_2023_PreTS1A = _base_2023.clone(
    L_int = 7.599473874,
    ctppsLHCInfo = _base_2023.ctppsLHCInfo.clone(
        xangleBetaStarHistogramObject = "2023/h2_betaStar_vs_xangle"
    ),
    ctppsRPAlignmentCorrectionsDataXML = _base_2023.ctppsRPAlignmentCorrectionsDataXML.clone(
        MisalignedFiles = ["Validation/CTPPS/alignment/alignment_2023_PreTS1_366403_367840.xml"],
        RealFiles = ["Validation/CTPPS/alignment/alignment_2023_PreTS1_366403_367840.xml"]
    ),
    ctppsDirectSimuData = _base_2023.ctppsDirectSimuData.clone(
        timeResolutionDiamonds45 = "0.200",
        timeResolutionDiamonds56 = "0.200"
    )
)

profile_2023_PreTS1B = _base_2023.clone(
    L_int = 8.815992386,
    ctppsLHCInfo = _base_2023.ctppsLHCInfo.clone(
        xangleBetaStarHistogramObject = "2023/h2_betaStar_vs_xangle"
    ),
    ctppsRPAlignmentCorrectionsDataXML = _base_2023.ctppsRPAlignmentCorrectionsDataXML.clone(
        MisalignedFiles = ["Validation/CTPPS/alignment/alignment_2023_PreTS1_367881_368765.xml"],
        RealFiles = ["Validation/CTPPS/alignment/alignment_2023_PreTS1_367881_368765.xml"]
    ),
    ctppsDirectSimuData = _base_2023.ctppsDirectSimuData.clone(
        timeResolutionDiamonds45 = "0.200",
        timeResolutionDiamonds56 = "0.200"
    )
)

profile_2023_PostTS1 = _base_2023.clone(
    L_int = 9.682846668,
    ctppsLHCInfo = _base_2023.ctppsLHCInfo.clone(
        xangleBetaStarHistogramObject = "2023/h2_betaStar_vs_xangle"
    ),
    ctppsRPAlignmentCorrectionsDataXML = _base_2023.ctppsRPAlignmentCorrectionsDataXML.clone(
        MisalignedFiles = ["Validation/CTPPS/alignment/alignment_2023_PostTS1.xml"],
        RealFiles = ["Validation/CTPPS/alignment/alignment_2023_PostTS1.xml"]
    ),
    ctppsDirectSimuData = _base_2023.ctppsDirectSimuData.clone(
        timeResolutionDiamonds45 = "0.200",
        timeResolutionDiamonds56 = "0.200"
    )
)
