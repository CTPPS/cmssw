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
        empiricalAperture45 = "(1.0509e-07*[xangle]+-2.8678e-05)+(([xi]<(0.00035092*[xangle]+0.10406))*(-6.0606e-07*[xangle]+-0.0042798)+([xi]>=(0.00035092*[xangle]+0.10406))*(5.93568e-05*[xangle]+-0.016018))*([xi]-(0.00035092*[xangle]+0.10406))",
        empiricalAperture56 = "(-1.8215e-07*[xangle]+1.5374e-05)+(([xi]<(0.00030184*[xangle]+0.10506))*(3.9308e-06*[xangle]+-0.0053814)+([xi]>=(0.00030184*[xangle]+0.10506))*(4.3242e-05*[xangle]+-0.013780))*([xi]-(0.00030184*[xangle]+0.10506))"
    )
)

profile_2023_PreTS1A = _base_2023.clone(
    L_int = 1.,
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
    L_int = 1.,
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
    L_int = 1.,
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
