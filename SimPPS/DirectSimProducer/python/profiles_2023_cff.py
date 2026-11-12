import FWCore.ParameterSet.Config as cms
from SimPPS.DirectSimProducer.profile_base_cff import profile_base as _base
from CalibPPS.ESProducers.ctppsOpticalFunctions_non_DB_cff import optics_2023 as _optics

# base profile settings for 2023
_base_2023 = _base.clone(
    ctppsOpticalFunctions = _base.ctppsOpticalFunctions.clone(
        opticalFunctions = _optics.opticalFunctions,
        scoringPlanes = _optics.scoringPlanes,
    ),
    ctppsDirectSimuData = _base.ctppsDirectSimuData.clone(
        # empiricalAperture45 = "1.e3*([xi] - 0.20)",
        # empiricalAperture56 = "1.e3*([xi] - 0.20)"
        empiricalAperture45 = "-((-5.8817288613671971e-08*[xangle]+-2.4318829173592709e-06)+(([xi]<(0.00038194884207164474*[xangle]+0.099158200288645887))*(1.0640688972741761e-05*[xangle]+-0.0060423360587366541)+([xi]>=(0.00038194884207164474*[xangle]+0.099158200288645887))*(1.958230380556284e-05*[xangle]+-0.0099670700279281058))*([xi]-(0.00038194884207164474*[xangle]+0.099158200288645887)))",
        empiricalAperture56 = "-((-9.2964189150945028e-07*[xangle]+0.00013012378721884344)+(([xi]<(0.00061198337229620826*[xangle]+0.056930491854876891))*(-7.7137585092833934e-06*[xangle]+-0.003593109939745568)+([xi]>=(0.00061198337229620826*[xangle]+0.056930491854876891))*(7.0535324118609869e-05*[xangle]+-0.017972678645228026))*([xi]-(0.00061198337229620826*[xangle]+0.056930491854876891)))"
    )
)

profile_2023_PreTS1A = _base_2023.clone(
    L_int = 7.599473874,
    ctppsLHCInfo = _base_2023.ctppsLHCInfo.clone(
        xangleBetaStarHistogramObject = "2023/h2_betaStar_vs_xangle"
    ),
    ctppsRPAlignmentCorrectionsDataXML = _base_2023.ctppsRPAlignmentCorrectionsDataXML.clone(
        MisalignedFiles = ["Validation/CTPPS/alignment/2023_PreTS1_366403_367840.xml"],
        RealFiles = ["Validation/CTPPS/alignment/2023_PreTS1_366403_367840.xml"]
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
        MisalignedFiles = ["Validation/CTPPS/alignment/2023_PreTS1_367881_368765.xml"],
        RealFiles = ["Validation/CTPPS/alignment/2023_PreTS1_367881_368765.xml"]
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
        MisalignedFiles = ["Validation/CTPPS/alignment/2023_PostTS1.xml"],
        RealFiles = ["Validation/CTPPS/alignment/2023_PostTS1.xml"]
    ),
    ctppsDirectSimuData = _base_2023.ctppsDirectSimuData.clone(
        timeResolutionDiamonds45 = "0.200",
        timeResolutionDiamonds56 = "0.200"
    )
)
