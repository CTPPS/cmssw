#include "CalibPPS/TimingCalibration/interface/DoublePeakCorrection.h"

#include "FWCore/Utilities/interface/EDMException.h"

#include "TF1.h"
#include "TFitResult.h"

void DoublePeakCorrection::extractLsAndTimeOffset(const std::string& tVsLsFilename,
                                                  const unsigned int run,
                                                  const std::vector<CTPPSDiamondDetId>& detIds) {
  if (!tVsLsFilename.empty()) {
    TFile tVsLsFile{tVsLsFilename.c_str()};
    if (tVsLsFile.IsOpen()) {
      const std::string runNumber{std::to_string(run)};
      for (const auto& detId : detIds) {
        const PlaneKey planeKey{detId.arm(), detId.station(), detId.plane()};
        fillLsAndTimeOffset(getTVsLs(tVsLsFile, runNumber, detId), planeKey);
      }
    } else {
      throw edm::Exception{edm::errors::FileOpenError} << "Can't open the file with t vs LS: " << tVsLsFilename << '.';
    }
  }
}

const TH2F* DoublePeakCorrection::getTVsLs(TFile& tVsLsFile,
                                           const std::string& runNumber,
                                           const CTPPSDiamondDetId& detId) {
  std::string planeName;
  detId.planeName(planeName);
  std::string tVsLsHistPath{"DQMData/Run " + runNumber + "/AlCaReco/Run summary/PPSTimingCalibrationPCL/tvsls_" +
                            planeName};
  const auto* tVsLs = tVsLsFile.Get<TH2F>(tVsLsHistPath.c_str());
  if (tVsLs) {
    return tVsLs;
  }
  throw edm::Exception{edm::errors::FileReadError} << "Can't open the t vs LS histogram: " << tVsLsHistPath << '.';
}

void DoublePeakCorrection::fillLsAndTimeOffset(const TH2F* tVsLs, const PlaneKey& planeKey) {
  if (!lsAndTimeOffsets_.contains(planeKey)) {
    const auto [doublePeakLs, firstPeakTWithMaxCount, secondPeakTWithMaxCount] = findLsAndTimePeaks(tVsLs, planeKey);
    if (doublePeakLs != 1) {
      lsAndTimeOffsets_[planeKey] = {doublePeakLs,
                                    findTimeOffset(tVsLs, firstPeakTWithMaxCount, secondPeakTWithMaxCount)};
    }
  }
}

std::tuple<unsigned int, double, double> DoublePeakCorrection::findLsAndTimePeaks(const TH2F* tVsLs,
                                                                                  const PlaneKey& planeKey) const {
  auto numOfLsBins = static_cast<unsigned int>(tVsLs->GetNbinsX());
  auto numOfTBins = static_cast<unsigned int>(tVsLs->GetNbinsY());
  double firstPeakTWithMaxCount{0.0};
  double secondPeakTWithMaxCount{0.0};
  for (unsigned int lsBin{1}; lsBin <= numOfLsBins; ++lsBin) {
    double tMaxCount{0};
    for (unsigned int tBin{1}; tBin <= numOfTBins; ++tBin) {
      const double tCount{tVsLs->GetBinContent(lsBin, tBin)};
      const double tBinCenter{tVsLs->GetYaxis()->GetBinCenter(tBin)};
      constexpr double minTCount{10.0};
      if (tCount >= minTCount && tCount > tMaxCount) {
        tMaxCount = tCount;
        secondPeakTWithMaxCount = tBinCenter;
      }
    }

    if (firstPeakTWithMaxCount != 0.0 && std::abs(secondPeakTWithMaxCount - firstPeakTWithMaxCount) > TMaxDiff_) {
      return {lsBin, firstPeakTWithMaxCount, secondPeakTWithMaxCount};
    }

    firstPeakTWithMaxCount = secondPeakTWithMaxCount;
  }

  return {1, firstPeakTWithMaxCount, secondPeakTWithMaxCount};
}

double DoublePeakCorrection::findTimeOffset(const TH2F* tVsLs,
                                            const double firstPeakEstimatedMean,
                                            const double secondPeakEstimatedMean) const {
  const std::unique_ptr<TH1D> tProjection{tVsLs->ProjectionY()};
  return findGaussianMean(tProjection, secondPeakEstimatedMean) - findGaussianMean(tProjection, firstPeakEstimatedMean);
}

double DoublePeakCorrection::findGaussianMean(const std::unique_ptr<TH1D>& tProjection,
                                              const double estimatedMean) const {
  constexpr unsigned int meanParamIndex{1};
  const double fitLeftBound{estimatedMean - TMaxDiff_};
  const double fitRightBound{estimatedMean + TMaxDiff_};
  TF1 fitFunction{"peak", "gaus"};
  fitFunction.SetParLimits(meanParamIndex, fitLeftBound, fitRightBound);
  fitFunction.SetParameter(meanParamIndex, estimatedMean);
  const TFitResultPtr& peakFit{tProjection->Fit(&fitFunction, "NS", "", fitLeftBound, fitRightBound)};
  if (peakFit->IsValid()) {
    return peakFit->Parameter(meanParamIndex);
  }
  throw edm::Exception{edm::errors::FatalRootError} << "Double peak Gaussian fit not valid.";
}

bool DoublePeakCorrection::isCorrectionNeeded(const PlaneKey& planeKey) const {
  return lsAndTimeOffsets_.contains(planeKey);
}

double DoublePeakCorrection::getCorrectedLeadingTime(const double leadingTime,
                                                     const unsigned int ls,
                                                     const PlaneKey& planeKey) const {
  if (auto it = lsAndTimeOffsets_.find(planeKey); it != std::end(lsAndTimeOffsets_)) {
    const auto [doublePeakLs, doublePeakTimeOffset] = it->second;
    if (ls >= doublePeakLs) {
      return leadingTime - doublePeakTimeOffset;
    }
  }
  return leadingTime;
}

double DoublePeakCorrection::getEncodedLsAndTimeOffset(const PlaneKey& planeKey) const {
  if (auto it = lsAndTimeOffsets_.find(planeKey); it != std::end(lsAndTimeOffsets_)) {
    constexpr double encodingMultiple = 100'000.0;
    const auto [doublePeakLs, doublePeakTimeOffset] = it->second;
    return doublePeakLs * encodingMultiple + doublePeakTimeOffset;
  }
  return 0.0;
}
