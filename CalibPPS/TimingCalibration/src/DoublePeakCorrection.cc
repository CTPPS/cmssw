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
      lsAndTimeOffsets_[planeKey] = {
          doublePeakLs, findTimeOffset(tVsLs, doublePeakLs, firstPeakTWithMaxCount, secondPeakTWithMaxCount)};
    }
  }
}

std::tuple<unsigned int, double, double> DoublePeakCorrection::findLsAndTimePeaks(const TH2F* tVsLs,
                                                                                  const PlaneKey& planeKey) const {
  auto numOfLsBins = static_cast<unsigned int>(tVsLs->GetNbinsX());
  auto numOfTBins = static_cast<unsigned int>(tVsLs->GetNbinsY());
  double firstPeakTWithMaxCount{0.0};
  double secondPeakTWithMaxCount{0.0};
  unsigned int tDiffCount{0};
  unsigned int doublePeakLs{1};
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

    constexpr double tMinDiff{1.0};
    if (firstPeakTWithMaxCount != 0.0 && std::abs(secondPeakTWithMaxCount - firstPeakTWithMaxCount) > tMinDiff) {
      constexpr unsigned int minTDiffCount{5};
      ++tDiffCount;
      if (tDiffCount == 1) {
        doublePeakLs = lsBin;
      } else if (tDiffCount == minTDiffCount) {
        return {doublePeakLs, firstPeakTWithMaxCount, secondPeakTWithMaxCount};
      }
    } else {
      tDiffCount = 0;
      firstPeakTWithMaxCount = secondPeakTWithMaxCount;
    }
  }

  return {1, firstPeakTWithMaxCount, secondPeakTWithMaxCount};
}

double DoublePeakCorrection::findTimeOffset(const TH2F* tVsLs,
                                            const unsigned int doublePeakLs,
                                            const double firstPeakEstimatedMean,
                                            const double secondPeakEstimatedMean) const {
  const std::unique_ptr<TH1D> firstPeak{tVsLs->ProjectionY("_py", 1, doublePeakLs - 1)};
  const std::unique_ptr<TH1D> secondPeak{tVsLs->ProjectionY("_py", doublePeakLs, -1)};
  return findGaussianMean(secondPeak, secondPeakEstimatedMean) - findGaussianMean(firstPeak, firstPeakEstimatedMean);
}

double DoublePeakCorrection::findGaussianMean(const std::unique_ptr<TH1D>& peak, const double estimatedMean) const {
  constexpr unsigned int meanParamIndex{1};
  constexpr double fitSigma{2.5};
  const double fitLeftBound{estimatedMean - fitSigma};
  const double fitRightBound{estimatedMean + fitSigma};
  TF1 fitFunction{"peak", "gaus"};
  fitFunction.SetParLimits(meanParamIndex, fitLeftBound, fitRightBound);
  fitFunction.SetParameter(meanParamIndex, estimatedMean);
  const TFitResultPtr& peakFit{peak->Fit(&fitFunction, "NS", "", fitLeftBound, fitRightBound)};
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
    const auto [doublePeakLs, doublePeakTimeOffset] = it->second;
    constexpr unsigned int lsEncodingMultiple{100};
    if (doublePeakTimeOffset >= 0.0) {
      return doublePeakLs * lsEncodingMultiple + doublePeakTimeOffset;
    }
    return -(doublePeakLs * lsEncodingMultiple - doublePeakTimeOffset);
  }
  return 0.0;
}

double DoublePeakCorrection::getCorrectedLeadingTime(const double leadingTime,
                                                     const unsigned int ls,
                                                     const double encodedLsAndTimeOffset) {
  const unsigned int doublePeakLs = std::abs(encodedLsAndTimeOffset) / LsEncodingMultiple_;
  if (ls >= doublePeakLs) {
    const double doublePeakTimeOffset = encodedLsAndTimeOffset >= 0
                                            ? encodedLsAndTimeOffset - doublePeakLs * LsEncodingMultiple_
                                            : encodedLsAndTimeOffset + doublePeakLs * LsEncodingMultiple_;
    return leadingTime - doublePeakTimeOffset;
  }
  return leadingTime;
}
