#ifndef CalibPPS_TimingCalibration_DoublePeakCorrection_h
#define CalibPPS_TimingCalibration_DoublePeakCorrection_h

// #include "CalibPPS/TimingCalibration/interface/PlaneMap.h"
#include "../interface/PlaneMap.h"

#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"

#include "TFile.h"
#include "TH2F.h"

#include <utility>

class DoublePeakCorrection {
public:
  void extractLsAndTimeOffset(const std::string&, const unsigned int, const std::vector<CTPPSDiamondDetId>&);
  void fillLsAndTimeOffset(const TH2F*, const PlaneKey&);
  bool isCorrectionNeeded(const PlaneKey&) const;
  double getCorrectedLeadingTime(const double, const unsigned int, const PlaneKey&) const;
  double getEncodedLsAndTimeOffset(const PlaneKey&) const;

private:
  const TH2F* getTVsLs(TFile&, const std::string&, const CTPPSDiamondDetId&);
  std::tuple<unsigned int, double, double> findLsAndTimePeaks(const TH2F*, const PlaneKey&) const;
  double findTimeOffset(const TH2F*, const unsigned int, const double, const double) const;
  double findGaussianMean(const std::unique_ptr<TH1D>&, const double) const;

  std::unordered_map<PlaneKey, std::pair<unsigned int, double>, PlaneKeyHash> lsAndTimeOffsets_;
};

#endif
