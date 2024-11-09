/****************************************************************************
 *
 * This is a part of PPS offline software.
 * Authors:
 *   Edoardo Bossini
 *   Piotr Maciej Cwiklicki
 *   Laurent Forthomme
 *   Tomasz Ostafin
 *
 ****************************************************************************/

#include "CalibPPS/TimingCalibration/interface/TimingCalibrationStruct.h"

#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"

#include "CondFormats/PPSObjects/interface/PPSTimingCalibration.h"

#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"

#include "DQMServices/Core/interface/DQMEDHarvester.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSGeometry.h"

#include "Math/MinimizerOptions.h"
#include "TFitResult.h"

//------------------------------------------------------------------------------

class PPSTimingCalibrationPCLHarvester : public DQMEDHarvester {
public:
  PPSTimingCalibrationPCLHarvester(const edm::ParameterSet&);

  void beginRun(const edm::Run&, const edm::EventSetup&) override;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void dqmEndJob(DQMStore::IBooker&, DQMStore::IGetter&) override;

  bool fetchWorkerHistograms(DQMStore::IGetter&,
                             TimingCalibrationHistograms&,
                             const CTPPSDiamondDetId&,
                             const uint32_t,
                             const std::string&) const;
  TProfile* createTVsTotProfile(DQMStore::IBooker&, dqm::reco::MonitorElement*, const std::string&) const;
  std::pair<double, double> findFitRange(dqm::reco::MonitorElement*, const double, const double) const;

  const std::string dqmDir_;
  const std::string formula_;
  std::vector<CTPPSDiamondDetId> detIds_;
  const edm::ESGetToken<CTPPSGeometry, VeryForwardRealGeometryRecord> geomEsToken_;
  const unsigned int minEntries_;

  static constexpr double FixedFitBoundIndication_{-1.0};
};

//------------------------------------------------------------------------------

PPSTimingCalibrationPCLHarvester::PPSTimingCalibrationPCLHarvester(const edm::ParameterSet& iConfig)
    : dqmDir_{iConfig.getParameter<std::string>("dqmDir")},
      formula_{iConfig.getParameter<std::string>("formula")},
      geomEsToken_{esConsumes<edm::Transition::BeginRun>()},
      minEntries_{iConfig.getParameter<unsigned int>("minEntries")} {
  // first ensure DB output service is available
  edm::Service<cond::service::PoolDBOutputService> poolDbService;
  if (!poolDbService.isAvailable()) {
    throw cms::Exception{"PPSTimingCalibrationPCLHarvester"} << "PoolDBService required";
  }
}

//------------------------------------------------------------------------------

void PPSTimingCalibrationPCLHarvester::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
  const auto& geom = iSetup.getData(geomEsToken_);
  for (auto it = geom.beginSensor(); it != geom.endSensor(); ++it) {
    if (CTPPSDiamondDetId::check(it->first)) {
      const CTPPSDiamondDetId detId{it->first};
      detIds_.push_back(detId);
    }
  }
}

//------------------------------------------------------------------------------

void PPSTimingCalibrationPCLHarvester::dqmEndJob(DQMStore::IBooker& iBooker, DQMStore::IGetter& iGetter) {
  TF1 interp{"interp", formula_.c_str()};
  interp.SetParLimits(0, 0.5, 5.0);
  interp.SetParLimits(1, 4.0, 15.0);
  interp.SetParLimits(2, 0.1, 4.0);
  interp.SetParLimits(3, 0.1, 15.0);

  // set a higher max function calls limit for the ROOT fit algorithm
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(5'000);

  // book the parameters containers
  PPSTimingCalibration::ParametersMap calibParams;
  PPSTimingCalibration::TimingMap calibTime;

  iGetter.cd();
  iGetter.setCurrentFolder(dqmDir_);

  constexpr double defaultFitSlope{0.8};
  constexpr double defaultOffset{0.0};
  constexpr double defaultResolution{0.1};
  constexpr std::array<double, 13> thresholds{
      {FixedFitBoundIndication_, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07}};

  // compute the fit parameters for all monitored channels
  TimingCalibrationHistograms hists;
  std::string channelName;
  for (const auto& detId : detIds_) {
    const uint32_t channelId{detId.rawId()};
    detId.channelName(channelName);
    if (fetchWorkerHistograms(iGetter, hists, detId, channelId, channelName)) {
      const PPSTimingCalibration::Key armKey{static_cast<int>(detId.arm()),
                                             static_cast<int>(detId.station()),
                                             static_cast<int>(detId.plane()),
                                             static_cast<int>(detId.channel())};

      TProfile* const tVsTotProfile{createTVsTotProfile(iBooker, hists.leadingTimeVsToT.at(channelId), channelName)};

      const double defaultUpperLowerAsymptotesDiff = hists.leadingTime[channelId]->getRMS();
      const double defaultCenterOfDistribution = hists.toT[channelId]->getMean();
      const double defaultLowerAsymptote = hists.leadingTime[channelId]->getMean() - defaultUpperLowerAsymptotesDiff;

      double bestChiSqDivNdf{std::numeric_limits<double>::max()};
      double bestLowerTotRange{0.0};
      double bestUpperTotRange{0.0};
      for (const double upperThresholdFractionOfMax : thresholds) {
        for (const double lowerThresholdFractionOfMax : thresholds) {
          interp.SetParameters(
              defaultUpperLowerAsymptotesDiff, defaultCenterOfDistribution, defaultFitSlope, defaultLowerAsymptote);
          const auto [lowerTotRange, upperTotRange] =
              findFitRange(hists.toT.at(channelId), lowerThresholdFractionOfMax, upperThresholdFractionOfMax);

          const TFitResultPtr& tVsTotProfileFitResult{
              tVsTotProfile->Fit(&interp, "BNS", "", lowerTotRange, upperTotRange)};
          if (tVsTotProfileFitResult->IsValid()) {
            const double chiSqDivNdf{tVsTotProfileFitResult->Chi2() / tVsTotProfileFitResult->Ndf()};
            if (chiSqDivNdf < bestChiSqDivNdf) {
              bestChiSqDivNdf = chiSqDivNdf;
              bestUpperTotRange = upperTotRange;
              bestLowerTotRange = lowerTotRange;
            }
          }
        }
      }

      calibParams[armKey] = {0.0, 0.0, 0.0, 0.0};
      calibTime[armKey] = {defaultOffset, defaultResolution};
      if (bestUpperTotRange != 0.0) {
        tVsTotProfile->Fit(&interp, "B", "", bestLowerTotRange, bestUpperTotRange);
        calibParams[armKey] = {
            interp.GetParameter(0), interp.GetParameter(1), interp.GetParameter(2), interp.GetParameter(3)};
        calibTime[armKey] = {defaultOffset,
                             defaultResolution};  // hardcoded offset/resolution placeholder for the time being
      } else {
        edm::LogWarning{"PPSTimingCalibrationPCLHarvester:dqmEndJob"} << "Fit did not converge for channel (" << detId
                                                                      << ").";
      }
    }

    // fill the DB object record
    PPSTimingCalibration calib{formula_, calibParams, calibTime};

    // write the object
    edm::Service<cond::service::PoolDBOutputService> poolDbService;
    poolDbService->writeOneIOV(calib, poolDbService->currentTime(), "PPSTimingCalibrationRcd_HPTDC");
  }
}

//------------------------------------------------------------------------------

bool PPSTimingCalibrationPCLHarvester::fetchWorkerHistograms(DQMStore::IGetter& iGetter,
                                                             TimingCalibrationHistograms& hists,
                                                             const CTPPSDiamondDetId& detId,
                                                             const uint32_t channelId,
                                                             const std::string& channelName) const {
  hists.leadingTime[channelId] = iGetter.get(dqmDir_ + "/t_" + channelName);
  if (!hists.leadingTime.at(channelId)) {
    edm::LogWarning{"PPSTimingCalibrationPCLHarvester:fetchWorkerHistograms"}
        << "Failed to retrieve leading time monitor for channel (" << detId << "). Skipping calibration.";
    return false;
  }

  hists.toT[channelId] = iGetter.get(dqmDir_ + "/tot_" + channelName);
  if (!hists.toT.at(channelId)) {
    edm::LogWarning{"PPSTimingCalibrationPCLHarvester:fetchWorkerHistograms"}
        << "Failed to retrieve time over threshold monitor for channel (" << detId << "). Skipping calibration.";
    return false;
  }

  hists.leadingTimeVsToT[channelId] = iGetter.get(dqmDir_ + "/tvstot_" + channelName);
  if (!hists.leadingTimeVsToT.at(channelId)) {
    edm::LogWarning{"PPSTimingCalibrationPCLHarvester:fetchWorkerHistograms"}
        << "Failed to retrieve leading time vs. time over threshold monitor for channel (" << detId
        << "). Skipping calibration.";
    return false;
  }

  const auto tVsTotEntries = static_cast<unsigned int>(hists.leadingTimeVsToT.at(channelId)->getEntries());
  if (tVsTotEntries < minEntries_) {
    edm::LogWarning{"PPSTimingCalibrationPCLHarvester:fetchWorkerHistograms"}
        << "Not enough entries for channel (" << detId << "): " << tVsTotEntries << " < " << minEntries_
        << ". Skipping calibration.";
    return false;
  }

  return true;
}

//------------------------------------------------------------------------------

TProfile* PPSTimingCalibrationPCLHarvester::createTVsTotProfile(DQMStore::IBooker& iBooker,
                                                                dqm::reco::MonitorElement* tVsTot,
                                                                const std::string& channelName) const {
  const MonitorElement* const monitorProfile{
      iBooker.bookProfile(channelName, channelName, 240, 0.0, 60.0, 450, -20.0, 25.0)};
  TProfile* const monitorTProfile{monitorProfile->getTProfile()};
  const std::unique_ptr<const TProfile> tVsTotProfile{tVsTot->getTH2F()->ProfileX()};
  *(monitorTProfile) = *(static_cast<TProfile*>(tVsTotProfile->Clone()));

  const char* const profileName = channelName.c_str();
  monitorTProfile->SetTitle(profileName);
  monitorTProfile->SetName(profileName);
  monitorTProfile->SetYTitle("Average t (ns)");

  return monitorTProfile;
}

//------------------------------------------------------------------------------

std::pair<double, double> PPSTimingCalibrationPCLHarvester::findFitRange(
    dqm::reco::MonitorElement* tot,
    const double lowerThresholdFractionOfMax,
    const double upperThresholdFractionOfMax) const {
  constexpr double totUpperLimitMaxSearch{20.0};
  int maxTotBin{1};
  const int numOfToTBins{tot->getNbinsX()};
  const TAxis* const totXAxis{tot->getTH1()->GetXaxis()};
  for (int i{2}; i <= numOfToTBins || totXAxis->GetBinCenter(i) <= totUpperLimitMaxSearch; ++i) {
    if (tot->getBinContent(i) > tot->getBinContent(maxTotBin)) {
      maxTotBin = i;
    }
  }

  constexpr double totLowerLimitRangeSearch{8.0};
  double lowerTotRange{8};
  if (lowerThresholdFractionOfMax != FixedFitBoundIndication_) {
    int lowerLimitPos{maxTotBin};
    const double lowerThreshold{lowerThresholdFractionOfMax * tot->getBinContent(maxTotBin)};
    while (tot->getBinContent(lowerLimitPos) >= lowerThreshold &&
           totXAxis->GetBinCenter(lowerLimitPos) > totLowerLimitRangeSearch) {
      --lowerLimitPos;
    }
    lowerTotRange = totXAxis->GetBinCenter(lowerLimitPos);
  }

  constexpr double totUpperLimitRangeSearch{20.0};
  double upperTotRange{15};
  if (upperThresholdFractionOfMax != FixedFitBoundIndication_) {
    int upperLimitPos{maxTotBin};
    const double upperThreshold{upperThresholdFractionOfMax * tot->getBinContent(maxTotBin)};
    while (tot->getBinContent(upperLimitPos) >= upperThreshold &&
           totXAxis->GetBinCenter(upperLimitPos) < totUpperLimitRangeSearch) {
      ++upperLimitPos;
    }
    upperTotRange = totXAxis->GetBinCenter(upperLimitPos);
  }

  return {lowerTotRange, upperTotRange};
}

//------------------------------------------------------------------------------

void PPSTimingCalibrationPCLHarvester::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("dqmDir", "AlCaReco/PPSTimingCalibrationPCL")
      ->setComment("input path for the various DQM plots");
  desc.add<std::string>("formula", "[0]/(exp((x-[1])/[2])+1)+[3]")
      ->setComment("interpolation formula for the time walk component");
  desc.add<unsigned int>("minEntries", 100)->setComment("minimal number of hits to extract calibration");
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(PPSTimingCalibrationPCLHarvester);
