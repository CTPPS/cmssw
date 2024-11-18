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

#include "CalibPPS/TimingCalibration/interface/DoublePeakCorrection.h"
#include "CalibPPS/TimingCalibration/interface/PlaneMap.h"
#include "CalibPPS/TimingCalibration/interface/TimingCalibrationData.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondRecHit.h"

#include "DQMServices/Core/interface/DQMGlobalEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSGeometry.h"

//------------------------------------------------------------------------------

class PPSTimingCalibrationPCLWorker : public DQMGlobalEDAnalyzer<TimingCalibrationData> {
public:
  explicit PPSTimingCalibrationPCLWorker(const edm::ParameterSet&);

  void dqmBeginRun(const edm::Run&, const edm::EventSetup&, TimingCalibrationData&) const override;
  void dqmAnalyze(const edm::Event&, const edm::EventSetup&, const TimingCalibrationData&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  using DiamondRecHitVector = edm::DetSetVector<CTPPSDiamondRecHit>;

  void bookHistograms(DQMStore::IBooker&,
                      const edm::Run&,
                      const edm::EventSetup&,
                      TimingCalibrationData&) const override;

  void searchForProduct(const edm::Event& iEvent,
                        const std::vector<edm::EDGetTokenT<DiamondRecHitVector>>& tokens,
                        const std::vector<edm::InputTag>& tags,
                        edm::Handle<DiamondRecHitVector>& handle) const;

  const std::string dqmDir_;
  const std::string tVsLsFilename_;
  std::vector<edm::EDGetTokenT<DiamondRecHitVector>> diamondRecHitTokens_;
  const std::vector<edm::InputTag> recHitTags_;
  const edm::ESGetToken<CTPPSGeometry, VeryForwardRealGeometryRecord> geomEsToken_;
};

//------------------------------------------------------------------------------

PPSTimingCalibrationPCLWorker::PPSTimingCalibrationPCLWorker(const edm::ParameterSet& iConfig)
    : dqmDir_{iConfig.getParameter<std::string>("dqmDir")},
      tVsLsFilename_{iConfig.getParameter<std::string>("tVsLsFilename")},
      recHitTags_{iConfig.getParameter<std::vector<edm::InputTag>>("diamondRecHitTags")},
      geomEsToken_{esConsumes<edm::Transition::BeginRun>()} {
  for (const auto& tag : recHitTags_) {
    diamondRecHitTokens_.push_back(consumes<DiamondRecHitVector>(tag));
  }
}

//------------------------------------------------------------------------------

void PPSTimingCalibrationPCLWorker::dqmBeginRun(const edm::Run& iRun,
                                                const edm::EventSetup& iSetup,
                                                TimingCalibrationData& timingCalibrationData) const {
  std::vector<CTPPSDiamondDetId> detIds;
  const auto& geom = iSetup.getData(geomEsToken_);
  for (auto it = geom.beginSensor(); it != geom.endSensor(); ++it) {
    if (CTPPSDiamondDetId::check(it->first)) {
      const CTPPSDiamondDetId detId{it->first};
      detIds.push_back(detId);
    }
  }
  timingCalibrationData.doublePeakCorrection.extractLsAndTimeOffset(tVsLsFilename_, iRun.run(), detIds);
}

//------------------------------------------------------------------------------

void PPSTimingCalibrationPCLWorker::bookHistograms(DQMStore::IBooker& iBooker,
                                                   const edm::Run& iRun,
                                                   const edm::EventSetup& iSetup,
                                                   TimingCalibrationData& timingCalibrationData) const {
  iBooker.cd();
  iBooker.setCurrentFolder(dqmDir_);

  std::string planeName;
  std::string channelName;
  const auto& geom = iSetup.getData(geomEsToken_);
  for (auto it = geom.beginSensor(); it != geom.endSensor(); ++it) {
    if (CTPPSDiamondDetId::check(it->first)) {
      const CTPPSDiamondDetId detId{it->first};
      const uint32_t channelId{detId.rawId()};
      detId.channelName(channelName);
      const PlaneKey planeKey{detId.arm(), detId.station(), detId.plane()};

      timingCalibrationData.leadingTime[channelId] =
          iBooker.book1D("t_" + channelName, channelName + ";t (ns);Entries", 1200, -60.0, 60.0);
      timingCalibrationData.toT[channelId] =
          iBooker.book1D("tot_" + channelName, channelName + ";ToT (ns);Entries", 160, -20.0, 20.0);
      timingCalibrationData.leadingTimeVsToT[channelId] =
          iBooker.book2D("tvstot_" + channelName, channelName + ";ToT (ns);t (ns)", 240, 0.0, 60.0, 450, -20.0, 25.0);
      if (!timingCalibrationData.leadingTimeVsLs.contains(planeKey)) {
        detId.planeName(planeName);
        timingCalibrationData.leadingTimeVsLs[planeKey] =
            iBooker.book2D("tvsls_" + planeName, planeName + ";LS;t (ns)", 3000, 1.0, 3001.0, 500, 0.0, 20.0);
      }
    }
  }
}

//------------------------------------------------------------------------------

void PPSTimingCalibrationPCLWorker::dqmAnalyze(const edm::Event& iEvent,
                                               const edm::EventSetup& iSetup,
                                               const TimingCalibrationData& timingCalibrationData) const {
  edm::Handle<DiamondRecHitVector> dsvRechits;
  // then extract the rechits information for later processing
  searchForProduct(iEvent, diamondRecHitTokens_, recHitTags_, dsvRechits);

  // ensure timing detectors rechits are found in the event content
  if (dsvRechits->empty()) {
    edm::LogWarning{"PPSTimingCalibrationPCLWorker:dqmAnalyze"}
        << "No rechits retrieved from the event content for event " << iEvent.id() << '.';
    return;
  }

  const unsigned int ls{iEvent.luminosityBlock()};
  std::string planeName;
  for (const auto& dsRechits : *dsvRechits) {
    const CTPPSDiamondDetId detId{dsRechits.detId()};
    const uint32_t channelId{detId.rawId()};
    if (!timingCalibrationData.leadingTimeVsToT.contains(channelId)) {
      edm::LogWarning{"PPSTimingCalibrationPCLWorker:dqmAnalyze"} << "Pad with Detector ID =" << detId
                                                                  << " is not set to be monitored.";
      continue;
    }

    detId.planeName(planeName);
    const PlaneKey planeKey{detId.arm(), detId.station(), detId.plane()};
    for (const auto& rechit : dsRechits) {
      // skip invalid rechits
      if (rechit.time() != 0.0 && rechit.toT() >= 0.0) {
        const double correctedLeadingTime{
            timingCalibrationData.doublePeakCorrection.getCorrectedLeadingTime(rechit.time(), ls, planeKey)};
        timingCalibrationData.leadingTime.at(channelId)->Fill(correctedLeadingTime);
        timingCalibrationData.toT.at(channelId)->Fill(rechit.toT());
        timingCalibrationData.leadingTimeVsToT.at(channelId)->Fill(rechit.toT(), correctedLeadingTime);
        timingCalibrationData.leadingTimeVsLs.at(planeKey)->Fill(ls, correctedLeadingTime);
      }
    }
  }
}

//------------------------------------------------------------------------------

void PPSTimingCalibrationPCLWorker::searchForProduct(const edm::Event& iEvent,
                                                     const std::vector<edm::EDGetTokenT<DiamondRecHitVector>>& tokens,
                                                     const std::vector<edm::InputTag>& tags,
                                                     edm::Handle<DiamondRecHitVector>& handle) const {
  bool foundProduct{false};
  for (unsigned int i{0}; i < tokens.size(); ++i)
    if (const auto h = iEvent.getHandle(tokens[i])) {
      handle = h;
      foundProduct = true;
      edm::LogInfo{"searchForProduct"} << "Found a product with " << tags[i] << '.';
      break;
    }

  if (!foundProduct) {
    throw edm::Exception{edm::errors::ProductNotFound} << "Could not find a product with any of the selected labels.";
  }
}

//------------------------------------------------------------------------------

void PPSTimingCalibrationPCLWorker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("dqmDir", "AlCaReco/PPSTimingCalibrationPCL")
      ->setComment("output path for the various DQM plots");
  desc.add<std::string>("tVsLsFilename", "")
      ->setComment("ROOT filename with t vs LS histogram for double peak correction");
  desc.add<std::vector<edm::InputTag>>("diamondRecHitTags", {edm::InputTag{"ctppsDiamondRecHits"}})
      ->setComment("input tag for the PPS diamond detectors rechits");

  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(PPSTimingCalibrationPCLWorker);
