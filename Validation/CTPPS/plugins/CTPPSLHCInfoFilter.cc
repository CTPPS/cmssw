/****************************************************************************
 * Authors:
 *   Jan Kašpar
 ****************************************************************************/

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/global/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "CondTools/RunInfo/interface/LHCInfoCombined.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

//----------------------------------------------------------------------------------------------------

class CTPPSLHCInfoFilter : public edm::global::EDFilter<> {
public:
  explicit CTPPSLHCInfoFilter(const edm::ParameterSet &);

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  bool filter(edm::StreamID, edm::Event &, const edm::EventSetup &) const override;

private:
  const edm::ESGetToken<LHCInfo, LHCInfoRcd> lhcInfoToken_;
  const edm::ESGetToken<LHCInfoPerLS, LHCInfoPerLSRcd> lhcInfoPerLSToken_;
  const edm::ESGetToken<LHCInfoPerFill, LHCInfoPerFillRcd> lhcInfoPerFillToken_;
  const bool useNewLHCInfo_;
  double minXangle_;
  double maxXangle_;
  double minBetaStarX_;
  double maxBetaStarX_;
  double minBetaStarY_;
  double maxBetaStarY_;
};

//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

CTPPSLHCInfoFilter::CTPPSLHCInfoFilter(const edm::ParameterSet &iConfig)
    : lhcInfoToken_(esConsumes(ESInputTag("", iConfig.getParameter<std::string>("lhcInfoLabel")))),
      lhcInfoPerLSToken_(esConsumes(ESInputTag("", iConfig.getParameter<std::string>("lhcInfoPerLSLabel")))),
      lhcInfoPerFillToken_(esConsumes(ESInputTag("", iConfig.getParameter<std::string>("lhcInfoPerFillLabel")))),
      useNewLHCInfo_(iConfig.getParameter<bool>("useNewLHCInfo")),
      minXangle_(iConfig.getParameter<double>("minCrossingAngle")),
      maxXangle_(iConfig.getParameter<double>("maxCrossingAngle")),
      minBetaStarX_(iConfig.getParameter<double>("minBetaStarX")),
      maxBetaStarX_(iConfig.getParameter<double>("maxBetaStarX")),
      minBetaStarY_(iConfig.getParameter<double>("minBetaStarY")),
      maxBetaStarY_(iConfig.getParameter<double>("maxBetaStarY")){}

void CTPPSLHCInfoFilter::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("lhcInfoLabel", "")->setComment("label of the LHCInfo record");
  desc.add<std::string>("lhcInfoPerLSLabel", "")->setComment("label of the LHCInfoPerLS record");
  desc.add<std::string>("lhcInfoPerFillLabel", "")->setComment("label of the LHCInfoPerFill record");
  desc.add<bool>("useNewLHCInfo", false)->setComment("flag whether to use new LHCInfoPer* records or old LHCInfo");
  desc.add<double>("minCrossingAngle", -1.)->setComment("minimum crossing angle (urad)");
  desc.add<double>("maxCrossingAngle", -1.)->setComment("maximum crossing angle (urad)");
  desc.add<double>("minBetaStarX", -1.)->setComment("minimum beta star X (m)");
  desc.add<double>("maxBetaStarX", -1.)->setComment("maximum beta star X (m)");
  desc.add<double>("minBetaStarY", -1.)->setComment("minimum beta star Y (m)");
  desc.add<double>("maxBetaStarY", -1.)->setComment("maximum beta star Y (m)");

  descriptions.add("ctppsLHCInfoFilterDefault", desc);
}

bool CTPPSLHCInfoFilter::filter(edm::StreamID, edm::Event &iEvent, const edm::EventSetup &iSetup) const {
  LHCInfoCombined lhcInfoCombined(iSetup, lhcInfoPerLSToken_, lhcInfoPerFillToken_, lhcInfoToken_, useNewLHCInfo_);

  double crossingAngle = lhcInfoCombined.crossingAngle();
  double betaStarX = lhcInfoCombined.betaStarX; 
  double betaStarY = lhcInfoCombined.betaStarY;

  if (minXangle_ >= 0. && crossingAngle < minXangle_) return false;
  if (maxXangle_ >= 0. && crossingAngle > maxXangle_) return false;
  if (minBetaStarX_ >= 0. && betaStarX < minBetaStarX_) return false;
  if (maxBetaStarX_ >= 0. && betaStarX > maxBetaStarX_) return false;
  if (minBetaStarY_ >= 0. && betaStarY < minBetaStarY_) return false;
  if (maxBetaStarY_ >= 0. && betaStarY > maxBetaStarY_) return false;

  return true;
}


DEFINE_FWK_MODULE(CTPPSLHCInfoFilter);
