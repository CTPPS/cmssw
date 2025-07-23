/****************************************************************************
 *
 * This is a part of PPS offline software.
 * Authors:
 *   Laurent Forthomme (laurent.forthomme@cern.ch)
 *   Nicola Minafra (nicola.minafra@cern.ch)
 *
 ****************************************************************************/

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Common/interface/DetSetVector.h"

#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondRecHit.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondLocalTrack.h"

#include "RecoPPS/Local/interface/CTPPSDiamondTrackRecognition.h"
#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSGeometry.h"

class CTPPSDiamondLocalTrackFitter : public edm::stream::EDProducer<> {
public:
  explicit CTPPSDiamondLocalTrackFitter(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<edm::DetSetVector<CTPPSDiamondRecHit> > recHitsToken_;
  const edm::ParameterSet trk_algo_params_;
  std::unordered_map<CTPPSDetId, std::unique_ptr<CTPPSDiamondTrackRecognition> > trk_algo_;
  edm::ESGetToken<CTPPSGeometry, VeryForwardRealGeometryRecord> ctppsGeometryEventToken_;
};

CTPPSDiamondLocalTrackFitter::CTPPSDiamondLocalTrackFitter(const edm::ParameterSet& iConfig)
    : recHitsToken_(
          consumes<edm::DetSetVector<CTPPSDiamondRecHit> >(iConfig.getParameter<edm::InputTag>("recHitsTag"))),
          trk_algo_params_(iConfig.getParameter<edm::ParameterSet>("trackingAlgorithmParams")),
          ctppsGeometryEventToken_(esConsumes<CTPPSGeometry, VeryForwardRealGeometryRecord>()){
  produces<edm::DetSetVector<CTPPSDiamondLocalTrack> >();
}

void CTPPSDiamondLocalTrackFitter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // prepare the output
  auto pOutLocal = std::make_unique<edm::DetSetVector<CTPPSDiamondLocalTrack> >();
  auto pOutGlobal = std::make_unique<edm::DetSetVector<CTPPSDiamondLocalTrack> >();

  edm::Handle<edm::DetSetVector<CTPPSDiamondRecHit> > recHits;
  iEvent.getByToken(recHitsToken_, recHits);

  const CTPPSGeometry* ctppsGeometry = &iSetup.getData(ctppsGeometryEventToken_);

  // clear all hits possibly inherited from previous event
  for (auto& algo_vs_id : trk_algo_)
    algo_vs_id.second->clear();

  // feed hits to the track producers
  for (const auto& vec : *recHits) {
    const CTPPSDiamondDetId raw_detid(vec.detId()), detid(raw_detid.arm(), raw_detid.station(), raw_detid.rp());
    // if algorithm is not found, build it
    if (trk_algo_.count(detid) == 0)
      trk_algo_[detid] = std::make_unique<CTPPSDiamondTrackRecognition>(trk_algo_params_);
    for (const auto& hit : vec)
      // skip hits without a leading edge
      if (hit.ootIndex() != CTPPSDiamondRecHit::TIMESLICE_WITHOUT_LEADING){

        // create a copy of the hit to transform it to the local coordinates in X and Y
        CTPPSDiamondRecHit hitCopy(hit);
        auto localVector = CTPPSGeometry::Vector(hit.x(), hit.y(), hit.z());
        const auto diam = ctppsGeometry->sensor(detid);
        // do the global to local transformation
        // Global = Rotation * Local + Translation
        // Local = Rotation^-1 * (Global - Translation)
        localVector -= diam->translation();
        localVector = diam->rotation().Inverse() * localVector;
        hitCopy.setX(localVector.x());
        hitCopy.setY(localVector.y());
        // print X,Y,Z coordinates before and after transformation
        // std::cout << "Hit before transformation: X=" << hit.x() << ", Y=" << hit.y() << ", Z=" << hit.z() << std::endl;
        // std::cout << "Hit after transformation: X=" << localVector.x() << std::endl;
       
        trk_algo_[detid]->addHit(hitCopy);
      }
        
  }

  // build the local tracks for all stations
  for (auto& algo_vs_id : trk_algo_) {
    auto& tracks = pOutLocal->find_or_insert(algo_vs_id.first);
    algo_vs_id.second->produceTracks(tracks);
  }

  for (const auto& vec : *pOutLocal) {
    const CTPPSDiamondDetId raw_detid(vec.detId()), detid(raw_detid.arm(), raw_detid.station(), raw_detid.rp());
    const auto diam = ctppsGeometry->sensor(detid);
    pOutGlobal->find_or_insert(detid);
    for (const auto& track : vec) {
      auto trackCopy = track; // make a copy of the track to transform it
      auto globalVector = CTPPSGeometry::Vector(track.x0(), track.y0(), track.z0());

      // do the global to local transformation, doing first the rotation and then the translation
      globalVector = diam->rotation() * globalVector;
      globalVector += diam->translation();
      auto globalPoint =  math::XYZPoint(globalVector.x(), globalVector.y(), track.z0());

      trackCopy.setPosition(globalPoint);
      pOutGlobal->operator[](detid).push_back(trackCopy);
    }
  }

  // print the local tracks
  std::cout << "Local tracks for event " << iEvent.id() << ":" << std::endl;
  for (const auto& vec : *pOutLocal) {
    const CTPPSDiamondDetId raw_detid(vec.detId()), detid(raw_detid.arm(), raw_detid.station(), raw_detid.rp());
    std::cout << "  DetId: " << detid << std::endl;
    for (const auto& track : vec) {
      std::cout << "    Track: x0=" << track.x0() << ", y0=" << track.y0() << ", z0=" << track.z0() << std::endl;
    }
  }
  // print the global tracks
  std::cout << "Global tracks for event " << iEvent.id() << ":" << std::endl;
  for (const auto& vec : *pOutGlobal) {
    const CTPPSDiamondDetId raw_detid(vec.detId()), detid(raw_detid.arm(), raw_detid.station(), raw_detid.rp());
    std::cout << "  DetId: " << detid << std::endl;
    for (const auto& track : vec) {
      std::cout << "    Track: x0=" << track.x0() << ", y0=" << track.y0() << ", z0=" << track.z0() << std::endl;
    }
  }

  iEvent.put(std::move(pOutGlobal));
}

void CTPPSDiamondLocalTrackFitter::fillDescriptions(edm::ConfigurationDescriptions& descr) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("recHitsTag", edm::InputTag("ctppsDiamondRecHits"))
      ->setComment("input rechits collection to retrieve");

  edm::ParameterSetDescription trackingAlgoParams;
  trackingAlgoParams.add<double>("threshold", 1.5)
      ->setComment("minimal number of rechits to be observed before launching the track recognition algorithm");
  trackingAlgoParams.add<double>("thresholdFromMaximum", 0.5);
  trackingAlgoParams.add<double>("resolution", 0.01 /* mm */)
      ->setComment("spatial resolution on the horizontal coordinate (in mm)");
  trackingAlgoParams.add<double>("sigma", 0.1);
  trackingAlgoParams.add<double>("startFromX", -0.5 /* mm */)
      ->setComment("starting horizontal coordinate of rechits for the track recognition");
  trackingAlgoParams.add<double>("stopAtX", 19.5 /* mm */)
      ->setComment("ending horizontal coordinate of rechits for the track recognition");
  trackingAlgoParams.add<double>("tolerance", 0.1 /* mm */)
      ->setComment("tolerance used for checking if the track contains certain hit");

  trackingAlgoParams.add<std::string>("pixelEfficiencyFunction", "(x>[0]-0.5*[1])*(x<[0]+0.5*[1])+0*[2]")
      ->setComment(
          "efficiency function for single pixel\n"
          "can be defined as:\n"
          " * Precise: (TMath::Erf((x-[0]+0.5*[1])/([2]/4)+2)+1)*TMath::Erfc((x-[0]-0.5*[1])/([2]/4)-2)/4\n"
          " * Fast: "
          "(x>[0]-0.5*[1])*(x<[0]+0.5*[1])+((x-[0]+0.5*[1]+[2])/"
          "[2])*(x>[0]-0.5*[1]-[2])*(x<[0]-0.5*[1])+(2-(x-[0]-0.5*[1]+[2])/[2])*(x>[0]+0.5*[1])*(x<[0]+0.5*[1]+[2])\n"
          " * Legacy: (1/(1+exp(-(x-[0]+0.5*[1])/[2])))*(1/(1+exp((x-[0]-0.5*[1])/[2])))\n"
          " * Default (sigma ignored): (x>[0]-0.5*[1])*(x<[0]+0.5*[1])+0*[2]\n"
          "with:\n"
          "  [0]: centre of pad\n"
          "  [1]: width of pad\n"
          "  [2]: sigma: distance between efficiency ~100 -> 0 outside width");

  trackingAlgoParams.add<double>("yPosition", 0.0)->setComment("vertical offset of the outcoming track centre");
  trackingAlgoParams.add<double>("yWidth", 0.0)->setComment("vertical track width");
  trackingAlgoParams.add<bool>("excludeSingleEdgeHits", true)
      ->setComment("exclude rechits with missing leading/trailing edge");

  desc.add<edm::ParameterSetDescription>("trackingAlgorithmParams", trackingAlgoParams)
      ->setComment("list of parameters associated to the track recognition algorithm");

  descr.add("ctppsDiamondLocalTracks", desc);
}

DEFINE_FWK_MODULE(CTPPSDiamondLocalTrackFitter);
