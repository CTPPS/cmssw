/**
 * PPSTrackDumper - Dumps track parameters to a simple flat ROOT TTree
 *
 * This analyzer reads CTPPSPixelLocalTracks and writes track parameters
 * (x0, y0, tx, ty) to a simple TTree that Python tools can easily read.
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSReco/interface/CTPPSPixelLocalTrack.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSPixelDetId.h"

#include "TTree.h"

class PPSTrackDumper : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit PPSTrackDumper(const edm::ParameterSet&);
  ~PPSTrackDumper() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<edm::DetSetVector<CTPPSPixelLocalTrack>> trackToken_;

  TTree* tree_;

  // Tree branches
  int eventNum_;
  int nTracks_;
  std::vector<float> x0_;
  std::vector<float> y0_;
  std::vector<float> tx_;
  std::vector<float> ty_;
  std::vector<float> chiSq_;
  std::vector<int> nHits_;
  std::vector<int> valid_;
};

PPSTrackDumper::PPSTrackDumper(const edm::ParameterSet& iConfig)
    : trackToken_(consumes<edm::DetSetVector<CTPPSPixelLocalTrack>>(
          iConfig.getParameter<edm::InputTag>("tracksTag"))) {
  usesResource("TFileService");
}

void PPSTrackDumper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracksTag", edm::InputTag("ctppsPixelLocalTracks"));
  descriptions.add("ppsTrackDumper", desc);
}

void PPSTrackDumper::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tracks", "PPS Reconstructed Tracks");

  tree_->Branch("event", &eventNum_, "event/I");
  tree_->Branch("nTracks", &nTracks_, "nTracks/I");
  tree_->Branch("x0", &x0_);
  tree_->Branch("y0", &y0_);
  tree_->Branch("tx", &tx_);
  tree_->Branch("ty", &ty_);
  tree_->Branch("chiSq", &chiSq_);
  tree_->Branch("nHits", &nHits_);
  tree_->Branch("valid", &valid_);
}

void PPSTrackDumper::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  // Clear vectors
  x0_.clear();
  y0_.clear();
  tx_.clear();
  ty_.clear();
  chiSq_.clear();
  nHits_.clear();
  valid_.clear();

  eventNum_ = iEvent.id().event();
  nTracks_ = 0;

  edm::Handle<edm::DetSetVector<CTPPSPixelLocalTrack>> tracks;
  iEvent.getByToken(trackToken_, tracks);

  if (tracks.isValid()) {
    for (const auto& ds : *tracks) {
      for (const auto& track : ds) {
        x0_.push_back(track.x0());
        y0_.push_back(track.y0());
        tx_.push_back(track.tx());
        ty_.push_back(track.ty());
        chiSq_.push_back(track.chiSquared());
        nHits_.push_back(track.ndf() + 4);  // ndf = nHits - 4 parameters
        valid_.push_back(track.isValid() ? 1 : 0);
        nTracks_++;
      }
    }
  }

  tree_->Fill();
}

DEFINE_FWK_MODULE(PPSTrackDumper);
