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
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

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
  edm::EDGetTokenT<edm::HepMCProduct> hepMCToken_;

  TTree* tree_;

  // Tree branches
  int eventNum_;
  int nTracks_;
  std::vector<int> arm_;     // PPS arm (0=sector45, 1=sector56)
  std::vector<int> station_; // station (0=210m, 2=220m)
  std::vector<float> x0_;
  std::vector<float> y0_;
  std::vector<float> z0_;  // track reference z (RP translation)
  std::vector<float> tx_;
  std::vector<float> ty_;
  std::vector<float> chiSq_;
  std::vector<int> nHits_;
  std::vector<int> valid_;
};

PPSTrackDumper::PPSTrackDumper(const edm::ParameterSet& iConfig)
    : trackToken_(consumes<edm::DetSetVector<CTPPSPixelLocalTrack>>(
          iConfig.getParameter<edm::InputTag>("tracksTag"))),
      hepMCToken_(consumes<edm::HepMCProduct>(edm::InputTag("ppsG4NtupleToRecHitProducer", "unsmeared"))) {
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
  tree_->Branch("arm", &arm_);
  tree_->Branch("station", &station_);
  tree_->Branch("x0", &x0_);
  tree_->Branch("y0", &y0_);
  tree_->Branch("z0", &z0_);  // track reference z (mm)
  tree_->Branch("tx", &tx_);
  tree_->Branch("ty", &ty_);
  tree_->Branch("chiSq", &chiSq_);
  tree_->Branch("nHits", &nHits_);
  tree_->Branch("valid", &valid_);
}

void PPSTrackDumper::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  // Clear vectors
  arm_.clear();
  station_.clear();
  x0_.clear();
  y0_.clear();
  z0_.clear();
  tx_.clear();
  ty_.clear();
  chiSq_.clear();
  nHits_.clear();
  valid_.clear();

  // Get G4 event ID from HepMC (set by PPSG4NtupleToRecHitProducer)
  edm::Handle<edm::HepMCProduct> hepMC;
  iEvent.getByToken(hepMCToken_, hepMC);
  if (hepMC.isValid() && hepMC->GetEvent()) {
    eventNum_ = hepMC->GetEvent()->event_number();
  } else {
    eventNum_ = iEvent.id().event();  // Fallback to CMSSW event number
  }
  nTracks_ = 0;

  edm::Handle<edm::DetSetVector<CTPPSPixelLocalTrack>> tracks;
  iEvent.getByToken(trackToken_, tracks);

  if (tracks.isValid()) {
    for (const auto& ds : *tracks) {
      // Extract arm and station from DetId
      CTPPSPixelDetId detId(ds.detId());
      int arm = detId.arm();
      int sta = detId.station();

      for (const auto& track : ds) {
        arm_.push_back(arm);
        station_.push_back(sta);
        x0_.push_back(track.x0());
        y0_.push_back(track.y0());
        z0_.push_back(track.z0());  // track reference z in mm
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
