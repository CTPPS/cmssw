/**
 * PPSProtonDumper - Dumps reconstructed proton kinematics to a flat ROOT TTree
 *
 * Reads reco::ForwardProton from ctppsProtons and writes the 5 kinematic variables:
 *   - vtx_x, vtx_y (vertex position, usually 0 for single-RP)
 *   - xi = Δp/p (fractional momentum loss)
 *   - th_x, th_y (scattering angles in rad)
 *
 * Additional variables:
 *   - t (four-momentum transfer squared, GeV^2)
 *   - p (momentum magnitude, GeV)
 *   - valid flag and reconstruction method
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/ProtonReco/interface/ForwardProton.h"
#include "DataFormats/ProtonReco/interface/ForwardProtonFwd.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "TTree.h"

#include <cmath>

class PPSProtonDumper : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit PPSProtonDumper(const edm::ParameterSet&);
  ~PPSProtonDumper() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  // Input tokens
  edm::EDGetTokenT<reco::ForwardProtonCollection> singleRPToken_;
  edm::EDGetTokenT<reco::ForwardProtonCollection> multiRPToken_;
  edm::EDGetTokenT<edm::HepMCProduct> hepMCToken_;
  bool useSingleRP_;
  bool useMultiRP_;

  TTree* tree_;

  // Tree branches - per event
  int eventNum_;
  int nProtons_;

  // Kinematic variables (vectors for multiple protons per event)
  std::vector<float> xi_;        // Δp/p (fractional momentum loss)
  std::vector<float> th_x_;      // horizontal scattering angle (rad)
  std::vector<float> th_y_;      // vertical scattering angle (rad)
  std::vector<float> vtx_x_;     // vertex x (mm)
  std::vector<float> vtx_y_;     // vertex y (mm)
  std::vector<float> t_;         // four-momentum transfer squared (GeV^2)
  std::vector<float> p_;         // momentum (GeV)
  std::vector<float> chi2_;      // fit chi2
  std::vector<int> ndf_;         // degrees of freedom
  std::vector<int> valid_;       // validity flag
  std::vector<int> method_;      // 0=invalid, 1=singleRP, 2=multiRP
  std::vector<int> arm_;         // 0=sector45, 1=sector56
  std::vector<int> station_;     // 0=210m (near), 2=220m (far)
};

PPSProtonDumper::PPSProtonDumper(const edm::ParameterSet& iConfig)
    : hepMCToken_(consumes<edm::HepMCProduct>(edm::InputTag("ppsG4NtupleToRecHitProducer", "unsmeared"))),
      useSingleRP_(iConfig.getParameter<bool>("useSingleRP")),
      useMultiRP_(iConfig.getParameter<bool>("useMultiRP")) {
  usesResource("TFileService");

  if (useSingleRP_) {
    singleRPToken_ = consumes<reco::ForwardProtonCollection>(
        iConfig.getParameter<edm::InputTag>("singleRPTag"));
  }
  if (useMultiRP_) {
    multiRPToken_ = consumes<reco::ForwardProtonCollection>(
        iConfig.getParameter<edm::InputTag>("multiRPTag"));
  }
}

void PPSProtonDumper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<bool>("useSingleRP", true)->setComment("Read single-RP protons");
  desc.add<bool>("useMultiRP", false)->setComment("Read multi-RP protons");
  desc.add<edm::InputTag>("singleRPTag", edm::InputTag("ctppsProtons", "singleRP"));
  desc.add<edm::InputTag>("multiRPTag", edm::InputTag("ctppsProtons", "multiRP"));
  descriptions.add("ppsProtonDumper", desc);
}

void PPSProtonDumper::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("protons", "PPS Reconstructed Protons");

  tree_->Branch("event", &eventNum_, "event/I");
  tree_->Branch("nProtons", &nProtons_, "nProtons/I");

  // Kinematic variables
  tree_->Branch("xi", &xi_);           // Δp/p
  tree_->Branch("th_x", &th_x_);       // θx (rad)
  tree_->Branch("th_y", &th_y_);       // θy (rad)
  tree_->Branch("vtx_x", &vtx_x_);     // vertex x (mm)
  tree_->Branch("vtx_y", &vtx_y_);     // vertex y (mm)
  tree_->Branch("t", &t_);             // -t (GeV^2)
  tree_->Branch("p", &p_);             // momentum (GeV)

  // Quality
  tree_->Branch("chi2", &chi2_);
  tree_->Branch("ndf", &ndf_);
  tree_->Branch("valid", &valid_);
  tree_->Branch("method", &method_);   // 1=singleRP, 2=multiRP
  tree_->Branch("arm", &arm_);         // 0=sector45, 1=sector56
  tree_->Branch("station", &station_); // 0=210m (near), 2=220m (far)
}

void PPSProtonDumper::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  // Clear vectors
  xi_.clear();
  th_x_.clear();
  th_y_.clear();
  vtx_x_.clear();
  vtx_y_.clear();
  t_.clear();
  p_.clear();
  chi2_.clear();
  ndf_.clear();
  valid_.clear();
  method_.clear();
  arm_.clear();
  station_.clear();

  // Get G4 event ID from HepMC (set by PPSG4NtupleToRecHitProducer)
  edm::Handle<edm::HepMCProduct> hepMC;
  iEvent.getByToken(hepMCToken_, hepMC);
  if (hepMC.isValid() && hepMC->GetEvent()) {
    eventNum_ = hepMC->GetEvent()->event_number();
  } else {
    eventNum_ = iEvent.id().event();  // Fallback to CMSSW event number
  }
  nProtons_ = 0;

  // Helper lambda to process a proton collection
  auto processProtons = [&](const reco::ForwardProtonCollection& protons, int methodFlag) {
    for (const auto& proton : protons) {
      xi_.push_back(proton.xi());
      th_x_.push_back(proton.thetaX());
      th_y_.push_back(proton.thetaY());
      vtx_x_.push_back(proton.vertex().x());
      vtx_y_.push_back(proton.vertex().y());
      t_.push_back(proton.t());
      // Momentum magnitude from the momentum vector
      const auto& mom = proton.momentum();
      p_.push_back(std::sqrt(mom.x()*mom.x() + mom.y()*mom.y() + mom.z()*mom.z()));
      chi2_.push_back(proton.chi2());
      ndf_.push_back(proton.ndof());
      valid_.push_back(proton.validFit() ? 1 : 0);
      method_.push_back(methodFlag);

      // Determine arm and station from contributing tracks
      // For single-RP, there's one track; for multi-RP, check first track
      int armVal = -1;
      int stationVal = -1;
      const auto& tracks = proton.contributingLocalTracks();
      if (!tracks.empty()) {
        const auto& trackRef = tracks[0];
        if (trackRef.isNonnull()) {
          // CTPPSDetId encodes arm, station, rp
          // arm = 0 for sector 45, arm = 1 for sector 56
          // station = 0 for 210m (near), station = 2 for 220m (far)
          CTPPSDetId detId(trackRef->rpId());
          armVal = detId.arm();
          stationVal = detId.station();
        }
      }
      arm_.push_back(armVal);
      station_.push_back(stationVal);

      nProtons_++;
    }
  };

  // Process single-RP protons
  if (useSingleRP_) {
    edm::Handle<reco::ForwardProtonCollection> singleRPProtons;
    iEvent.getByToken(singleRPToken_, singleRPProtons);
    if (singleRPProtons.isValid()) {
      processProtons(*singleRPProtons, 1);  // method=1 for single-RP
    }
  }

  // Process multi-RP protons
  if (useMultiRP_) {
    edm::Handle<reco::ForwardProtonCollection> multiRPProtons;
    iEvent.getByToken(multiRPToken_, multiRPProtons);
    if (multiRPProtons.isValid()) {
      processProtons(*multiRPProtons, 2);  // method=2 for multi-RP
    }
  }

  tree_->Fill();
}

DEFINE_FWK_MODULE(PPSProtonDumper);
