/****************************************************************************
 * PPSG4NtupleToRecHitProducer
 *
 * Bridge between standalone Geant4 PPS simulation ROOT ntuples and the
 * CMSSW PPS reconstruction chain. Reads per-plane rechit trees and scoring
 * plane trees from the G4 simulation, produces:
 *   - edm::DetSetVector<CTPPSPixelRecHit> for local track fitting
 *   - edm::HepMCProduct("unsmeared") for MC truth correlation
 *
 * Designed to feed into recoDirectSimPPS (RecoPPS/Configuration/python/recoCTPPS_cff.py)
 * which starts from RecHits and runs local track fitting + proton reconstruction.
 *
 * Authors: PPS Collaboration / G4 Simulation Team
 ****************************************************************************/

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSReco/interface/CTPPSPixelRecHit.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSPixelDetId.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"

#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <cmath>

class PPSG4NtupleToRecHitProducer : public edm::one::EDProducer<> {
public:
  explicit PPSG4NtupleToRecHitProducer(const edm::ParameterSet&);
  ~PPSG4NtupleToRecHitProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  void beginJob() override;
  void endJob() override;

  // Build CTPPSPixelDetId from plane index (1-12 in G4 convention for dual station)
  // Planes 1-6 -> station 0 (210m near), Planes 7-12 -> station 2 (220m far)
  CTPPSPixelDetId buildDetId(int g4PlaneIndex) const;

  // Configuration parameters
  std::string inputFileName_;
  unsigned int arm_;      // 0 = sector 45, 1 = sector 56
  unsigned int station_;  // 0 = 210m, 2 = 220m (used only in single-station mode)
  unsigned int rp_;       // Roman Pot index within station
  bool dualStation_;      // true for multi-RP mode (both 210m and 220m stations)

  // Hit resolution (sigma, not variance) - matches PPSDirectProtonSimulation convention
  double sigmaHitX_;  // mm
  double sigmaHitY_;  // mm

  // Proton mass for HepMC
  static constexpr double protonMass_ = 0.93827208816;  // GeV/c^2

  // Number of pixel planes per RP (6 per station, 12 total for dual station)
  static constexpr int kPlanesPerStation = 6;
  static constexpr int kMaxPlanes = 12;

  // ROOT file and trees
  std::unique_ptr<TFile> rootFile_;

  // RecHit trees (one per plane, planes numbered 1-6)
  std::vector<TTree*> rechitTrees_;

  // Scoring plane tree (plane0 has the truth kinematics)
  TTree* scoringTree_;

  // Current entry counter (not event number - we read sequentially)
  Long64_t currentEntry_;
  Long64_t totalEntries_;

  // Branch variables for rechit trees
  Int_t rechit_eventID_;
  Int_t rechit_row_;
  Int_t rechit_col_;
  Double_t rechit_x_;
  Double_t rechit_y_;
  Double_t rechit_edep_;
  Int_t rechit_nhits_;

  // Branch variables for scoring plane tree (truth info)
  Double_t truth_x_;       // mm, position at scoring plane
  Double_t truth_y_;
  Double_t truth_theta_x_; // mrad, angle
  Double_t truth_theta_y_;
  Double_t truth_p_;       // GeV/c, momentum magnitude
  Double_t truth_energy_;  // GeV, primary energy
  Int_t truth_eventID_;

  // Cached data structures for event-based reading
  struct RecHitData {
    int plane;  // 1-6
    int row;
    int col;
    double x;   // mm, local
    double y;   // mm, local
    double edep;
    int nhits;
  };

  struct TruthData {
    double x, y;           // mm
    double theta_x, theta_y;  // mrad
    double p;              // GeV
    double energy;         // GeV
    int eventID;
  };

  // Maps from eventID to data
  std::map<int, std::vector<RecHitData>> rechitsByEvent_;
  std::map<int, TruthData> truthByEvent_;
  std::vector<int> eventIDList_;  // Ordered list of event IDs
  size_t currentEventIndex_;

  // Load all data from ROOT file into memory
  void loadAllData();

  unsigned int verbosity_;
};

//----------------------------------------------------------------------------------------------------

PPSG4NtupleToRecHitProducer::PPSG4NtupleToRecHitProducer(const edm::ParameterSet& iConfig)
    : inputFileName_(iConfig.getParameter<std::string>("inputFileName")),
      arm_(iConfig.getParameter<unsigned int>("arm")),
      station_(iConfig.getParameter<unsigned int>("station")),
      rp_(iConfig.getParameter<unsigned int>("rp")),
      dualStation_(iConfig.getParameter<bool>("dualStation")),
      sigmaHitX_(iConfig.getParameter<double>("sigmaHitX")),
      sigmaHitY_(iConfig.getParameter<double>("sigmaHitY")),
      scoringTree_(nullptr),
      currentEntry_(0),
      totalEntries_(0),
      currentEventIndex_(0),
      verbosity_(iConfig.getUntrackedParameter<unsigned int>("verbosity", 0)) {
  // Register products
  produces<edm::DetSetVector<CTPPSPixelRecHit>>();
  produces<edm::HepMCProduct>("unsmeared");

  // Initialize rechit tree pointers (6 for single station, 12 for dual station)
  int numPlanes = dualStation_ ? kMaxPlanes : kPlanesPerStation;
  rechitTrees_.resize(numPlanes, nullptr);
}

//----------------------------------------------------------------------------------------------------

PPSG4NtupleToRecHitProducer::~PPSG4NtupleToRecHitProducer() {}

//----------------------------------------------------------------------------------------------------

void PPSG4NtupleToRecHitProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("inputFileName", "g4_output.root")
      ->setComment("Path to the G4 simulation ROOT ntuple");

  desc.add<unsigned int>("arm", 0)->setComment("PPS arm: 0 = sector 45, 1 = sector 56");
  desc.add<unsigned int>("station", 0)->setComment("Station: 0 = 210m, 2 = 220m (ignored if dualStation=true)");
  desc.add<unsigned int>("rp", 3)->setComment("Roman Pot index within station");
  desc.add<bool>("dualStation", false)
      ->setComment("Enable dual-station (multi-RP) mode: planes 1-6 -> station 0, planes 7-12 -> station 2");

  // Default hit resolution: pitch / sqrt(12)
  // Pixel pitch: 100 um (X) x 150 um (Y) -> sigma ~ 29 um (X), 43 um (Y)
  desc.add<double>("sigmaHitX", 0.029)->setComment("Hit X resolution in mm (sigma, not variance)");
  desc.add<double>("sigmaHitY", 0.043)->setComment("Hit Y resolution in mm (sigma, not variance)");

  desc.addUntracked<unsigned int>("verbosity", 0);

  descriptions.add("ppsG4NtupleToRecHitProducer", desc);
}

//----------------------------------------------------------------------------------------------------

void PPSG4NtupleToRecHitProducer::beginJob() {
  // Open ROOT file
  rootFile_ = std::make_unique<TFile>(inputFileName_.c_str(), "READ");
  if (!rootFile_ || rootFile_->IsZombie()) {
    throw cms::Exception("FileOpenError")
        << "Could not open G4 ntuple file: " << inputFileName_;
  }

  edm::LogInfo("PPSG4NtupleToRecHitProducer") << "Opened G4 ntuple: " << inputFileName_;

  // Get scoring plane tree (plane0 contains truth info)
  scoringTree_ = dynamic_cast<TTree*>(rootFile_->Get("plane0"));
  if (!scoringTree_) {
    throw cms::Exception("TreeNotFound") << "Scoring plane tree 'plane0' not found in " << inputFileName_;
  }

  // Set up scoring tree branches
  scoringTree_->SetBranchAddress("x", &truth_x_);
  scoringTree_->SetBranchAddress("y", &truth_y_);
  scoringTree_->SetBranchAddress("theta_x", &truth_theta_x_);
  scoringTree_->SetBranchAddress("theta_y", &truth_theta_y_);
  scoringTree_->SetBranchAddress("p", &truth_p_);
  scoringTree_->SetBranchAddress("primary_energy", &truth_energy_);
  scoringTree_->SetBranchAddress("eventID", &truth_eventID_);

  // Get rechit trees (planes numbered 1-6 for single station, 1-12 for dual station)
  int numPlanes = dualStation_ ? kMaxPlanes : kPlanesPerStation;
  for (int i = 0; i < numPlanes; ++i) {
    std::string treeName = "rechits_plane" + std::to_string(i + 1);
    rechitTrees_[i] = dynamic_cast<TTree*>(rootFile_->Get(treeName.c_str()));
    if (!rechitTrees_[i]) {
      edm::LogWarning("PPSG4NtupleToRecHitProducer")
          << "RecHit tree '" << treeName << "' not found - plane will have no hits";
    } else {
      rechitTrees_[i]->SetBranchAddress("eventID", &rechit_eventID_);
      rechitTrees_[i]->SetBranchAddress("row", &rechit_row_);
      rechitTrees_[i]->SetBranchAddress("col", &rechit_col_);
      rechitTrees_[i]->SetBranchAddress("x", &rechit_x_);
      rechitTrees_[i]->SetBranchAddress("y", &rechit_y_);
      rechitTrees_[i]->SetBranchAddress("edep", &rechit_edep_);
      rechitTrees_[i]->SetBranchAddress("nhits", &rechit_nhits_);
    }
  }

  if (dualStation_) {
    edm::LogInfo("PPSG4NtupleToRecHitProducer")
        << "Dual-station (multi-RP) mode enabled: arm=" << arm_
        << ", planes 1-6 -> station 0 (210m), planes 7-12 -> station 2 (220m)";
  }

  // Load all data into memory for event-based access
  loadAllData();

  edm::LogInfo("PPSG4NtupleToRecHitProducer")
      << "Loaded " << eventIDList_.size() << " events from G4 ntuple";
}

//----------------------------------------------------------------------------------------------------

void PPSG4NtupleToRecHitProducer::loadAllData() {
  // Load truth data from scoring plane
  Long64_t nTruth = scoringTree_->GetEntries();
  for (Long64_t i = 0; i < nTruth; ++i) {
    scoringTree_->GetEntry(i);
    TruthData td;
    td.x = truth_x_;
    td.y = truth_y_;
    td.theta_x = truth_theta_x_;
    td.theta_y = truth_theta_y_;
    td.p = truth_p_;
    td.energy = truth_energy_;
    td.eventID = truth_eventID_;
    truthByEvent_[truth_eventID_] = td;
  }

  // Load rechit data from all planes (6 for single station, 12 for dual station)
  int numPlanes = dualStation_ ? kMaxPlanes : kPlanesPerStation;
  for (int planeIdx = 0; planeIdx < numPlanes; ++planeIdx) {
    if (!rechitTrees_[planeIdx])
      continue;

    Long64_t nRechits = rechitTrees_[planeIdx]->GetEntries();
    for (Long64_t i = 0; i < nRechits; ++i) {
      rechitTrees_[planeIdx]->GetEntry(i);
      RecHitData rd;
      rd.plane = planeIdx + 1;  // G4 planes are 1-6 (station 0) or 7-12 (station 2)
      rd.row = rechit_row_;
      rd.col = rechit_col_;
      rd.x = rechit_x_;
      rd.y = rechit_y_;
      rd.edep = rechit_edep_;
      rd.nhits = rechit_nhits_;
      rechitsByEvent_[rechit_eventID_].push_back(rd);
    }
  }

  // Build ordered list of event IDs (union of truth and rechit events)
  std::set<int> allEvents;
  for (const auto& kv : truthByEvent_)
    allEvents.insert(kv.first);
  for (const auto& kv : rechitsByEvent_)
    allEvents.insert(kv.first);

  eventIDList_.assign(allEvents.begin(), allEvents.end());
  totalEntries_ = eventIDList_.size();

}

//----------------------------------------------------------------------------------------------------

void PPSG4NtupleToRecHitProducer::endJob() {
  if (rootFile_) {
    rootFile_->Close();
    rootFile_.reset();
  }
}

//----------------------------------------------------------------------------------------------------

CTPPSPixelDetId PPSG4NtupleToRecHitProducer::buildDetId(int g4PlaneIndex) const {
  if (dualStation_) {
    // Dual station (multi-RP) mode:
    //   Planes 1-6  -> station 0 (210m near), CMSSW plane 0-5
    //   Planes 7-12 -> station 2 (220m far),  CMSSW plane 0-5
    if (g4PlaneIndex <= kPlanesPerStation) {
      unsigned int plane = static_cast<unsigned int>(g4PlaneIndex - 1);
      return CTPPSPixelDetId(arm_, 0, rp_, plane);
    } else {
      unsigned int plane = static_cast<unsigned int>(g4PlaneIndex - kPlanesPerStation - 1);
      return CTPPSPixelDetId(arm_, 2, rp_, plane);
    }
  } else {
    // Single station mode: use configured station parameter
    // g4PlaneIndex: 1-6 from the G4 simulation -> CMSSW plane index: 0-5
    unsigned int plane = static_cast<unsigned int>(g4PlaneIndex - 1);
    return CTPPSPixelDetId(arm_, station_, rp_, plane);
  }
}

//----------------------------------------------------------------------------------------------------

void PPSG4NtupleToRecHitProducer::produce(edm::Event& iEvent, const edm::EventSetup&) {
  // Check if we have more events to process
  if (currentEventIndex_ >= eventIDList_.size()) {
    throw cms::Exception("EndOfFile")
        << "Ran out of events in G4 ntuple after " << currentEventIndex_ << " events";
  }

  int eventID = eventIDList_[currentEventIndex_];
  ++currentEventIndex_;

  if (verbosity_ > 0) {
    edm::LogInfo("PPSG4NtupleToRecHitProducer")
        << "Processing G4 eventID " << eventID << " (index " << (currentEventIndex_ - 1) << ")";
  }

  // Create output collections
  auto pixelRecHits = std::make_unique<edm::DetSetVector<CTPPSPixelRecHit>>();

  // Fill rechits for this event
  auto rechitIt = rechitsByEvent_.find(eventID);
  if (rechitIt != rechitsByEvent_.end()) {
    for (const auto& hit : rechitIt->second) {
      CTPPSPixelDetId detId = buildDetId(hit.plane);

      // Apply coordinate corrections to align G4 geometry with CMSSW geometry.
      //
      // G4 rechit coordinates are "wafer-centered, world-aligned":
      //   hit.x = dispersion direction (varies with ξ)
      //   hit.y = vertical direction (~0 when θ_y=0)
      //
      // CMSSW local frame has a 90° rotation relative to G4:
      //   local_x → track_y0 (vertical)
      //   local_y → track_x0 (dispersion, used for ξ reconstruction)
      //
      // Therefore we must swap: local_x = hit.y, local_y = -hit.x
      //
      // Additionally, station 0 (210m) has an 8° physical tilt in CMSSW geometry
      // that requires a rotation correction.

      double lp_x, lp_y;

      // CMSSW geometry: local_y → track_x0, local_x → track_y0
      // We need: rechit_x (dispersion) → track_x0 for single-RP ξ reco
      //
      // Station 0 (210m) has 8° physical tilt encoded in CMSSW geometry.
      // When CMSSW transforms local → global, it applies +8° rotation.
      // We must PRE-apply -8° rotation to compensate, so that:
      //   final track_y0 = hit.y (vertical, not mixed with dispersion)
      //   final track_x0 ∝ hit.x (dispersion)
      //
      // Station 2 (220m) has no tilt - just apply the 90° swap.

      if (hit.plane <= kPlanesPerStation) {
        // Station 0 (210m): pre-rotate by -8° to counteract CMSSW's +8° rotation
        constexpr double cos8 = 0.99027;  // cos(8°)
        constexpr double sin8 = 0.13917;  // sin(8°)
        constexpr double st0_y_offset = +1.40;  // mm, align y0 with station 2

        // After swapping: local_x = hit.y (vert), local_y = hit.x (disp)
        // CMSSW applies +8°: track_y0 = cos8*lp_x - sin8*lp_y
        //                    track_x0 = sin8*lp_x + cos8*lp_y
        //
        // We want: track_y0 = hit.y, track_x0 ∝ hit.x
        // Pre-apply -8° rotation: lp_x' = cos8*lp_x + sin8*lp_y
        //                         lp_y' = -sin8*lp_x + cos8*lp_y
        //
        // After CMSSW's +8°: track_y0 = cos8*(cos8*lp_x + sin8*lp_y) - sin8*(-sin8*lp_x + cos8*lp_y)
        //                             = (cos²8 + sin²8)*lp_x = lp_x = hit.y ✓

        double swapped_x = hit.y;  // vertical → local_x before pre-rotation
        double swapped_y = hit.x;  // dispersion → local_y before pre-rotation

        // Pre-rotate by -8° and add y0 alignment offset
        lp_x =  cos8 * swapped_x + sin8 * swapped_y + st0_y_offset;
        lp_y = -sin8 * swapped_x + cos8 * swapped_y;
      } else {
        // Station 2 (220m): no tilt, just swap
        lp_x = hit.y;   // vertical → local_x → track_y0
        lp_y = hit.x;   // dispersion → local_y → track_x0
      }

      LocalPoint lp(lp_x, lp_y, 0.);
      LocalError le(sigmaHitX_, 0., sigmaHitY_);  // (sigmaX, 0, sigmaY) - same as PPSDirectProtonSimulation

      // Get or create DetSet for this detector ID
      edm::DetSet<CTPPSPixelRecHit>& ds = pixelRecHits->find_or_insert(detId);
      ds.emplace_back(lp, le);

      if (verbosity_ > 1) {
        edm::LogInfo("PPSG4NtupleToRecHitProducer")
            << "  RecHit: plane=" << hit.plane << " detId=" << detId.rawId()
            << " x=" << hit.x << " y=" << hit.y
            << " row=" << hit.row << " col=" << hit.col;
      }
    }
  }

  // Create HepMC product with truth information
  auto hepMCProduct = std::make_unique<edm::HepMCProduct>();

  auto truthIt = truthByEvent_.find(eventID);
  if (truthIt != truthByEvent_.end()) {
    const TruthData& truth = truthIt->second;

    // Create HepMC event (following IOMC/ParticleGuns pattern)
    auto* genEvent = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::MM);
    genEvent->set_event_number(eventID);

    // Create vertex at origin (the scoring plane position is encoded in the particle momentum)
    // In a more sophisticated version, we could use the actual vertex position
    auto* genVertex = new HepMC::GenVertex(HepMC::FourVector(0., 0., 0., 0.));
    genEvent->add_vertex(genVertex);

    // Calculate momentum components from angles and magnitude
    // theta_x, theta_y are in mrad, p is in GeV
    // For small angles: px = p * theta_x, py = p * theta_y
    double theta_x_rad = truth.theta_x * 1e-3;  // mrad -> rad
    double theta_y_rad = truth.theta_y * 1e-3;
    double px = truth.p * theta_x_rad;
    double py = truth.p * theta_y_rad;
    // pz from normalization: p^2 = px^2 + py^2 + pz^2
    double pz_sq = truth.p * truth.p - px * px - py * py;
    double pz = (pz_sq > 0) ? std::sqrt(pz_sq) : 0.;

    // Determine z-direction based on arm
    // arm 0 (sector 45): proton travels in -z direction
    // arm 1 (sector 56): proton travels in +z direction
    if (arm_ == 0)
      pz = -pz;

    // Energy from momentum and proton mass
    double E = std::sqrt(truth.p * truth.p + protonMass_ * protonMass_);

    // Create proton (pdgId = 2212, status = 1 for stable)
    auto* genParticle = new HepMC::GenParticle(HepMC::FourVector(px, py, pz, E), 2212, 1);
    genVertex->add_particle_out(genParticle);

    hepMCProduct->addHepMCData(genEvent);

    if (verbosity_ > 0) {
      edm::LogInfo("PPSG4NtupleToRecHitProducer")
          << "  Truth: p=" << truth.p << " GeV, theta_x=" << truth.theta_x
          << " mrad, theta_y=" << truth.theta_y << " mrad"
          << " -> px=" << px << " py=" << py << " pz=" << pz;
    }
  } else {
    // No truth for this event - create empty HepMC event
    auto* genEvent = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::MM);
    genEvent->set_event_number(eventID);
    hepMCProduct->addHepMCData(genEvent);
  }

  // Put products into event
  iEvent.put(std::move(pixelRecHits));
  iEvent.put(std::move(hepMCProduct), "unsmeared");
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(PPSG4NtupleToRecHitProducer);
