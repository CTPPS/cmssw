/****************************************************************************
 * PPSG4NtupleValidator
 *
 * Simple analyzer to validate the G4 ntuple -> RecHit -> LocalTrack chain.
 * Compares fitted track parameters against the MC truth from HepMCProduct.
 ****************************************************************************/

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSReco/interface/CTPPSPixelRecHit.h"
#include "DataFormats/CTPPSReco/interface/CTPPSPixelLocalTrack.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSPixelDetId.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "TH1D.h"
#include "TH2D.h"

#include <cmath>

class PPSG4NtupleValidator : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit PPSG4NtupleValidator(const edm::ParameterSet&);
  ~PPSG4NtupleValidator() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // Input tokens
  edm::EDGetTokenT<edm::DetSetVector<CTPPSPixelRecHit>> rechitToken_;
  edm::EDGetTokenT<edm::DetSetVector<CTPPSPixelLocalTrack>> localTrackToken_;
  edm::EDGetTokenT<std::vector<CTPPSLocalTrackLite>> trackLiteToken_;
  edm::EDGetTokenT<edm::HepMCProduct> hepMCToken_;

  // Histograms
  TH1D* h_nRecHits_;
  TH1D* h_nTracks_;
  TH1D* h_nPlanesWithHits_;

  TH1D* h_rechit_x_;
  TH1D* h_rechit_y_;

  TH1D* h_track_x0_;
  TH1D* h_track_y0_;
  TH1D* h_track_tx_;
  TH1D* h_track_ty_;

  TH1D* h_truth_px_;
  TH1D* h_truth_py_;
  TH1D* h_truth_pz_;

  // Resolution histograms (reco - truth)
  TH1D* h_delta_tx_;
  TH1D* h_delta_ty_;

  // 2D correlations
  TH2D* h2_reco_vs_truth_tx_;
  TH2D* h2_reco_vs_truth_ty_;

  unsigned int verbosity_;
  int nEventsProcessed_;
  int nEventsWithTracks_;
  int nEventsWithTruth_;
};

//----------------------------------------------------------------------------------------------------

PPSG4NtupleValidator::PPSG4NtupleValidator(const edm::ParameterSet& iConfig)
    : rechitToken_(consumes<edm::DetSetVector<CTPPSPixelRecHit>>(
          iConfig.getParameter<edm::InputTag>("rechitTag"))),
      localTrackToken_(consumes<edm::DetSetVector<CTPPSPixelLocalTrack>>(
          iConfig.getParameter<edm::InputTag>("localTrackTag"))),
      trackLiteToken_(consumes<std::vector<CTPPSLocalTrackLite>>(
          iConfig.getParameter<edm::InputTag>("trackLiteTag"))),
      hepMCToken_(consumes<edm::HepMCProduct>(
          iConfig.getParameter<edm::InputTag>("hepMCTag"))),
      verbosity_(iConfig.getUntrackedParameter<unsigned int>("verbosity", 0)),
      nEventsProcessed_(0),
      nEventsWithTracks_(0),
      nEventsWithTruth_(0) {
  usesResource("TFileService");
}

//----------------------------------------------------------------------------------------------------

void PPSG4NtupleValidator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("rechitTag", edm::InputTag("ppsG4NtupleToRecHitProducer"));
  desc.add<edm::InputTag>("localTrackTag", edm::InputTag("ctppsPixelLocalTracks"));
  desc.add<edm::InputTag>("trackLiteTag", edm::InputTag("ctppsLocalTrackLiteProducer"));
  desc.add<edm::InputTag>("hepMCTag", edm::InputTag("ppsG4NtupleToRecHitProducer", "unsmeared"));

  desc.addUntracked<unsigned int>("verbosity", 0);

  descriptions.add("ppsG4NtupleValidator", desc);
}

//----------------------------------------------------------------------------------------------------

void PPSG4NtupleValidator::beginJob() {
  edm::Service<TFileService> fs;

  // Multiplicity histograms
  h_nRecHits_ = fs->make<TH1D>("nRecHits", "Number of RecHits per event;N_{rechits};Events", 50, 0, 50);
  h_nTracks_ = fs->make<TH1D>("nTracks", "Number of tracks per event;N_{tracks};Events", 10, 0, 10);
  h_nPlanesWithHits_ = fs->make<TH1D>("nPlanesWithHits", "Number of planes with hits;N_{planes};Events", 8, 0, 8);

  // RecHit positions
  h_rechit_x_ = fs->make<TH1D>("rechit_x", "RecHit X position;x (mm);Entries", 100, -10, 10);
  h_rechit_y_ = fs->make<TH1D>("rechit_y", "RecHit Y position;y (mm);Entries", 100, -10, 10);

  // Track parameters
  h_track_x0_ = fs->make<TH1D>("track_x0", "Track intercept X;x_0 (mm);Tracks", 100, -10, 10);
  h_track_y0_ = fs->make<TH1D>("track_y0", "Track intercept Y;y_0 (mm);Tracks", 100, -10, 10);
  h_track_tx_ = fs->make<TH1D>("track_tx", "Track slope X;t_x;Tracks", 100, -0.02, 0.02);
  h_track_ty_ = fs->make<TH1D>("track_ty", "Track slope Y;t_y;Tracks", 100, -0.02, 0.02);

  // Truth momentum
  h_truth_px_ = fs->make<TH1D>("truth_px", "Truth p_x;p_x (GeV);Events", 100, -100, 100);
  h_truth_py_ = fs->make<TH1D>("truth_py", "Truth p_y;p_y (GeV);Events", 100, -100, 100);
  h_truth_pz_ = fs->make<TH1D>("truth_pz", "Truth p_z;p_z (GeV);Events", 100, -7000, 7000);

  // Resolution
  h_delta_tx_ = fs->make<TH1D>("delta_tx", "Track t_x resolution;t_x^{reco} - t_x^{truth};Entries", 100, -0.01, 0.01);
  h_delta_ty_ = fs->make<TH1D>("delta_ty", "Track t_y resolution;t_y^{reco} - t_y^{truth};Entries", 100, -0.01, 0.01);

  // 2D correlations
  h2_reco_vs_truth_tx_ = fs->make<TH2D>("reco_vs_truth_tx",
      "Track t_x: reco vs truth;t_x^{truth};t_x^{reco}", 50, -0.02, 0.02, 50, -0.02, 0.02);
  h2_reco_vs_truth_ty_ = fs->make<TH2D>("reco_vs_truth_ty",
      "Track t_y: reco vs truth;t_y^{truth};t_y^{reco}", 50, -0.02, 0.02, 50, -0.02, 0.02);
}

//----------------------------------------------------------------------------------------------------

void PPSG4NtupleValidator::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  ++nEventsProcessed_;

  // Get RecHits
  edm::Handle<edm::DetSetVector<CTPPSPixelRecHit>> rechits;
  iEvent.getByToken(rechitToken_, rechits);

  int nRecHits = 0;
  std::set<unsigned int> planesWithHits;
  if (rechits.isValid()) {
    for (const auto& ds : *rechits) {
      CTPPSPixelDetId detId(ds.detId());
      planesWithHits.insert(detId.plane());
      for (const auto& hit : ds) {
        ++nRecHits;
        h_rechit_x_->Fill(hit.point().x());
        h_rechit_y_->Fill(hit.point().y());
      }
    }
  }
  h_nRecHits_->Fill(nRecHits);
  h_nPlanesWithHits_->Fill(planesWithHits.size());

  // Get local tracks
  edm::Handle<edm::DetSetVector<CTPPSPixelLocalTrack>> localTracks;
  iEvent.getByToken(localTrackToken_, localTracks);

  int nTracks = 0;
  double reco_tx = 0, reco_ty = 0;
  bool hasTrack = false;

  if (localTracks.isValid()) {
    for (const auto& ds : *localTracks) {
      for (const auto& track : ds) {
        if (!track.isValid())
          continue;
        ++nTracks;
        hasTrack = true;

        h_track_x0_->Fill(track.x0());
        h_track_y0_->Fill(track.y0());
        h_track_tx_->Fill(track.tx());
        h_track_ty_->Fill(track.ty());

        reco_tx = track.tx();
        reco_ty = track.ty();

        if (verbosity_ > 1) {
          edm::LogInfo("PPSG4NtupleValidator")
              << "Track: x0=" << track.x0() << " y0=" << track.y0()
              << " tx=" << track.tx() << " ty=" << track.ty()
              << " chi2=" << track.chiSquared();
        }
      }
    }
  }
  h_nTracks_->Fill(nTracks);
  if (nTracks > 0)
    ++nEventsWithTracks_;

  // Get HepMC truth
  edm::Handle<edm::HepMCProduct> hepMC;
  iEvent.getByToken(hepMCToken_, hepMC);

  double truth_tx = 0, truth_ty = 0;
  bool hasTruth = false;

  if (hepMC.isValid()) {
    const HepMC::GenEvent* evt = hepMC->GetEvent();
    if (evt) {
      for (auto it = evt->particles_begin(); it != evt->particles_end(); ++it) {
        const HepMC::GenParticle* p = *it;
        if (p && p->pdg_id() == 2212 && p->status() == 1) {
          double px = p->momentum().px();
          double py = p->momentum().py();
          double pz = p->momentum().pz();

          h_truth_px_->Fill(px);
          h_truth_py_->Fill(py);
          h_truth_pz_->Fill(pz);

          // Calculate truth angles
          if (std::abs(pz) > 1e-6) {
            truth_tx = px / std::abs(pz);
            truth_ty = py / std::abs(pz);
            hasTruth = true;
          }

          if (verbosity_ > 1) {
            edm::LogInfo("PPSG4NtupleValidator")
                << "Truth: px=" << px << " py=" << py << " pz=" << pz
                << " -> tx=" << truth_tx << " ty=" << truth_ty;
          }
        }
      }
    }
  }
  if (hasTruth)
    ++nEventsWithTruth_;

  // Fill resolution histograms if we have both track and truth
  if (hasTrack && hasTruth) {
    h_delta_tx_->Fill(reco_tx - truth_tx);
    h_delta_ty_->Fill(reco_ty - truth_ty);
    h2_reco_vs_truth_tx_->Fill(truth_tx, reco_tx);
    h2_reco_vs_truth_ty_->Fill(truth_ty, reco_ty);
  }
}

//----------------------------------------------------------------------------------------------------

void PPSG4NtupleValidator::endJob() {
  edm::LogInfo("PPSG4NtupleValidator")
      << "Summary:\n"
      << "  Events processed: " << nEventsProcessed_ << "\n"
      << "  Events with tracks: " << nEventsWithTracks_ << "\n"
      << "  Events with truth: " << nEventsWithTruth_;
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(PPSG4NtupleValidator);
