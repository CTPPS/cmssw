#include <cstdlib>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>
#include <TVirtualPad.h>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/join.hpp>
#include <boost/algorithm/string/join.hpp>

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondRecHit.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"

void calculateAndSaveHistograms(int maxEvents_,
                                unsigned int outputEvery_,
                                double totCut_,
                                const std::string& outputFile_,
                                int mode_,
                                const std::vector<std::string>& inFiles_,
                                int useOnlyPixelOn_ = 0,
                                float maxXDiff = 20.0,
                                std::optional<std::vector<std::pair<int, int>>> goodLumisections_ = std::nullopt,
                                const std::set<int>& pickedBunches_ = {}) {
  // ----------------------------------------------------------------------
  // First Part:
  //
  //  * enable FWLite
  //  * book the histograms of interest
  //  * open the input file
  // ----------------------------------------------------------------------

  // book a set of histograms
  fwlite::TFileService fs = fwlite::TFileService(outputFile_);
  TFileDirectory dir = fs.mkdir("diamondHistograms");

  // Pixel-diamond x-correlation to define matching cuts
  TH2F* xtpix21045_ = dir.make<TH2F>("xtpix21045", "xtpix21045", 40, -10, 30, 40, -10, 30);
  TH2F* xtpix21056_ = dir.make<TH2F>("xtpix21056", "xtpix21056", 40, -10, 30, 40, -10, 30);

  // Denominator for all efficiencies
  TH2F* deffden45_ = dir.make<TH2F>("deffden45", "deffden45", 220, -2, 20, 240, -8, 16);
  TH2F* deffden56_ = dir.make<TH2F>("deffden56", "deffden56", 220, -2, 20, 240, -8, 16);

  // Numerator for per-arm/track efficiencies
  TH2F* deffnum45_ = dir.make<TH2F>("deffnum45", "deffnum45", 220, -2, 20, 240, -8, 16);
  TH2F* deffnum56_ = dir.make<TH2F>("deffnum56", "deffnum56", 220, -2, 20, 240, -8, 16);

  TH2F* dboxeffnum45_ = dir.make<TH2F>("dboxeffnum45", "dboxeffnum45", 220, -2, 20, 240, -8, 16);
  TH2F* dboxeffnum56_ = dir.make<TH2F>("dboxeffnum56", "dboxeffnum56", 220, -2, 20, 240, -8, 16);

  // Numerator for per-plane efficiencies, cylindrical RPs
  TH2F* deffnum45plane0_ = dir.make<TH2F>("deffnum45plane0", "deffnum45plane0", 220, -2, 20, 240, -8, 16);
  TH2F* deffnum56plane0_ = dir.make<TH2F>("deffnum56plane0", "deffnum56plane0", 220, -2, 20, 240, -8, 16);
  TH2F* deffnum45plane1_ = dir.make<TH2F>("deffnum45plane1", "deffnum45plane1", 220, -2, 20, 240, -8, 16);
  TH2F* deffnum56plane1_ = dir.make<TH2F>("deffnum56plane1", "deffnum56plane1", 220, -2, 20, 240, -8, 16);
  TH2F* deffnum45plane2_ = dir.make<TH2F>("deffnum45plane2", "deffnum45plane2", 220, -2, 20, 240, -8, 16);
  TH2F* deffnum56plane2_ = dir.make<TH2F>("deffnum56plane2", "deffnum56plane2", 220, -2, 20, 240, -8, 16);
  TH2F* deffnum45plane3_ = dir.make<TH2F>("deffnum45plane3", "deffnum45plane3", 220, -2, 20, 240, -8, 16);
  TH2F* deffnum56plane3_ = dir.make<TH2F>("deffnum56plane3", "deffnum56plane3", 220, -2, 20, 240, -8, 16);

  // Numerator for per-plane efficiencies, box RPs
  TH2F* deffnum45boxplane0_ = dir.make<TH2F>("deffnum45boxplane0", "deffnum45boxplane0", 220, -2, 20, 240, -8, 16);
  TH2F* deffnum56boxplane0_ = dir.make<TH2F>("deffnum56boxplane0", "deffnum56boxplane0", 220, -2, 20, 240, -8, 16);
  TH2F* deffnum45boxplane1_ = dir.make<TH2F>("deffnum45boxplane1", "deffnum45boxplane1", 220, -2, 20, 240, -8, 16);
  TH2F* deffnum56boxplane1_ = dir.make<TH2F>("deffnum56boxplane1", "deffnum56boxplane1", 220, -2, 20, 240, -8, 16);
  TH2F* deffnum45boxplane2_ = dir.make<TH2F>("deffnum45boxplane2", "deffnum45boxplane2", 220, -2, 20, 240, -8, 16);
  TH2F* deffnum56boxplane2_ = dir.make<TH2F>("deffnum56boxplane2", "deffnum56boxplane2", 220, -2, 20, 240, -8, 16);
  TH2F* deffnum45boxplane3_ = dir.make<TH2F>("deffnum45boxplane3", "deffnum45boxplane3", 220, -2, 20, 240, -8, 16);
  TH2F* deffnum56boxplane3_ = dir.make<TH2F>("deffnum56boxplane3", "deffnum56boxplane3", 220, -2, 20, 240, -8, 16);

  // Average ToT, cylindrical RPs
  TH2F* tot45plane0_ = dir.make<TH2F>("tot45plane0", "tot45plane0", 220, -2, 20, 240, -8, 16);
  TH2F* tot56plane0_ = dir.make<TH2F>("tot56plane0", "tot56plane0", 220, -2, 20, 240, -8, 16);
  TH2F* tot45plane1_ = dir.make<TH2F>("tot45plane1", "tot45plane1", 220, -2, 20, 240, -8, 16);
  TH2F* tot56plane1_ = dir.make<TH2F>("tot56plane1", "tot56plane1", 220, -2, 20, 240, -8, 16);
  TH2F* tot45plane2_ = dir.make<TH2F>("tot45plane2", "tot45plane2", 220, -2, 20, 240, -8, 16);
  TH2F* tot56plane2_ = dir.make<TH2F>("tot56plane2", "tot56plane2", 220, -2, 20, 240, -8, 16);
  TH2F* tot45plane3_ = dir.make<TH2F>("tot45plane3", "tot45plane3", 220, -2, 20, 240, -8, 16);
  TH2F* tot56plane3_ = dir.make<TH2F>("tot56plane3", "tot56plane3", 220, -2, 20, 240, -8, 16);

  // Average ToT, box RPs
  TH2F* tot45boxplane0_ = dir.make<TH2F>("tot45boxplane0", "tot45boxplane0", 220, -2, 20, 240, -8, 16);
  TH2F* tot56boxplane0_ = dir.make<TH2F>("tot56boxplane0", "tot56boxplane0", 220, -2, 20, 240, -8, 16);
  TH2F* tot45boxplane1_ = dir.make<TH2F>("tot45boxplane1", "tot45boxplane1", 220, -2, 20, 240, -8, 16);
  TH2F* tot56boxplane1_ = dir.make<TH2F>("tot56boxplane1", "tot56boxplane1", 220, -2, 20, 240, -8, 16);
  TH2F* tot45boxplane2_ = dir.make<TH2F>("tot45boxplane2", "tot45boxplane2", 220, -2, 20, 240, -8, 16);
  TH2F* tot56boxplane2_ = dir.make<TH2F>("tot56boxplane2", "tot56boxplane2", 220, -2, 20, 240, -8, 16);
  TH2F* tot45boxplane3_ = dir.make<TH2F>("tot45boxplane3", "tot45boxplane3", 220, -2, 20, 240, -8, 16);
  TH2F* tot56boxplane3_ = dir.make<TH2F>("tot56boxplane3", "tot56boxplane3", 220, -2, 20, 240, -8, 16);

  TH2F* numvsls_ = dir.make<TH2F>("numvsls_", "numvsls_", 1000, 0, 1000, 18, -2, 16);

  // Radiographies and anti-radiographies
  TH2F* drad45_ = dir.make<TH2F>("drad45", "drad45", 220, -2, 20, 240, -8, 16);
  TH2F* drad56_ = dir.make<TH2F>("drad56", "drad56", 220, -2, 20, 240, -8, 16);
  TH2F* dantirad45_ = dir.make<TH2F>("dantirad45", "dantirad45", 220, -2, 20, 240, -8, 16);
  TH2F* dantirad56_ = dir.make<TH2F>("dantirad56", "dantirad56", 220, -2, 20, 240, -8, 16);
  TH2F* dboxrad45_ = dir.make<TH2F>("dboxrad45", "dboxrad45", 220, -2, 20, 240, -8, 16);
  TH2F* dboxrad56_ = dir.make<TH2F>("dboxrad56", "dboxrad56", 220, -2, 20, 240, -8, 16);
  TH2F* dboxantirad45_ = dir.make<TH2F>("dboxantirad45", "dboxantirad45", 220, -2, 20, 240, -8, 16);
  TH2F* dboxantirad56_ = dir.make<TH2F>("dboxantirad56", "dboxantirad56", 220, -2, 20, 240, -8, 16);

  TH1F* ls_ = dir.make<TH1F>("ls", "ls", 2000, 0, 2000);

  // Considered bunches
  TH1F* bunchNumbers_ = dir.make<TH1F>("bunchPresence", "Bunch Presence;BX;Presence", 3564, 0.5, 3564.5);

  // loop the events
  int ievt = 0;

  int lumiblock_ = -99;

  for (unsigned int iFile = 0; iFile < inFiles_.size(); ++iFile) {
    // open input file (can be located on castor)
    TFile* inFile = TFile::Open(inFiles_[iFile].c_str());
    if (!inFile || inFile->IsZombie()) {
      throw edm::Exception{edm::errors::FileOpenError} << "Could not open file!: " << inFiles_[iFile] << std::endl;
    } else {
      // ----------------------------------------------------------------------
      // Second Part:
      //
      //  * loop the events in the input file
      //  * receive the collections of interest via fwlite::Handle
      //  * fill the histograms
      //  * after the loop close the input file
      // ----------------------------------------------------------------------
      fwlite::Event ev(inFile);

      for (ev.toBegin(); !ev.atEnd(); ++ev, ++ievt) {
        edm::EventBase const& event = ev;

        // if pickedBunchesCSV provided, filter-out all that arent on the list
        if (pickedBunches_.size() && pickedBunches_.count(event.bunchCrossing()) == 0)
          continue;

        // break loop if maximal number of events is reached
        if (maxEvents_ > 0 ? ievt + 1 > maxEvents_ : false)
          break;

        // simple event counter
        if (outputEvery_ > 0 && ievt > 0 && ievt % outputEvery_ == 0) {
          edm::LogWarning("TimingEfficiencyRadiography") << "processing event: " << ievt;
        }

        // LumiSection
        lumiblock_ = ev.luminosityBlock();

        if (goodLumisections_) {
          bool found = false;
          for (auto [startLumi, endLumi] : goodLumisections_.value()) {
            if (startLumi <= lumiblock_ && lumiblock_ <= endLumi) {
              found = true;
              break;
            }
          }
          if (!found) {
            continue;
          }
        }

        bunchNumbers_->Fill(event.bunchCrossing());

        ls_->Fill(lumiblock_);

        // Long RP IDs for Lite tracks
        // Pixels
        // 2031616000
        // 2040004608
        // 2014838784
        // 2023227392
        //
        // Diamonds cylindrical
        // 2070937600
        // 2054160384

        // Diamonds box
        // 2073034752
        // 2056257536

        // If 220 pixels are needed
        float x56 = -999.;
        float y56 = -999.;
        float x45 = -999.;
        float y45 = -999.;

        float xtimetrack45 = -999.;
        float xtimetrack56 = -999.;
        float xboxtimetrack45 = -999.;
        float xboxtimetrack56 = -999.;

        int n45210 = 0;
        int n56210 = 0;
        int ntimetrack56 = 0;
        int ntimetrack45 = 0;
        int n45220 = 0;
        int n56220 = 0;
        int nboxtimetrack45 = 0;
        int nboxtimetrack56 = 0;

        bool takePixelTrack45 = false;
        bool takePixelTrack56 = false;

        // Handle to the collection of lite tracks
        edm::Handle<std::vector<CTPPSLocalTrackLite>> ppstracks;

        // Switch to run on AlCaPPS or standard Physics AOD
        if (mode_ == 1)
          event.getByLabel(std::string("ctppsLocalTrackLiteProducerAlCaRecoProducer"), ppstracks);
        else if (mode_ == 2)
          event.getByLabel(std::string("ctppsLocalTrackLiteProducer"), ppstracks);

        /*
         * Loop on tracks to get pixel tracks for the efficiency denominator
         */
        for (std::vector<CTPPSLocalTrackLite>::const_iterator track0 = ppstracks->begin(); track0 != ppstracks->end();
             ++track0) {
          // 210
          if (track0->rpId() == 2014838784) {
            n45210++;
            if (useOnlyPixelOn_ == 210) {
              x45 = track0->x();
              y45 = track0->y();
            }
          }
          if (track0->rpId() == 2031616000) {
            n56210++;
            if (useOnlyPixelOn_ == 210) {
              x56 = track0->x();
              y56 = track0->y();
            }
          }
          // 220
          if (track0->rpId() == 2023227392) {
            n45220++;
            if (useOnlyPixelOn_ != 210) {
              x45 = track0->x();
              y45 = track0->y();
            }
          }
          if (track0->rpId() == 2040004608) {
            n56220++;
            if (useOnlyPixelOn_ != 210) {
              x56 = track0->x();
              y56 = track0->y();
            }
          }
        }

        if (useOnlyPixelOn_ == 220) {
          takePixelTrack45 = (n45220 == 1);
          takePixelTrack56 = (n56220 == 1);
        } else if (useOnlyPixelOn_ == 210) {
          takePixelTrack45 = (n45210 == 1);
          takePixelTrack56 = (n56210 == 1);
        } else {
          takePixelTrack45 = (n45210 == 1 && n45220 == 1);
          takePixelTrack56 = (n56210 == 1 && n56220 == 1);
        }

        // Skip everything unless the event has exactly 1 pixel track on at least one arm
        if (takePixelTrack45 || takePixelTrack56) {
          // Denominator for efficiency: events with exactly 1 pixel track in the 45-210 pixels
          if (takePixelTrack45)
            deffden45_->Fill(x45, y45);

          // Denominator for efficiency: events with exactly 1 pixel track in the 56-210 pixels
          if (takePixelTrack56)
            deffden56_->Fill(x56, y56);

          //Loop again to check time-tracks for the per-arm efficiency numerator
          bool deffnum45_changed_ = false, deffnum56_changed_ = false, dboxeffnum45_changed_ = false,
               dboxeffnum56_changed = false;
          for (std::vector<CTPPSLocalTrackLite>::const_iterator track1 = ppstracks->begin(); track1 != ppstracks->end();
               ++track1) {
            if (track1->rpId() == 2054160384) {
              xtimetrack45 = track1->x();

              if (takePixelTrack45 && fabs(x45 - xtimetrack45) <= maxXDiff) {
                // Pixel-diamond correlation to define matching cuts
                xtpix21045_->Fill(x45, xtimetrack45);

                // Histogram for radiography
                drad45_->Fill(x45, y45);

                if (!deffnum45_changed_) {  // x-matching between pixel+time tracks for eff. numerator
                  deffnum45_->Fill(x45, y45);
                  deffnum45_changed_ = true;
                }
                ntimetrack45++;
              }
            }

            if (track1->rpId() == 2070937600) {
              xtimetrack56 = track1->x();

              if (takePixelTrack56 && fabs(x56 - xtimetrack56) <= maxXDiff) {
                /* Place where runs should be filtered out based on xtimetrack56 >= 8 if runs are from the end of 2024 */
                // Pixel-diamond correlation to define matching cuts
                xtpix21056_->Fill(x56, xtimetrack56);

                // Histogram for radiography
                drad56_->Fill(x56, y56);

                if (!deffnum56_changed_) {  // x-matching between pixel+time tracks for eff. numerator
                  deffnum56_->Fill(x56, y56);
                  deffnum56_changed_ = true;
                }
                ntimetrack56++;
              }
            }

            if (track1->rpId() == 2056257536) {
              xboxtimetrack45 = track1->x();

              if (takePixelTrack45 && fabs(x45 - xboxtimetrack45) < maxXDiff) {
                // Pixel-diamond correlation to define matching cuts
                xtpix21045_->Fill(x45, xboxtimetrack45);

                // Histogram for radiography
                dboxrad45_->Fill(x45, y45);

                if (!dboxeffnum45_changed_) {  // x-matching between pixel+time tracks for eff. numerator
                  dboxeffnum45_->Fill(x45, y45);
                  dboxeffnum45_changed_ = true;
                }
                nboxtimetrack45++;
              }
            }

            if (track1->rpId() == 2073034752) {
              xboxtimetrack56 = track1->x();

              if (takePixelTrack56 && fabs(x56 - xboxtimetrack56) < maxXDiff) {
                /* Place where runs should be filtered out based on xboxtimetrack56 >= 8 if runs are from the end of 2024 */
                // Pixel-diamond correlation to define matching cuts
                xtpix21056_->Fill(x56, xboxtimetrack56);

                // Histogram for radiography
                dboxrad56_->Fill(x56, y56);

                if (!dboxeffnum56_changed) {  // x-matching between pixel+time tracks for eff. numerator
                  dboxeffnum56_->Fill(x56, y56);
                  dboxeffnum56_changed = true;
                }
                nboxtimetrack56++;
              }
            }
          }

          // Histograms for anti-radiography: events with 1 track in pixels and none in the diamonds
          if (takePixelTrack45 && ntimetrack45 == 0)
            dantirad45_->Fill(x45, y45);
          if (takePixelTrack56 && ntimetrack56 == 0)
            dantirad56_->Fill(x56, y56);

          if (takePixelTrack45 && nboxtimetrack45 == 0)
            dboxantirad45_->Fill(x45, y45);
          if (takePixelTrack56 && nboxtimetrack56 == 0)
            dboxantirad56_->Fill(x56, y56);

          // Now loop on Diamond rechits to do plane-by-plane efficiencies
          edm::Handle<edm::DetSetVector<CTPPSDiamondRecHit>> diamondRecHits;

          // Switch to run on AlCaPPS or standard Physics AOD
          if (mode_ == 1)
            event.getByLabel(std::string("ctppsDiamondRecHitsAlCaRecoProducer"), diamondRecHits);
          else if (mode_ == 2)
            event.getByLabel(std::string("ctppsDiamondRecHits"), diamondRecHits);

          for (const auto& rechits_ds : *diamondRecHits) {
            const CTPPSDiamondDetId detidforrh(rechits_ds.detId());
            for (const auto& rechit : rechits_ds) {
              int arm = detidforrh.arm();
              int plane = detidforrh.plane();
              float xrh = rechit.x();

              int station = detidforrh.station();

              // Not currently used
              float tot = rechit.toT();

              // Sector 45
              if (station == 1 && plane == 0 && arm == 0 && takePixelTrack45 && tot >= totCut_) {
                numvsls_->Fill(lumiblock_, plane + (arm * 4));
                if (fabs(x45 - xrh) <= maxXDiff)  //  x-matching between pixel+diamond rechits for eff. numerator
                {
                  deffnum45plane0_->Fill(x45, y45);
                  tot45plane0_->Fill(x45, y45, tot);
                }
              }
              if (station == 1 && plane == 1 && arm == 0 && takePixelTrack45 && tot >= totCut_) {
                numvsls_->Fill(lumiblock_, plane + (arm * 4));
                if (fabs(x45 - xrh) <= maxXDiff)  //  x-matching between pixel+diamond rechits for eff. numerator
                {
                  deffnum45plane1_->Fill(x45, y45);
                  tot45plane1_->Fill(x45, y45, tot);
                }
              }
              if (station == 1 && plane == 2 && arm == 0 && takePixelTrack45 && tot >= totCut_) {
                numvsls_->Fill(lumiblock_, plane + (arm * 4));
                if (fabs(x45 - xrh) <= maxXDiff)  //  x-matching between pixel+diamond rechits for eff. numerator
                {
                  deffnum45plane2_->Fill(x45, y45);
                  tot45plane2_->Fill(x45, y45, tot);
                }
              }
              if (station == 1 && plane == 3 && arm == 0 && takePixelTrack45 && tot >= totCut_) {
                numvsls_->Fill(lumiblock_, plane + (arm * 4));
                if (fabs(x45 - xrh) <= maxXDiff)  //  x-matching between pixel+diamond rechits for eff. numerator
                {
                  deffnum45plane3_->Fill(x45, y45);
                  tot45plane3_->Fill(x45, y45, tot);
                }
              }

              // Sector 56
              if (station == 1 && plane == 0 && arm == 1 && takePixelTrack56 && tot >= totCut_) {
                numvsls_->Fill(lumiblock_, plane + (arm * 4));
                if (fabs(x56 - xrh) <= maxXDiff)  //  x-matching between pixel+diamond rechits for eff. numerator
                {
                  deffnum56plane0_->Fill(x56, y56);
                  tot56plane0_->Fill(x56, y56, tot);
                }
              }
              if (station == 1 && plane == 1 && arm == 1 && takePixelTrack56 && tot >= totCut_) {
                numvsls_->Fill(lumiblock_, plane + (arm * 4));
                if (fabs(x56 - xrh) <= maxXDiff)  //  x-matching between pixel+diamond rechits for eff. numerator
                {
                  deffnum56plane1_->Fill(x56, y56);
                  tot56plane1_->Fill(x56, y56, tot);
                }
              }
              if (station == 1 && plane == 2 && arm == 1 && takePixelTrack56 && tot >= totCut_) {
                numvsls_->Fill(lumiblock_, plane + (arm * 4));
                if (fabs(x56 - xrh) <= maxXDiff)  //  x-matching between pixel+diamond rechits for eff. numerator
                {
                  deffnum56plane2_->Fill(x56, y56);
                  tot56plane2_->Fill(x56, y56, tot);
                }
              }
              if (station == 1 && plane == 3 && arm == 1 && takePixelTrack56 && tot >= totCut_) {
                numvsls_->Fill(lumiblock_, plane + (arm * 4));
                if (fabs(x56 - xrh) <= maxXDiff)  //  x-matching between pixel+diamond rechits for eff. numerator
                {
                  deffnum56plane3_->Fill(x56, y56);
                  tot56plane3_->Fill(x56, y56, tot);
                }
              }

              // Sector 45 box
              if (station == 2 && plane == 0 && arm == 0 && takePixelTrack45 && tot >= totCut_) {
                numvsls_->Fill(lumiblock_, plane + (arm * 4) + 8);
                if (fabs(x45 - xrh) <= maxXDiff)  //  x-matching between pixel+diamond rechits for eff. numerator
                {
                  deffnum45boxplane0_->Fill(x45, y45);
                  tot45boxplane0_->Fill(x45, y45, tot);
                }
              }
              if (station == 2 && plane == 1 && arm == 0 && takePixelTrack45 && tot >= totCut_) {
                numvsls_->Fill(lumiblock_, plane + (arm * 4) + 8);
                if (fabs(x45 - xrh) <= maxXDiff)  //  x-matching between pixel+diamond rechits for eff. numerator
                {
                  deffnum45boxplane1_->Fill(x45, y45);
                  tot45boxplane1_->Fill(x45, y45, tot);
                }
              }
              if (station == 2 && plane == 2 && arm == 0 && takePixelTrack45 && tot >= totCut_) {
                numvsls_->Fill(lumiblock_, plane + (arm * 4) + 8);
                if (fabs(x45 - xrh) <= maxXDiff)  //  x-matching between pixel+diamond rechits for eff. numerator
                {
                  deffnum45boxplane2_->Fill(x45, y45);
                  tot45boxplane2_->Fill(x45, y45, tot);
                }
              }
              if (station == 2 && plane == 3 && arm == 0 && takePixelTrack45 && tot >= totCut_) {
                numvsls_->Fill(lumiblock_, plane + (arm * 4) + 8);
                if (fabs(x45 - xrh) <= maxXDiff)  //  x-matching between pixel+diamond rechits for eff. numerator
                {
                  deffnum45boxplane3_->Fill(x45, y45);
                  tot45boxplane3_->Fill(x45, y45, tot);
                }
              }

              // Sector 56 box
              if (station == 2 && plane == 0 && arm == 1 && takePixelTrack56 && tot >= totCut_) {
                numvsls_->Fill(lumiblock_, plane + (arm * 4) + 8);
                if (fabs(x56 - xrh) <= maxXDiff)  //  x-matching between pixel+diamond rechits for eff. numerator
                {
                  deffnum56boxplane0_->Fill(x56, y56);
                  tot56boxplane0_->Fill(x56, y56, tot);
                }
              }
              if (station == 2 && plane == 1 && arm == 1 && takePixelTrack56 && tot >= totCut_) {
                numvsls_->Fill(lumiblock_, plane + (arm * 4) + 8);
                if (fabs(x56 - xrh) <= maxXDiff)  //  x-matching between pixel+diamond rechits for eff. numerator
                {
                  deffnum56boxplane1_->Fill(x56, y56);
                  tot56boxplane1_->Fill(x56, y56, tot);
                }
              }
              if (station == 2 && plane == 2 && arm == 1 && takePixelTrack56 && tot >= totCut_) {
                numvsls_->Fill(lumiblock_, plane + (arm * 4) + 8);
                if (fabs(x56 - xrh) <= maxXDiff)  //  x-matching between pixel+diamond rechits for eff. numerator
                {
                  deffnum56boxplane2_->Fill(x56, y56);
                  tot56boxplane2_->Fill(x56, y56, tot);
                }
              }
              if (station == 2 && plane == 3 && arm == 1 && takePixelTrack56 && tot >= totCut_) {
                numvsls_->Fill(lumiblock_, plane + (arm * 4) + 8);
                if (fabs(x56 - xrh) <= maxXDiff)  //  x-matching between pixel+diamond rechits for eff. numerator
                {
                  deffnum56boxplane3_->Fill(x56, y56);
                  tot56boxplane3_->Fill(x56, y56, tot);
                }
              }
            }
          }
        }
      }

      // close input file
      inFile->Close();
    }
    // break loop if maximal number of events is reached:
    // this has to be done twice to stop the file loop as well
    if (maxEvents_ > 0 ? ievt + 1 > maxEvents_ : false)
      break;
  }
}

int main(int argc, char* argv[]) {
  // define what muon you are using; this is necessary as FWLite is not
  // capable of reading edm::Views
  using optutl::CommandLineParser;
  using reco::Muon;

  // load framework libraries
  gSystem->Load("libFWCoreFWLite");
  FWLiteEnabler::enable();

  // initialize command line parser
  optutl::CommandLineParser parser("Analyze FWLite Histograms");

  parser.integerValue("maxEvents") = 120000000;
  parser.integerValue("outputEvery") = 10000;
  parser.addOption("outputFileAllBunches", CommandLineParser::kString, "output file for all bunches", "");
  parser.addOption("outputFilePickedBunches", CommandLineParser::kString, "output file for picked bunches", "");
  parser.addOption("certificationJSONPath",
                   CommandLineParser::kString,
                   "path to certification JSON (used to filter out unwanted lumisections)",
                   "");
  parser.addOption("inputPathsCSV", CommandLineParser::kString, "Comma-separated list of input root files", "");
  parser.addOption("pickedBunchesCSV", CommandLineParser::kString, "Comma-separated list of picked bunches", "");
  parser.addOption("maxXDiff",
                   CommandLineParser::kDouble,
                   "Max x coord diff between tracks in pixel and diamond to be paired together",
                   20.0);
  parser.addOption("minimumToT", CommandLineParser::kDouble, "minimum ToT for rechits", -999.0);
  parser.addOption("mode", CommandLineParser::kInteger, "use AlCaPPS or PromptReco", 1);
  parser.addOption("useOnlyPixelOn",
                   CommandLineParser::kInteger,
                   "If some pixel has radiation damage, use only given to calculate denominators; useOnlyPixelOn=210 "
                   "if you want to use only 210, useOnlyPixelOn=220 if you only want to use 220. Providing nothing or "
                   "any other value uses their logical AND.",
                   0);

  // parse arguments
  parser.parseArguments(argc, argv);
  int maxEvents_ = parser.integerValue("maxEvents");
  unsigned int outputEvery_ = parser.integerValue("outputEvery");
  double totCut_ = parser.doubleValue("minimumToT");
  double maxXDiff = parser.doubleValue("maxXDiff");
  std::string outputFileAllBunches_ = parser.stringValue("outputFileAllBunches");
  std::string outputFilePickedBunches_ = parser.stringValue("outputFilePickedBunches");
  std::string certificationJSONPath_ = parser.stringValue("certificationJSONPath");
  std::string inputPathsCSV = parser.stringValue("inputPathsCSV");
  std::string pickedBunchesCSV = parser.stringValue("pickedBunchesCSV");
  int mode_ = parser.integerValue("mode");
  int useOnlyPixelOn_ = parser.integerValue("useOnlyPixelOn");

  if (useOnlyPixelOn_ == 220) {
    edm::LogWarning("TimingEfficiencyRadiography") << "Tracks will be suggested by only pixel 220!";
  } else if (useOnlyPixelOn_ == 210) {
    edm::LogWarning("TimingEfficiencyRadiography") << "Tracks will be suggested by only pixel 210!";
  }

  if (outputFileAllBunches_.empty() && outputFilePickedBunches_.empty()) {
    throw edm::Exception{edm::errors::NotFound}
        << "At least one output path has to be provided (outputFileAllBunches, outputFilePickedBunches).\n";
  }

  edm::LogWarning("TimingEfficiencyRadiography") << "Max X diff between diamond and pixel tracks: " << maxXDiff;

  // AOD input files
  std::vector<std::string> inFiles_;
  {
    std::vector<std::string> tokens;
    boost::split(tokens, inputPathsCSV, boost::is_any_of(","));
    for (const std::string& token : tokens) {
      inFiles_.push_back(boost::algorithm::trim_copy(token));
    }
  }

  std::optional<std::vector<std::pair<int, int>>> goodLumisections_ = std::nullopt;
  {
    auto firstNonemptyFileIt =
        std::find_if(inFiles_.begin(), inFiles_.end(), [](const std::string& s) { return !s.empty(); });
    if (firstNonemptyFileIt == inFiles_.end()) {
      throw edm::Exception{edm::errors::FileOpenError} << "No input file provided: " << inFiles_.front() << std::endl;
    }

    if (!certificationJSONPath_.empty()) {
      boost::property_tree::ptree root;

      try {
        boost::property_tree::read_json(certificationJSONPath_, root);
      } catch (const boost::property_tree::json_parser_error& e) {
        throw edm::Exception(edm::errors::FileReadError)
            << "Failed to parse JSON file: " << certificationJSONPath_ << "\n"
            << e.what() << "\n";
      }

      TFile* file = TFile::Open(firstNonemptyFileIt->c_str());
      if (!file || file->IsZombie()) {
        throw edm::Exception{edm::errors::FileOpenError} << "Failed to open ROOT file: " << inFiles_.front()
                                                         << std::endl;
      }

      fwlite::Event event(file);
      event.to(0);
      int runNumber = event.eventAuxiliary().run();
      file->Close();

      for (const auto& item : root) {
        int key = std::stoi(item.first);
        if (key != runNumber) {
          continue;
        }

        std::vector<std::pair<int, int>> pairs_vec;

        for (const auto& inner_array : item.second) {
          auto it = inner_array.second.begin();

          int first_val = it->second.get_value<int>();
          ++it;
          int second_val = it->second.get_value<int>();

          pairs_vec.emplace_back(first_val, second_val);
        }

        goodLumisections_ = std::move(pairs_vec);
        break;
      }

      if (goodLumisections_) {
        auto pair_to_string = [](const std::pair<int, int>& x) {
          std::ostringstream oss;
          oss << "[" << x.first << ", " << x.second << "]";
          return oss.str();
        };

        edm::LogWarning("TimingEfficiencyRadiography")
            << "Good lumisection ranges: ["
            << boost::algorithm::join(goodLumisections_.value() | boost::adaptors::transformed(pair_to_string), ", ")
            << "]";
      } else {
        edm::LogWarning("TimingEfficiencyRadiography")
            << "Currently processed run (" << runNumber
            << ") not found in the provided certification JSON! Processing ALL lumisections...";
      }
    }
  }

  std::set<int> pickedBunches_;
  if (!pickedBunchesCSV.empty()) {
    std::vector<std::string> tokens;
    boost::split(tokens, pickedBunchesCSV, boost::is_any_of(","));

    for (const std::string& token : tokens) {
      try {
        int bunch = boost::lexical_cast<int>(boost::algorithm::trim_copy(token));
        pickedBunches_.insert(bunch);
      } catch (const boost::bad_lexical_cast& e) {
        throw edm::Exception{edm::errors::NotFound} << "Invalid bunch number: '" << token << "'\n";
      }
    }
  }

  if (!outputFilePickedBunches_.empty() && pickedBunches_.size() == 0) {
    edm::LogWarning("TimingEfficiencyRadiography")
        << "outputFilePickedBunches argument provided but no bunches given in pickedBunchesCSV. Program will NOT "
           "run for picked bunches.";
  }

  if (!outputFilePickedBunches_.empty() && pickedBunches_.size() > 0) {
    edm::LogWarning("TimingEfficiencyRadiography")
        << "outputFilePickedBunches and pickedBunchesCSV arguments provided. \n"
        << "Calculating histograms for picked bunches: ["
        << boost::algorithm::join(pickedBunches_ | boost::adaptors::transformed(boost::lexical_cast<std::string, int>),
                                  ", ")
        << "]";

    calculateAndSaveHistograms(maxEvents_,
                               outputEvery_,
                               totCut_,
                               outputFilePickedBunches_,
                               mode_,
                               inFiles_,
                               useOnlyPixelOn_,
                               maxXDiff,
                               goodLumisections_,
                               pickedBunches_);
    edm::LogWarning("TimingEfficiencyRadiography") << "DONE Calculating histograms for picked bunches.";
  }

  if (!outputFileAllBunches_.empty()) {
    edm::LogWarning("TimingEfficiencyRadiography")
        << "outputFileAllBunches arguments provided. \nCalculating histograms for all bunches.";

    calculateAndSaveHistograms(maxEvents_,
                               outputEvery_,
                               totCut_,
                               outputFileAllBunches_,
                               mode_,
                               inFiles_,
                               useOnlyPixelOn_,
                               maxXDiff,
                               goodLumisections_);
    edm::LogWarning("TimingEfficiencyRadiography") << "DONE Calculating histograms for all bunches.";
  }

  return 0;
}
