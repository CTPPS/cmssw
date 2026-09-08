/**
 * PPSOpticsDumper - Dumps optical functions from CMSSW conditions DB
 *
 * Usage: cmsRun dump_optics_cfg.py
 *
 * Extracts optical functions for all RPs at all crossing angles
 * from the Global Tag conditions database.
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "CondFormats/PPSObjects/interface/LHCOpticalFunctionsSetCollection.h"
#include "CondFormats/DataRecord/interface/CTPPSOpticsRcd.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"

#include <iostream>
#include <iomanip>
#include <fstream>

class PPSOpticsDumper : public edm::one::EDAnalyzer<> {
public:
  explicit PPSOpticsDumper(const edm::ParameterSet&);
  ~PPSOpticsDumper() override = default;

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  edm::ESGetToken<LHCOpticalFunctionsSetCollection, CTPPSOpticsRcd> opticsToken_;
  bool dumped_;
};

PPSOpticsDumper::PPSOpticsDumper(const edm::ParameterSet& iConfig)
    : opticsToken_(esConsumes<LHCOpticalFunctionsSetCollection, CTPPSOpticsRcd>()),
      dumped_(false) {}

void PPSOpticsDumper::analyze(const edm::Event&, const edm::EventSetup& iSetup) {
  if (dumped_) return;
  dumped_ = true;

  std::cout << "\n";
  std::cout << "============================================================\n";
  std::cout << "    PPS OPTICAL FUNCTIONS FROM CMSSW CONDITIONS DB\n";
  std::cout << "============================================================\n\n";

  // Open CSV output file
  std::ofstream csv("cmssw_optics.csv");
  csv << "xangle,rpId,arm,station,rp,z_scoring,xi,v_x,L_x,x_D,v_y,L_y,y_D\n";

  try {
    const auto& optics = iSetup.getData(opticsToken_);

    std::cout << "Number of xangle entries: " << optics.size() << "\n\n";

    // LHCOpticalFunctionsSetCollection is: map<xangle, map<rpId, LHCOpticalFunctionsSet>>
    for (const auto& [xangle, rpMap] : optics) {
      std::cout << "========================================\n";
      std::cout << "Crossing angle: " << xangle << " urad\n";
      std::cout << "Number of RPs: " << rpMap.size() << "\n\n";

      for (const auto& [rpId, opticsSet] : rpMap) {
        std::cout << "----------------------------------------\n";
        std::cout << "RP ID (rawId): " << rpId << " (0x" << std::hex << rpId << std::dec << ")\n";

        // Decode RP ID using CTPPSDetId
        // The rawId format for CTPPS is defined in CTPPSDetId
        CTPPSDetId detId(rpId);
        unsigned int arm = detId.arm();
        unsigned int station = detId.station();
        unsigned int rp = detId.rp();

        std::cout << "  Decoded: arm=" << arm << ", station=" << station << ", rp=" << rp << "\n";

        // Get the scoring plane z
        double scoringPlaneZ = opticsSet.getScoringPlaneZ();
        std::cout << "  Scoring plane Z: " << scoringPlaneZ << " m\n";

        // Get optical functions - directly read from the stored arrays
        // LHCOpticalFunctionsSet::getFcnValues() returns vector<vector<double>>
        // Testing both orientations to determine actual structure
        const auto& xiVals = opticsSet.getXiValues();
        const auto& fcnVals = opticsSet.getFcnValues();
        std::cout << "  xi points: " << xiVals.size() << ", fcnVals: " << fcnVals.size();
        if (!fcnVals.empty()) std::cout << " x " << fcnVals[0].size();
        std::cout << "\n";

        // Access ALL xi points: fcnVals[fcn_index][xi_index]
        // fcn indices: evx=0, eLx=1, e14=2, exd=3, evpx=4, eLpx=5, e24=6, expd=7
        //              e32=8, evy=9, eLy=10, eyd=11, e42=12, evpy=13, eLpy=14, eypd=15
        if (fcnVals.size() >= 12 && !xiVals.empty()) {
          // Print summary for console (every 10th point to keep output manageable)
          std::cout << "  Dumping " << xiVals.size() << " xi points to CSV...\n";

          // Output ALL points to CSV
          for (size_t xi_idx = 0; xi_idx < xiVals.size(); ++xi_idx) {
            if (xi_idx >= fcnVals[0].size()) continue;
            double xi = xiVals[xi_idx];

            // Transposed access: fcnVals[fcn][xi]
            double v_x = fcnVals[0][xi_idx];   // evx
            double L_x = fcnVals[1][xi_idx];   // eLx (m)
            double x_D = fcnVals[3][xi_idx];   // exd (m)
            double v_y = fcnVals[9][xi_idx];   // evy
            double L_y = fcnVals[10][xi_idx];  // eLy (m)
            double y_D = fcnVals[11][xi_idx];  // eyd (m)

            // Print every 10th point to console for verification
            // Note: CMSSW stores optical functions in cm (per LHCOpticalFunctionsSet.h)
            if (xi_idx % 10 == 0) {
              std::cout << "  xi=" << std::fixed << std::setprecision(4) << xi << ": ";
              std::cout << "x_D=" << std::setw(8) << std::setprecision(4) << x_D << " cm, ";
              std::cout << "L_x=" << std::setw(9) << std::setprecision(1) << L_x << " cm, ";
              std::cout << "L_y=" << std::setw(9) << std::setprecision(1) << L_y << " cm, ";
              std::cout << "v_y=" << std::setw(7) << std::setprecision(3) << v_y << "\n";
            }

            // Write ALL points to CSV (units: cm, same as stored in CMSSW conditions)
            csv << xangle << "," << rpId << "," << arm << "," << station << "," << rp << ","
                << std::fixed << std::setprecision(6) << scoringPlaneZ << ","
                << std::setprecision(6) << xi << ","
                << std::setprecision(6) << v_x << ","
                << std::setprecision(6) << L_x << ","
                << std::setprecision(6) << x_D << ","
                << std::setprecision(6) << v_y << ","
                << std::setprecision(6) << L_y << ","
                << std::setprecision(6) << y_D << "\n";
          }
        }

        std::cout << "  Xi range: [" << xiVals.front() << ", " << xiVals.back() << "]\n";

        std::cout << "\n";
      }
    }

  } catch (const std::exception& e) {
    std::cout << "ERROR accessing optics: " << e.what() << "\n";
    std::cout << "Optics may not be available in this GlobalTag.\n";
  }

  csv.close();
  std::cout << "Wrote optics to cmssw_optics.csv\n";
  std::cout << "============================================================\n\n";
}

DEFINE_FWK_MODULE(PPSOpticsDumper);
