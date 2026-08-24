/**
 * PPSAllRPsDumper - Dumps Z positions for all PPS Roman Pots
 *
 * Queries the CMSSW conditions database and prints the global Z position
 * for each Roman Pot sensor plane.
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSGeometry.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSPixelDetId.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"

#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include <cmath>

class PPSAllRPsDumper : public edm::one::EDAnalyzer<> {
public:
  explicit PPSAllRPsDumper(const edm::ParameterSet&);
  ~PPSAllRPsDumper() override = default;

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  edm::ESGetToken<CTPPSGeometry, VeryForwardRealGeometryRecord> geomToken_;
  bool dumped_;
};

PPSAllRPsDumper::PPSAllRPsDumper(const edm::ParameterSet& iConfig)
    : geomToken_(esConsumes<CTPPSGeometry, VeryForwardRealGeometryRecord>()),
      dumped_(false) {}

void PPSAllRPsDumper::analyze(const edm::Event&, const edm::EventSetup& iSetup) {
  if (dumped_) return;
  dumped_ = true;

  const auto& geom = iSetup.getData(geomToken_);

  std::cout << "\n";
  std::cout << "============================================================\n";
  std::cout << "    PPS ROMAN POT Z POSITIONS FROM CMSSW CONDITIONS DB\n";
  std::cout << "============================================================\n\n";

  // Map to store RP z positions: (arm, station, rp) -> z
  std::map<std::tuple<int,int,int>, double> rpPositions;

  // Sector 45 (arm 0) and Sector 56 (arm 1)
  // Station: 0 = 210m near, 2 = 220m near
  // RP: 3 = pixel horizontal, 23 = pixel far station
  //
  // Full RP IDs for pixel detectors:
  //   Sector 45: 3 (210m), 23 (220m)
  //   Sector 56: 103 (210m), 123 (220m)

  std::vector<std::tuple<unsigned int, unsigned int, unsigned int, std::string>> rps = {
    // arm, station, rp, name
    // Sector 45 (arm 0)
    {0, 0, 3, "sector45_210m_pixel (XRPH_D6L5_B2)"},
    {0, 2, 3, "sector45_220m_pixel (XRPH_B6L5_B2)"},
    // Sector 56 (arm 1)
    {1, 0, 3, "sector56_210m_pixel (XRPH_D6R5_B1)"},
    {1, 2, 3, "sector56_220m_pixel (XRPH_A6R5_B1)"},
  };

  std::cout << std::fixed << std::setprecision(3);
  std::cout << "Roman Pot               | Arm | Stn | RP |   Z (mm)    |  Z (m)\n";
  std::cout << "------------------------|-----|-----|----|-------------|--------\n";

  std::map<std::string, double> zValues;

  for (const auto& [arm, station, rp, name] : rps) {
    // Get plane 0 position (representative of the RP)
    CTPPSPixelDetId detId(arm, station, rp, 0);

    try {
      const auto* sensor = geom.sensor(detId);
      if (sensor) {
        auto trans = sensor->translation();
        double z_mm = trans.z();
        double z_m = z_mm / 1000.0;

        std::cout << std::left << std::setw(24) << name
                  << "| " << arm << "   | " << station << "   | " << rp
                  << "  | " << std::right << std::setw(11) << z_mm
                  << " | " << std::setw(7) << z_m << "\n";

        zValues[name] = z_mm;
      } else {
        std::cout << std::left << std::setw(24) << name
                  << "| " << arm << "   | " << station << "   | " << rp
                  << "  | NOT FOUND\n";
      }
    } catch (const std::exception& e) {
      std::cout << std::left << std::setw(24) << name
                << "| " << arm << "   | " << station << "   | " << rp
                << "  | ERROR: " << e.what() << "\n";
    }
  }

  // Calculate separations for multi-RP
  std::cout << "\n";
  std::cout << "============================================================\n";
  std::cout << "    STATION SEPARATIONS FOR MULTI-RP RECONSTRUCTION\n";
  std::cout << "============================================================\n\n";

  // Sector 45
  if (zValues.count("sector45_210m_pixel (XRPH_D6L5_B2)") &&
      zValues.count("sector45_220m_pixel (XRPH_B6L5_B2)")) {
    double z210 = zValues["sector45_210m_pixel (XRPH_D6L5_B2)"];
    double z220 = zValues["sector45_220m_pixel (XRPH_B6L5_B2)"];
    double sep = std::abs(z220 - z210);
    std::cout << "Sector 45 (arm 0):\n";
    std::cout << "  210m station: " << z210 << " mm = " << z210/1000.0 << " m\n";
    std::cout << "  220m station: " << z220 << " mm = " << z220/1000.0 << " m\n";
    std::cout << "  Separation:   " << sep << " mm = " << sep/1000.0 << " m\n\n";
  }

  // Sector 56
  if (zValues.count("sector56_210m_pixel (XRPH_D6R5_B1)") &&
      zValues.count("sector56_220m_pixel (XRPH_A6R5_B1)")) {
    double z210 = zValues["sector56_210m_pixel (XRPH_D6R5_B1)"];
    double z220 = zValues["sector56_220m_pixel (XRPH_A6R5_B1)"];
    double sep = std::abs(z220 - z210);
    std::cout << "Sector 56 (arm 1):\n";
    std::cout << "  210m station: " << z210 << " mm = " << z210/1000.0 << " m\n";
    std::cout << "  220m station: " << z220 << " mm = " << z220/1000.0 << " m\n";
    std::cout << "  Separation:   " << sep << " mm = " << sep/1000.0 << " m\n\n";
  }

  std::cout << "============================================================\n\n";
}

DEFINE_FWK_MODULE(PPSAllRPsDumper);
