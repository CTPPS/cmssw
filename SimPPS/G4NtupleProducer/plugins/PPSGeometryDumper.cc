/**
 * PPSGeometryDumper - Dumps PPS detector geometry information
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

#include <iostream>

class PPSGeometryDumper : public edm::one::EDAnalyzer<> {
public:
  explicit PPSGeometryDumper(const edm::ParameterSet&);
  ~PPSGeometryDumper() override = default;

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  edm::ESGetToken<CTPPSGeometry, VeryForwardRealGeometryRecord> geomToken_;
  unsigned int arm_;
  unsigned int station_;
  unsigned int rp_;
  bool dumped_;
};

PPSGeometryDumper::PPSGeometryDumper(const edm::ParameterSet& iConfig)
    : geomToken_(esConsumes<CTPPSGeometry, VeryForwardRealGeometryRecord>()),
      arm_(iConfig.getParameter<unsigned int>("arm")),
      station_(iConfig.getParameter<unsigned int>("station")),
      rp_(iConfig.getParameter<unsigned int>("rp")),
      dumped_(false) {}

void PPSGeometryDumper::analyze(const edm::Event&, const edm::EventSetup& iSetup) {
  if (dumped_) return;
  dumped_ = true;

  const auto& geom = iSetup.getData(geomToken_);

  std::cout << "\n========== PPS GEOMETRY DUMP ==========\n";
  std::cout << "Detector: arm=" << arm_ << " station=" << station_ << " rp=" << rp_ << "\n\n";

  // Dump all 6 planes
  for (unsigned int plane = 0; plane < 6; ++plane) {
    CTPPSPixelDetId detId(arm_, station_, rp_, plane);

    try {
      const auto* sensor = geom.sensor(detId);
      if (sensor) {
        auto trans = sensor->translation();
        auto rot = sensor->rotation();

        std::cout << "Plane " << plane << " (detId=" << detId.rawId() << "):\n";
        std::cout << "  Translation: x=" << trans.x() << " mm, y=" << trans.y() << " mm, z=" << trans.z() << " mm\n";

        // Extract rotation angles
        double xx, xy, xz, yx, yy, yz, zx, zy, zz;
        rot.GetComponents(xx, xy, xz, yx, yy, yz, zx, zy, zz);
        std::cout << "  Rotation matrix:\n";
        std::cout << "    [" << xx << ", " << xy << ", " << xz << "]\n";
        std::cout << "    [" << yx << ", " << yy << ", " << yz << "]\n";
        std::cout << "    [" << zx << ", " << zy << ", " << zz << "]\n";

        // Test localToGlobal for origin
        CTPPSGeometry::Vector localOrigin(0., 0., 0.);
        auto globalOrigin = geom.localToGlobal(detId, localOrigin);
        std::cout << "  Local (0,0,0) -> Global: (" << globalOrigin.x() << ", "
                  << globalOrigin.y() << ", " << globalOrigin.z() << ")\n";

        // Test localToGlobal for unit vectors
        CTPPSGeometry::Vector localX(1., 0., 0.);
        CTPPSGeometry::Vector localY(0., 1., 0.);
        auto globalX = geom.localToGlobal(detId, localX);
        auto globalY = geom.localToGlobal(detId, localY);
        std::cout << "  Local (1,0,0) -> Global: (" << globalX.x() << ", "
                  << globalX.y() << ", " << globalX.z() << ")\n";
        std::cout << "  Local (0,1,0) -> Global: (" << globalY.x() << ", "
                  << globalY.y() << ", " << globalY.z() << ")\n";
        std::cout << "\n";
      } else {
        std::cout << "Plane " << plane << ": NOT FOUND in geometry\n\n";
      }
    } catch (const std::exception& e) {
      std::cout << "Plane " << plane << ": ERROR - " << e.what() << "\n\n";
    }
  }

  std::cout << "========================================\n\n";
}

DEFINE_FWK_MODULE(PPSGeometryDumper);
