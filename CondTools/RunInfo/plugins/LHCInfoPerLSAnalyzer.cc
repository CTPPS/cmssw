#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"

#include "CondFormats/RunInfo/interface/LHCInfoPerLS.h"
#include "CondFormats/DataRecord/interface/LHCInfoPerLSRcd.h"

#include <memory>
#include <iostream>

class LHCInfoPerLSAnalyzer : public edm::one::EDAnalyzer<> {
public:
  explicit LHCInfoPerLSAnalyzer(const edm::ParameterSet&)
      : tokenInfoPerLS_(esConsumes<LHCInfoPerLS, LHCInfoPerLSRcd>()) {}

private:
  void beginJob() override {}
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {}

  edm::ESWatcher<LHCInfoPerLSRcd> infoPerLSWatcher_;

  edm::ESGetToken<LHCInfoPerLS, LHCInfoPerLSRcd> tokenInfoPerLS_;
};

void LHCInfoPerLSAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // get InfoPerLS
  if (infoPerLSWatcher_.check(iSetup)) {
    edm::LogInfo("LHCInfoPerLSAnalyzer") << "LHCInfoPerLS retrieved:\n" << iSetup.getData(tokenInfoPerLS_);
  }
}

DEFINE_FWK_MODULE(LHCInfoPerLSAnalyzer);