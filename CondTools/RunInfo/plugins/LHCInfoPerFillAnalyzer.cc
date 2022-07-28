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

#include "CondFormats/RunInfo/interface/LHCInfoPerFill.h"
#include "CondFormats/DataRecord/interface/LHCInfoPerFillRcd.h"

#include <memory>
#include <iostream>

class LHCInfoPerFillAnalyzer : public edm::one::EDAnalyzer<> {
public:
  explicit LHCInfoPerFillAnalyzer(const edm::ParameterSet&)
      : tokenInfoPerFill_(esConsumes<LHCInfoPerFill, LHCInfoPerFillRcd>()) {}

private:
  void beginJob() override {}
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {}

  edm::ESWatcher<LHCInfoPerFillRcd> infoPerFillWatcher_;

  edm::ESGetToken<LHCInfoPerFill, LHCInfoPerFillRcd> tokenInfoPerFill_;
};

void LHCInfoPerFillAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // get InfoPerFill
  if (infoPerFillWatcher_.check(iSetup)) {
    edm::LogInfo("LHCInfoPerFillAnalyzer") << "LHCInfoPerFill retrieved:\n" << iSetup.getData(tokenInfoPerFill_);
  }
}

DEFINE_FWK_MODULE(LHCInfoPerFillAnalyzer);