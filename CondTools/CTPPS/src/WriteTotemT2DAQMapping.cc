/****************************************************************************
*
* Offline analyzer for writing Totem T2 DAQ Mapping sqlite file 
* Arkadiusz Cwikla
* based on code from H. Malbouisson,
* which was based on TOTEM code from  Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CondCore/CondDB/interface/Time.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"

#include "CondFormats/DataRecord/interface/TotemT2DAQMappingRcd.h"
// #include "CondFormats/DataRecord/interface/CTPPSPixelAnalysisMaskRcd.h"
#include "CondFormats/PPSObjects/interface/TotemDAQMapping.h"
// #include "CondFormats/PPSObjects/interface/TotemAnalysisMask.h"
#include "FWCore/Utilities/interface/ESInputTag.h"
#include "CondFormats/DataRecord/interface/TotemReadoutRcd.h"
#include "CondFormats/PPSObjects/interface/TotemDAQMapping.h"

#include <cstdint>

//----------------------------------------------------------------------------------------------------
class WriteTotemT2DAQMapping : public edm::one::EDAnalyzer<> {
public:
  WriteTotemT2DAQMapping(const edm::ParameterSet &ps);
  ~WriteTotemT2DAQMapping() override {}

private:
  void analyze(const edm::Event &e, const edm::EventSetup &es) override;
  cond::Time_t daqmappingiov_;
  std::string record_;
  std::string label_;

  edm::ESGetToken<TotemDAQMapping, TotemReadoutRcd> tokenMapping_;
};

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

WriteTotemT2DAQMapping::WriteTotemT2DAQMapping(const edm::ParameterSet &ps)
    : daqmappingiov_(ps.getParameter<unsigned long long>("daqmappingiov")),
      record_(ps.getParameter<string>("record")),
      label_(ps.getParameter<string>("label")),
      tokenMapping_(esConsumes<TotemDAQMapping, TotemReadoutRcd>(edm::ESInputTag("", label_))) {}

void WriteTotemT2DAQMapping::analyze(const edm::Event &, edm::EventSetup const &es) {
  // print mapping
  

  // Write DAQ Mapping to sqlite file:

  const auto &mapping = es.getData(tokenMapping_);

  printf("* VFAT mapping\n");
  for (const auto &p : mapping.VFATMapping)
    cout << "    " << p.first << " -> " << p.second << endl;
  printf("* Channel mapping\n");
  for (const auto &p : mapping.totemTimingChannelMap)
    cout << "    " << p.first << " plane " << p.second.plane << " channel " << p.second.channel << endl;
  edm::Service<cond::service::PoolDBOutputService> poolDbService;

  cout << "is pool db service available: " << poolDbService.isAvailable() << endl;
  if (poolDbService.isAvailable()) {
    poolDbService->writeOneIOV(mapping, daqmappingiov_, record_);
  }
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(WriteTotemT2DAQMapping);
