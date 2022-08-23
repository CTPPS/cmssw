// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/DetSetVector.h"

#include "DataFormats/CTPPSDetId/interface/TotemT2DetId.h"
#include "DataFormats/TotemReco/interface/TotemT2Digi.h"

class TotemT2DigiAnalyzer : public edm::one::EDAnalyzer<> {
public:
  explicit TotemT2DigiAnalyzer(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override; 
  
  edm::EDGetTokenT<edm::DetSetVector<TotemT2Digi>> digiToken_;
  std::vector<std::vector<unsigned int>> testCasesVec;
  unsigned int numberOfAllEvents;

  std::unordered_map<unsigned int, unsigned int> event_counters_;
  unsigned int allEventCounter_ = 0;

};

TotemT2DigiAnalyzer::TotemT2DigiAnalyzer(const edm::ParameterSet& iConfig) : 
        digiToken_(consumes<edm::DetSetVector<TotemT2Digi>>(iConfig.getParameter<edm::InputTag>("digisTag"))),
        numberOfAllEvents(iConfig.getParameter<unsigned int>("numberOfAllEvents")) {

    for (const auto& idsAndLimits : iConfig.getParameter<std::vector<edm::ParameterSet>>("testCasesSet")){
        std::vector<unsigned int> ids = idsAndLimits.getParameter<std::vector<unsigned int>>("detId");
        unsigned int limit = idsAndLimits.getParameter<unsigned int>("eventLimit");

        testCasesVec.push_back(ids);
        for(auto id : ids){
            if(event_counters_.find(id) == event_counters_.end()){
                event_counters_[id] = 0;
            }
            event_counters_[id] += limit;
        }
    }

}

void TotemT2DigiAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    allEventCounter_++;
    for (const auto& ds_digis : iEvent.get(digiToken_)) {
        const TotemT2DetId detid(ds_digis.detId());
        for (const auto& digi : ds_digis) {
            event_counters_[detid]--;
            (void)digi;
        }
    }
}

// ------------ method called once each job just before starting event loop  ------------
void TotemT2DigiAnalyzer::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void TotemT2DigiAnalyzer::endJob() {
    edm::LogWarning("End job summary: events")
            << "total number of events: " << allEventCounter_ << " expected: " << numberOfAllEvents;
    for (const auto kv : event_counters_) {
        const auto id = kv.first;
        const TotemT2DetId detid(id);
        unsigned int totalDigisPerId = event_counters_[id];
        edm::LogWarning("End job summary: digis")
            << "for det id: " << id << " - " << detid << "\n"
            << "digi counter is : " << totalDigisPerId << " expected: " << 0;
        if(totalDigisPerId != 0){
            throw cms::Exception("TotemT2DigiAnalyzer") 
            << "total number of digis and expected do not match \n"
            << "for det id: " << id << " - " << detid << "\n"
            << "digi counter is : " << totalDigisPerId << " expected: " << 0;
        }
    }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TotemT2DigiAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("digisTag")->setComment("digis tag");
    desc.add<unsigned int>("numberOfAllEvents")->setComment("how much events should be analyzed");

    edm::ParameterSetDescription idmap_valid;
    idmap_valid.add<std::vector<unsigned int>>("detId")->setComment("mapped TotemT2DetId's for some test case");
    idmap_valid.add<unsigned int>("eventLimit")->setComment("event limit for some test case");

    desc.addVPSet("testCasesSet", idmap_valid);

    descriptions.add("TotemT2DigiAnalyzer", desc);
}

DEFINE_FWK_MODULE(TotemT2DigiAnalyzer);
