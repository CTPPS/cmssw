/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Author:
 *   Arkadiusz Cwikla
 *
 ****************************************************************************/

#include <memory>
#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Common/interface/DetSetVector.h"

#include "DataFormats/CTPPSDetId/interface/TotemT2DetId.h"
#include "DataFormats/TotemReco/interface/TotemT2Digi.h"

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"

class TotemT2DigiProducer : public edm::stream::EDProducer<>{
public:
    explicit TotemT2DigiProducer(const edm::ParameterSet&);

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    void produce(edm::Event&, const edm::EventSetup&) override;

    std::string t2DataFile;
    std::vector<std::vector<unsigned int>> testCasesVec;
    std::unique_ptr<TChain> inputTree_;

    // event counter for each id
    std::unordered_map<unsigned int, unsigned int> event_counters_;
    
    // variables to read data from branches    
    int eventNo;
    int numberOfChannels;
    int channelNo;
    double maxAmplitude;
    double leadingEdge;
    double trailingEdge;
    double tot;

    // CONSTANTS
    // only channel 4 contains meaningful data
    const short channelNumber_ = 4;
    // dummy chars just to satisfy TotemT2Digi constructor
    const unsigned char dummyGeo_ = 'g';
    const unsigned char dummyId_ = 'i';
    const unsigned char dummyMarker_ = 'm';
    // conversion constant between time slice and absolute time (in ns) - but inverted
    const double timeSliceNsInverted_ = 1024.0 / 25.0;
    // conversion constant between seconds and nanoseconds
    const double sToNs_ = std::pow(10, 9);
};

TotemT2DigiProducer::TotemT2DigiProducer(const edm::ParameterSet& iConfig)
    : t2DataFile(iConfig.getParameter<std::string>("t2DataFile")){
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

    inputTree_ = std::make_unique<TChain>("tree");

    inputTree_->Add(t2DataFile.c_str());

    inputTree_->SetBranchAddress("eventNo", &eventNo);
    inputTree_->SetBranchAddress("numberOfChannels", &numberOfChannels);
    inputTree_->SetBranchAddress("channelNo", &channelNo);
    inputTree_->SetBranchAddress("MaxAmplitude", &maxAmplitude);
    inputTree_->SetBranchAddress("LeadingEdge", &leadingEdge);
    inputTree_->SetBranchAddress("TrailingEdge", &trailingEdge);
    inputTree_->SetBranchAddress("TOT", &tot);

    produces<edm::DetSetVector<TotemT2Digi>>("TotemT2");
}

void TotemT2DigiProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
    int eventNum = iEvent.id().event();
    auto digi = std::make_unique<edm::DetSetVector<TotemT2Digi>>();

    inputTree_->GetEntry(eventNum);
    
    if(channelNo == channelNumber_){
        unsigned short le = (unsigned short) floor(leadingEdge * sToNs_ * timeSliceNsInverted_);
        unsigned short te = (unsigned short) floor(trailingEdge * sToNs_ * timeSliceNsInverted_);
        TotemT2Digi digiTmp(dummyGeo_, dummyId_, dummyMarker_, le, te); 

        for(auto vec : testCasesVec){
            for (const auto& id : vec) {
                // check if limit was reached
                if(event_counters_[id] > 0){
                    TotemT2DetId detId(id);
                    edm::DetSet<TotemT2Digi>& digis_for_detid = digi->find_or_insert(detId);
                    digis_for_detid.push_back(digiTmp);
                    
                    event_counters_[id]--;
                }
            }
        }
    }

    iEvent.put(std::move(digi), "TotemT2");
}

void TotemT2DigiProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<std::string>("t2DataFile")->setComment("path to fake T2 root data");

    edm::ParameterSetDescription idmap_valid;
    idmap_valid.add<std::vector<unsigned int>>("detId")->setComment("mapped TotemT2DetId's for some test case");
    idmap_valid.add<unsigned int>("eventLimit")->setComment("event limit for some test case");

    desc.addVPSet("testCasesSet", idmap_valid);

    descriptions.add("TotemT2DigiProducer", desc);
}

DEFINE_FWK_MODULE(TotemT2DigiProducer);