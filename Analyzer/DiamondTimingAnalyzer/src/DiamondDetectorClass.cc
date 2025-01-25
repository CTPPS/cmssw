#include "Analyzer/DiamondTimingAnalyzer/interface/DiamondDetectorClass.h"

#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSGeometry.h"

namespace pt = boost::property_tree;

DiamondDetectorClass::DiamondDetectorClass(int validOOT,
    const CTPPSGeometry& geom,
    const edm::DetSetVector<CTPPSDiamondRecHit>& recHits,
    const edm::DetSetVector<CTPPSDiamondLocalTrack>& localTracks,
    const DiamondTimingCalibration& calib)
    :
    recHits(recHits),
    localTracks(localTracks),
    valid_OOT_(validOOT){

    //setup calibration
    for(auto it = geom.beginSensor(); it != geom.endSensor(); ++it){
        if (!CTPPSDiamondDetId::check(it->first))
            continue;

        const CTPPSDiamondDetId detid(it->first);
        const ChannelKey key(detid);

        SPC_map_[key] = calib.calibPar(key);
    }

    this->ExtractData();
}

DiamondDetectorClass::~DiamondDetectorClass(){}

void DiamondDetectorClass::ExtractData(){
    this->localTracks = localTracks;
    this->recHits = recHits;

    RecHit_map_.clear();
    LocalTrack_map_.clear();
    Mux_map_.clear();
    Mux_validT_map_.clear();
    Mux_inTrack_map_.clear();
    TrackMuxInStation_map_.clear();

    //extract timing  tracks
    for (const auto& locTracks : localTracks){ //rechits = array of hits in one channel
        for (const auto& locTrack : locTracks){
            if((valid_OOT_==2) && (locTrack.ootIndex() != valid_OOT_) && (locTrack.ootIndex() != valid_OOT_+1)) continue;
            if (valid_OOT_!=2 && (locTrack.ootIndex() != valid_OOT_) &&  valid_OOT_!=-1 ) continue;
            //TODO: confirm if LocalTrackMap is filled corectly
            LocalTrack_map_.insert(std::make_pair(locTrack, std::vector<std::pair<ChannelKey,CTPPSDiamondRecHit>>()));
        }
    }

    // extract reco hit
    for (const auto& _recHits : recHits){ //rechits = array of hits in one channel
        const CTPPSDiamondDetId detId(_recHits.detId());
        const ChannelKey key(detId);

        // retrieve and order all events in map.
        for (const auto& recHit : _recHits){ //rechit
            // if(((recHit.ootIndex() != 0) // TODO: To be used with newest tag
            // edm::LogWarning("Offset") <<"Offset:" << SPC_map_[key].offset;
            if(((recHit.ootIndex() != (int)(SPC_map_[key].offset/25)) // TODO: To be used with newest tag
                &&  valid_OOT_!=-1) ||  recHit.multipleHits()) continue;


            Mux_map_[key.planeKey]++;

            //put in hit map ("select hit with valid leading time and TOT")
            if(recHit.time()!=0.0 && recHit.toT()> 0.0) {
                RecHit_map_[key].push_back(recHit);
                Mux_validT_map_[key.planeKey]++;

                bool counted = false;

                for (auto& locTrack_mapIter : LocalTrack_map_) {
                    if (locTrack_mapIter.first.containsHit(recHit, 0.1)){
                        locTrack_mapIter.second.push_back(std::make_pair(key,recHit));

                        if (!counted){
                            Mux_inTrack_map_[key.planeKey]++;
                            counted=true;
                        }
                    }
                }
            }
        }
    }

    for (const auto& locTrack_mapIter : LocalTrack_map_)
    {
        if (locTrack_mapIter.second.size() > 0)
        {
            int station = locTrack_mapIter.second.at(0).first.planeKey.station;
            int sector = locTrack_mapIter.second.at(0).first.planeKey.sector;
            TrackMuxInStation_map_[std::make_pair(sector,station)]++;
        }
    }
}
