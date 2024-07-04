#include "Analyzer/DiamondTimingAnalyzer/interface/DiamondTimingCalibration.h"

std::ostream& operator<<(std::ostream& os, const DiamondTimingCalibration& data){
    os<<data.calib;
    return os;
}

std::ostream& operator<<(std::ostream& os, const ChannelKey& key){
    os<<key.to_string();
    return os;
}

std::ostream& operator<<(std::ostream& os, const PlaneKey& planeKey){
    os<<planeKey.to_string();
    return os;
}