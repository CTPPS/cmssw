#ifndef CalibPPS_TimingCalibration_TimingCalibrationStruct_h
#define CalibPPS_TimingCalibration_TimingCalibrationStruct_h

#include "DQMServices/Core/interface/DQMStore.h"
#include <cstdint>
#include <tuple>
#include <unordered_map>

namespace std {

template<>
struct hash<std::tuple<uint8_t, uint8_t, uint8_t>> {
    size_t operator()(const std::tuple<uint8_t, uint8_t, uint8_t>& t) const {
        return std::hash<uint8_t>()(std::get<0>(t)) ^ std::hash<uint8_t>()(std::get<1>(t)) ^ std::hash<uint8_t>()(std::get<2>(t));
    }
};

} // namespace std

struct TimingCalibrationHistograms {
public:
  TimingCalibrationHistograms() = default;

  using MonitorMap = std::unordered_map<uint32_t, dqm::reco::MonitorElement*>;

  MonitorMap leadingTime, toT;
  MonitorMap leadingTimeVsToT;
  std::unordered_map<std::tuple<uint8_t, uint8_t, uint8_t>, dqm::reco::MonitorElement*> leadingTimeVsLs;
};

#endif
