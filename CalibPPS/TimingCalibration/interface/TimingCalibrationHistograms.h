#ifndef CalibPPS_TimingCalibration_TimingCalibrationHistograms_h
#define CalibPPS_TimingCalibration_TimingCalibrationHistograms_h

// #include "CalibPPS/TimingCalibration/interface/PlaneMap.h"
#include "PlaneMap.h"

#include "DQMServices/Core/interface/MonitorElement.h"

struct TimingCalibrationHistograms {
  using PlaneMonitorMap = std::unordered_map<PlaneKey, dqm::reco::MonitorElement*, PlaneKeyHash>;
  using ChannelMonitorMap = std::unordered_map<uint32_t, dqm::reco::MonitorElement*>;

  ChannelMonitorMap leadingTime;
  ChannelMonitorMap toT;
  ChannelMonitorMap leadingTimeVsToT;
  PlaneMonitorMap leadingTimeVsLs;
};

#endif
