#ifndef CalibPPS_TimingCalibration_TimingCalibrationData_h
#define CalibPPS_TimingCalibration_TimingCalibrationData_h

#include "CalibPPS/TimingCalibration/interface/DoublePeakCorrection.h"
#include "CalibPPS/TimingCalibration/interface/PlaneMap.h"

#include "DQMServices/Core/interface/MonitorElement.h"

struct TimingCalibrationData {
  using PlaneMonitorMap = std::unordered_map<PlaneKey, dqm::reco::MonitorElement*, PlaneKeyHash>;
  using ChannelMonitorMap = std::unordered_map<uint32_t, dqm::reco::MonitorElement*>;

  ChannelMonitorMap leadingTime;
  ChannelMonitorMap toT;
  ChannelMonitorMap leadingTimeVsToT;
  PlaneMonitorMap leadingTimeVsLs;

  DoublePeakCorrection doublePeakCorrection;
};

#endif
