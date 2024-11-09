#ifndef CalibPPS_TimingCalibration_TimingCalibrationData_h
#define CalibPPS_TimingCalibration_TimingCalibrationData_h

// #include "CalibPPS/TimingCalibration/interface/DoublePeakCorrection.h"
#include "DoublePeakCorrection.h"
// #include "CalibPPS/TimingCalibration/interface/TimingCalibrationHistograms.h"
#include "TimingCalibrationHistograms.h"

struct TimingCalibrationData : TimingCalibrationHistograms {
  DoublePeakCorrection doublePeakCorrection;
};

#endif
