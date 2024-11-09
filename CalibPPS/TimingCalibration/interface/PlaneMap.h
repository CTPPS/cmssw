#ifndef CalibPPS_TimingCalibration_PlaneMap_h
#define CalibPPS_TimingCalibration_PlaneMap_h

#include <cstdint>
#include <functional>
#include <tuple>

using PlaneKey = std::tuple<uint32_t, uint32_t, uint32_t>;

struct PlaneKeyHash {
  std::size_t operator()(const PlaneKey& planeKey) const noexcept {
    return std::hash<uint32_t>()(std::get<0>(planeKey)) ^ std::hash<uint32_t>()(std::get<1>(planeKey)) ^
           std::hash<uint32_t>()(std::get<2>(planeKey));
  }
};

#endif
