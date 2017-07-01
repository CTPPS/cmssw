/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#ifndef _AlignmentGeometry_h_
#define _AlignmentGeometry_h_

#include "FWCore/Utilities/interface/Exception.h"

#include <map>
#include <set>
#include <string>

//----------------------------------------------------------------------------------------------------

/**
 *\brief A structure to hold relevant geometrical information about one detector/sensor.
 **/
struct DetGeometry
{
  double z;                     ///< z postion at detector centre; mm

  double sx, sy;                ///< detector nominal shift = detector center in global coordinates; in mm

  struct DirectionData
  {
    double dx, dy, dz;          ///< x, y and z components of the direction unit vector in global coordinates
    double s;                   ///< projection of (sx, sy) to (dx, dy)
  };

  std::map<unsigned int, DirectionData> directionData;

  bool isU;                     ///< only relevant for strips: true for U detectors, false for V detectors
                                ///< global U, V frame is used - that matches with u, v frame of the 1200 detector

  DetGeometry(double _z = 0., double _sx = 0., double _sy = 0., bool _isU = false) :
      z(_z), sx(_sx), sy(_sy), isU(_isU)
  {
  }

  void SetDirection(unsigned int idx, double dx, double dy, double dz)
  {
    directionData[idx] = { dx, dy, dz, dx*sx + dy*sy};
  }

  const DirectionData& GetDirectionData(unsigned int idx) const
  {
    auto it = directionData.find(idx);
    if (it == directionData.end())
      throw cms::Exception("DetGeometry") << "direction index " << idx << " not in the mapping.";

    return it->second;
  }
};

//----------------------------------------------------------------------------------------------------

/**
 *\brief A collection of geometrical information.
 * Map: (decimal) detector ID --> DetGeometry
 **/
class AlignmentGeometry : public std::map<unsigned int, DetGeometry>
{
  protected:
    std::set<unsigned int> rps;

  public:
    /// a characteristic z in mm
    double z0;

    /// puts an element to the map
    void Insert(unsigned int id, const DetGeometry &g);

    /// returns the number of detectors in the collection
    unsigned int GetNumberOfDetectors() const
    {
      return size();
    }

    /// returns the number of RPs in the collection
    unsigned int GetNumberOfRPs() const
    {
      return rps.size();
    }
    
    /// check whether the sensor Id is valid (present in the map)
    bool IsValidSensorId(unsigned int id) const
    {
      return (find(id) != end());
    }

    /// check whether the RP Id is valid (present in the set)
    bool IsValidRPId(unsigned int id) const
    {
      return (rps.find(id) != rps.end());
    }

    /// Prints the geometry.
    void Print() const;
};

#endif
