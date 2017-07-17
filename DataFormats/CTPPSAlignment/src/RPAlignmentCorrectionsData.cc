/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *	Jan Kašpar (jan.kaspar@gmail.com)
 *
 ****************************************************************************/

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/typelookup.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/CTPPSAlignment/interface/RPAlignmentCorrectionsData.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"

#include <set>

using namespace std;


//----------------------------------------------------------------------------------------------------

RPAlignmentCorrectionData& RPAlignmentCorrectionsData::GetRPCorrection(unsigned int id)
{
  return rps[id];
}

//----------------------------------------------------------------------------------------------------

RPAlignmentCorrectionData RPAlignmentCorrectionsData::GetRPCorrection(unsigned int id) const
{
  RPAlignmentCorrectionData a;
  mapType::const_iterator it = rps.find(id);
  if (it != rps.end())
	  a = it->second;
  return a;
} 

//----------------------------------------------------------------------------------------------------

RPAlignmentCorrectionData& RPAlignmentCorrectionsData::GetSensorCorrection(unsigned int id)
{
  return sensors[id];
}

RPAlignmentCorrectionData RPAlignmentCorrectionsData::GetSensorCorrection(unsigned int id) const
{
  RPAlignmentCorrectionData a;
  mapType::const_iterator it = sensors.find(id);
  if (it != sensors.end())
	  a = it->second;
  return a;
}

//----------------------------------------------------------------------------------------------------

RPAlignmentCorrectionData RPAlignmentCorrectionsData::GetFullSensorCorrection(unsigned int id,
  bool useRPErrors) const
{
  RPAlignmentCorrectionData c;

  mapType::const_iterator it = rps.find(CTPPSDetId(id).getRPId());
  if (it != rps.end())
    c = it->second;

  it = sensors.find(id);
  if (it != sensors.end())
    c.add(it->second, useRPErrors);

  return c;
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentCorrectionsData::SetRPCorrection(unsigned int id, const RPAlignmentCorrectionData& ac)
{
  rps[id] = ac;
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentCorrectionsData::SetSensorCorrection(unsigned int id, const RPAlignmentCorrectionData& ac)
{
  sensors[id] = ac;
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentCorrectionsData::AddRPCorrection(unsigned int id, const RPAlignmentCorrectionData &a,
  bool sumErrors)
{
  auto it = rps.find(id);
  if (it == rps.end())
    rps.insert({id, a});
  else
    it->second.add(a, sumErrors);
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentCorrectionsData::AddSensorCorrection(unsigned int id, const RPAlignmentCorrectionData &a,
  bool sumErrors)
{
  auto it = sensors.find(id);
  if (it == sensors.end())
    sensors.insert({id, a});
  else
    it->second.add(a, sumErrors);
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentCorrectionsData::AddCorrections(const RPAlignmentCorrectionsData &nac, bool sumErrors)
{
  for (auto it = nac.rps.begin(); it != nac.rps.end(); ++it)
    AddRPCorrection(it->first, it->second, sumErrors);
  
  for (auto it = nac.sensors.begin(); it != nac.sensors.end(); ++it)
    AddSensorCorrection(it->first, it->second, sumErrors);
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentCorrectionsData::Clear()
{
  rps.clear();
  sensors.clear();
}

//----------------------------------------------------------------------------------------------------

TYPELOOKUP_DATA_REG (RPAlignmentCorrectionsData);
