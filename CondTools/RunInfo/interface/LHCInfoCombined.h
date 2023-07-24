#ifndef CondTools_RunInfo_LHCInfoCombined_H
#define CondTools_RunInfo_LHCInfoCombined_H

#include "FWCore/Framework/interface/EventSetup.h"

#include "CondFormats/RunInfo/interface/LHCInfo.h"
#include "CondFormats/RunInfo/interface/LHCInfoPerLS.h"
#include "CondFormats/RunInfo/interface/LHCInfoPerFill.h"

#include "CondFormats/DataRecord/interface/LHCInfoPerLSRcd.h"
#include "CondFormats/DataRecord/interface/LHCInfoPerFillRcd.h"
#include "CondFormats/DataRecord/interface/LHCInfoRcd.h"

#include <bitset>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

class LHCInfoCombined {
public:
  LHCInfoCombined() = default;

  LHCInfoCombined(const LHCInfo& lhcInfo);
  LHCInfoCombined(const LHCInfoPerLS& infoPerLS, const LHCInfoPerFill& infoPerFill);
  LHCInfoCombined(const edm::EventSetup& iSetup,
                  const edm::ESGetToken<LHCInfoPerLS, LHCInfoPerLSRcd>& tokenInfoPerLS,
                  const edm::ESGetToken<LHCInfoPerFill, LHCInfoPerFillRcd>& tokenInfoPerFill,
                  const edm::ESGetToken<LHCInfo, LHCInfoRcd>& tokenInfo, bool useNewLHCInfo);

  void setFromLHCInfo(const LHCInfo& lhcInfo);
  void setFromPerLS(const LHCInfoPerLS& infoPerLS);
  void setFromPerFill(const LHCInfoPerFill& infoPerFill);

  float crossingAngleX;
  float crossingAngleY;
  float betaStarX;
  float betaStarY;
  float energy;

  void print(std::ostream& os) const;
};

std::ostream& operator<<(std::ostream& os, LHCInfoCombined beamInfo);

#endif  // CondTools_RunInfo_LHCInfoCombined_H
