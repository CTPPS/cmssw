
#include "CondFormats/CTPPSObjects/interface/CTPPSRPPosition.h"

CTPPSRPPosition::CTPPSRPPosition(double dist,double offset,
				 double LVDT,double MOTOR,double RESOLVER){
  m_distBPCenter = dist;
  m_offset       = offset;
  m_lvdt         = LVDT;
  m_motor        = MOTOR;
  m_resolver     = RESOLVER;
};
