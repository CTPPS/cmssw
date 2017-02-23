
#include "CondFormats/CTPPSObjects/interface/CTPPSRPPosition.h"

CTPPSRPPosition::CTPPSRPPosition(double dist,double offset,
				 double LVD,double MOTOR,double RESOLVER){
  m_distBPCenter = dist;
  m_offset       = offset;
  m_lvd          = LVD;
  m_motor        = MOTOR;
  m_resolver     = RESOLVER;
};
