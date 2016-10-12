#ifndef _CTPPSRPPosition_h_
#define _CTPPSRPPosition_h_



#include "CondFormats/Serialization/interface/Serializable.h"

class CTPPSRPPosition
{
 public:
  //   CTPPSRPPosition(){}
  virtual ~CTPPSRPPosition(){}

  CTPPSRPPosition(double dist=0.0,double offset=0.0, double LVDT=-1.0, double MOTOR=-1.0, double RESOLVER=-1.0 );


  double getDistanceToBeamPipeCenter() const  {return m_distBPCenter;}
  double getOffset() const        {return m_offset;}
  double getRawLVDT() const        {return m_lvdt;}
  double getRawMotor() const      {return m_motor;}
  double getRawResolver() const   {return m_resolver;}

 private:
  double m_distBPCenter;
  double m_offset;
  double m_lvdt;
  double m_motor;
  double m_resolver;

  COND_SERIALIZABLE;

};

#endif


