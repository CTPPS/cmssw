#ifndef _CTPPSRPPositions_h_
#define _CTPPSRPPositions_h_

#include "CondFormats/CTPPSObjects/interface/CTPPSRPPosition.h"
#include "CondFormats/Serialization/interface/Serializable.h"

#include <map>
#include <vector>

using namespace std;

class CTPPSRPPositions{
 public:
  typedef map<unsigned int,CTPPSRPPosition>  posmap;

  CTPPSRPPositions(){}
  virtual ~CTPPSRPPositions(){}


  void setRPPosition(unsigned int RPid, const CTPPSRPPosition & RPPos);
  void  setRPPosition(unsigned int RPid, double distBPCenter=0.0);
  void  setRPPositions(const posmap & RPpositions);
  void  setRPPositions(const vector<unsigned int>& RPids,const vector<double>& distsBPCent);

  const posmap& getPosmap () const  {return m_positions;}

  CTPPSRPPosition   getRPPosition(unsigned int RPid) const;
  CTPPSRPPosition & getRPPosition(unsigned int RPid);
  
  double getRPDistBPCenter(unsigned int RPid) const;

 private:

  posmap m_positions;


  COND_SERIALIZABLE;

};

#endif
