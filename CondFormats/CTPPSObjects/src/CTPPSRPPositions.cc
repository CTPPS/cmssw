
#include "CondFormats/CTPPSObjects/interface/CTPPSRPPositions.h"
#include "FWCore/Utilities/interface/Exception.h" 
#include <iostream>


void CTPPSRPPositions::setRPPosition(unsigned int RPid, const CTPPSRPPosition & RPPos) {
  m_positions[RPid]=RPPos;
}


void CTPPSRPPositions::setRPPosition(unsigned int RPid, double dist) {
  m_positions[RPid]=CTPPSRPPosition(dist);
}

void CTPPSRPPositions::setRPPositions(const posmap & RPPositions){
  m_positions=RPPositions;
}

void CTPPSRPPositions::setRPPositions(const vector<unsigned int> & RPids,const vector<double> & distsBPCent){
  int n = RPids.size();
  int nd = distsBPCent.size();
  if (n != nd){
    //throw exception here
    //cout<<"The vectors have different sizes !! "<<endl;
    throw cms::Exception("CTPPSRPPositions")<<"The ids and positions vectors have different sizes !!\n";
    return;
  }

  for (int i=0; i<n ; ++i){
    m_positions[RPids.at(i)]=distsBPCent.at(i);
  }
}


CTPPSRPPosition CTPPSRPPositions::getRPPosition(unsigned int id) const {
  posmap::const_iterator it = m_positions.find(id);

  if (it != m_positions.end())
    return it->second;

  else
    //    std::cout<<"No positions defined for RP id "<< id <<std::endl;
    throw cms::Exception("CTPPSRPPositions")<< "No positions defined for RP id " << id << "\n";

  CTPPSRPPosition a;
  return a;
}

CTPPSRPPosition& CTPPSRPPositions::getRPPosition(unsigned int id) {
  return m_positions[id];
}

double CTPPSRPPositions::getRPDistanceToBeamPipeCenter(unsigned int id) const {
  CTPPSRPPosition a = getRPPosition(id);
  double d= a.getDistanceToBeamPipeCenter();
  return d;
}


double CTPPSRPPositions::getRPOffset(unsigned int id) const {
  CTPPSRPPosition a = getRPPosition(id);
  double o= a.getOffset();
  return o;
}

double CTPPSRPPositions::getRPRawLVDT(unsigned int id) const {
  CTPPSRPPosition a = getRPPosition(id);
  double lvdt= a.getRawLVDT();
  return lvdt;
}

double CTPPSRPPositions::getRPRawMotor(unsigned int id) const {
  CTPPSRPPosition a = getRPPosition(id);
  double motor= a.getRawMotor();
  return motor;
}

double CTPPSRPPositions::getRPRawResolver(unsigned int id) const {
  CTPPSRPPosition a = getRPPosition(id);
  double res= a.getRawResolver();
  return res;
}
