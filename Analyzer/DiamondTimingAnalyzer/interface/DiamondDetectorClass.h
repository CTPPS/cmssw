#ifndef BOSS__DIAMONDDETECTORCLASS__LIBGUARD
#define BOSS__DIAMONDDETECTORCLASS__LIBGUARD


// system include files
#include <memory> 

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
// user include files
//#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSGeometry.h"
#include "DataFormats/CTPPSDigi/interface/CTPPSDiamondDigi.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondRecHit.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondLocalTrack.h"

#include "TF1.h"
 
#include <map>
#include <math.h>

#include "Analyzer/DiamondTimingAnalyzer/interface/DiamondTimingCalibration.h"

typedef std::map<CTPPSDiamondLocalTrack, std::vector<std::pair<ChannelKey, CTPPSDiamondRecHit>>> LocalTrack_map_type;
		
class DiamondDetectorClass{
public:
    explicit DiamondDetectorClass(int,
		const CTPPSGeometry&, 
		const edm::DetSetVector<CTPPSDiamondRecHit>&,
		const edm::DetSetVector<CTPPSDiamondLocalTrack>&,
		const DiamondTimingCalibration&);
	
    ~ DiamondDetectorClass();
	
	void ExtractData();
	
	inline int GetMux(const PlaneKey& planeKey) 			{return Mux_map_[planeKey];}
	inline int GetMuxValidT(const PlaneKey& planeKey) 	{return Mux_validT_map_[planeKey];}
	inline int GetMuxInTrack(const PlaneKey& planeKey) 	{return Mux_inTrack_map_[planeKey];}
	inline std::map<PlaneKey, int> GetMuxInTrackMap() 	{return Mux_inTrack_map_;}
	
	inline bool PadActive(const ChannelKey& key)   			
	{ return (RecHit_map_.find(key) != RecHit_map_.end());}
	
	inline double GetTime(const ChannelKey& key)
	{return RecHit_map_[key].at(0).time();}

	inline std::vector<CTPPSDiamondRecHit> getRecHitVector(const ChannelKey& key)
	{return RecHit_map_[key];}

	inline std::map<ChannelKey, std::vector<CTPPSDiamondRecHit>> getRecHitMap()
	{return RecHit_map_;}

	inline bool isRecHitEmpty(const ChannelKey& key)
	{return RecHit_map_[key].size() == 0;}
	
	inline double GetPadPrecision(const ChannelKey& key)
	{return SPC_map_[key].precision;}
	
	inline double GetPadWeight(const ChannelKey& key)
	{return (pow(SPC_map_[key].precision, -2));}
	
	inline double GetToT(const ChannelKey& key)
	{return RecHit_map_[key].at(0).toT();}	
	
	inline LocalTrack_map_type GetDiamondTrack_map()
	{return LocalTrack_map_;}	  
	
	int GetTrackMuxInSector(int sector);

private:
	std::map<ChannelKey, Calib_par> SPC_map_;

	// ---------- objects to retrieve ---------------------------
	edm::DetSetVector<CTPPSDiamondRecHit> recHits;
	edm::DetSetVector<CTPPSDiamondLocalTrack> localTracks;

	// ---------- data extracted---------------------------

	std::map<ChannelKey, std::vector<CTPPSDiamondRecHit>> RecHit_map_;
	LocalTrack_map_type LocalTrack_map_;
	
	
	// ---------- mux map ---------------------------
	std::map<PlaneKey, int> Mux_map_;   //arm, station. plane
	std::map<PlaneKey, int> Mux_validT_map_;  //arm, station, plane
	std::map<PlaneKey, int> Mux_inTrack_map_;  //arm, station, plane
	
	int valid_OOT_;
};
#endif
