#ifndef BOSS__JSONPRODUCER__LIBGUARD
#define BOSS__JSONPRODUCER__LIBGUARD

#include <iostream>
#include <map>
#include <vector>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TF1.h"

#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSGeometry.h"
#include "Analyzer/DiamondTimingAnalyzer/interface/DiamondTimingCalibration.h"

class JSON{
public:
	static void save_calib(const CTPPSGeometry& geom,
		const DiamondTimingCalibration& calib,
		std::map<ChannelKey, double>& res,
		std::string output_file_name);

	static PPSTimingCalibration read_calib(const std::string& path);

	static std::map<ChannelKey, int> read_planes_config(const std::string& path);
};

#endif