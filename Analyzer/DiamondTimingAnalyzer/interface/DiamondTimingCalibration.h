#ifndef BOSS__DIAMONDTIMINGCALIBRATION__LIBGUARD
#define BOSS__DIAMONDTIMINGCALIBRATION__LIBGUARD

#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"
#include "CondFormats/PPSObjects/interface/PPSTimingCalibration.h"

enum SECTOR{
	_45_ID,
	_56_ID
};

struct PlaneKey{
	/**
	* Identifies speciified plane.
	* Used to create object per plane (for example histogram for each plane).
	**/
    int sector, station, plane;
	
	bool operator < (const PlaneKey &rhs) const {
		return std::tie(this->sector, this->station, this->plane) < std::tie(rhs.sector, rhs.station, rhs.plane);
	};
	
	PlaneKey(const int sector=-1, const int station=-1, const int plane=-1)
	: sector(sector), station(station), plane(plane)
	{};

	PlaneKey(const CTPPSDiamondDetId& detid)
	: sector(detid.arm()), station(detid.station()), plane(detid.plane())
	{};

	std::string to_string() const {
		return "Plane{sector=" + std::to_string(this->sector) + ", station=" + std::to_string(this->station) + ", plane=" + std::to_string(this->plane) + "}";
	}

	friend std::ostream& operator<<(std::ostream& os, const PlaneKey& planeKey);
};

struct ChannelKey{
	/**
	* Same as PlaneKey, but for each channel.
	**/
    PlaneKey planeKey;
	int channel;
	
	bool operator < (const ChannelKey &rhs) const {
		return std::tie(planeKey, this->channel) < std::tie(rhs.planeKey, rhs.channel);
	}
	
	ChannelKey(const int sector=-1, const int station=-1, const int plane=-1, const int channel=-1)
	: planeKey(sector, station, plane), channel(channel)
	{};

	ChannelKey(const CTPPSDiamondDetId& detid)
	: planeKey(detid), channel(detid.channel())
	{};

	std::string to_string() const {
		return "Channel{planeKey=" + this->planeKey.to_string() + ", channel=" + std::to_string(this->channel) + "}";
	}

	friend std::ostream& operator<<(std::ostream& os, const ChannelKey& key);
};
	
struct Calib_par{
	/**
	* Calibration parameters each channel contains.
	**/
	double offset, precision;
	std::vector<double> params;
	
	Calib_par(const double offset=0, const double precision=200, const std::vector<double> params={0.0,0.0,0.0,0.0})
	: offset(offset), precision(precision), params(params)
	{};
};

class DiamondTimingCalibration{
	/**
	* Interface to PPSTimingCalibration, which uses different naming.
	**/
public:
	DiamondTimingCalibration() {};
    DiamondTimingCalibration(const PPSTimingCalibration& calib) : calib(calib) {};

    std::vector<double> parameters(const ChannelKey& key) const {
        return calib.parameters(key.planeKey.sector, key.planeKey.station, key.planeKey.plane, key.channel);
    };

    inline const std::string& formula() const { 
        return calib.formula(); 
    };

    double timeOffset(const ChannelKey& key) const {
        return calib.timeOffset(key.planeKey.sector, key.planeKey.station, key.planeKey.plane, key.channel);
    };

    double timePrecision(const ChannelKey& key) const {
        return calib.timePrecision(key.planeKey.sector, key.planeKey.station, key.planeKey.plane, key.channel);
    };

    Calib_par calibPar(const ChannelKey& key) const {
        return Calib_par(timeOffset(key), timePrecision(key), parameters(key));
    }

    friend std::ostream& operator<<(std::ostream& os, const DiamondTimingCalibration& data);

private:
    PPSTimingCalibration calib;
};

#endif