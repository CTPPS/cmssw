// -*- C++ -*-
//
// Package:    Analyzer/DiamondTimingHarvester
// Class:      DiamondTimingHarvester
//
/**\class DiamondTimingHarvester DiamondTimingHarvester.cc Analyzer/DiamondTimingHarvester/plugins/DiamondTimingHarvester.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Arkadiusz Wolk
//         Created:  Mon, 19 Jul 2021 13:22:33 GMT
//
//

#include <string>
#include <bitset>

#include <iostream>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DQMServices/Core/interface/DQMEDHarvester.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSGeometry.h"
#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Analyzer/DiamondTimingAnalyzer/interface/DiamondDetectorClass.h"

#include "CondFormats/PPSObjects/interface/PPSTimingCalibration.h"
#include "CondFormats/DataRecord/interface/PPSTimingCalibrationRcd.h"

#include "Analyzer/DiamondTimingAnalyzer/interface/DiamondTimingCalibration.h"
#include "Analyzer/DiamondTimingAnalyzer/interface/JSON.h"

#include "TFile.h"
#include "TFitResult.h"
#include "TGraph.h"

//
// class declaration
//

class DiamondTimingHarvester : public DQMEDHarvester{
public:
  explicit DiamondTimingHarvester(const edm::ParameterSet&);
  ~DiamondTimingHarvester() = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void dqmEndJob(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter) override;
  void dqmEndRun(DQMStore::IBooker &iBooker,
                 DQMStore::IGetter &iGetter,
                 edm::Run const &iRun,
                 edm::EventSetup const &iSetup) override;

  edm::ESGetToken<CTPPSGeometry, VeryForwardRealGeometryRecord> geomEsToken_;
  edm::ESGetToken<PPSTimingCalibration, PPSTimingCalibrationRcd> calibEsToken_;

  std::string output_file;
  std::vector<std::string> calib_files; //files from previous iterations, where:
  std::vector<DiamondTimingCalibration> calibs;
  //[0] - input to first iterstion: DB or JSON
  //[1] - output from first iteration and input to second
  //[loop_index] - input to this iteration and output from previous

  int loop_index;

  double mean_max;
  double rms_max;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DiamondTimingHarvester::DiamondTimingHarvester(const edm::ParameterSet& iConfig)
    :
    geomEsToken_(esConsumes<CTPPSGeometry, VeryForwardRealGeometryRecord, edm::Transition::EndRun>()),
    // calibEsToken_(esConsumes<PPSTimingCalibration, PPSTimingCalibrationRcd, edm::Transition::EndRun>(edm::ESInputTag("PoolDBESSource:PPSTestCalibration"))),
    output_file(iConfig.getParameter<std::string>("calib_json_output")),
    calib_files(iConfig.getParameter<std::vector<std::string>>("calibFiles")),
    loop_index(iConfig.getParameter<int>("loopIndex")),
    mean_max(iConfig.getParameter<double>("meanMax")),
    rms_max(iConfig.getParameter<double>("rmsMax"))
    {
    if(calib_files.size() != (long unsigned int)loop_index){
        edm::LogError("DiamondTimingHarvester")<<"Not enough calibration files";
    }

    calibEsToken_ = esConsumes<PPSTimingCalibration, PPSTimingCalibrationRcd,  edm::Transition::EndRun>(
    edm::ESInputTag(iConfig.getParameter<std::string>("timingCalibrationTag")));

    calibs.push_back(DiamondTimingCalibration()); //empty to be replaced by input calibration to first iteration

    for(auto& file : calib_files){
        edm::LogWarning("DiamondTimingHarvester")<<"Opening file "<<file;
        calibs.push_back(DiamondTimingCalibration(JSON::read_calib(file)));
    }
}

//
// member functions
//

// ------------ method called for each event  ------------

void DiamondTimingHarvester::dqmEndJob(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter) {
  std::cout << "LOOP INDEX! = " << loop_index << std::endl;
  std::cout << "CALIB FILES! = ";
  for (long unsigned int i = 0; i < calib_files.size(); i++)
  {
    std::cout << calib_files[i] << ",";
  }
  std::cout << std::endl;
  std::cout<<"######## EndJob ########"<<std::endl;
}

void DiamondTimingHarvester::dqmEndRun(DQMStore::IBooker &iBooker,
               DQMStore::IGetter &iGetter,
               edm::Run const &iRun,
               edm::EventSetup const &iSetup) {

    std::map<ChannelKey, double> Resolution_L2_map_;

    ///////////////////////////////////////////
    // deriving full track based L2 resolution
    ///////////////////////////////////////////
    std::cout<<"####### EndRun #######"<<std::endl;

    //get data
    const auto& geom = iSetup.getData(geomEsToken_);
    edm::ESHandle<PPSTimingCalibration> calibEsTokenHandle_ = iSetup.getHandle(calibEsToken_);
    calibs[0] = DiamondTimingCalibration(*calibEsTokenHandle_);

    edm::LogWarning("Calibs") << "calibs[0]" << calibs[0];
    const auto& calib = calibs[loop_index];
    edm::LogWarning("Calibs") << "calibs["<< loop_index << "]" << calibs[loop_index];

    std::map<int, MonitorElement*> sector_track_time;
    std::map<std::pair<int, int>, MonitorElement*> station_track_time;

    std::string ch_name, ch_path;
    for (auto it = geom.beginSensor(); it != geom.endSensor(); ++it) {
        if (!CTPPSDiamondDetId::check(it->first))
            continue;

        //channel id and path
        const CTPPSDiamondDetId detid(it->first);
        ChannelKey histo_key(detid);
        ch_name = histo_key.to_string();
        detid.channelName(ch_path, CTPPSDiamondDetId::nPath);

        auto* t = iGetter.get(ch_path + "/" + "Valid Time distribution " + ch_name);

        //warnings
        auto m = t->getMean();
        auto rms = t->getRMS();
        if(std::abs(m) >= mean_max)
            edm::LogWarning("DiamondTimingHarvester::dqmEndRun")<<"[WARNING] Mean value for "<<histo_key<<" is bigger than "<<mean_max<<std::endl;
        if(rms >= rms_max)
            edm::LogWarning("DiamondTimingHarvester::dqmEndRun")<<"[WARNING] RMS value for "<<histo_key<<" is bigger than "<<rms_max<<std::endl;

        //RESOLUTION
        auto* l2_res = iGetter.get(ch_path + "/" + "l2_res_" + ch_name);
        auto* expected_trk_time = iGetter.get(ch_path + "/" + "Expected track time resolution distribution " + ch_name);
        edm::LogWarning("l2_res")<<"l2 res number: "<<  l2_res->getEntries() << " for path" <<ch_path << std::endl;
        if(l2_res->getEntries() > 100){
			l2_res->getTH1F()->Fit("gaus","+Q","",-10,10);

			if(l2_res->getTH1F()->GetFunction("gaus") != NULL){
				double ResL2_mean = l2_res->getTH1F()->GetFunction("gaus")->GetParameter(1);
				double ResL2_sigma = l2_res->getTH1F()->GetFunction("gaus")->GetParameter(2);
				l2_res->getTH1F()->Fit("gaus","","",ResL2_mean-(2.2*ResL2_sigma),ResL2_mean+(2.2*ResL2_sigma));
				ResL2_sigma = l2_res->getTH1F()->GetFunction("gaus")->GetParameter(2);

				double Exp_sigma = expected_trk_time->getTH1F()->GetMean();
				if (ResL2_sigma > Exp_sigma)
					Resolution_L2_map_[histo_key] = pow(pow(ResL2_sigma,2)-pow(Exp_sigma,2), 0.5);
                    // Resolution_L2_map_[histo_key] = 1.111;
				else
					Resolution_L2_map_[histo_key] = 0.050;
			}else{
                Resolution_L2_map_[histo_key] = 0.400;
            }
		}else{
            edm::LogWarning("DiamondTimingHarvester::dqmEndRun")<<"[WARNING] Not enough entries for "<<histo_key;
            Resolution_L2_map_[histo_key] = calib.timePrecision(histo_key);
        }

        iBooker.setCurrentFolder(ch_path);

        //current resolution
        auto* resHist = iBooker.book1D("curr_res_" + ch_name, "Resolution in current step;Resolution (ns);Entries", 1200, 0, 0.5);
        resHist->Fill(Resolution_L2_map_[histo_key]);

        //res vs step
        //auto* resStepHist = iBooker.book2D("res_vs_step_" + ch_name, "Resolution vs Step;Step;Resolution (ns)", 1200, -1, loop_index+2, 1200, 0, 1);
        auto* resStepHist = iBooker.book1D("res_vs_step_" + ch_name, "Resolution vs Step;Step;Resolution (ns)", 1200, -1, loop_index+2);
        resStepHist->getTH1F()->SetMarkerStyle(20);



        for(int i=0; i<=loop_index; i++){
            resStepHist->Fill(i, calibs[i].timePrecision(histo_key));
        }
        resStepHist->Fill(loop_index+1, Resolution_L2_map_[histo_key]);

        //res diff
        if(loop_index > 0){
            auto* mEl = iBooker.book1D("diff_res_" + ch_name, "Resolution difference;Difference (ns);Entries", 1200, -1, 1);

            double diff = std::abs(Resolution_L2_map_[histo_key] - calibs[loop_index].timePrecision(histo_key));
            mEl->Fill(diff);
        }

        const auto sector = detid.arm();
        if (!sector_track_time.contains(sector)) {
            std::string sector_path;
            detid.armName(sector_path, CTPPSDiamondDetId::nPath);
            std::string sector_one_digit{std::to_string(sector)};
            sector_track_time[sector] = iGetter.get(sector_path + "/" + "Timing track time sector " + sector_one_digit);
            const auto& sector_track_time_fit = sector_track_time[sector]->getTH1F()->Fit("gaus", "S");
            if (!sector_track_time_fit->IsValid()) {
                throw edm::Exception{edm::errors::FatalRootError} << "Can't fit the track time for sector " + sector_one_digit;
            }

        }

        const auto station = detid.station();
        const auto station_key = std::pair{sector, station};
        if (!station_track_time.contains(station_key)) {
            std::string station_path;
            detid.stationName(station_path, CTPPSDiamondDetId::nPath);
            const std::string sector_one_digit{std::to_string(sector)};
            const std::string station_one_digit{std::to_string(station)};
            station_track_time[station_key] = iGetter.get(station_path + "/" + "Timing track time SPC sector " + sector_one_digit + " station " + station_one_digit);
            station_track_time[station_key]->getTH1F()->Fit("gaus");
            const auto& station_track_time_fit = sector_track_time[sector]->getTH1F()->Fit("gaus", "S");
            if (!station_track_time_fit->IsValid()) {
                throw edm::Exception{edm::errors::FatalRootError} << "Can't fit the track time for sector " + sector_one_digit + "station " + station_one_digit;
            }
        }
    }

    //create JSON file
    JSON::save_calib(geom, calib, Resolution_L2_map_, output_file);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DiamondTimingHarvester::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    // The following says we do not know what parameters are allowed so do no
    // validation
    // Please change this to state exactly what you do use, even if it is no
    // parameters
    edm::ParameterSetDescription desc;
    desc.setAllowAnything();
    // desc.add<std::string>("folder", "MY_FOLDER");
    // descriptions.add("DiamondTimingHarvester", desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(DiamondTimingHarvester);
