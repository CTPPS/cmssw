// -*- C++ -*-
//
// Package:    Analyzer/DiamondTimingWorker
// Class:      DiamondTimingWorker
//
/**\class DiamondTimingWorker DiamondTimingWorker.cc Analyzer/DiamondTimingWorker/plugins/DiamondTimingWorker.cc

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
#include <array>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSGeometry.h"
#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"

#include "DataFormats/CTPPSDigi/interface/CTPPSDiamondDigi.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondRecHit.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondLocalTrack.h"

#include "DataFormats/CTPPSReco/interface/CTPPSPixelLocalTrack.h"

#include "CondFormats/DataRecord/interface/PPSTimingCalibrationRcd.h"

#include "Analyzer/DiamondTimingAnalyzer/interface/DiamondDetectorClass.h"
#include "Analyzer/DiamondTimingAnalyzer/interface/JSON.h"

#include "TFitResult.h"

//
// class declaration
//
class DiamondTimingWorker : public DQMEDAnalyzer {
public:
    explicit DiamondTimingWorker(const edm::ParameterSet&);
    ~DiamondTimingWorker() = default;

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    void bookHistograms(DQMStore::IBooker& iBooker, edm::Run const&, edm::EventSetup const& iSetup) override;
    void dqmBeginRun(edm::Run const&, edm::EventSetup const&) override;
    void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

    // ---------- objects to retrieve ---------------------------
    edm::EDGetTokenT<edm::DetSetVector<CTPPSDiamondDigi>> tokenDigi_;
    edm::EDGetTokenT<edm::DetSetVector<CTPPSDiamondRecHit>> tokenRecHit_;
    edm::EDGetTokenT<edm::DetSetVector<CTPPSDiamondLocalTrack>> tokenLocalTrack_;
    edm::EDGetTokenT<edm::DetSetVector<CTPPSPixelLocalTrack>> tokenPixelLocalTrack_;
    edm::ESGetToken<CTPPSGeometry, VeryForwardRealGeometryRecord> geomEsToken_;
    edm::ESGetToken<PPSTimingCalibration, PPSTimingCalibrationRcd> calibEsToken_;

    // ------------ member data ------------
    struct Histograms_DiamondTiming {
        std::map<ChannelKey, MonitorElement*> t;
        std::map<ChannelKey, MonitorElement*> valid_t;
        std::map<ChannelKey, MonitorElement*> tot;
        std::map<ChannelKey, MonitorElement*> valid_tot;
        std::map<ChannelKey, MonitorElement*> t_vs_tot;

        std::map<ChannelKey, MonitorElement*> l2_res;
        std::map<ChannelKey, MonitorElement*> trk_time;
        std::map<ChannelKey, MonitorElement*> expected_trk_time_res;

        //pair == (sector, station)
        std::map<std::pair<int, int>, MonitorElement*> trk_time_SPC;
        std::map<std::pair<int, int>, MonitorElement*> trk_res;
        std::map<std::pair<int, int>, MonitorElement*> trk_time_SPC_vs_BX;
        std::map<std::pair<int, int>, MonitorElement*> trk_time_SPC_vs_LS;

        std::map<int, MonitorElement*> track_x_box_vs_cyl;
        std::map<int, MonitorElement*> track_time_box_vs_cyl;
        std::map<int, MonitorElement*> track_dt_vs_dx;
        std::map<int, MonitorElement*> track_time;
        std::map<int, MonitorElement*> track_resolution;
        std::map<int, MonitorElement*> timing_tracks_population;
    };
    Histograms_DiamondTiming histos;

    edm::ESHandle<CTPPSGeometry> geom;

    int validOOT;
    std::map<std::pair<int, int>, std::pair<int, int>> Ntracks_cuts_map_;  //arm, station ,, Lcut,Ucut
    std::map<ChannelKey, int> planes_config;
    int required_active_planes;
    std::string track_plots_filename;
    std::array<double, 2> track_x_max_dev{{-50, -50}};
    std::array<double, 2> track_time_max_dev{{-20, -20}};
};

//
// constants, enums and typedefs
//

enum STATION { _210_M_ID, _TIMING_ID, _220_M_ID };

const int PLANES_X_DETECTOR = 4;

//
// constructors and destructor
//
DiamondTimingWorker::DiamondTimingWorker(const edm::ParameterSet& iConfig)
    :
    tokenDigi_(
        consumes<edm::DetSetVector<CTPPSDiamondDigi>>(iConfig.getParameter<edm::InputTag>("tagDigi"))),
    tokenRecHit_(
        consumes<edm::DetSetVector<CTPPSDiamondRecHit>>(iConfig.getParameter<edm::InputTag>("tagRecHit"))),
    tokenLocalTrack_(
        consumes<edm::DetSetVector<CTPPSDiamondLocalTrack>>(iConfig.getParameter<edm::InputTag>("tagLocalTrack"))),
    tokenPixelLocalTrack_(
        consumes<edm::DetSetVector<CTPPSPixelLocalTrack>>(iConfig.getParameter<edm::InputTag>("tagPixelLocalTrack"))),
    geomEsToken_(esConsumes<CTPPSGeometry, VeryForwardRealGeometryRecord, edm::Transition::BeginRun>()),
    // calibEsToken_(esConsumes<PPSTimingCalibration, PPSTimingCalibrationRcd>()),
    validOOT(iConfig.getParameter<int>("tagValidOOT")),
    required_active_planes(iConfig.getParameter<int>("requiredActivePlanes")),
    track_plots_filename(iConfig.getParameter<std::string>("trackPlotsFilename")) {
        //TODO check if this tag is provde or no
    calibEsToken_ = esConsumes<PPSTimingCalibration, PPSTimingCalibrationRcd>(
        edm::ESInputTag(iConfig.getParameter<std::string>("timingCalibrationTag")));
    Ntracks_cuts_map_[std::make_pair(SECTOR::_45_ID, STATION::_210_M_ID)] =
        std::make_pair(iConfig.getParameter<std::vector<int>>("Ntracks_Lcuts")[0],
                       iConfig.getParameter<std::vector<int>>("Ntracks_Ucuts")[0]);

    Ntracks_cuts_map_[std::make_pair(SECTOR::_45_ID, STATION::_220_M_ID)] =
        std::make_pair(iConfig.getParameter<std::vector<int>>("Ntracks_Lcuts")[1],
                       iConfig.getParameter<std::vector<int>>("Ntracks_Ucuts")[1]);

    Ntracks_cuts_map_[std::make_pair(SECTOR::_56_ID, STATION::_210_M_ID)] =
        std::make_pair(iConfig.getParameter<std::vector<int>>("Ntracks_Lcuts")[2],
                       iConfig.getParameter<std::vector<int>>("Ntracks_Ucuts")[2]);

    Ntracks_cuts_map_[std::make_pair(SECTOR::_56_ID, STATION::_220_M_ID)] =
        std::make_pair(iConfig.getParameter<std::vector<int>>("Ntracks_Lcuts")[3],
                       iConfig.getParameter<std::vector<int>>("Ntracks_Ucuts")[3]);

    //read planes config
    planes_config = JSON::read_planes_config(iConfig.getParameter<std::string>("planesConfig"));
}

//
// member functions
//

void DiamondTimingWorker::dqmBeginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
    if (!track_plots_filename.empty()) {
        TFile tracks_plot_file{track_plots_filename.c_str()};
        if (tracks_plot_file.IsOpen()) {
            const std::string runNumber{std::to_string(iRun.run())};
            const std::array<std::pair<int, std::string>, 2> sectors{{{0, "45"}, {1, "56"}}};
            constexpr unsigned int std_dev_param_index{2};
            for (const auto& [sector_one_digit, sector_two_digits] : sectors) {
                std::string track_dt_vs_dx_plot_path{"DQMData/Run " + runNumber + "/CTPPS/Run summary/TimingDiamond/sector " + sector_two_digits + "/Track dt vs dx cyl - box sector " + std::to_string(sector_one_digit)};
                const auto* track_dt_vs_dx_plot = tracks_plot_file.Get<TH2F>(track_dt_vs_dx_plot_path.c_str());
                if (track_dt_vs_dx_plot) {
                    const std::unique_ptr<TH1D> track_x{track_dt_vs_dx_plot->ProjectionX("_py1")};
                    const TFitResultPtr& track_x_fit{track_x->Fit("gaus", "SN", "", -3, 3)};
                    const std::unique_ptr<TH1D> track_time{track_dt_vs_dx_plot->ProjectionY("_py2")};
                    const TFitResultPtr& track_time_fit{track_time->Fit("gaus", "SN")};
                    if (track_x_fit->IsValid() && track_time_fit->IsValid()) {
                        track_x_max_dev[sector_one_digit] = 3 * track_x_fit->Parameter(std_dev_param_index);
                        track_time_max_dev[sector_one_digit] = 3 * track_time_fit->Parameter(std_dev_param_index);
                    } else {
                        throw edm::Exception{edm::errors::FatalRootError} << "Can't fit the track dt vs dx projections.";
                    }
                }
            }
        }
    }
}

// ------------ method called for each event  ------------

void DiamondTimingWorker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    // edm::LogWarning("MyInfoLog") << "Entering analyze method";
    //retrieve data
    edm::Handle<edm::DetSetVector<CTPPSDiamondRecHit>> recHits; // TODO
    edm::Handle<edm::DetSetVector<CTPPSDiamondLocalTrack>> localTracks;
    edm::Handle<edm::DetSetVector<CTPPSPixelLocalTrack>> pixelLocalTracks;
    edm::ESHandle<PPSTimingCalibration> calib;

    iEvent.getByToken(tokenRecHit_, recHits);
    iEvent.getByToken(tokenLocalTrack_, localTracks);
    iEvent.getByToken(tokenPixelLocalTrack_, pixelLocalTracks);
    calib = iSetup.getHandle(calibEsToken_);
    // edm::LogWarning("MyInfoLog") << "analyze method read tokens complete";

    DiamondDetectorClass DiamondDet(validOOT, *geom, *recHits, *localTracks, DiamondTimingCalibration(*calib));

    ////////////////////////////////////////////////////////////////
    //
    //		EXTRACT PIXELS TRACK NUMBER
    //      Will be used for sector independent event selection
    //
    /////////////////////////////////////////////////////////////////
    std::map<std::pair<int, int>, int> Pixel_Mux_map_;  //arm, station
    Pixel_Mux_map_.clear();
    std::vector<bool> Sector_TBA(2, true);

    for (const auto& RP_trks : *pixelLocalTracks) {  //array of tracks
        const CTPPSDetId detId(RP_trks.detId());
        // edm::LogWarning("TrackInArm") << "Tracks in arm " << detId.arm() << ", station " << detId.station() << ", rp " << detId.rp();

        for (const auto& trk : RP_trks) {
            if (!trk.isValid())
                continue;
            Pixel_Mux_map_[std::make_pair(detId.arm(), detId.station())]++;
        }
    }
    // edm::LogWarning("MyInfoLog") << "analyze method: Sector_TBA filled complete";

    for (const auto& Ntracks_cuts_iter_ : Ntracks_cuts_map_) {
        if ((Ntracks_cuts_iter_.second.first < 0) || (Ntracks_cuts_iter_.second.second < 0))
            continue;  // don't care condition
        if ((Pixel_Mux_map_[Ntracks_cuts_iter_.first] < Ntracks_cuts_iter_.second.first) ||
            (Pixel_Mux_map_[Ntracks_cuts_iter_.first] > Ntracks_cuts_iter_.second.second))  //condition violated
            Sector_TBA[Ntracks_cuts_iter_.first.first] = false;
    }
    // edm::LogWarning("MyInfoLog") << "analyze method: Pixel MAP complete";


    if (!(Sector_TBA[0] || Sector_TBA[1]))
        return;

    ////////////////////////////////////////////////////////////////
    //
    //		control over PCL calibration quality
    //
    /////////////////////////////////////////////////////////////////
    // edm::LogWarning("MyInfoLog") << "Analyze method: begin of control over PCL calibration quality";

    for (const auto& _recHits : *recHits) {  //rechits = array of hits in one channel
        const CTPPSDiamondDetId detid(_recHits.detId());
        const ChannelKey key(detid);
        if (!(Sector_TBA[detid.arm()]))
            continue;

        // Perform channel histogram
        for (const auto& recHit : _recHits) {  //rechit
            if (((recHit.ootIndex() != 0) && validOOT != -1) || recHit.multipleHits())
                continue;

            //T,TOT and OOT for all hits, important for monitoring the calibration
            histos.t[key]->Fill(recHit.time());
            histos.tot[key]->Fill(recHit.toT());

            // T,TOT and OOT complete hits (T and TOT available), important for monitoring the calibration
            if (DiamondDet.PadActive(key)) {
                histos.valid_tot[key]->Fill(DiamondDet.GetToT(key));
                histos.t_vs_tot[key]->Fill(DiamondDet.GetToT(key), DiamondDet.GetTime(key));
                histos.valid_t[key]->Fill(DiamondDet.GetTime(key));
            }
        }
    }

    ////////////////////////////////////////////////////////////////
    //
    //		RESOLUTION STUDIES
    //
    /////////////////////////////////////////////////////////////////
    // edm::LogWarning("MyInfoLog") << "Analyze method: begin of resolution studies";

    //which planes are active
    std::array<bool, 4> active_plane{ {false, false, false, false} };
    // for (auto& item: DiamondDet.GetMuxInTrackMap()){
               // edm::LogWarning("MuxInTrackMap") << "Mux In Track Map item1:" << item.first << " item2: " << item.second;
    // }

    for (const auto& LocalTrack_mapIter : DiamondDet.GetDiamondTrack_map()) {  // loop on predigested tracks
        // edm::LogWarning("DiamondTrackMapAKALocal") << "DiamondTrackMapSector:  z (sector): " << LocalTrack_mapIter.first.z0();

        // TOTEM z convention
        int sector = LocalTrack_mapIter.first.z0() < 0.0 ? SECTOR::_45_ID : SECTOR::_56_ID; //TODO

        // if (!(Sector_TBA[sector]))
        //     continue;

        if (LocalTrack_mapIter.second.size() < 2) continue;

        //station id
        int station = LocalTrack_mapIter.second.at(0).first.planeKey.station;
        auto stationKey = std::pair<int, int>{sector, station};
        // edm::LogWarning("SectorStationNumber") << "Sector Station:  sector: " << sector << " station: " << station;
        //which planes are active
        active_plane[0] = DiamondDet.GetMuxInTrack(PlaneKey(sector, station, 0)) == 1;
        active_plane[1] = DiamondDet.GetMuxInTrack(PlaneKey(sector, station, 1)) == 1;
        active_plane[2] = DiamondDet.GetMuxInTrack(PlaneKey(sector, station, 2)) == 1;
        active_plane[3] = DiamondDet.GetMuxInTrack(PlaneKey(sector, station, 3)) == 1;
        // edm::LogWarning("GetMuxInTrackVal") <<"sector:" << sector << "station: " << station << "plane: " << 0 << "Get Mux In Track value:" << DiamondDet.GetMuxInTrack(PlaneKey(sector, station, 0));
        // edm::LogWarning("GetMuxInTrackVal") <<"sector:" << sector << "station: " << station << "plane: " << 1 << "Get Mux In Track value:" << DiamondDet.GetMuxInTrack(PlaneKey(sector, station, 1));
        // edm::LogWarning("GetMuxInTrackVal") <<"sector:" << sector << "station: " << station << "plane: " << 2 << "Get Mux In Track value:" << DiamondDet.GetMuxInTrack(PlaneKey(sector, station, 2));
        // edm::LogWarning("GetMuxInTrackVal") <<"sector:" << sector << "station: " << station << "plane: " << 3 << "Get Mux In Track value:" << DiamondDet.GetMuxInTrack(PlaneKey(sector, station, 3));

        //number of active planes
        int active_num = std::count_if(active_plane.begin(), active_plane.end(), [](bool it) -> bool{return it;});
        // edm::LogWarning("ActivePlaneNumber") << "Active Plane Number: " << active_num;

        //we don't check active planes here, because each channel might require different number of them
        // edm::LogWarning("GetTrackMuxInSector") << "GetTrackMuxInSector: " << DiamondDet.GetTrackMuxInSector(sector);

        bool mark_tag = DiamondDet.GetTrackMuxInStation(stationKey) == 1 && active_num >= required_active_planes;

        std::vector<ChannelKey> hit_selected(PLANES_X_DETECTOR);
        std::vector<std::pair<ChannelKey, CTPPSDiamondRecHit>>::const_iterator hit_iter;

        double Track_time_SPC = 100.0;
        double Track_precision_SPC = 100.0;

        // hits in track loop LocalTrack_mapIter.second = std::vector<std::pair<ChannelKey,CTPPSDiamondRecHit>>
        for (hit_iter = LocalTrack_mapIter.second.begin(); hit_iter < LocalTrack_mapIter.second.end(); hit_iter++) {
            auto& key = (*hit_iter).first;

            double hit_time_SPC = DiamondDet.GetTime(key);
            // double hit_prec_SPC = DiamondDet.GetPadPrecision(key); //TODO: Unused variable
            double hit_weig_SPC = DiamondDet.GetPadWeight(key);

            if (mark_tag && active_plane[key.planeKey.plane]) {
                // edm::LogWarning("HitSelectedIsFilling")  << "mark_tag true. " << "Key "<< key << "key.planeKey.plane "<< key.planeKey.plane;
                hit_selected[key.planeKey.plane] = key;  // save for resolution reco
            } else {
                // edm::LogWarning("HitSelectedIsNotFilling") << "mark_tag false. " << "Key "<< key;
            }
            Track_time_SPC = (Track_time_SPC * pow(Track_precision_SPC, -2) + hit_time_SPC * hit_weig_SPC) /
                             (pow(Track_precision_SPC, -2) + hit_weig_SPC);
            Track_precision_SPC = pow((pow(Track_precision_SPC, -2) + hit_weig_SPC), -0.5);
        }

        histos.trk_time_SPC[stationKey]->Fill(Track_time_SPC);
        histos.trk_res[stationKey]->Fill(Track_precision_SPC);

        histos.trk_time_SPC_vs_BX[stationKey]->Fill(iEvent.bunchCrossing(), Track_time_SPC);
        histos.trk_time_SPC_vs_LS[stationKey]->Fill(iEvent.luminosityBlock(), Track_time_SPC);

        bool is_any_track_in_other_station{false};
        for (const auto& LocalTrackOther_mapIter : DiamondDet.GetDiamondTrack_map()) {
            int sector_other = LocalTrackOther_mapIter.first.z0() < 0.0 ? SECTOR::_45_ID : SECTOR::_56_ID;
            if (LocalTrackOther_mapIter.second.size() == 0) {
                continue;
            }

            int station_other = LocalTrackOther_mapIter.second.at(0).first.planeKey.station;
            if (sector == sector_other && station == 1 && station_other == 2) {
                is_any_track_in_other_station = true;
                const auto track_x_cyl = LocalTrack_mapIter.first.x0();
                const auto track_x_box = LocalTrackOther_mapIter.first.x0();
                const auto track_time_cyl = LocalTrack_mapIter.first.time();
                const auto track_time_box = LocalTrackOther_mapIter.first.time();
                histos.track_x_box_vs_cyl[sector]->Fill(track_x_cyl, track_x_box);
                histos.track_time_box_vs_cyl[sector]->Fill(track_time_cyl, track_time_box);
                histos.track_dt_vs_dx[sector]->Fill(track_x_cyl - track_x_box, track_time_cyl - track_time_box);
                const bool are_compatible_in_x{std::abs(track_x_cyl - track_x_box) < track_x_max_dev[sector]};
                const bool are_compatible_in_time{std::abs(track_time_cyl - track_time_box) < track_time_max_dev[sector]};
                if (!are_compatible_in_x) {
                    histos.timing_tracks_population[sector]->Fill(2);
                }
                if (!are_compatible_in_time) {
                    histos.timing_tracks_population[sector]->Fill(3);
                }
                if (are_compatible_in_x && are_compatible_in_time) {
                    histos.timing_tracks_population[sector]->Fill(4);
                    std::vector<std::pair<ChannelKey, CTPPSDiamondRecHit>>::const_iterator other_hit_iter;
                    double Other_Track_time_SPC = 100.0;
                    double Other_Track_precision_SPC = 100.0;
                    for (other_hit_iter = LocalTrackOther_mapIter.second.begin(); other_hit_iter < LocalTrackOther_mapIter.second.end(); other_hit_iter++) {
                        auto& other_key = (*other_hit_iter).first;

                        double hit_time_SPC = DiamondDet.GetTime(other_key);
                        double hit_weig_SPC = DiamondDet.GetPadWeight(other_key);

                        Other_Track_time_SPC = (Other_Track_time_SPC * pow(Other_Track_precision_SPC, -2) + hit_time_SPC * hit_weig_SPC) /
                                            (pow(Other_Track_precision_SPC, -2) + hit_weig_SPC);
                        Other_Track_precision_SPC = pow((pow(Other_Track_precision_SPC, -2) + hit_weig_SPC), -0.5);
                    }
                    const double track_weight = pow(Track_precision_SPC, -2);
                    const double other_track_weight = pow(Other_Track_precision_SPC, -2);
                    const double sector_track_time = (Track_time_SPC * track_weight + Other_Track_time_SPC * other_track_weight) / (track_weight + other_track_weight);
                    const double sector_track_resolution = pow(track_weight + other_track_weight, -0.5);
                    histos.track_time[sector]->Fill(sector_track_time);
                    histos.track_resolution[sector]->Fill(sector_track_resolution);
                }
            }
        }
        if (!is_any_track_in_other_station) {
            if (station == 1) {
                histos.timing_tracks_population[sector]->Fill(0);
            } else {
                histos.timing_tracks_population[sector]->Fill(1);
            }
        }

        if (mark_tag) {
            // edm::LogWarning("RecHitMap") << "Analyze method: RecHitMap Size: " << DiamondDet.getRecHitMap().size();

            // edm::LogWarning("HitSelectedVector") << "Analyze method: hit_selected vector:  ";
            // for (auto item_hs: hit_selected) {
            //     // edm::LogWarning("HitSelectedVector") << "hit_selected vector item:  " << item_hs;
            // }
            // MarkLoop
            for (int pl_mark = 0; pl_mark < PLANES_X_DETECTOR; pl_mark++) {
                if(!active_plane[pl_mark]) {
                    // if(pl_mark + 1 == PLANES_X_DETECTOR)   edm::LogWarning("MarkLoopNoPlaneActive") << "No Plane active";
                    continue;
                }
                // edm::LogWarning("MarkLoopBegin") << "****************** BEGIN OF THE LOOP ********************* ";
                // edm::LogWarning("MarkLoopPlMarkInfo") << "pl_mark: " << pl_mark;


                double Marked_track_time = 12.5;
                double Marked_track_precision = 25.0;


                // // edm::LogWarning("MyLogInfoGetTimeDebug") << "Analyze method: get item hit_selected for pl_mark: " << pl_mark  << " value " << hit_selected[pl_mark];

                // edm::LogWarning("MyLogInfoGetTimeDebug") << "Analyze method: recHit vector size " << DiamondDet.getRecHitVector(hit_selected[pl_mark]).size();
                // for(const CTPPSDiamondRecHit& item: DiamondDet.getRecHitVector(hit_selected[pl_mark])) {
                //     // edm::LogWarning("MyLogInfoGetTimeDebug") << "Analyze method: recHit vector element: " << item.tPrecision();
                // }

                if(DiamondDet.isRecHitEmpty(hit_selected[pl_mark])) {
                    // edm::LogWarning("MarkLoopContinueGtErr") << "REC VECTOR IS EMPTY FOR PL_MARK: " << pl_mark;
                    // edm::LogWarning("Info") << "-------------- - --------------------------------------------";
                    continue;
                }


                double Marked_hit_time = DiamondDet.GetTime(hit_selected[pl_mark]); //TODO: error is thrown here
                // edm::LogWarning("MarkLoopGetTimeWorked") << "Get time worked for pl mark: " << pl_mark;

                int Marked_hit_channel = hit_selected[pl_mark].channel;

                ChannelKey key(sector, station, pl_mark, Marked_hit_channel);

                //we continue calculations for this channel only if the number of active planes is the same as in the configuration
                if(active_num != planes_config[key]) continue;

                for (int pl_loop = 0; pl_loop < PLANES_X_DETECTOR; pl_loop++) {
                    if(!active_plane[pl_loop]) continue;
                    if (pl_loop == pl_mark) continue;
                    if(DiamondDet.isRecHitEmpty(hit_selected[pl_loop])) {
                            // // edm::LogWarning("GetTimeWorked") << hit_selected[pl_loop];
			                // edm::LogWarning("MarkLoopPlLoopContinueGtErr") << "Get Time is Empty ";
			                continue;
		            }
                    double Others_hit_time = DiamondDet.GetTime(hit_selected[pl_loop]);
                    double Others_hit_weig = DiamondDet.GetPadWeight(hit_selected[pl_loop]);
                    // edm::LogWarning("MarkLoopPlLoopComplete") << "After the get time protection";

                    Marked_track_time =
                        (Marked_track_time * pow(Marked_track_precision, -2) + Others_hit_time * Others_hit_weig) /
                        (pow(Marked_track_precision, -2) + Others_hit_weig);
                    Marked_track_precision = pow((pow(Marked_track_precision, -2) + Others_hit_weig), -0.5);
                }

                //std::cout<<"T = "<<Marked_hit_time<<"; Trk = "<<Marked_track_time<<std::endl;
                double Marked_hit_difference = Marked_hit_time - Marked_track_time;

                histos.l2_res[key]->Fill(Marked_hit_difference);
                histos.trk_time[key]->Fill(Marked_track_time);
                histos.expected_trk_time_res[key]->Fill(Marked_track_precision);
                // edm::LogWarning("MarkLoopComplete") << "Mark loop complete key: " << key.to_string();
            }
        }
    }
    // edm::LogWarning("AnalyzeComplete") << "Analyze method complete";

}

void DiamondTimingWorker::bookHistograms(DQMStore::IBooker& iBooker,
                                         edm::Run const& run,
                                         edm::EventSetup const& iSetup) {
    // edm::LogWarning("MyInfoLog") << "Entering bookHistograms method";
    std::string ch_path;
    geom = iSetup.getHandle(geomEsToken_);
    for (auto it = (*geom).beginSensor(); it != (*geom).endSensor(); ++it) {
        if (!CTPPSDiamondDetId::check(it->first))
            continue;

        const CTPPSDiamondDetId detid(it->first);
        ChannelKey key(detid);

        // Booking summary histograms for sector
        if (histos.track_x_box_vs_cyl.count(detid.arm()) == 0) {
            std::string sector_path;

            detid.armName(sector_path, CTPPSDiamondDetId::nPath);

            iBooker.setCurrentFolder(sector_path);

            const auto sector = detid.arm();
            const auto arm_name = std::to_string(detid.arm());

            std::string name = "Track x box vs cyl sector " + arm_name;
            histos.track_x_box_vs_cyl[sector] = iBooker.book2D(name.c_str(), name + ";x cyl (mm);x box (mm)", 25, 0, 25, 25, 0, 25);

            name = "Track time box vs cyl sector " + arm_name;
            histos.track_time_box_vs_cyl[sector] = iBooker.book2D(name.c_str(), name + ";t cyl (ns);t box (ns)", 2000, -10, 10, 2000, -10, 10);

            name = "Track dt vs dx cyl - box sector " + arm_name;
            histos.track_dt_vs_dx[sector] = iBooker.book2D(name.c_str(), name + ";dx (mm);dt (ns)", 100, -50, 50, 2000, -20, 20);

            name = "Timing track time sector " + std::to_string(detid.arm());
            histos.track_time[sector] = iBooker.book1D(name.c_str(), name.c_str(), 1200, -60, 60);

            name = "Timing track resolution sector " + std::to_string(detid.arm());
            histos.track_resolution[sector] = iBooker.book1D(name.c_str(), name.c_str(), 1000, 0, 1);

            name = "Timing tracks population sector " + std::to_string(detid.arm());
            histos.timing_tracks_population[sector] = iBooker.book1D(name.c_str(), name.c_str(), 5, 0, 5);
            histos.timing_tracks_population[sector]->setBinLabel(1, "Only cyl");
            histos.timing_tracks_population[sector]->setBinLabel(2, "Only box");
            histos.timing_tracks_population[sector]->setBinLabel(3, "Cyl and box not compatible in x");
            histos.timing_tracks_population[sector]->setBinLabel(4, "Cyl and box not compatible in time");
            histos.timing_tracks_population[sector]->setBinLabel(5, "Cyl and box compatible in x and time");
        }

        // Booking summary histograms for station
        if (histos.trk_time_SPC.count({detid.arm(), detid.station()}) == 0) {
            std::string station_path;
            // edm::LogWarning("bookHistograms") << "Before stationName method";

            detid.stationName(station_path, CTPPSDiamondDetId::nPath);

            iBooker.setCurrentFolder(station_path);

            auto stationKey = std::make_pair<int, int>(detid.arm(), detid.station());

            std::string name = "Timing track time SPC sector " + std::to_string(detid.arm()) + " station " + std::to_string(detid.station());
            histos.trk_time_SPC[stationKey] = iBooker.book1D(name.c_str(), name.c_str(), 1200, -60, 60);

            name = "Timing track resolution sector " + std::to_string(detid.arm()) + " station " + std::to_string(detid.station());
            histos.trk_res[stationKey] = iBooker.book1D(name.c_str(), name.c_str(), 1000, 0, 1);

            name = "Timing track time SPC Vs BX sector " + std::to_string(detid.arm()) + " station " + std::to_string(detid.station());
            histos.trk_time_SPC_vs_BX[stationKey] = iBooker.book2D(name.c_str(), name.c_str(), 4000, 0, 4000, 500, -5, 5);

            name = "Timing track time SPC Vs LS sector " + std::to_string(detid.arm()) + " station " + std::to_string(detid.station());
            histos.trk_time_SPC_vs_LS[stationKey] = iBooker.book2D(name.c_str(), name.c_str(), 4000, 0, 4000, 500, -5, 5);
        }

        auto ch_name = key.to_string();

        detid.channelName(ch_path, CTPPSDiamondDetId::nPath);
        iBooker.setCurrentFolder(ch_path);

        std::string name;

        // Booking histograms for channels
        name = "Time distribution " + ch_name;
        histos.t[key] = iBooker.book1D(name.c_str(), (name + ";t (ns);Entries").c_str(), 1200, -60., 60.);

        name = "Valid Time distribution " + ch_name;
        histos.valid_t[key] = iBooker.book1D(name.c_str(), (name + ";t (ns);Entries").c_str(), 1200, -60., 60.);

        name = "TOT distribution " + ch_name;
        histos.tot[key] = iBooker.book1D(name.c_str(), (name + ";TOT (ns);Entries").c_str(), 100, -20., 20.);

        name = "Valid TOT distribution " + ch_name;
        histos.valid_tot[key] = iBooker.book1D(name.c_str(), (name + ";TOT (ns);Entries").c_str(), 100, -20., 20.);

        name = "Time vs TOT distribution " + ch_name;
        histos.t_vs_tot[key] = iBooker.book2D(name.c_str(), (name + ";TOT (ns);t (ns)").c_str(), 240, 0., 60., 450, -20., 25.);

        name = "l2_res_" + ch_name;
        histos.l2_res[key] = iBooker.book1D(name.c_str(), ("Difference between track time and hit time " + ch_name + ";t (ns);Entries").c_str(), 1200, -60, 60);

        name = "Track time distribution " + ch_name;
        histos.trk_time[key] = iBooker.book1D(name.c_str(), (name + ";t (ns);Entries").c_str(), 1200, -60, 60);

        name = "Expected track time resolution distribution " + ch_name;
        histos.expected_trk_time_res[key] =
            iBooker.book1D(name.c_str(), (name + ";t (ns);Entries").c_str(), 1000, 0, 2);
    }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DiamondTimingWorker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    // The following says we do not know what parameters are allowed so do no
    // validation
    // Please change this to state exactly what you do use, even if it is no
    // parameters
    edm::ParameterSetDescription desc;
    desc.setAllowAnything();
    // desc.add<std::string>("folder", "MY_FOLDER");
    // descriptions.add("DiamondTimingWorker", desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(DiamondTimingWorker);
