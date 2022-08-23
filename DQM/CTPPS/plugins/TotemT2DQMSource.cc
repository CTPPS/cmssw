/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Author:
 *   Laurent Forthomme
 *
 ****************************************************************************/

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"

#include "DataFormats/CTPPSDetId/interface/TotemT2DetId.h"
#include "DataFormats/TotemReco/interface/TotemT2Digi.h"
#include "DataFormats/TotemReco/interface/TotemT2RecHit.h"

#include "Geometry/Records/interface/TotemGeometryRcd.h"
#include "DQM/CTPPS/interface/TotemT2Segmentation.h"

#include <string>

class TotemT2DQMSource : public DQMEDAnalyzer {
public:
  TotemT2DQMSource(const edm::ParameterSet&);
  ~TotemT2DQMSource() override;

protected:
  void dqmBeginRun(const edm::Run&, const edm::EventSetup&) override;
  void bookHistograms(DQMStore::IBooker&, const edm::Run&, const edm::EventSetup&) override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  edm::ESGetToken<TotemGeometry, TotemGeometryRcd> geometryToken_;
  edm::EDGetTokenT<edm::DetSetVector<TotemT2Digi>> digiToken_;
  edm::EDGetTokenT<edmNew::DetSetVector<TotemT2RecHit>> rechitToken_;

  std::unique_ptr<TotemT2Segmentation> segm_;

  void fillActivePlanes(std::unordered_map<unsigned int, std::set<unsigned int>>&, const TotemT2DetId);

  void fillTriggerBitset(const TotemT2DetId detid);
  void clearTriggerBitset();
  bool areChannelsTriggered(const TotemT2DetId detid);

  double HPTDC_BIN_WIDTH_NS = 25. / 1024;
  MonitorElement* HPTDCErrorFlags_2D = nullptr;
  void bookErrorFlagsHistogram(DQMStore::IBooker&);
  void fillErrorFlagsHistogram(TotemT2Digi digi);

  void fillEdges(TotemT2Digi digi, TotemT2DetId detid);

  struct SectorPlots{
    MonitorElement* activePlanes = nullptr;
    MonitorElement* activePlanesCount = nullptr;

    MonitorElement* triggerEmulator = nullptr;
    std::bitset<(TotemT2DetId::maxPlane+1) * (TotemT2DetId::maxChannel+1)> hitTilesArray;
    static const unsigned int MINIMAL_TRIGGER = 3;

    MonitorElement *leadingEdge = nullptr, *trailingEdge = nullptr, *timeOverTreshold = nullptr;
  
    SectorPlots() = default;
    SectorPlots(DQMStore::IBooker& ibooker, unsigned int id);
  };

  struct PlanePlots{
    MonitorElement* digisMultiplicity = nullptr;
    MonitorElement* rechitMultiplicity = nullptr;
    
    PlanePlots() = default;
    PlanePlots(DQMStore::IBooker& ibooker, unsigned int id);
  };

  // struct ChannelPlots{
  //   ChannelPlots() = default;
  //   ChannelPlots(DQMStore::IBooker& ibooker, unsigned int id);
  // };

  std::unordered_map<unsigned int, SectorPlots> sectorPlots_;
  std::unordered_map<unsigned int, PlanePlots> planePlots_;

};

TotemT2DQMSource::SectorPlots::SectorPlots(DQMStore::IBooker& ibooker, unsigned int id){
  std::string title, path;

  TotemT2DetId(id).armName(path, TotemT2DetId::nPath);
  ibooker.setCurrentFolder(path);

  TotemT2DetId(id).armName(title, TotemT2DetId::nFull);

  activePlanes = ibooker.book1D("active planes",
                                  title + " which planes are active;plane number",
                                  8,
                                  -0.5,
                                  7.5);
                                  
  activePlanesCount = ibooker.book1D("number of active planes",
                                  title + " how many planes are active;number of active planes",
                                  8,
                                  -0.5,
                                  7.5);

  const size_t summary_nbinsx = 25, summary_nbinsy = 25; // TODO make it a parameter
  triggerEmulator = ibooker.book2DD("trigger emulator",
                                  title + " trigger emulator;arbitrary unit;arbitrary unit",
                                  summary_nbinsx,
                                  0.,
                                  summary_nbinsx,
                                  summary_nbinsy,
                                  0.,
                                  summary_nbinsy);
  int windowsNum = 20; // TODO make it a parameter
  leadingEdge = ibooker.book1D("leading edge",
                                  title + " leading edge (DIGIs); leading edge (ns)",
                                  25 * windowsNum,
                                  0,
                                  25 * windowsNum);
  trailingEdge = ibooker.book1D("trailing edge",
                                  title + " trailing edge (DIGIs); trailing edge (ns)",
                                  25 * windowsNum,
                                  0,
                                  25 * windowsNum);

  timeOverTreshold = ibooker.book1D("time over threshold",
                                  title + " time over threshold;time over threshold (ns)",
                                  250, 
                                  -25, 
                                  100);                                         

}

TotemT2DQMSource::PlanePlots::PlanePlots(DQMStore::IBooker& ibooker, unsigned int id){
  std::string title, path;
  TotemT2DetId(id).planeName(title, TotemT2DetId::nFull);
  TotemT2DetId(id).planeName(path, TotemT2DetId::nPath);
  ibooker.setCurrentFolder(path);

  const size_t summary_nbinsx = 25, summary_nbinsy = 25; // TODO make it a parameter
  digisMultiplicity = ibooker.book2DD("digis multiplicity",
                                        title + " digis multiplicity;x;y", // TODO change x and y 
                                        summary_nbinsx,
                                        0.,
                                        summary_nbinsx,
                                        summary_nbinsy,
                                        0.,
                                        summary_nbinsy);
  rechitMultiplicity = ibooker.book2DD("rechits multiplicity",
                                        title + " rechits multiplicity;x;y",
                                        summary_nbinsx,
                                        0.,
                                        summary_nbinsx,
                                        summary_nbinsy,
                                        0.,
                                        summary_nbinsy);
}

TotemT2DQMSource::TotemT2DQMSource(const edm::ParameterSet& iConfig)
    : geometryToken_(esConsumes<TotemGeometry, TotemGeometryRcd, edm::Transition::BeginRun>()),
      digiToken_(consumes<edm::DetSetVector<TotemT2Digi>>(iConfig.getParameter<edm::InputTag>("digisTag"))),
      rechitToken_(consumes<edmNew::DetSetVector<TotemT2RecHit>>(iConfig.getParameter<edm::InputTag>("rechitsTag"))) {}

TotemT2DQMSource::~TotemT2DQMSource() {}

void TotemT2DQMSource::dqmBeginRun(const edm::Run&, const edm::EventSetup&) {}

void TotemT2DQMSource::bookHistograms(DQMStore::IBooker& ibooker, const edm::Run&, const edm::EventSetup& iSetup) {
  ibooker.cd();
  ibooker.setCurrentFolder("TotemT2");

  bookErrorFlagsHistogram(ibooker);

  const size_t summary_nbinsx = 25, summary_nbinsy = 25;

  for (unsigned int arm = 0; arm <= CTPPSDetId::maxArm; ++arm){
    for (unsigned int pl = 0; pl <= TotemT2DetId::maxPlane; ++pl) {
      const TotemT2DetId detid(arm, pl, 0);
      const TotemT2DetId planeId(detid.planeId());
      planePlots_[planeId] = PlanePlots(ibooker, planeId);
    }
    const TotemT2DetId detid(arm, 0, 0);
    const TotemT2DetId secId(detid.armId());
    sectorPlots_[secId] = SectorPlots(ibooker, secId);
  }

  // build a segmentation helper for the size of histograms previously booked
  segm_ = std::make_unique<TotemT2Segmentation>(iSetup.getData(geometryToken_), summary_nbinsx, summary_nbinsy);
}

void TotemT2DQMSource::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // fill digis information
  for (const auto& ds_digis : iEvent.get(digiToken_)) {
    const TotemT2DetId detid(ds_digis.detId());
    const TotemT2DetId planeId(detid.planeId());
    for (const auto& digi : ds_digis) {
      segm_->fill(planePlots_[planeId].digisMultiplicity->getTH2D(), detid);
      // (void)digi;  //FIXME make use of them

      fillTriggerBitset(detid);
      fillErrorFlagsHistogram(digi);
      fillEdges(digi, detid);
    }
  }

  // fill rechits information
  std::unordered_map<unsigned int, std::set<unsigned int>> planes;
  for (const auto& ds_rechits : iEvent.get(rechitToken_)) {
    const TotemT2DetId detid(ds_rechits.detId());
    const TotemT2DetId planeId(detid.planeId());
    for (const auto& rechit : ds_rechits) {
      segm_->fill(planePlots_[planeId].rechitMultiplicity->getTH2D(), detid);
      (void)rechit;  //FIXME make use of them

      fillActivePlanes(planes, detid);
    }
  }

  for (const auto& plt : sectorPlots_) {
    plt.second.activePlanesCount->Fill(planes[plt.first].size());
  }

  for (unsigned short arm = 0; arm <= CTPPSDetId::maxArm; ++arm)
    for (unsigned short plane = 0; plane <= 1; ++plane)
      for (unsigned short id = 0; id <= TotemT2DetId::maxChannel; ++id) {
        const TotemT2DetId detid(arm, plane, id);
        if(areChannelsTriggered(detid)){
          const TotemT2DetId secId(detid.armId());
          segm_->fill(sectorPlots_[secId].triggerEmulator->getTH2D(), detid);
        }
      }

  clearTriggerBitset();
}

void TotemT2DQMSource::fillActivePlanes(std::unordered_map<unsigned int, std::set<unsigned int>>& planes, const TotemT2DetId detid){
  const TotemT2DetId secId(detid.armId());
  unsigned short pl = detid.plane();

  planes[secId].insert(pl);    
  // this metod is not SRP
  sectorPlots_[secId].activePlanes->Fill(pl); // move it to another method
}

void TotemT2DQMSource::fillTriggerBitset(const TotemT2DetId detid){
  const TotemT2DetId secId(detid.armId());
  unsigned short pl = detid.plane();
  unsigned short ch = detid.channel();
  sectorPlots_[secId].hitTilesArray[4*pl + ch] = 1;
}

void TotemT2DQMSource::clearTriggerBitset(){
  for(auto& sectorPlot : sectorPlots_){
    sectorPlot.second.hitTilesArray.reset();
  }
}

bool TotemT2DQMSource::areChannelsTriggered(const TotemT2DetId detid){
  unsigned int channel = detid.channel();

  // prepare mask
  std::bitset<(TotemT2DetId::maxPlane+1) * (TotemT2DetId::maxChannel+1)> mask; 
  // check if plane is even or not
  unsigned int pl = detid.plane() % 2 == 0 ? 0 : 1;
  // set only even or only odd plane bits for this channel
  for(; pl <= TotemT2DetId::maxPlane; pl+=2){
    mask[4*pl+channel] = 1;  
  }
  const TotemT2DetId secId(detid.armId());
  // check how many masked channels were hit
  unsigned int triggeredChannelsNumber = (mask & sectorPlots_[secId].hitTilesArray).count();
  
  return triggeredChannelsNumber >= SectorPlots::MINIMAL_TRIGGER;
}

void TotemT2DQMSource::bookErrorFlagsHistogram(DQMStore::IBooker& ibooker){
  HPTDCErrorFlags_2D = ibooker.book2D("HPTDC Errors?", " HPTDC Errors?", 16, -0.5, 16.5, 9, -0.5, 8.5);
  for (unsigned short error_index = 1; error_index < 16; ++error_index)
    HPTDCErrorFlags_2D->setBinLabel(error_index, "Good EC");
  HPTDCErrorFlags_2D->setBinLabel(16, "Wrong EC");

  int tmpIndex = 0;
  HPTDCErrorFlags_2D->setBinLabel(++tmpIndex, "DB 0 TDC 18", /* axis */ 2);
  HPTDCErrorFlags_2D->setBinLabel(++tmpIndex, "DB 0 TDC 17", /* axis */ 2);
  HPTDCErrorFlags_2D->setBinLabel(++tmpIndex, "DB 0 TDC 16", /* axis */ 2);
  HPTDCErrorFlags_2D->setBinLabel(++tmpIndex, "DB 0 TDC 15", /* axis */ 2);
  HPTDCErrorFlags_2D->setBinLabel(++tmpIndex, "DB 1 TDC 18", /* axis */ 2);
  HPTDCErrorFlags_2D->setBinLabel(++tmpIndex, "DB 1 TDC 17", /* axis */ 2);
  HPTDCErrorFlags_2D->setBinLabel(++tmpIndex, "DB 1 TDC 16", /* axis */ 2);
  HPTDCErrorFlags_2D->setBinLabel(++tmpIndex, "DB 1 TDC 15", /* axis */ 2);
}

void TotemT2DQMSource::fillErrorFlagsHistogram(TotemT2Digi digi){
  // dummy filling, just to test histogram
  int hptdcErrorIndex = std::rand() % 16;
  int verticalIndex = std::rand() % 8;
  HPTDCErrorFlags_2D->Fill(hptdcErrorIndex, verticalIndex);
  (void)digi;
}

void TotemT2DQMSource::fillEdges(TotemT2Digi digi, TotemT2DetId detid){
  const TotemT2DetId secId(detid.armId());
  sectorPlots_[secId].leadingEdge->Fill(HPTDC_BIN_WIDTH_NS * digi.leadingEdge());
  sectorPlots_[secId].trailingEdge->Fill(HPTDC_BIN_WIDTH_NS * digi.trailingEdge());
}

DEFINE_FWK_MODULE(TotemT2DQMSource);
