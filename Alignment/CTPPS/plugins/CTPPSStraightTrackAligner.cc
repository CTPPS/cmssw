/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*  Cristian Baldenegro (crisx.baldenegro@gmail.com)
****************************************************************************/

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESWatcher.h"

#include "Alignment/CTPPS/interface/StraightTrackAlignment.h"

#include "Geometry/VeryForwardGeometryBuilder/interface/TotemRPGeometry.h"
#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSReco/interface/TotemRPRecHit.h"
#include "DataFormats/CTPPSReco/interface/TotemRPUVPattern.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondRecHit.h"
#include "DataFormats/CTPPSReco/interface/CTPPSPixelRecHit.h"

/**
 *\brief An EDAnalyzer that runs StraightTrackAlignment.
 **/
class CTPPSStraightTrackAligner : public edm::EDAnalyzer
{
  public:
    CTPPSStraightTrackAligner(const edm::ParameterSet &ps); 
    ~CTPPSStraightTrackAligner() {}

  private:
    unsigned int verbosity;

    edm::EDGetTokenT<edm::DetSetVector<TotemRPUVPattern>> tokenUVPatternsStrip;
    edm::EDGetTokenT<edm::DetSetVector<CTPPSDiamondRecHit>> tokenDiamondHits;
    edm::EDGetTokenT<edm::DetSetVector<CTPPSPixelRecHit>> tokenPixelHits;

    bool worker_initialized;
    StraightTrackAlignment worker;

    edm::ESWatcher<VeryForwardRealGeometryRecord> geometryWatcher;

    virtual void beginJob() override {}

    virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;

    virtual void analyze(const edm::Event &e, const edm::EventSetup &es) override;

    virtual void endJob() override;
};

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

CTPPSStraightTrackAligner::CTPPSStraightTrackAligner(const ParameterSet &ps) : 
  verbosity(ps.getUntrackedParameter<unsigned int>("verbosity", 0)),

  tokenUVPatternsStrip(consumes<DetSetVector<TotemRPUVPattern>>(ps.getParameter<edm::InputTag>("tagUVPatternsStrip"))),
  tokenDiamondHits(consumes<DetSetVector<CTPPSDiamondRecHit>>(ps.getParameter<edm::InputTag>("tagDiamondHits"))),
  tokenPixelHits(consumes<DetSetVector<CTPPSPixelRecHit>>(ps.getParameter<edm::InputTag>("tagPixelHits"))),

  worker_initialized(false),
  worker(ps)
{
}

//----------------------------------------------------------------------------------------------------

void CTPPSStraightTrackAligner::beginRun(edm::Run const&, edm::EventSetup const& es)
{
}

//----------------------------------------------------------------------------------------------------

void CTPPSStraightTrackAligner::analyze(const edm::Event &event, const edm::EventSetup &es)
{
  // check if geometry hasn't changed
  if (geometryWatcher.check(es))
  {
    if (worker_initialized)
      throw cms::Exception("CTPPSStraightTrackAligner") <<
        "CTPPSStraightTrackAligner can't cope with changing geometry - change in event " << event.id() << endl;
  }

  // check if worker already initialised
  if (!worker_initialized)
  {
    worker.Begin(es);
    worker_initialized = true;
  }

  // get input
  Handle<DetSetVector<TotemRPUVPattern>> inputUVPatternsStrip;
  event.getByToken(tokenUVPatternsStrip, inputUVPatternsStrip);

  Handle<DetSetVector<CTPPSDiamondRecHit>> inputDiamondHits;
  event.getByToken(tokenDiamondHits, inputDiamondHits);

  Handle<DetSetVector<CTPPSPixelRecHit>> inputPixelHits;
  event.getByToken(tokenPixelHits, inputPixelHits);

  // feed worker
  worker.ProcessEvent(*inputUVPatternsStrip, *inputDiamondHits, *inputPixelHits);
}

//----------------------------------------------------------------------------------------------------

void CTPPSStraightTrackAligner::endJob()
{
  worker.Finish();
}

DEFINE_FWK_MODULE(CTPPSStraightTrackAligner);
