/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*  
****************************************************************************/

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandExponential.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSReco/interface/TotemRPRecHit.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondRecHit.h"
#include "DataFormats/CTPPSReco/interface/CTPPSPixelRecHit.h"

#include "Geometry/VeryForwardRPTopology/interface/RPTopology.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/TotemRPGeometry.h"
#include "Geometry/Records/interface/VeryForwardMisalignedGeometryRecord.h"

#include "Alignment/CTPPS/interface/utilities.h"

#include "TMath.h"
#include "TMatrixD.h"
#include "TVectorD.h"

/**
 *\brief Fast (no G4) proton simulation in within one station.
 * Uses misaligned geometry.
 */
class CTPPSFastLocalSimulation : public edm::EDProducer
{
  public:
	CTPPSFastLocalSimulation(const edm::ParameterSet &);
	virtual ~CTPPSFastLocalSimulation();

  protected:
    /// verbosity level
	unsigned int verbosity;

    /// whether a HepMC description of the proton shall be saved in the event
    bool makeHepMC;

    /// whether the hits of the proton shall be calculated and saved
    bool makeHits;

    /// the list of RPs to simulate
    std::vector<unsigned int> RPs;

    /// number of particles to generate per event
    unsigned int particlesPerEvent;

    /// particle energy and momentum
    double particle_E, particle_p;

    /// the "origin" of tracks, in mm
    double z0;                                  

    /// whether measurement values shall be rounded to the nearest strip
    bool roundToPitch;

    /// in mm
    double pitchStrips, pitchDiamonds, pitchPixels;

    /// size of insensitive margin at sensor's edge facing the beam, in mm
    double insensitiveMarginStrips;

    struct Distribution
    {
      enum Type { dtBox, dtGauss, dtGaussLimit } type;
      double x_mean, x_width, x_min, x_max;
      double y_mean, y_width, y_min, y_max;

      Distribution(const edm::ParameterSet &);

      void Generate(CLHEP::HepRandomEngine &rndEng, double &x, double &y);
    };

    /// position parameters in mm
    Distribution position_dist;

    /// angular parameters in rad
    Distribution angular_dist;

    //---------- internal parameters ----------

    /// v position of strip 0, in mm
    double stripZeroPosition;

    edm::ESHandle<TotemRPGeometry> geometry;

    void GenerateTrack(unsigned int pi, CLHEP::HepRandomEngine &rndEng,
        HepMC::GenEvent* gEv, std::unique_ptr<edm::DetSetVector<TotemRPRecHit>> &stripHitColl,
        std::unique_ptr<edm::DetSetVector<CTPPSDiamondRecHit>> &diamondHitColl,
        std::unique_ptr<edm::DetSetVector<CTPPSPixelRecHit>> &pixelHitColl
      );

    //---------- framework methods ----------

    virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;

	virtual void produce(edm::Event &, const edm::EventSetup&) override;
};

//----------------------------------------------------------------------------------------------------

using namespace edm;
using namespace std;
using namespace CLHEP;
using namespace HepMC;

//----------------------------------------------------------------------------------------------------

CTPPSFastLocalSimulation::CTPPSFastLocalSimulation(const edm::ParameterSet &ps) :
  verbosity(ps.getUntrackedParameter<unsigned int>("verbosity", 0)),

  makeHepMC(ps.getParameter<bool>("makeHepMC")),
  makeHits(ps.getParameter<bool>("makeHits")),

  RPs(ps.getParameter< vector<unsigned int> >("RPs")),

  particlesPerEvent(ps.getParameter<unsigned int>("particlesPerEvent")),
  particle_E(ps.getParameter<double>("particle_E")),
  particle_p(ps.getParameter<double>("particle_p")),
  z0(ps.getParameter<double>("z0")),

  roundToPitch(ps.getParameter<bool>("roundToPitch")),
  pitchStrips(ps.getParameter<double>("pitchStrips")),
  pitchDiamonds(ps.getParameter<double>("pitchDiamonds")),
  pitchPixels(ps.getParameter<double>("pitchPixels")),

  insensitiveMarginStrips(ps.getParameter<double>("insensitiveMarginStrips")),

  position_dist(ps.getParameterSet("position_distribution")),
  angular_dist(ps.getParameterSet("angular_distribution"))
{
  // v position of strip 0
  stripZeroPosition = RPTopology::last_strip_to_border_dist_ + (RPTopology::no_of_strips_-1)*RPTopology::pitch_
    - RPTopology::y_width_/2.;

  // register the output
  if (makeHepMC)
    produces<HepMCProduct>();

  if (makeHits)
  {
    produces<DetSetVector<TotemRPRecHit>>();
    produces<DetSetVector<CTPPSDiamondRecHit>>();
    produces<DetSetVector<CTPPSPixelRecHit>>();
  }
}

//----------------------------------------------------------------------------------------------------

CTPPSFastLocalSimulation::Distribution::Distribution(const edm::ParameterSet &ps)
{
  // get type
  string typeName = ps.getParameter<string>("type");
  if (!typeName.compare("box"))
    type = dtBox;
  else if (!typeName.compare("gauss"))
      type = dtGauss;
    else if (!typeName.compare("gauss-limit"))
        type = dtGaussLimit;
      else
        throw cms::Exception("CTPPSFastLocalSimulation") << "Unknown distribution type `" << typeName << "'.";

  x_mean = ps.getParameter<double>("x_mean");
  x_width = ps.getParameter<double>("x_width");
  x_min = ps.getParameter<double>("x_min");
  x_max = ps.getParameter<double>("x_max");

  y_mean = ps.getParameter<double>("y_mean");
  y_width = ps.getParameter<double>("y_width");
  y_min = ps.getParameter<double>("y_min");
  y_max = ps.getParameter<double>("y_max");
}

//----------------------------------------------------------------------------------------------------

CTPPSFastLocalSimulation::~CTPPSFastLocalSimulation()
{
}

//----------------------------------------------------------------------------------------------------

void CTPPSFastLocalSimulation::beginRun(edm::Run const&, edm::EventSetup const& es)
{
  // get geometry
  es.get<VeryForwardMisalignedGeometryRecord>().get(geometry);
}

//----------------------------------------------------------------------------------------------------

void CTPPSFastLocalSimulation::Distribution::Generate(CLHEP::HepRandomEngine &rndEng, double &x, double &y)
{
  switch (type)
  {
    case dtBox:
      x = x_mean + x_width * (rndEng.flat() - 0.5);
      y = y_mean + y_width * (rndEng.flat() - 0.5);
      break;

    case dtGauss:
      x = x_mean + RandGauss::shoot(&rndEng) * x_width;
      y = y_mean + RandGauss::shoot(&rndEng) * y_width;
      break;

    case dtGaussLimit:
      {
        double u_x = rndEng.flat(), u_y = rndEng.flat();

        double cdf_x_min = (1. + TMath::Erf((x_min - x_mean) / x_width / sqrt(2.))) / 2.;
        double cdf_x_max = (1. + TMath::Erf((x_max - x_mean) / x_width / sqrt(2.))) / 2.;
        double a_x = cdf_x_max - cdf_x_min, b_x = cdf_x_min;

        double cdf_y_min = (1. + TMath::Erf((y_min - y_mean) / y_width / sqrt(2.))) / 2.;
        double cdf_y_max = (1. + TMath::Erf((y_max - y_mean) / y_width / sqrt(2.))) / 2.;
        double a_y = cdf_y_max - cdf_y_min, b_y = cdf_y_min;

        x = x_mean + x_width * sqrt(2.) * TMath::ErfInverse(2.*(a_x * u_x + b_x) - 1.); 
        y = y_mean + y_width * sqrt(2.) * TMath::ErfInverse(2.*(a_y * u_y + b_y) - 1.); 
      }

      break;

    default:
      x = y = 0.;
  }
}

//----------------------------------------------------------------------------------------------------

void CTPPSFastLocalSimulation::GenerateTrack(unsigned int idx, CLHEP::HepRandomEngine &rndEng,
  HepMC::GenEvent* gEv, unique_ptr<edm::DetSetVector<TotemRPRecHit>> &stripHitColl,
  unique_ptr<edm::DetSetVector<CTPPSDiamondRecHit>> &diamondHitColl,
  unique_ptr<edm::DetSetVector<CTPPSPixelRecHit>> &pixelHitColl)
{
  // generate track
  double bx = 0., by = 0., ax = 0., ay = 0.;
  position_dist.Generate(rndEng, bx, by);
  angular_dist.Generate(rndEng, ax, ay);

  if (verbosity > 5)
    printf("\tax = %.3f mrad, bx = %.3f mm, ay = %.3f mrad, by = %.3f mm, z0 = %.3f m\n", ax*1E3, bx, ay*1E3, by, z0*1E-3);

  // add HepMC track description
  if (makeHepMC)
  { 
    GenVertex* gVx = new GenVertex(HepMC::FourVector(bx, by, z0, 0.));
    gEv->add_vertex(gVx);

    GenParticle* gPe; 
    double az = sqrt(1. - ax*ax - ay*ay);
    gPe = new GenParticle(HepMC::FourVector(particle_p*ax, particle_p*ay, particle_p*az, particle_E), 2212, 1); // add a proton in final state
    gPe->suggest_barcode(idx + 1);
    gVx->add_particle_out(gPe);
  }

  if (makeHits)
  {
    // check all sensors known to geometry
    for (TotemRPGeometry::mapType::const_iterator it = geometry->beginDet(); it != geometry->endDet(); ++it)
    {
      // get RP decimal id
      CTPPSDetId detId(it->first);
      unsigned int decRPId = detId.arm()*100 + detId.station()*10 + detId.rp();

      // stop if the RP is not selected
      if (find(RPs.begin(), RPs.end(), decRPId) == RPs.end())
        continue;

      // keep only 1 diamond channel to represent 1 plane
      if (detId.subdetId() == CTPPSDetId::sdTimingDiamond)
      {
        CTPPSDiamondDetId channelId(it->first);
        if (channelId.channel() != 0)
          continue;
      }

      if (verbosity > 5)
      {
        printf("        ");
        PrintId(it->first);
        printf(": ");
      }

      // determine the track impact point (in global coordinates)
      // !! this assumes that local axes (1, 0, 0) and (0, 1, 0) describe the sensor surface
      Hep3Vector gl_o = geometry->LocalToGlobal(detId, Hep3Vector(0, 0, 0));
      Hep3Vector gl_a1 = geometry->LocalToGlobal(detId, Hep3Vector(1, 0, 0)) - gl_o;
      Hep3Vector gl_a2 = geometry->LocalToGlobal(detId, Hep3Vector(0, 1, 0)) - gl_o;

      TMatrixD A(3, 3);
      TVectorD B(3);
      A(0, 0) = ax; A(0, 1) = -gl_a1.x(); A(0, 2) = -gl_a2.x(); B(0) = gl_o.x() - bx;
      A(1, 0) = ay; A(1, 1) = -gl_a1.y(); A(1, 2) = -gl_a2.y(); B(1) = gl_o.y() - by;
      A(2, 0) = 1.; A(2, 1) = -gl_a1.z(); A(2, 2) = -gl_a2.z(); B(2) = gl_o.z() - z0;
      TMatrixD Ai(3, 3);
      Ai = A.Invert();
      TVectorD P(3);
      P = Ai * B;
      //printf("            de z = %f, p1 = %f, p2 = %f\n", P(0), P(1), P(2));

      double de_z = P(0);
      Hep3Vector h_glo(ax * de_z + bx, ay * de_z + by, de_z + z0);
      h_glo.setX(ax * de_z + bx);
      h_glo.setY(ay * de_z + by);
      h_glo.setZ(de_z + z0);
      //printf("            h_glo: x = %f, y = %f, z = %f\n", h_glo.x(), h_glo.y(), h_glo.z());

      // hit in local coordinates
      Hep3Vector h_loc = geometry->GlobalToLocal(detId, h_glo);
      //printf("            h_loc: c1 = %f, c2 = %f, c3 = %f\n", h_loc.x(), h_loc.y(), h_loc.z());

      // strips
      if (detId.subdetId() == CTPPSDetId::sdTrackingStrip)
      {
        double u = h_loc.x();
        double v = h_loc.y();

        if (verbosity > 5)
          printf("            u=%+8.4f, v=%+8.4f", u, v);

        // is it within detector?
        if (!RPTopology::IsHit(u, v, insensitiveMarginStrips))
        {
          if (verbosity > 5)
            printf(" | no hit\n");
          continue;
        } 

        // round the measurement
        if (roundToPitch)
        {
          double m = stripZeroPosition - v;
          signed int strip = (int) floor(m / pitchStrips + 0.5);

          v = stripZeroPosition - pitchStrips * strip;

          if (verbosity > 5)
            printf(" | strip=%+4i", strip);
        }

        double sigma = pitchStrips / sqrt(12.);

        if (verbosity > 5)
          printf(" | m=%+8.4f, sigma=%+8.4f\n", v, sigma);

        DetSet<TotemRPRecHit> &hits = stripHitColl->find_or_insert(detId);
        hits.push_back(TotemRPRecHit(v, sigma));
      }

      // diamonds
      if (detId.subdetId() == CTPPSDetId::sdTimingDiamond)
      {
        if (roundToPitch)
        {
          h_loc.setX( pitchDiamonds * floor(h_loc.x()/pitchDiamonds + 0.5) );
        }

        if (verbosity > 5)
          printf("            m = %.3f\n", h_loc.x());

        const double width = pitchDiamonds;

        DetSet<CTPPSDiamondRecHit> &hits = diamondHitColl->find_or_insert(detId);
        hits.push_back(CTPPSDiamondRecHit(h_loc.x(), width, 0., 0., 0., 0., 0, 0, false));
      }

      // pixels
      if (detId.subdetId() == CTPPSDetId::sdTrackingPixel)
      {
        if (roundToPitch)
        {
          h_loc.setX( pitchPixels * floor(h_loc.x()/pitchPixels + 0.5) );
          h_loc.setY( pitchPixels * floor(h_loc.y()/pitchPixels + 0.5) );
        }

        if (verbosity > 5)
          printf("            m1 = %.3f, m2 = %.3f\n", h_loc.x(), h_loc.y());

        const double sigma = pitchPixels / sqrt(12.);

        const LocalPoint lp(h_loc.x(), h_loc.y(), h_loc.z());
        const LocalError le(sigma, 0., sigma);

        DetSet<CTPPSPixelRecHit> &hits = pixelHitColl->find_or_insert(detId);
        hits.push_back(CTPPSPixelRecHit(lp, le));
      }
    }
  } 
}

//----------------------------------------------------------------------------------------------------

void CTPPSFastLocalSimulation::produce(edm::Event &event, const edm::EventSetup &es)
{
  if (verbosity > 2)
    printf(">> CTPPSFastLocalSimulation::produce > event %llu\n", event.id().event());

  Service<edm::RandomNumberGenerator> rng;
  HepRandomEngine &rndEng = rng->getEngine(event.streamID());

  if (verbosity > 5)
    printf("\tseed = %li\n", rndEng.getSeed());

  // initialize products
  GenEvent* gEv = new GenEvent();
  gEv->set_event_number(event.id().event());

  unique_ptr<DetSetVector<TotemRPRecHit>> stripHitColl(new DetSetVector<TotemRPRecHit>());
  unique_ptr<DetSetVector<CTPPSDiamondRecHit>> diamondHitColl(new DetSetVector<CTPPSDiamondRecHit>());
  unique_ptr<DetSetVector<CTPPSPixelRecHit>> pixelHitColl(new DetSetVector<CTPPSPixelRecHit>());

  // run particle loop
  for (unsigned int pi = 0; pi < particlesPerEvent; pi++)
  {
    if (verbosity > 5)
      printf("    generating track %u\n", pi);

    GenerateTrack(pi, rndEng, gEv, stripHitColl, diamondHitColl, pixelHitColl);
  }

  // save products
  if (makeHepMC)
  { 
    unique_ptr<HepMCProduct> hepMCoutput(new HepMCProduct());
    hepMCoutput->addHepMCData(gEv);
    event.put(move(hepMCoutput));
  }

  if (makeHits)
  {
    event.put(move(stripHitColl));
    event.put(move(diamondHitColl));
    event.put(move(pixelHitColl));
  }
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(CTPPSFastLocalSimulation);
