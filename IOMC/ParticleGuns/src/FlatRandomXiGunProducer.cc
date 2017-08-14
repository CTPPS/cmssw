#include "IOMC/ParticleGuns/interface/FlatRandomXiGunProducer.h"

namespace edm
{
  FlatRandomXiGunProducer::FlatRandomXiGunProducer( const edm::ParameterSet& iConfig ) :
    partGunParams_( iConfig.getParameter<edm::ParameterSet>( "PGunParameters" ) ),
    partIds_( partGunParams_.getParameter< std::vector<int> >( "PartID" ) ),
    sqrtS_  ( partGunParams_.getParameter<double>( "SqrtS" ) ),
    minXi_  ( partGunParams_.getParameter<double>( "MinXi" ) ),
    maxXi_  ( partGunParams_.getParameter<double>( "MaxXi" ) ),
    minT_   ( partGunParams_.getParameter<double>( "MinT" ) ),
    maxT_   ( partGunParams_.getParameter<double>( "MaxT" ) ),
    minPhi_ ( partGunParams_.getUntrackedParameter<double>( "MinPhi", -M_PI ) ),
    maxPhi_ ( partGunParams_.getUntrackedParameter<double>( "MaxPhi", +M_PI ) )
  {
    produces<edm::HepMCProduct>( "unsmeared" );
    produces<GenEventInfoProduct>();
  }

  FlatRandomXiGunProducer::~FlatRandomXiGunProducer()
  {}

  void
  FlatRandomXiGunProducer::produce( edm::Event& iEvent, const edm::EventSetup& iSetup )
  {
    edm::Service<edm::RandomNumberGenerator> rng;
    if ( !rng.isAvailable() ) {
      throw cms::Exception("Configuration")
        << "The module inheriting from FlatRandomXiGunProducer requires the\n"
           "RandomNumberGeneratorService, which appears to be absent. Please\n"
           "add that service to your configuration or remove the modules that"
           "require it.";
    }

    CLHEP::HepRandomEngine* rnd = &( rng->getEngine( iEvent.streamID() ) );

    std::unique_ptr<edm::HepMCProduct> pOut( new edm::HepMCProduct() );

    // generate event
    auto pEvt = new HepMC::GenEvent(); // cleanup in HepMCProduct

    // generate vertex
    auto pVertex = new HepMC::GenVertex(); // cleanup in HepMCProduct
    pVertex->set_position( HepMC::FourVector() );

    unsigned short barcode = 0;
    for ( const auto& part : partIds_ ) {
      auto part_data = pdgTable_->particle( HepPDT::ParticleID( abs( part ) ) );
      const double mass = part_data->mass().value();

      int dir = ( CLHEP::RandFlat::shoot( rnd )<0.5 ) ? -1 : 1;

      auto p = new HepMC::GenParticle( shoot( rnd, mass, dir ), part, 1 ); // cleanup in HepMCProduct
      p->suggest_barcode( barcode );
      pVertex->add_particle_out( p );

      barcode++;
    }

    pEvt->add_vertex( pVertex );
    pEvt->set_event_number( iEvent.id().event() );
    pOut->addHepMCData( pEvt );

    iEvent.put( std::move( pOut ), "unsmeared" );

    std::unique_ptr<GenEventInfoProduct> pGenEventInfo( new GenEventInfoProduct( pEvt ) );
    iEvent.put( std::move( pGenEventInfo ) );
  }

  //----------------------------------------------------------------------------------------------------

  HepMC::FourVector
  FlatRandomXiGunProducer::shoot( CLHEP::HepRandomEngine* rnd, double mass, int z_direction )
  {
    // generate xi
    const double xi = CLHEP::RandFlat::shoot( rnd, minXi_, maxXi_ );
    const double e_0 = sqrtS_*0.5, e_part = e_0 * ( 1.-xi );
    const double p_0 = sqrt( e_0*e_0-mass*mass ), p = sqrt( e_part*e_part-mass*mass );

    // generate phi
    const double phi = CLHEP::RandFlat::shoot( rnd, minPhi_, maxPhi_ );

    // generate t
    const double min_t = std::max( minT_, -2. * ( sqrt( e_0*e_0-mass*mass ) * p - e_0*e_part + mass*mass ) );
    const double t = CLHEP::RandFlat::shoot( rnd, min_t, maxT_ );

    double theta = acos( ( -0.5*t - mass*mass + e_part*e_0 ) / ( p*p_0 ) );
    if ( z_direction < 0 ) theta = M_PI - theta;

    const double px = p * cos( phi ) * sin( theta ) * z_direction;
    const double py = p * sin( phi ) * sin( theta );
    const double pz = p * cos( theta );
    return HepMC::FourVector( px, py, pz, e_part );
  }

  void
  FlatRandomXiGunProducer::beginRun( const edm::Run&, const edm::EventSetup& iSetup )
  {
    iSetup.getData( pdgTable_ );
  }

  void
  FlatRandomXiGunProducer::endRun( const edm::Run&, const edm::EventSetup& )
  {}

  void
  FlatRandomXiGunProducer::fillDescriptions( edm::ConfigurationDescriptions& descriptions )
  {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault( desc );
  }
}
