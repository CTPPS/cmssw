/****************************************************************************
 *
 * This is a part of CTPPS offline software.
 * Authors:
 *   Jan Kašpar
 *   Laurent Forthomme
 *
 ****************************************************************************/

#ifndef DataFormats_ProtonReco_ProtonTrackFwd_h
#define DataFormats_ProtonReco_ProtonTrackFwd_h

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/RefVector.h"

#include <vector>

namespace reco
{
  class ProtonTrack;
  /// Collection of ProtonTrack objects
  typedef std::vector<ProtonTrack> ProtonTrackCollection;
  /// Persistent reference to a ProtonTrack
  typedef edm::Ref<ProtonTrackCollection> ProtonTrackRef;
  /// Reference to a ProtonTrack collection
  typedef edm::RefProd<ProtonTrackCollection> ProtonTrackRefProd;
  /// Vector of references to ProtonTrack in the same collection
  typedef edm::RefVector<ProtonTrackCollection> ProtonTrackRefVector;
}

#endif

