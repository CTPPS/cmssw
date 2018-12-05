/****************************************************************************
 *
 * This is a part of CTPPS offline software.
 * Authors:
 *   Jan Kašpar
 *   Laurent Forthomme
 *
 ****************************************************************************/

#ifndef DataFormats_ProtonReco_ProtonTrackExtraFwd_h
#define DataFormats_ProtonReco_ProtonTrackExtraFwd_h

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/RefVector.h"

#include <vector>

namespace reco
{
  class ProtonTrackExtra;
  /// Collection of ProtonTrackExtra objects
  typedef std::vector<ProtonTrackExtra> ProtonTrackExtraCollection;
  /// Persistent reference to a ProtonTrackExtra
  typedef edm::Ref<ProtonTrackExtraCollection> ProtonTrackExtraRef;
  /// Reference to a ProtonTrackExtra collection
  typedef edm::RefProd<ProtonTrackExtraCollection> ProtonTrackExtraRefProd;
  /// Vector of references to ProtonTrackExtra in the same collection
  typedef edm::RefVector<ProtonTrackExtraCollection> ProtonTrackExtraRefVector;
}

#endif

