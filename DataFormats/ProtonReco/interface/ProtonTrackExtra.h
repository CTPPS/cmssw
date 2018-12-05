/****************************************************************************
 *
 * This is a part of CTPPS offline software.
 * Authors:
 *   Jan Kašpar
 *   Laurent Forthomme
 *
 ****************************************************************************/

#ifndef DataFormats_ProtonReco_ProtonTrackExtra_h
#define DataFormats_ProtonReco_ProtonTrackExtra_h

#include <set>
#include "DataFormats/Common/interface/Ref.h"

namespace reco
{
  class ProtonTrackExtra
  {
    public:
      typedef std::set<unsigned int> RPList;
      /// Type of reconstruction for this track
      enum class ReconstructionMethod { invalid = -1, singleRP, multiRP };
      /// LHC sector for this track
      enum class LHCSector { invalid = -1, sector45, sector56 };

      /// Empty (invalid track extra info) constructor
      ProtonTrackExtra();
      /// Default constructor
      ProtonTrackExtra( bool valid, const ReconstructionMethod& method, const LHCSector& sector, const RPList& rp_list );

      /// Set the flag for the fit validity
      void setValidFit( bool valid = true ) { valid_fit_ = valid; }
      /// Flag for the fit validity
      bool validFit() const { return valid_fit_; }
      /// Set the reconstruction method for this track
      void setMethod( const ReconstructionMethod& method ) { method_ = method; }
      /// Reconstruction method for this track
      ReconstructionMethod method() const { return method_; }
      /// Set the LHC sector (0=sector45, 1=sector56)
      void setSector( const LHCSector& sector ) { sector_ = sector; }
      /// LHC sector for this track (0=sector45, 1=sector56)
      LHCSector sector() const { return sector_; }
      /// Store the list of RP tracks that contributed to this global track
      void setContributingRPs( const RPList& list ) { contributing_rp_ids_ = list; }
      /// List of RP tracks that contributed to this global track
      const RPList& contributingRPs() const { return contributing_rp_ids_; }

    private:
      bool valid_fit_;
      ReconstructionMethod method_;
      LHCSector sector_;
      RPList contributing_rp_ids_;
  };
}

#endif

