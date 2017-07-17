/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*	Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/CTPPSAlignment/interface/RPAlignmentCorrectionData.h"

#include <Math/RotationZYX.h>

using namespace std;

//----------------------------------------------------------------------------------------------------
RPAlignmentCorrectionData::RPAlignmentCorrectionData(double sh_r, double sh_r_e, double sh_x, double sh_x_e,
  double sh_y, double sh_y_e, double sh_z, double sh_z_e, double rot_x, double rot_x_e, double rot_y,
  double rot_y_e, double rot_z, double rot_z_e) :
  translation(sh_x, sh_y, sh_z), translation_error(sh_x_e, sh_y_e, sh_z_e),
  translation_r1(sh_r), translation_r1_error(sh_r_e),
  rotation_x(rot_x), rotation_y(rot_y), rotation_z(rot_z), 
  rotation_x_error(rot_x_e), rotation_y_error(rot_y_e), rotation_z_error(rot_z_e)
{
}

//----------------------------------------------------------------------------------------------------

RPAlignmentCorrectionData::RPAlignmentCorrectionData(double sh_r1, double sh_r1_e, double sh_r2, double sh_r2_e,
    double sh_x, double sh_x_e, double sh_y, double sh_y_e, double sh_z, double sh_z_e, double rot_z, double rot_z_e) :
  translation(sh_x, sh_y, sh_z), translation_error(sh_x_e, sh_y_e, sh_z_e),
  translation_r1(sh_r1), translation_r2(sh_r2),
  translation_r1_error(sh_r1_e), translation_r2_error(sh_r2_e),
  rotation_x(0.), rotation_y(0.), rotation_z(rot_z), 
  rotation_x_error(0.), rotation_y_error(0.), rotation_z_error(rot_z_e)
{
}

//----------------------------------------------------------------------------------------------------

RPAlignmentCorrectionData::RPAlignmentCorrectionData(double sh_x, double sh_y, double sh_z, double rot_z) : 
  translation(sh_x, sh_y, sh_z), translation_error(0., 0., 0.),
  translation_r1(0.), translation_r2(0.),
  translation_r1_error(0.), translation_r2_error(0.),
  rotation_x(0.), rotation_y(0.), rotation_z(rot_z), 
  rotation_x_error(0.), rotation_y_error(0.), rotation_z_error(0.)
{
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentCorrectionData::add(const RPAlignmentCorrectionData &a, bool sumErrors)
{
  // add translations
  translation.SetX(a.translation.X() + translation.X());
  translation.SetY(a.translation.Y() + translation.Y());
  translation.SetZ(a.translation.Z() + translation.Z());

  if (sumErrors)
  {
    translation_error.SetX(sqrt(a.translation_error.X()*a.translation_error.X() + translation_error.X()*translation_error.X()));
    translation_error.SetY(sqrt(a.translation_error.Y()*a.translation_error.Y() + translation_error.Y()*translation_error.Y()));
    translation_error.SetZ(sqrt(a.translation_error.Z()*a.translation_error.Z() + translation_error.Z()*translation_error.Z()));
  } else {
    translation_error.SetX(a.translation_error.X());
    translation_error.SetY(a.translation_error.Y());
    translation_error.SetZ(a.translation_error.Z());
  }
  
  // add rotations
  rotation_x = a.rotation_x + rotation_x;
  rotation_y = a.rotation_y + rotation_y;
  rotation_z = a.rotation_z + rotation_z;
  if (sumErrors)
  {
    rotation_x_error = sqrt(a.rotation_x_error*a.rotation_x_error + rotation_x_error*rotation_x_error);
    rotation_y_error = sqrt(a.rotation_y_error*a.rotation_y_error + rotation_y_error*rotation_y_error);
    rotation_z_error = sqrt(a.rotation_z_error*a.rotation_z_error + rotation_z_error*rotation_z_error);
  } else {
    rotation_x_error = a.rotation_x_error;
    rotation_y_error = a.rotation_y_error;
    rotation_z_error = a.rotation_z_error;
  }
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentCorrectionData::readoutTranslationsToXY(double d1x, double d1y, double d2x, double d2y)
{
  const double tr_z = translation.z();
  const double D = d1x*d2y - d1y*d2x;
  const double tr_x = (d2y * translation_r1 - d1y * translation_r2) / D;
  const double tr_y = (-d2x * translation_r1 + d1x * translation_r2) / D;

  translation.SetXYZ(tr_x, tr_y, tr_z);

  // TODO: error propagation
  translation_error.SetXYZ(0., 0., 0.);
  /*
  tr_z = translation_error.z();
  translation_error.SetXYZ(translation_r1_error*dx, translation_r1_error*dy, tr_z);
  */
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentCorrectionData::xyTranslationToReadout(double d1x, double d1y, double d2x, double d2y)
{
  double dot1 = d1x*translation.x() + d1y*translation.y();
  translation_r1 = dot1;

  // TODO: why this?
  //translation.SetXYZ(dot*dx, dot*dy, translation.z());

  // there is a very high correlation between x and y components of translation_error
  //double dot_error = sqrt(dx*dx * translation_error.x()*translation_error.x() + dy*dy * translation_error.y()*translation_error.y());
  //double dot_error = sqrt(translation_error.x()*translation_error.x() + translation_error.y()*translation_error.y());
  //translation_r1_error = dot_error;
  translation_r1_error = 0.;

  // TODO: why this?
  //translation_error.SetXYZ(dot_error*dx, dot_error*dy, translation_error.z());

  double dot2 = d2x*translation.x() + d2y*translation.y();
  translation_r2 = dot2;

  //translation_r2_error = dot_error;
  translation_r2_error = 0.;
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentCorrectionData::setTranslationR1(double sh_r1, double sh_r1_e)
{
  translation_r1 = sh_r1;
  translation_r1_error = sh_r1_e;
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentCorrectionData::setTranslationR2(double sh_r2, double sh_r2_e)
{
  translation_r2 = sh_r2;
  translation_r2_error = sh_r2_e;
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentCorrectionData::setTranslationZ(double sh_z, double sh_z_e)
{
  translation.SetZ(sh_z);
  translation_error.SetZ(sh_z_e);
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentCorrectionData::setRotationZ(double rot_z, double rot_z_e)
{
  rotation_z = rot_z;
  rotation_z_error = rot_z_e;
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentCorrectionData::normalizeRotationZ()
{
  rotation_z -= floor( (rotation_z + M_PI) / 2. / M_PI ) * 2. * M_PI;
}


//----------------------------------------------------------------------------------------------------

void RPAlignmentCorrectionData::print() const
{
  printf("shift: r1=%.1f, r2=%.1f, x=%.1f, y=%.1f, z=%.1f, rotation: z=%.1f\n", sh_r1()*1E3, sh_r2()*1E3,
    getTranslation().x()*1E3, getTranslation().y()*1E3, getTranslation().z()*1E3, rot_z()*1E3);
}


