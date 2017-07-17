/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#include "Alignment/CTPPS/interface/AlignmentTask.h"
#include "Alignment/CTPPS/interface/AlignmentConstraint.h"

#include "Geometry/VeryForwardGeometryBuilder/interface/TotemRPGeometry.h"

#include <algorithm>

//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;
using namespace CLHEP;

//----------------------------------------------------------------------------------------------------


AlignmentTask::AlignmentTask(const ParameterSet& ps) :
  resolveShR(ps.getParameter<bool>("resolveShR")),
  resolveShZ(ps.getParameter<bool>("resolveShZ")),
  resolveRotZ(ps.getParameter<bool>("resolveRotZ")),

  useExtendedRotZConstraint(ps.getParameter<bool>("useExtendedRotZConstraint")),
  useZeroThetaRotZConstraint(ps.getParameter<bool>("useZeroThetaRotZConstraint")),
  useExtendedShZConstraints(ps.getParameter<bool>("useExtendedShZConstraints")),
  useEqualMeanUMeanVRotZConstraint(ps.getParameter<bool>("useEqualMeanUMeanVRotZConstraint")),
  oneRotZPerPot(ps.getParameter<bool>("oneRotZPerPot")),

  homogeneousConstraints(ps.getParameterSet("homogeneousConstraints")),
  fixedDetectorsConstraints(ps.getParameterSet("fixedDetectorsConstraints")),
  standardConstraints(ps.getParameterSet("standardConstraints"))
{
  if (resolveShR)
  {
    quantityClasses.push_back(qcShR1);
    quantityClasses.push_back(qcShR2);
  }

  if (resolveShZ)
  {
    quantityClasses.push_back(qcShZ);
  }

  if (resolveRotZ)
  {
    quantityClasses.push_back(qcRotZ);
  }
}

//----------------------------------------------------------------------------------------------------

void AlignmentTask::BuildGeometry(const vector<unsigned int> &rpDecIds,
  const vector<unsigned int> &excludedSensors, const TotemRPGeometry *input, double z0, AlignmentGeometry &geometry)
{
  geometry.clear();
  geometry.z0 = z0; 

  // traverse full known geometry
  for (auto it = input->beginDet(); it != input->endDet(); ++it)
  {
    // skip excluded sensors
    if (find(excludedSensors.begin(), excludedSensors.end(), it->first) != excludedSensors.end())
      continue;

    // is RP selected?
    const CTPPSDetId detId(it->first);
    const unsigned int rpDecId = 100*detId.arm() + 10*detId.station() + detId.rp();
    if (find(rpDecIds.begin(), rpDecIds.end(), rpDecId) == rpDecIds.end())
      continue;

    // extract geometry data
    Hep3Vector c = input->LocalToGlobal(detId, Hep3Vector(0., 0., 0.));
    Hep3Vector d1 = input->LocalToGlobal(detId, Hep3Vector(1., 0., 0.)) - c;
    Hep3Vector d2 = input->LocalToGlobal(detId, Hep3Vector(0., 1., 0.)) - c;

    // for strips: is it U plane?
    bool isU = false;
    if (detId.subdetId() == CTPPSDetId::sdTrackingStrip)
    {
      TotemRPDetId stripDetId(it->first);
      unsigned int rpNum = stripDetId.rp();
      unsigned int plNum = stripDetId.plane();
      isU = (plNum % 2 != 0);
      if (rpNum == 2 || rpNum == 3)
        isU = !isU;
    }

    DetGeometry dg(c.z() - z0, c.x(), c.y(), isU);
    dg.SetDirection(1, d1.x(), d1.y(), d1.z());
    dg.SetDirection(2, d2.x(), d2.y(), d2.z());
    geometry.Insert(it->first, dg);
  }
}

//----------------------------------------------------------------------------------------------------

void AlignmentTask::BuildIndexMaps()
{
  // remove old mapping
  mapMeasurementIndeces.clear();
  mapQuantityIndeces.clear();

  // loop over all classes
  for (const auto &qcl : quantityClasses)
  {
    // create entry for this class
    mapMeasurementIndeces[qcl];

    // loop over all sensors
    unsigned int idxMeas = 0;
    unsigned int idxQuan = 0;
    for (const auto &git : geometry)
    {
      const unsigned int detId = git.first;
      const unsigned int subdetId = CTPPSDetId(git.first).subdetId();

      // update measurement map
      if (qcl == qcShR1)
      {
        if (subdetId == CTPPSDetId::sdTimingDiamond) mapMeasurementIndeces[qcl][{detId, 1}] = idxMeas++;
        if (subdetId == CTPPSDetId::sdTrackingPixel) mapMeasurementIndeces[qcl][{detId, 1}] = idxMeas++;
      }

      if (qcl == qcShR2)
      {
        if (subdetId == CTPPSDetId::sdTrackingStrip) mapMeasurementIndeces[qcl][{detId, 2}] = idxMeas++;
        if (subdetId == CTPPSDetId::sdTrackingPixel) mapMeasurementIndeces[qcl][{detId, 2}] = idxMeas++;
      }

      if (qcl == qcShZ)
      {
        if (subdetId == CTPPSDetId::sdTrackingStrip) mapMeasurementIndeces[qcl][{detId, 2}] = idxMeas++;
        if (subdetId == CTPPSDetId::sdTrackingPixel) mapMeasurementIndeces[qcl][{detId, 1}] = idxMeas++;
        if (subdetId == CTPPSDetId::sdTrackingPixel) mapMeasurementIndeces[qcl][{detId, 2}] = idxMeas++;
        if (subdetId == CTPPSDetId::sdTimingDiamond) mapMeasurementIndeces[qcl][{detId, 1}] = idxMeas++;
      }

      if (qcl == qcRotZ)
      {
        if (subdetId == CTPPSDetId::sdTrackingStrip) mapMeasurementIndeces[qcl][{detId, 2}] = idxMeas++;
        if (subdetId == CTPPSDetId::sdTrackingPixel) mapMeasurementIndeces[qcl][{detId, 1}] = idxMeas++;
        if (subdetId == CTPPSDetId::sdTrackingPixel) mapMeasurementIndeces[qcl][{detId, 2}] = idxMeas++;
        if (subdetId == CTPPSDetId::sdTimingDiamond) mapMeasurementIndeces[qcl][{detId, 1}] = idxMeas++;
      }

      // update quantity map
      if (qcl == qcShR1)
      {
        if (subdetId == CTPPSDetId::sdTimingDiamond) mapQuantityIndeces[qcl][detId] = idxQuan++;
        if (subdetId == CTPPSDetId::sdTrackingPixel) mapQuantityIndeces[qcl][detId] = idxQuan++;
      }

      if (qcl == qcShR2)
      {
        if (subdetId == CTPPSDetId::sdTrackingStrip) mapQuantityIndeces[qcl][detId] = idxQuan++;
        if (subdetId == CTPPSDetId::sdTrackingPixel) mapQuantityIndeces[qcl][detId] = idxQuan++;
      }

      if (qcl == qcShZ)
      {
        if (subdetId == CTPPSDetId::sdTrackingStrip) mapQuantityIndeces[qcl][detId] = idxQuan++;
        if (subdetId == CTPPSDetId::sdTimingDiamond) mapQuantityIndeces[qcl][detId] = idxQuan++;
        if (subdetId == CTPPSDetId::sdTrackingPixel) mapQuantityIndeces[qcl][detId] = idxQuan++;
      }

      if (qcl == qcRotZ)
      {
        if (subdetId == CTPPSDetId::sdTrackingStrip) mapQuantityIndeces[qcl][detId] = idxQuan++;
        if (subdetId == CTPPSDetId::sdTimingDiamond) mapQuantityIndeces[qcl][detId] = idxQuan++;
        if (subdetId == CTPPSDetId::sdTrackingPixel) mapQuantityIndeces[qcl][detId] = idxQuan++;
      }
    }
  }

  {
    printf(">> AlignmentTask::BuildIndexMaps > measurement indeces\n");

    for (const auto &cli : mapMeasurementIndeces)
    {
      printf("class %u\n", cli.first);
    
      for (const auto &it : cli.second)
      {
        printf("    detId=%u, dirIdx=%u --> matrixIndex=%u\n", it.first.detId, it.first.dirIdx, it.second);
      }
    }

    printf(">> AlignmentTask::BuildIndexMaps > quantity indeces\n");

    for (const auto &cli : mapQuantityIndeces)
    {
      printf("class %u\n", cli.first);
    
      for (const auto &it : cli.second)
      {
        printf("    detId=%u --> matrixIndex=%u\n", it.first, it.second);
      }
    }
  }
}

//----------------------------------------------------------------------------------------------------

signed int AlignmentTask::GetMeasurementIndex(QuantityClass cl, unsigned int detId, unsigned int dirIdx) const
{
  auto clit = mapMeasurementIndeces.find(cl);
  if (clit == mapMeasurementIndeces.end())
    return -1;

  auto it = clit->second.find({detId, dirIdx});
  if (it == clit->second.end())
    return -1;

  return it->second;
}

//----------------------------------------------------------------------------------------------------

signed int AlignmentTask::GetQuantityIndex(QuantityClass cl, unsigned int detId) const
{
  auto clit = mapQuantityIndeces.find(cl);
  if (clit == mapQuantityIndeces.end())
    return -1;

  auto it = clit->second.find(detId);
  if (it == clit->second.end())
    return -1;

  return it->second;
}

//----------------------------------------------------------------------------------------------------

string AlignmentTask::QuantityClassTag(QuantityClass qc) const
{
  switch (qc)
  {
    case qcShR1: return "ShR1";
    case qcShR2: return "ShR2";
    case qcShZ: return "ShZ";
    case qcRotZ: return "RotZ";
  }

  throw cms::Exception("AlignmentTask::QuantityClassTag") << "Unknown quantity class " << qc << ".";
}

//----------------------------------------------------------------------------------------------------

unsigned int AlignmentTask::MeasurementsOfClass(QuantityClass qc) const
{
  auto it = mapMeasurementIndeces.find(qc);
  return it->second.size();
}

//----------------------------------------------------------------------------------------------------

unsigned int AlignmentTask::QuantitiesOfClass(QuantityClass qc) const
{
  auto it = mapQuantityIndeces.find(qc);
  return it->second.size();
}

//----------------------------------------------------------------------------------------------------

void AlignmentTask::BuildHomogeneousConstraints(vector<AlignmentConstraint> &constraints) const
{
  printf(">> AlignmentTask::BuildHomogeneousConstraints > not yet ported, sorry!\n");
  throw 1;

#if 0
  for (unsigned int cl = 0; cl < quantityClasses.size(); cl++)
  {
    unsigned int size = ConstraintsForClass(quantityClasses[cl]);
    const string &tag = QuantityClassTag(quantityClasses[cl]);
    
    // just one basic constraint
    if (oneRotZPerPot && quantityClasses[cl] == qcRotZ)
      size = 1;

    // get constraint values
    char buf[20];
    sprintf(buf, "%s_values", tag.c_str());
    vector<double> values(homogeneousConstraints.getParameter< vector<double> >(buf));
    if (values.size() < size)
      throw cms::Exception("AlignmentTask::BuildHomogeneousConstraints") <<
        "Invalid number of constraint values for " << tag << ". Given " << values.size() <<
        ", expected " << size << ".";

    for (unsigned int j = 0; j < size; j++)
    {
      // prepare a constraint with coefficient vectors
      AlignmentConstraint ac;
      ac.forClass = quantityClasses[cl];
      for (unsigned int i = 0; i < quantityClasses.size(); i++)
      {
        ac.coef[quantityClasses[i]].ResizeTo(QuantitiesOfClass(quantityClasses[i]));
        ac.coef[quantityClasses[i]].Zero();
      }
      ac.val = values[j];
      ac.extended = false;

      unsigned int indeces = QuantitiesOfClass(quantityClasses[cl]);

      for (unsigned int idx = 0; idx < indeces; ++idx)
      {
        double &coef = ac.coef[quantityClasses[cl]][idx];
        const DetGeometry &dt = geometry.FindByMatrixIndex(idx)->second;
        double sc = -dt.d2x*dt.sy + dt.d2y*dt.sx;

        switch (quantityClasses[cl])
        {
          case qcShR:
            switch (j)
            {
              case 0: ac.name = "ShR: z*dx"; coef = dt.z * dt.d2x; break;
              case 1: ac.name = "ShR: dx"; coef = dt.d2x; break;
              case 2: ac.name = "ShR: z*dy"; coef = dt.z * dt.d2y; break;
              case 3: ac.name = "ShR: dy"; coef = dt.d2y; break;
            }
            break;

          case qcShZ:
            switch (j)
            {
              case 0: ac.name = "ShZ: z"; coef = dt.z; break;
              case 1: ac.name = "ShZ: 1"; coef = 1.; break;
              case 2: ac.name = "ShZ: z for V-det"; coef = (dt.isU) ? 0. : dt.z; ac.extended = true; break;
              case 3: ac.name = "ShZ: 1 for V-det"; coef = (dt.isU) ? 0. : 1.; ac.extended = true; break;
            }
            break;
        
          case qcRPShZ:
            switch (j)
            {
              case 0: ac.name = "RPShZ: 1"; coef = 1.; break;
              case 1: ac.name = "RPShZ: z"; coef = dt.z; ac.extended = true; break;
            }
            break;
          
          case qcRotZ:
            unsigned int je = j;
            if (!useExtendedRotZConstraint)
              je *= 2;
            switch (je)
            {
              case 0: ac.name = "RotZ: 1 all det"; coef = 1.; break;
              case 1: ac.name = "RotZ: 1 V-det"; coef = (dt.isU) ? 0. : 1.; ac.extended = true; break;
              case 2: ac.name = "RotZ: z all det"; coef = dt.z; ac.extended = true; break;
              case 3: ac.name = "RotZ: z V-det"; coef = (dt.isU) ? 0. : dt.z; ac.extended = true; break;
            }
            ac.coef[qcShR][idx] = sc * ac.coef[qcRotZ][idx];
            break;
        }
      }

      constraints.push_back(ac);
    } 
    
    // only 1 rot_z per RP
    if (oneRotZPerPot && quantityClasses[cl] == qcRotZ)
      BuildOneRotZPerPotConstraints(constraints);
  }
#endif
}

//----------------------------------------------------------------------------------------------------

void AlignmentTask::BuildFixedDetectorsConstraints(vector<AlignmentConstraint> &constraints) const
{
  for (auto &quantityClass : quantityClasses)
  {
    // get input
    const string &tag = QuantityClassTag(quantityClass);

    const ParameterSet &classSettings = fixedDetectorsConstraints.getParameterSet(tag.c_str());
    vector<unsigned int> ids(classSettings.getParameter<vector<unsigned int>>("ids"));
    vector<double> values(classSettings.getParameter< vector<double> >("values"));

    if (ids.size() != values.size())
      throw cms::Exception("AlignmentTask::BuildFixedDetectorsConstraints") << 
        "Different number of constraint ids and values for " << tag << ".";
    
    // determine number of constraints
    unsigned int size = ids.size();
    
    // just one basic constraint
    if (oneRotZPerPot && quantityClass == qcRotZ)
    {
      if (size > 1)
        size = 1;
    }
    
    // build constraints
    for (unsigned int j = 0; j < size; j++)
    {
      // prepare empty constraint
      AlignmentConstraint ac;
      ac.forClass = quantityClass;
      ac.extended = false;

      for (auto &qcit : quantityClasses)
      {
        ac.coef[qcit].ResizeTo(QuantitiesOfClass(qcit));
        ac.coef[qcit].Zero();
      }

      // set constraint name
      char buf[40];
      sprintf(buf, "%s: fixed plane %4u", tag.c_str(), ids[j]);
      ac.name = buf;

      // get quantity index
      signed int qIndex = GetQuantityIndex(quantityClass, ids[j]);
      if (qIndex < 0)
        throw cms::Exception("AlignmentTask::BuildFixedDetectorsConstraints") <<
          "Quantity index for class " << quantityClass << " and id " << ids[j] << " is " << qIndex;
      
      // set constraint coefficient and value
      ac.coef[quantityClass][qIndex] = 1.;
      ac.val = values[j] * 1E-3;

      // save constraint
      constraints.push_back(ac);
    }
    
    if (oneRotZPerPot && quantityClass == qcRotZ)
        BuildOneRotZPerPotConstraints(constraints);
  }
}

//----------------------------------------------------------------------------------------------------

void AlignmentTask::BuildOneRotZPerPotConstraints(std::vector<AlignmentConstraint> &constraints) const
{
  printf(">> AlignmentTask::BuildOneRotZPerPotConstraints > not yet ported\n");
  throw 1;

#if 0
  // geometry is sorted by the detector number
  unsigned int prev_rp = 12345;
  for (AlignmentGeometry::iterator it = geometry.begin(); it != geometry.end(); ++it)
  {
    // do not mix different pots
    if (it->first / 10 != prev_rp)
    {
      prev_rp = it->first / 10;
      continue;
    }

    AlignmentConstraint ac;
    ac.forClass = qcRotZ;
    for (unsigned int i = 0; i < quantityClasses.size(); i++)
    {
      ac.coef[quantityClasses[i]].ResizeTo(QuantitiesOfClass(quantityClasses[i]));
      ac.coef[quantityClasses[i]].Zero();
    }
    ac.val = 0;
    ac.extended = true;
    ac.coef[qcRotZ][geometry[it->first].matrixIndex] = 1.;
    ac.coef[qcRotZ][geometry[it->first-1].matrixIndex] = -1.;

    char buf[40];
    sprintf(buf, "RotZ: %u = %u", it->first, it->first-1);
    ac.name = buf;

    constraints.push_back(ac);
  }
#endif
}

//----------------------------------------------------------------------------------------------------

void AlignmentTask::BuildStandardConstraints(vector<AlignmentConstraint> &constraints) const
{
  const vector<unsigned int> &decUnitIds = standardConstraints.getParameter<vector<unsigned int>>("units");

// TODO needed?
#if 0
  // count planes in RPs
  struct PlaneCount { unsigned int all, U, V; };
  map<unsigned int, PlaneCount> planesPerPot;
  set<unsigned int> rpsToConstrainShR, rpsToConstrainRotZ;
  for (const auto it : geometry)
  {
    unsigned int rp = it.first / 10;
    planesPerPot[rp].all++;
    if (it.second.isU)
      planesPerPot[rp].U++;
    else
      planesPerPot[rp].V++;
  }
#endif

  // ShR constraints
  if (resolveShR)
  {
    for (const auto &decUnitId : decUnitIds)
    {
      // prepare empty constraints
      AlignmentConstraint ac_X;
      ac_X.forClass = qcShR1; // TODO: it is also for ShR2...
      for (auto &qcit : quantityClasses)
      {
        ac_X.coef[qcit].ResizeTo(QuantitiesOfClass(qcit));
        ac_X.coef[qcit].Zero();
      }
      ac_X.val = 0;
      ac_X.extended = false;
    
      AlignmentConstraint ac_Y(ac_X);
  
      // set constraint names
      char buf[50];
      sprintf(buf, "ShR: unit %u, MeanX=0", decUnitId);
      ac_X.name = buf;
      sprintf(buf, "ShR: unit %u, MeanY=0", decUnitId);
      ac_Y.name = buf;

      // traverse geometry
      for (const auto git : geometry)
      {
        // stop is sensor not in the selected arm
        CTPPSDetId senId(git.first);
        unsigned int senDecUnit = senId.arm()*100 + senId.station()*10;
        if (senId.rp() > 2)
          senDecUnit += 1;

        if (senDecUnit != decUnitId)
          continue;

        // fill constraint for strip sensors
        if (senId.subdetId() == CTPPSDetId::sdTrackingStrip)
        {
          signed int qIndex = GetQuantityIndex(qcShR2, git.first);
          if (qIndex < 0)
            throw cms::Exception("AlignmentTask::BuildStandardConstraints") <<
              "Cannot get quantity index for class " << qcShR2 << " and sensor id " << git.first << ".";

          // set constraint coefficients
          ac_X.coef[qcShR2][qIndex] = git.second.GetDirectionData(2).dx;
          ac_Y.coef[qcShR2][qIndex] = git.second.GetDirectionData(2).dy;
        }

        // fill constraint for strip sensors
        if (senId.subdetId() == CTPPSDetId::sdTrackingPixel)
        {
          // get quantity indeces
          signed int qIndex1 = GetQuantityIndex(qcShR1, git.first);
          if (qIndex1 < 0)
            throw cms::Exception("AlignmentTask::BuildStandardConstraints") <<
              "Cannot get quantity index for class " << qcShR1 << " and sensor id " << git.first << ".";

          signed int qIndex2 = GetQuantityIndex(qcShR2, git.first);
          if (qIndex2 < 0)
            throw cms::Exception("AlignmentTask::BuildStandardConstraints") <<
              "Cannot get quantity index for class " << qcShR2 << " and sensor id " << git.first << ".";

          // get geometry
          const double d1x = git.second.GetDirectionData(1).dx;
          const double d1y = git.second.GetDirectionData(1).dy;
          const double d2x = git.second.GetDirectionData(2).dx;
          const double d2y = git.second.GetDirectionData(2).dy;

          // calculate coefficients
          const double D = d1x*d2y - d1y*d2x;
          const double coef_x_s1 = + d2y / D;
          const double coef_y_s1 = - d2x / D;
          const double coef_x_s2 = - d1y / D;
          const double coef_y_s2 = + d1x / D;

          // set constraint coefficients
          ac_X.coef[qcShR1][qIndex1] = coef_x_s1;
          ac_Y.coef[qcShR1][qIndex1] = coef_y_s1;

          ac_X.coef[qcShR2][qIndex2] = coef_x_s2;
          ac_Y.coef[qcShR2][qIndex2] = coef_y_s2;
        }
      }

      // add constraints
      constraints.push_back(ac_X);
      constraints.push_back(ac_Y);
    }

#if 0
    for (const auto &coUnit : units)
    {
      unsigned int indeces = QuantitiesOfClass(qcShR);
      for (unsigned int idx = 0; idx < indeces; ++idx)
      {
        const auto &entry = geometry.FindByMatrixIndex(idx);
        unsigned int plId = entry->first;
        const DetGeometry &dt = entry->second;
  
        unsigned int rpId = plId / 10;
        unsigned int rpNum = rpId % 10;
        unsigned int stNum = (rpId / 10) % 10;
        unsigned int unNum = (rpNum > 2) ? 1 : 0;     
        unsigned int unitIdx = 10*stNum + unNum;
        bool hor = (rpNum == 2 || rpNum == 3);
  
        if (unitIdx == coUnit && !hor)
        {
          const auto &pl = planesPerPot[rpId];
          ac_X.coef[qcShR][idx] = dt.d2x / pl.all;
          ac_Y.coef[qcShR][idx] = dt.d2y / pl.all;
        }
      }
  
      constraints.push_back(ac_X);
      constraints.push_back(ac_Y);
    }
#endif
  }

  // RotZ constraints
  if (resolveRotZ)
  {
    throw cms::Exception("AlignmentTask::BuildStandardConstraints") << "RotZ constriants not yet ported.";

    // TODO
#if 0
    for (const auto &coUnit : units)
    {
      // prepare empty constraints
      AlignmentConstraint ac_U;
      ac_U.forClass = qcRotZ;
      for (unsigned int i = 0; i < quantityClasses.size(); i++)
      {
        ac_U.coef[quantityClasses[i]].ResizeTo(QuantitiesOfClass(quantityClasses[i]));
        ac_U.coef[quantityClasses[i]].Zero();
      }
      ac_U.val = 0;
      ac_U.extended = true;
    
      AlignmentConstraint ac_V(ac_U);
  
      char buf[50];
      sprintf(buf, "RotZ: unit %u, SmU=0", coUnit);
      ac_U.name = buf;
      sprintf(buf, "RotZ: unit %u, SmV=0", coUnit);
      ac_V.name = buf;
  
      unsigned int indeces = QuantitiesOfClass(qcRotZ);
      for (unsigned int idx = 0; idx < indeces; ++idx)
      {
        const auto &entry = geometry.FindByMatrixIndex(idx);
        unsigned int plId = entry->first;
        const DetGeometry &dt = entry->second;
  
        unsigned int rpId = plId / 10;
        unsigned int rpNum = rpId % 10;
        unsigned int stNum = (rpId / 10) % 10;
        unsigned int unNum = (rpNum > 2) ? 1 : 0;     
        unsigned int unitIdx = 10*stNum + unNum;
  
        if (unitIdx == coUnit)
        {
          const auto &pl = planesPerPot[rpId];
          ac_U.coef[qcRotZ][idx] = (dt.isU) ? 1. / pl.U : 0.;
          ac_V.coef[qcRotZ][idx] = (!dt.isU) ? 1. / pl.V : 0.;
        }
      }
  
      constraints.push_back(ac_U);
  
      if (!useEqualMeanUMeanVRotZConstraint)
        constraints.push_back(ac_V);
    }
#endif
  }

  // mean U = mean V RotZ constraints
  if (resolveRotZ && useEqualMeanUMeanVRotZConstraint)
  {
    throw cms::Exception("AlignmentTask::BuildStandardConstraints") << "MeanUMeanVRotZ constriants not yet ported.";

    // TODO
#if 0
    map<unsigned int, pair<unsigned int, unsigned int>> planesPerPot;
    for (const auto it : geometry)
    {
      unsigned int rp = it.first / 10;
      if (it.second.isU)
        planesPerPot[rp].first++;
      else
        planesPerPot[rp].second++;
    }

    for (const auto it : planesPerPot)
    {
      const unsigned int &rp = it.first;
      const unsigned int &rp_u_planes = it.second.first;
      const unsigned int &rp_v_planes = it.second.second;

      // prepare empty constraint
      AlignmentConstraint ac;
      ac.forClass = qcRotZ;
      for (unsigned int i = 0; i < quantityClasses.size(); i++)
      {
        ac.coef[quantityClasses[i]].ResizeTo(QuantitiesOfClass(quantityClasses[i]));
        ac.coef[quantityClasses[i]].Zero();
      }
      ac.val = 0;
      ac.extended = true;
  
      char buf[40];
      sprintf(buf, "RotZ: RP%u, mU=mV", rp);
      ac.name = buf;

      // set the coefficients
      const unsigned int indeces = QuantitiesOfClass(qcRotZ);
      for (unsigned int idx = 0; idx < indeces; idx++)
      {
        double &coef = ac.coef[qcRotZ][idx];
        const auto &entry = geometry.FindByMatrixIndex(idx);
        unsigned int idx_det = entry->first;
        unsigned int idx_rp = idx_det / 10;
        bool idx_isU = entry->second.isU;

        if (rp == idx_rp)
          coef = (idx_isU) ? +1./rp_u_planes : -1./rp_v_planes;
      }
  
      // push the constraint
      constraints.push_back(ac);
    }
#endif
  }

  // TODO: remove once debugging is over
#if 0
  // print constraints
  for (const auto &c : constraints)
  {
    printf("* name=%s, value=%.3f, forClass=%i\n", c.name.c_str(), c.val, c.forClass);

    for (const auto &ccl : c.coef)
    {
      printf("    class %u (%s)\n", ccl.first, QuantityClassTag(QuantityClass(ccl.first)).c_str());

      for (int i = 0; i < ccl.second.GetNrows(); i++)
        printf("    %u -> %+.4f\n", i, ccl.second(i));
    }  
  }
#endif
 
  // TODO: when oneRotZPerPot=true, make use of BuildOneRotZPerPotConstraints(constraints) ??
}
