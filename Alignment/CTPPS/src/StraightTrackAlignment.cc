/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/


#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "CondFormats/AlignmentRecord/interface/RPRealAlignmentRecord.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/TotemRPGeometry.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/RPAlignmentCorrectionsMethods.h"
#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"

#include "Alignment/CTPPS/interface/IdealResult.h"
#include "Alignment/CTPPS/interface/JanAlignmentAlgorithm.h"
#include "Alignment/CTPPS/interface/StraightTrackAlignment.h"
#include "Alignment/CTPPS/interface/utilities.h"
#include "Alignment/CTPPS/interface/CommonMethods.h"

#include <set>
#include <unordered_set>
#include <vector>
#include <string>

#include "TDecompLU.h"
#include "TH1D.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
  
//#define DEBUG

using namespace edm;
using namespace std;

//----------------------------------------------------------------------------------------------------

TH1D* StraightTrackAlignment::NewResiduaHist(const char *name)
{
  return new TH1D(name, ";residual   (mm)", 1000, -0.2, +0.2); // in mm
}

//----------------------------------------------------------------------------------------------------

TGraph* NewGraph(const string &name, const string &title)
{
  TGraph *g = new TGraph();
  g->SetName(name.c_str());
  g->SetTitle(title.c_str());
  return g;
}

//----------------------------------------------------------------------------------------------------

StraightTrackAlignment::RPSetPlots::RPSetPlots(const string &_name) : name(_name)
{
  chisqn_lin_fitted = new TH1D("chi^2 norm, lin, fitted", ";#chi^{2}/ndf;", 5000, 0., 500.);
  chisqn_lin_selected = new TH1D("chi^2 norm, lin, selected", ";#chi^{2}/ndf;", 5000, 0., 500.);
  chisqn_log_fitted = new TH1D("chi^2 norm, log, fitted", ";log_{10}(#chi^{2}/ndf);", 700, -1., 6.);
  chisqn_log_selected = new TH1D("chi^2 norm, log, selected", ";log_{10}(#chi^{2}/ndf);", 700, -1., 6.);

  fitAxVsAyGraph_fitted = NewGraph("ax vs. ay, fitted", ";a_{x}   (rad);a_{y}   (rad)");
  fitAxVsAyGraph_selected = NewGraph("ax vs. ay, selected", ";a_{x}   (rad);a_{y}   (rad)");
  fitBxVsByGraph_fitted = NewGraph("bx vs. by, fitted", ";b_{x}   (mm);b_{y}   (mm)");
  fitBxVsByGraph_selected = NewGraph("bx vs. by, selected", ";b_{x}   (mm);b_{y}   (mm)");
}
//----------------------------------------------------------------------------------------------------

void StraightTrackAlignment::RPSetPlots::Free()
{
  delete chisqn_lin_fitted;
  delete chisqn_lin_selected;
  delete chisqn_log_fitted;
  delete chisqn_log_selected;

  delete fitAxVsAyGraph_fitted;
  delete fitAxVsAyGraph_selected;
  delete fitBxVsByGraph_fitted;
  delete fitBxVsByGraph_selected;
}

//----------------------------------------------------------------------------------------------------

void StraightTrackAlignment::RPSetPlots::Write() const
{
  chisqn_lin_fitted->Write();
  chisqn_lin_selected->Write();
  chisqn_log_fitted->Write();
  chisqn_log_selected->Write();

  fitAxVsAyGraph_fitted->Write();
  fitAxVsAyGraph_selected->Write();
  fitBxVsByGraph_fitted->Write();
  fitBxVsByGraph_selected->Write();
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

StraightTrackAlignment::StraightTrackAlignment(const ParameterSet& ps) :
  verbosity(ps.getUntrackedParameter<unsigned int>("verbosity", 0)),
  factorizationVerbosity(ps.getUntrackedParameter<unsigned int>("factorizationVerbosity", 0)),

  rpIds(ps.getParameter< vector<unsigned int> >("rpIds")),
  excludePlanes(ps.getParameter< vector<unsigned int> >("excludePlanes")),
  z0(ps.getParameter<double>("z0")),

  maxEvents(ps.getParameter<signed int>("maxEvents")),

  removeImpossible(ps.getParameter<bool>("removeImpossible")),
  requireNumberOfUnits(ps.getParameter<unsigned int>("requireNumberOfUnits")),
  requireAtLeast3PotsInOverlap(ps.getParameter<bool>("requireAtLeast3PotsInOverlap")),
  requireOverlap(ps.getParameter<bool>("requireOverlap")),
  cutOnChiSqPerNdf(ps.getParameter<bool>("cutOnChiSqPerNdf")),
  chiSqPerNdfCut(ps.getParameter<double>("chiSqPerNdfCut")),
  maxTrackAx(ps.getParameter<double>("maxTrackAx")),
  maxTrackAy(ps.getParameter<double>("maxTrackAy")),

  fileNamePrefix(ps.getParameter<string>("fileNamePrefix")),
  expandedFileNamePrefix(ps.getParameter<string>("expandedFileNamePrefix")),
  factoredFileNamePrefix(ps.getParameter<string>("factoredFileNamePrefix")),
  preciseXMLFormat(ps.getParameter<bool>("preciseXMLFormat")),

  saveIntermediateResults(ps.getParameter<bool>("saveIntermediateResults")),
  taskDataFileName(ps.getParameter<string>("taskDataFileName")),
  taskDataFile(NULL),

  task(ps),
  fitter(ps),

  buildDiagnosticPlots(ps.getParameter<bool>("buildDiagnosticPlots")),
  diagnosticsFile(ps.getParameter<string>("diagnosticsFile")),
  fitNdfHist_fitted(new TH1D("ndf_fitted", ";ndf;", 41, -4.5, 36.5)),
  fitNdfHist_selected(new TH1D("ndf_selected", ";ndf;", 41, -4.5, 36.5)),
  fitPHist_fitted(new TH1D("p_fitted", ";p value;", 100, 0., 1.)),
  fitPHist_selected(new TH1D("p_selected", ";p value;", 100, 0., 1.)),
  fitAxHist_fitted(new TH1D("ax_fitted", ";a_{x}   (rad);", 10000, -0.1, 0.1)),
  fitAxHist_selected(new TH1D("ax_selected", ";a_{x}   (rad);", 10000, -0.1, 0.1)),
  fitAyHist_fitted(new TH1D("ay_fitted", ";a_{y}   (rad);", 10000, -0.1, 0.1)),
  fitAyHist_selected(new TH1D("ay_selected", ";a_{y}   (rad);", 10000, -0.1, 0.1)),
  fitBxHist_fitted(new TH1D("bx_fitted", ";b_{x}   (mm);", 500, -30., 30.)),
  fitBxHist_selected(new TH1D("bx_selected", ";b_{x}   (mm);", 500, -30., 30.)),
  fitByHist_fitted(new TH1D("by_fitted", ";b_{y}   (mm);", 500, -30., 30.)),
  fitByHist_selected(new TH1D("by_selected", ";b_{y}   (mm);", 500, -30., 30.)),

  globalPlots("global")
{
  // open task data file
  if (!taskDataFileName.empty())
    taskDataFile = new TFile(taskDataFileName.c_str(), "recreate");

  // instantiate algorithm objects
  // (and save them)
  vector<string> alNames(ps.getParameter< vector<string> >("algorithms"));
  for (unsigned int i = 0; i < alNames.size(); i++)
  {
    AlignmentAlgorithm *a = NULL;

    if (alNames[i].compare("Ideal") == 0)
      a = new IdealResult(ps, &task);

    if (alNames[i].compare("Jan") == 0)
    {
      JanAlignmentAlgorithm *jaa = new JanAlignmentAlgorithm(ps, &task);
      a = jaa;
    }

    if (a)
      algorithms.push_back(a);
    else
      throw cms::Exception("StraightTrackAlignment") << "Unknown alignment algorithm `" << alNames[i] << "'.";
  }

  // get constraints type
  string ct = ps.getParameter<string>("constraintsType");
  if (ct.compare("homogeneous") == 0) constraintsType = ctHomogeneous;
  else
    if (ct.compare("fixedDetectors") == 0) constraintsType = ctFixedDetectors;
    else
      if (ct.compare("standard") == 0) constraintsType = ctStandard;
      else
          throw cms::Exception("StraightTrackAlignment") << "Unknown constraints type `" << ct << "'.";

  // parse additional accepted RP sets
  string aars_str = ps.getParameter<string>("additionalAcceptedRPSets");

  size_t idx_b = 0, idx_e = string::npos;
  while (idx_b != string::npos)
  {
    // get one block - portion between successive ";"
    idx_e = aars_str.find(';', idx_b);
    size_t len = (idx_e == string::npos) ? string::npos : idx_e - idx_b;
    string block = aars_str.substr(idx_b, len);

    // process the block
    if (!block.empty())
    {
      set<unsigned int> rpSet;

      // isolate bits (= RP ids)
      size_t bi_b = 0, bi_e = string::npos;
      while (bi_b != string::npos)
      {
        bi_e = block.find(',', bi_b);
        size_t bit_len = (bi_e == string::npos) ? string::npos : bi_e - bi_b;
        const string &bit = aars_str.substr(bi_b, bit_len);

        unsigned int rp = atoi(bit.c_str());
        rpSet.insert(rp);       

        bi_b = (bi_e == string::npos) ? string::npos : bi_e + 1;
      }

      additionalAcceptedRPSets.push_back(rpSet);
    }

    // move to next block
    idx_b = (idx_e == string::npos) ? string::npos : idx_e + 1;
  }
}

//----------------------------------------------------------------------------------------------------

StraightTrackAlignment::~StraightTrackAlignment()
{
  if (taskDataFile);
    delete taskDataFile;

  for (vector<AlignmentAlgorithm *>::iterator it = algorithms.begin(); it != algorithms.end(); ++it)
    delete (*it);

  delete fitNdfHist_fitted;
  delete fitNdfHist_selected;
  delete fitPHist_fitted;
  delete fitPHist_selected;
  delete fitAxHist_fitted;
  delete fitAyHist_fitted;
  delete fitAxHist_selected;
  delete fitAyHist_selected;
  delete fitBxHist_fitted;
  delete fitByHist_fitted;
  delete fitBxHist_selected;
  delete fitByHist_selected;

  globalPlots.Free();

  for (auto &p : rpSetPlots)
    p.second.Free();

  for (map<unsigned int, ResiduaHistogramSet>::iterator it = residuaHistograms.begin(); it != residuaHistograms.end(); ++it)
  {
    delete it->second.total_fitted;
    delete it->second.total_selected;
    delete it->second.selected_vs_chiSq;
    for (map< set<unsigned int>, TH1D* >::iterator sit = it->second.perRPSet_fitted.begin();
          sit != it->second.perRPSet_fitted.end(); ++sit)
      delete sit->second;
    for (map< set<unsigned int>, TH1D* >::iterator sit = it->second.perRPSet_selected.begin();
          sit != it->second.perRPSet_selected.end(); ++sit)
      delete sit->second;
  }
}

//----------------------------------------------------------------------------------------------------

void StraightTrackAlignment::Begin(const EventSetup &es)
{
  printf(">> StraightTrackAlignment::Begin\n");

  // reset counters
  eventsTotal = 0;
  eventsFitted = 0;
  eventsSelected = 0;
  fittedTracksPerRPSet.clear();
  selectedTracksPerRPSet.clear();
  
  // prepare geometry (in fact, this should be done whenever es gets changed)
  ESHandle<TotemRPGeometry> hGeometry;
  es.get<VeryForwardRealGeometryRecord>().get(hGeometry);
  task.BuildGeometry(rpIds, excludePlanes, hGeometry.product(), z0, task.geometry);  

  // build matrix index mappings
  task.BuildIndexMaps();

  // print geometry info
  if (verbosity > 1)
    task.geometry.Print();
  
  // save task (including geometry) and fitter
  if (taskDataFile)
  {
    taskDataFile->WriteObject(&task, "task");
    taskDataFile->WriteObject(&fitter, "fitter");
  }

  // initiate the algorithms
  for (vector<AlignmentAlgorithm *>::iterator it = algorithms.begin(); it != algorithms.end(); ++it)
    (*it)->Begin(es);
  
  // get initial alignments
  try
  {
    ESHandle<RPAlignmentCorrectionsData> h;
    es.get<RPRealAlignmentRecord>().get(h);
    initialAlignments = *h;
  }
  catch (...) {}
}

//----------------------------------------------------------------------------------------------------

void StraightTrackAlignment::ProcessEvent(const DetSetVector<TotemRPUVPattern> &uvPatternsStrip, const DetSetVector<CTPPSDiamondRecHit> &hitsDiamond, const DetSetVector<CTPPSPixelRecHit> &hitsPixel)
{
  eventsTotal++;

  if (verbosity > 9)
    printf("\n---------- StraightTrackAlignment::ProcessEvent\n");

  // -------------------- STEP 1: get hits from selected RPs
  
  HitCollection hitSelection;

// strips

  for (const auto &pv : uvPatternsStrip)
  {
    const CTPPSDetId detId(pv.detId());
    const unsigned int rpDecId = detId.arm()*100 + detId.station()*10 + detId.rp();

    if (find(rpIds.begin(), rpIds.end(), rpDecId) == rpIds.end())
      continue;

    // is U-V association unique?
    unsigned int n_U=0, n_V=0;
    unsigned int idx_U=0, idx_V=0;
    for (unsigned int pi = 0; pi < pv.size(); pi++)
    {
      const TotemRPUVPattern &pattern = pv[pi];

       switch (pattern.getProjection())
      {
        case TotemRPUVPattern::projU:
          n_U++;
          idx_U=pi;
          break;
     
        case TotemRPUVPattern::projV:
          n_V++;
          idx_V=pi;
          break;

        default:
          break;
       }  
    }

    if (n_U != 1 || n_V != 1)
    {
	  std::cout<< ">> StraightTrackAlignment::ProcessEvent > Impossible to combine U and V patterns in RP " << detId << " (n_U=" << n_U << ", n_V=" << n_V << ")." << std::endl;
         continue;
    }
    
    if (!pv[idx_U].getFittable() || !pv[idx_V].getFittable())
	 continue;

    // combine U and V hits
    DetSetVector<TotemRPRecHit> hits;
    
    // U pattern
    for (auto &ids : pv[idx_U].getHits())
    {	
	const double &z = task.geometry[ids.detId()].z;

		for (auto &h : ids)
		{
	  		hitSelection.emplace_back(Hit(pv.detId(), 2, h.getPosition(), h.getSigma(), z));
		}  
    }

	// V pattern
    for (auto &ids : pv[idx_V].getHits())
    {	
	const double &z = task.geometry[ids.detId()].z;

		for (auto &h : ids)
		{
	  		hitSelection.emplace_back(Hit(pv.detId(), 2, h.getPosition(), h.getSigma(), z));
		}  
    }
  }

  // diamonds
  for (const auto &pv : hitsDiamond)
  {
    // skip if RP not selected
    const CTPPSDetId detId(pv.detId());
    const unsigned int rpDecId = detId.arm()*100 + detId.station()*10 + detId.rp();

    if (find(rpIds.begin(), rpIds.end(), rpDecId) == rpIds.end())
      continue;

    const double &z = task.geometry[pv.detId()].z;

    for (const auto &h : pv)
      hitSelection.emplace_back(Hit(pv.detId(), 1, h.getX(), h.getXWidth() / sqrt(12.), z));
  }

  // pixels
  for (const auto &pv : hitsPixel)
  {
    // skip if RP not selected
    const CTPPSDetId detId(pv.detId());
    const unsigned int rpDecId = detId.arm()*100 + detId.station()*10 + detId.rp();

    if (find(rpIds.begin(), rpIds.end(), rpDecId) == rpIds.end())
      continue;

    for (const auto &h : pv)
    {
      const auto &dg = task.geometry[pv.detId()];
      const double dz1 = dg.GetDirectionData(1).dz;
      const double dz2 = dg.GetDirectionData(2).dz;
      const double z = dg.z + h.getPoint().x() * dz1 + h.getPoint().y() * dz2;

      hitSelection.emplace_back(Hit(pv.detId(), 1, h.getPoint().x(), sqrt(h.getError().xx()), z));
      hitSelection.emplace_back(Hit(pv.detId(), 2, h.getPoint().y(), sqrt(h.getError().yy()), z));
    }
  }

  if (hitSelection.empty())
    return;

  // -------------------- STEP 2: fit + outlier rejection

  LocalTrackFit trackFit;
  if (! fitter.Fit(hitSelection, task.geometry, trackFit))
     return;

  set<unsigned int> selectedRPs;
  for (const auto &hit : hitSelection)
  {
    CTPPSDetId detId(hit.id);
    const unsigned int decRPId = detId.arm()*100 + detId.station()*10 + detId.rp();
    selectedRPs.insert(decRPId);
  }

  eventsFitted++;
  fittedTracksPerRPSet[selectedRPs]++;

  // -------------------- STEP 3: quality checks

  bool top = false, bottom = false, horizontal = false;
  unordered_set<unsigned int> units;
  for (const auto &rp : selectedRPs)
  {
    unsigned int rpIdx = rp % 10;
    unsigned int stId = rp / 10;
    unsigned int unitId = stId * 10;
    if (rpIdx > 2)
      unitId++;

    if (rpIdx == 0 || rpIdx == 4) top = true;
    if (rpIdx == 1 || rpIdx == 5) bottom = true;
    if (rpIdx == 2 || rpIdx == 3) horizontal = true;

    units.insert(unitId);
  }

  bool overlap = (top && horizontal) || (bottom && horizontal);

  bool rp_set_accepted = true;

  // impossible signature
  if (removeImpossible && top && bottom)
    rp_set_accepted = false;

  // cleanliness cuts
  if (units.size() < requireNumberOfUnits)
    rp_set_accepted = false;
  
  if (requireOverlap && !overlap)
    rp_set_accepted = false;

  if (requireAtLeast3PotsInOverlap && overlap && selectedRPs.size() < 3)
    rp_set_accepted = false;

  // is it an additional accepted RP set?
  if (find(additionalAcceptedRPSets.begin(), additionalAcceptedRPSets.end(), selectedRPs) != additionalAcceptedRPSets.end())
    rp_set_accepted = true;

  bool selected = rp_set_accepted;


  // too bad chisq
  if (cutOnChiSqPerNdf && trackFit.ChiSqPerNdf() > chiSqPerNdfCut)
    selected = false;

  // parallelity cut
  if (fabs(trackFit.ax) > maxTrackAx || fabs(trackFit.ay) > maxTrackAy)
    selected = false;

  UpdateDiagnosticHistograms(hitSelection, selectedRPs, trackFit, selected);

  if (verbosity > 5)
    printf("* SELECTED: %u\n", selected);

  if (!selected)
    return;
  
  eventsSelected++;
  selectedTracksPerRPSet[selectedRPs]++;
  
  // -------------------- STEP 4: FEED ALGORITHMS

  for (auto &a : algorithms)
    a->Feed(hitSelection, trackFit);

  // -------------------- STEP 5: ENOUGH TRACKS?

  if (eventsSelected == maxEvents)
      throw "StraightTrackAlignment: Number of tracks processed reached maximum";
}

//----------------------------------------------------------------------------------------------------

void StraightTrackAlignment::UpdateDiagnosticHistograms(const HitCollection &selection, 
      const set<unsigned int> &selectedRPs, const LocalTrackFit &trackFit, bool trackSelected)
{
  if (!buildDiagnosticPlots)
    return;

  fitNdfHist_fitted->Fill(trackFit.ndf);
  fitPHist_fitted->Fill(trackFit.PValue());
  fitAxHist_fitted->Fill(trackFit.ax);
  fitAyHist_fitted->Fill(trackFit.ay);
  fitBxHist_fitted->Fill(trackFit.bx);
  fitByHist_fitted->Fill(trackFit.by);
  
  globalPlots.chisqn_lin_fitted->Fill(trackFit.ChiSqPerNdf());
  globalPlots.chisqn_log_fitted->Fill(log10(trackFit.ChiSqPerNdf()));
  globalPlots.fitAxVsAyGraph_fitted->SetPoint(globalPlots.fitAxVsAyGraph_fitted->GetN(), trackFit.ax, trackFit.ay);
  globalPlots.fitBxVsByGraph_fitted->SetPoint(globalPlots.fitBxVsByGraph_fitted->GetN(), trackFit.bx, trackFit.by);

  if (trackSelected)
  {
    fitNdfHist_selected->Fill(trackFit.ndf);
    fitPHist_selected->Fill(trackFit.PValue());
    fitAxHist_selected->Fill(trackFit.ax);
    fitAyHist_selected->Fill(trackFit.ay);
    fitBxHist_selected->Fill(trackFit.bx);
    fitByHist_selected->Fill(trackFit.by);

    globalPlots.chisqn_lin_selected->Fill(trackFit.ChiSqPerNdf());
    globalPlots.chisqn_log_selected->Fill(log10(trackFit.ChiSqPerNdf()));
    globalPlots.fitAxVsAyGraph_selected->SetPoint(globalPlots.fitAxVsAyGraph_selected->GetN(), trackFit.ax, trackFit.ay);
    globalPlots.fitBxVsByGraph_selected->SetPoint(globalPlots.fitBxVsByGraph_selected->GetN(), trackFit.bx, trackFit.by);
  }

  map< set<unsigned int>, RPSetPlots >::iterator it = rpSetPlots.find(selectedRPs);
  if (it == rpSetPlots.end())
    it = rpSetPlots.insert( { selectedRPs, RPSetPlots(SetToString(selectedRPs)) } ).first;

  it->second.chisqn_lin_fitted->Fill(trackFit.ChiSqPerNdf());
  it->second.chisqn_log_fitted->Fill(log10(trackFit.ChiSqPerNdf()));
  it->second.fitAxVsAyGraph_fitted->SetPoint(it->second.fitAxVsAyGraph_fitted->GetN(), trackFit.ax, trackFit.ay);
  it->second.fitBxVsByGraph_fitted->SetPoint(it->second.fitBxVsByGraph_fitted->GetN(), trackFit.bx, trackFit.by);

  if (trackSelected)
  { 
    it->second.chisqn_lin_selected->Fill(trackFit.ChiSqPerNdf());
    it->second.chisqn_log_selected->Fill(log10(trackFit.ChiSqPerNdf()));
    it->second.fitAxVsAyGraph_selected->SetPoint(it->second.fitAxVsAyGraph_selected->GetN(), trackFit.ax, trackFit.ay);
    it->second.fitBxVsByGraph_selected->SetPoint(it->second.fitBxVsByGraph_selected->GetN(), trackFit.bx, trackFit.by);
  }
  
  for (HitCollection::const_iterator hitCollectionIterator = selection.begin(); hitCollectionIterator != selection.end(); ++hitCollectionIterator)
  {
    unsigned int id = hitCollectionIterator->id;

    AlignmentGeometry::iterator dit = task.geometry.find(id);
    if (dit == task.geometry.end())
      continue;

    DetGeometry &d = dit->second;
    const auto dirData = d.GetDirectionData(hitCollectionIterator->dirIdx);

    double m = hitCollectionIterator->position + dirData.s;
    double x = trackFit.ax * hitCollectionIterator->z + trackFit.bx;
    double y = trackFit.ay * hitCollectionIterator->z + trackFit.by;
    double f = x*dirData.dx + y*dirData.dy;
    double R = m - f;

    map<unsigned int, ResiduaHistogramSet>::iterator it = residuaHistograms.find(id);
    if (it == residuaHistograms.end())
    {
      it = residuaHistograms.insert(pair<unsigned int, ResiduaHistogramSet>(id, ResiduaHistogramSet())).first;
      char buf[30];
      sprintf(buf, "%u: total_fitted", id); it->second.total_fitted = NewResiduaHist(buf);
      sprintf(buf, "%u: total_selected", id); it->second.total_selected = NewResiduaHist(buf);
      it->second.selected_vs_chiSq = new TGraph();
      sprintf(buf, "%u: selected_vs_chiSq", id);
      it->second.selected_vs_chiSq->SetName(buf);
    }

    it->second.total_fitted->Fill(R);
    if (trackSelected)
    {
      it->second.total_selected->Fill(R);
      it->second.selected_vs_chiSq->SetPoint(it->second.selected_vs_chiSq->GetN(), trackFit.ChiSqPerNdf(), R);
    }

    map< set<unsigned int>, TH1D* >::iterator sit = it->second.perRPSet_fitted.find(selectedRPs);
    if (sit == it->second.perRPSet_fitted.end())
    {
      char buf[10];
      sprintf(buf, "%u: ", id);
      string label = buf;
      label += SetToString(selectedRPs);
      sit = it->second.perRPSet_fitted.insert(pair< set<unsigned int>, TH1D* >(selectedRPs, NewResiduaHist(label.c_str()))).first;
    }
    sit->second->Fill(R);
    
    if (trackSelected)
    {
      sit = it->second.perRPSet_selected.find(selectedRPs);
      if (sit == it->second.perRPSet_selected.end())
      {
        char buf[10];
        sprintf(buf, "%u: ", id);
        string label = buf;
        label += SetToString(selectedRPs);
        sit = it->second.perRPSet_selected.insert(pair< set<unsigned int>, TH1D* >(selectedRPs, NewResiduaHist(label.c_str()))).first;
      }
      sit->second->Fill(R);
    }
  }
}

//----------------------------------------------------------------------------------------------------

void StraightTrackAlignment::BuildConstraints(vector<AlignmentConstraint> &constraints)
{
  constraints.clear();

  switch (constraintsType)
  {
    case ctHomogeneous:
      task.BuildHomogeneousConstraints(constraints);
      return;

    case ctFixedDetectors:
      task.BuildFixedDetectorsConstraints(constraints);
      return;

    case ctStandard:
      task.BuildStandardConstraints(constraints);
      return;
  }
}

//----------------------------------------------------------------------------------------------------

void StraightTrackAlignment::Finish()
{
  // print statistics
  if (verbosity)
  {
    printf("----------------------------------------------------------------------------------------------------\n");
    printf("\n>> StraightTrackAlignment::Finish\n");
    printf("\tevents total = %i\n", eventsTotal);
    printf("\tevents fitted = %i\n", eventsFitted);
    printf("\tevents selected = %i\n", eventsSelected);
    printf("\n%30s  %10s%10s\n", "set of RPs", "fitted", "selected");
    for (map< set<unsigned int>, unsigned long >::iterator it = fittedTracksPerRPSet.begin();
        it != fittedTracksPerRPSet.end(); ++it)
    {
      const string &label = SetToString(it->first);

      map< set<unsigned int>, unsigned long >::iterator sit = selectedTracksPerRPSet.find(it->first);
      unsigned long sv = (sit == selectedTracksPerRPSet.end()) ? 0 : sit->second;

      printf("%30s :%10lu%10lu\n", label.c_str(), it->second, sv);
    }
  }

  // write diagnostics plots
  SaveDiagnostics();

  // run analysis
  for (auto a : algorithms)
    a->Analyze();

  // build constraints
  vector<AlignmentConstraint> constraints;
  BuildConstraints(constraints);

  // save constraints
  if (taskDataFile)
    taskDataFile->WriteObject(&constraints, "constraints");  
  
  printf("\n>> StraightTrackAlignment::Finish > %lu constraints built\n", constraints.size());
  for (unsigned int i = 0; i < constraints.size(); i++)
  {
    printf("\t%25s, qc = %i, extended = %i\n", constraints[i].name.c_str(), constraints[i].forClass, constraints[i].extended);
  }

  // solve
  vector<RPAlignmentCorrectionsData> results;
  for (auto algorithm : algorithms)
  {
    TDirectory *dir = NULL;
    if (taskDataFile && saveIntermediateResults)
      dir = taskDataFile->mkdir((algorithm->GetName() + "_data").c_str());

    results.resize(results.size() + 1);
    unsigned int rf = algorithm->Solve(constraints, results.back(), dir);

    if (rf)
      throw cms::Exception("StraightTrackAlignment") << "The Solve method of `" << algorithm->GetName() 
        << "' algorithm has failed (return value " << rf << ").";
  }

  // print results
  printf("\n>> StraightTrackAlignment::Finish > Print\n");

  PrintLineSeparator(results);
  PrintQuantitiesLine(results);
  PrintAlgorithmsLine(results);

  signed int prevRPId = -1;

  for (const auto &dit : task.geometry)
  {
    //  ═ ║ 

    signed int rpId = CTPPSDetId(dit.first).getRPId();
    if (rpId != prevRPId)
      PrintLineSeparator(results);
    prevRPId = rpId;

    PrintId(dit.first);
    printf(" ║");

    for (unsigned int q = 0; q < task.quantityClasses.size(); q++)
    {
      for (unsigned int a = 0; a < results.size(); a++)
      {
        RPAlignmentCorrectionsData::mapType::const_iterator it = results[a].sensors.find(dit.first);
        if (it == results[a].sensors.end())
        {
          if (algorithms[a]->HasErrorEstimate())
            printf("%21s", "----║");
          else
            printf("%10s", "----║");
          continue;
        }

        const RPAlignmentCorrectionData &ac = it->second;
        double v = 0., e = 0.;
        switch (task.quantityClasses[q])
        {
          case AlignmentTask::qcShR1: v = ac.sh_r1(); e = ac.sh_r1_e(); break;
          case AlignmentTask::qcShR2: v = ac.sh_r2(); e = ac.sh_r2_e(); break;
          case AlignmentTask::qcShZ: v = ac.sh_z(); e = ac.sh_z_e(); break;
          case AlignmentTask::qcRotZ: v = ac.rot_z(); e = ac.rot_z_e(); break;
        }

        if (algorithms[a]->HasErrorEstimate())
          printf("%+8.1f ± %7.1f", v*1E3, e*1E3);
        else
          printf("%+8.1f", v*1E3);

        if (a + 1 == results.size())
          printf("║");
        else
          printf("│");
      }
    }

    printf("\n");
  }
  
  PrintLineSeparator(results);
  PrintAlgorithmsLine(results);
  PrintQuantitiesLine(results);
  PrintLineSeparator(results);

  // save results
  for (unsigned int a = 0; a < results.size(); a++)
  {
    // convert readout corrections to X and Y
    for (RPAlignmentCorrectionsData::mapType::iterator it = results[a].sensors.begin(); it != results[a].sensors.end(); ++it)
    {
      DetGeometry &d = task.geometry[it->first];
      const auto dirData1 = d.GetDirectionData(1);
      const auto dirData2 = d.GetDirectionData(2);
      it->second.readoutTranslationsToXY(dirData1.dx, dirData1.dy, dirData2.dx, dirData2.dy);
    }

    // write non-cumulative results
    if (!fileNamePrefix.empty())
    {
      RPAlignmentCorrectionsMethods::WriteXMLFile(results[a], fileNamePrefix + algorithms[a]->GetName() + ".xml",
        preciseXMLFormat, algorithms[a]->HasErrorEstimate(), true, false, false, true);
    }

    // convert shr to xy, normalize z rotations
    RPAlignmentCorrectionsData convertedAlignments = results[a];
    for (auto it = convertedAlignments.sensors.begin(); it != convertedAlignments.sensors.end(); ++it)
    {
      const DetGeometry &g = task.geometry[it->first];
      auto d1 = g.GetDirectionData(1);
      auto d2 = g.GetDirectionData(2);

      it->second.readoutTranslationsToXY(d1.dx, d1.dy, d2.dx, d2.dy);

      it->second.normalizeRotationZ();
    }

    // merge alignments
    RPAlignmentCorrectionsData cumulativeAlignments;
    cumulativeAlignments.AddCorrections(initialAlignments, false);
    cumulativeAlignments.AddCorrections(convertedAlignments, false);

    // write expanded and factored results
    if (!expandedFileNamePrefix.empty() || !factoredFileNamePrefix.empty())
    {
      RPAlignmentCorrectionsData expandedAlignments;
      RPAlignmentCorrectionsData factoredAlignments;

      if (factorizationVerbosity)
        printf(">> Factorizing results of %s algorithm\n", algorithms[a]->GetName().c_str());
      
      const bool equalWeights = false;
      CommonMethods::FactorRPFromSensorCorrections(cumulativeAlignments, expandedAlignments, factoredAlignments,
        task.geometry, equalWeights, factorizationVerbosity);

      if (!expandedFileNamePrefix.empty())
      {
        RPAlignmentCorrectionsMethods::WriteXMLFile(expandedAlignments, expandedFileNamePrefix + algorithms[a]->GetName() + ".xml",
          preciseXMLFormat, false, false, true, false, true);
      }

      if (!factoredFileNamePrefix.empty())
      {
        RPAlignmentCorrectionsMethods::WriteXMLFile(factoredAlignments, factoredFileNamePrefix + algorithms[a]->GetName() + ".xml",
          preciseXMLFormat, false, false, true, false, true);
      }
    }
  }
  
  // prepare algorithms for destructions
  for (const auto &algorithm : algorithms)
    algorithm->End();
}

//----------------------------------------------------------------------------------------------------

string StraightTrackAlignment::SetToString(const set<unsigned int> &s)
{
  unsigned int N = s.size();
  if (N == 0)
    return "empty";

  string str;
  char buf[10];
  unsigned int i = 0;
  for (set<unsigned int>::iterator it = s.begin(); it != s.end(); ++it, ++i)
  {
    sprintf(buf, "%u", *it);
    str += buf;
    if (i < N - 1)
      str += ", ";
  }

  return str;
}

//----------------------------------------------------------------------------------------------------

void StraightTrackAlignment::PrintN(const char *str, unsigned int N)
{
  for (unsigned int i = 0; i < N; i++)
    printf("%s", str);
}

//----------------------------------------------------------------------------------------------------

void StraightTrackAlignment::PrintLineSeparator(const std::vector<RPAlignmentCorrectionsData> &results)
{
  printf("═════════════════════════╬");
  for (unsigned int q = 0; q < task.quantityClasses.size(); q++)
  {
    for (unsigned int a = 0; a < results.size(); a++)
    {
      PrintN("═", algorithms[a]->HasErrorEstimate() ? 18 : 8);
      if (a + 1 != results.size())
        printf("═");
    }
    printf("╬");
  }
  printf("\n");
}

//----------------------------------------------------------------------------------------------------

void StraightTrackAlignment::PrintQuantitiesLine(const std::vector<RPAlignmentCorrectionsData> &results)
{
  printf("                         ║");

  for (unsigned int q = 0; q < task.quantityClasses.size(); q++)
  {
    unsigned int size = 0;
    for (unsigned int a = 0; a < results.size(); a++)
      size += (algorithms[a]->HasErrorEstimate()) ? 18 : 8;
    size += algorithms.size() - 1; 

    const string &tag = task.QuantityClassTag(task.quantityClasses[q]);
    unsigned int space = (size - tag.size())/2;
    PrintN(" ", space);
    printf("%s", tag.c_str());
    PrintN(" ", size - space - tag.size());
    printf("║");
  }
  printf("\n");
}

//----------------------------------------------------------------------------------------------------

void StraightTrackAlignment::PrintAlgorithmsLine(const std::vector<RPAlignmentCorrectionsData> &results)
{
  printf("                         ║");

  for (unsigned int q = 0; q < task.quantityClasses.size(); q++)
  {
    for (unsigned int a = 0; a < results.size(); a++)
    {
      printf((algorithms[a]->HasErrorEstimate()) ? "%18s" : "%8s", algorithms[a]->GetName().substr(0, 8).c_str());

      if (a + 1 == results.size())
        printf("║");
      else
        printf("│");
    }
  }
  printf("\n");
}

//----------------------------------------------------------------------------------------------------

void StraightTrackAlignment::SaveDiagnostics() const
{
  if (diagnosticsFile.empty())
    return;

  TFile *df = new TFile(diagnosticsFile.c_str(), "recreate");
  if (df->IsZombie())
    throw cms::Exception("StraightTrackAlignment::SaveDiagnostics") << "Cannot open file `" << 
      diagnosticsFile << "' for writing.";

  if (buildDiagnosticPlots)
  {
    TDirectory *commonDir = df->mkdir("common");
    gDirectory = commonDir;

    fitNdfHist_fitted->Write();
    fitNdfHist_selected->Write();
    fitAxHist_fitted->Write();
    fitAyHist_fitted->Write();
    fitAxHist_selected->Write();
    fitAyHist_selected->Write();
    fitBxHist_fitted->Write();
    fitByHist_fitted->Write();
    fitBxHist_selected->Write();
    fitByHist_selected->Write();
    fitPHist_fitted->Write();
    fitPHist_selected->Write();

    gDirectory = commonDir->mkdir("plots global");
    globalPlots.Write();

    TDirectory *ppsDir = commonDir->mkdir("plots per RP set");
    for (map< set<unsigned int>, RPSetPlots >::const_iterator it = rpSetPlots.begin(); it != rpSetPlots.end(); ++it)
    {
      gDirectory = ppsDir->mkdir(SetToString(it->first).c_str());

      it->second.Write();
    } 

    TDirectory *resDir = commonDir->mkdir("residuals");
    for (map<unsigned int, ResiduaHistogramSet>::const_iterator it = residuaHistograms.begin(); it != residuaHistograms.end(); ++it)
    {
      char buf[10];
      sprintf(buf, "%u", it->first);
      gDirectory = resDir->mkdir(buf);
      it->second.total_fitted->Write();
      it->second.total_selected->Write();
      it->second.selected_vs_chiSq->Write();

/*
      gDirectory = gDirectory->mkdir("fitted per RP set");
      for (map< set<unsigned int>, TH1D* >::iterator sit = it->second.perRPSet_fitted.begin();
          sit != it->second.perRPSet_fitted.end(); ++sit)
        sit->second->Write();
      gDirectory->cd("..");
*/

      gDirectory = gDirectory->mkdir("selected per RP set");
      TCanvas *c = new TCanvas; c->SetName("alltogether");
      unsigned int idx = 0;
      for (map< set<unsigned int>, TH1D* >::const_iterator sit = it->second.perRPSet_selected.begin();
          sit != it->second.perRPSet_selected.end(); ++sit, ++idx) {
        sit->second->SetLineColor(idx+1);
        sit->second->Draw((idx == 0) ? "" : "same");
        sit->second->Write();
      }
      c->Write();
    }
  }

  // save diagnostics of algorithms
  for (vector<AlignmentAlgorithm *>::const_iterator it = algorithms.begin(); it != algorithms.end(); ++it)
  {
    TDirectory *algDir = df->mkdir((*it)->GetName().c_str());
    (*it)->SaveDiagnostics(algDir);
  }
  
  delete df;
}

