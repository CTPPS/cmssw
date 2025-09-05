#include <Rtypes.h>
#include <RtypesCore.h>
#include <TError.h>
#include <TMath.h>
#include <TString.h>
#include "TVector2.h"
#include "TROOT.h"
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <set>
#include <utility>
#include <nlohmann/json.hpp>

#include <iostream>

// Points manipulation
double Cross2D(const TVector2 &u, const TVector2 &v) { return u.X() * v.Y() - u.Y() * v.X(); }

TVector2 findCentroid(const std::vector<TVector2> &pts) {
  double cx = 0, cy = 0;
  for (auto &p : pts) {
    cx += p.X();
    cy += p.Y();
  }
  cx /= pts.size();
  cy /= pts.size();
  return TVector2(cx, cy);
}

std::vector<TVector2> orderCCW(std::vector<TVector2> pts) {
  TVector2 centroid = findCentroid(pts);

  std::sort(pts.begin(), pts.end(), [&](const TVector2 &a, const TVector2 &b) {
    return atan2(a.Y() - centroid.Y(), a.X() - centroid.X()) < atan2(b.Y() - centroid.Y(), b.X() - centroid.X());
  });
  return pts;
}

bool isConvex(const std::vector<TVector2> &poly) {  // points have to be ccw or cw
  int n = poly.size();
  if (n < 3)
    return false;

  double lastCross = 0;

  for (int i = 0; i < n; ++i) {
    TVector2 A = poly[i];
    TVector2 B = poly[(i + 1) % n];
    TVector2 C = poly[(i + 2) % n];

    TVector2 AB = B - A;
    TVector2 BC = C - B;

    double cross = Cross2D(AB, BC);

    if (i == 0)
      lastCross = cross;
    else {
      if (cross * lastCross < 0)
        return false;
    }

    if (std::abs(cross) > 1e-12)
      lastCross = cross;
  }

  return true;
}

bool pointInConvexPolygon(const std::vector<TVector2> &poly, const TVector2 &P) {
  bool hasPos = false, hasNeg = false;
  for (int i = 0; i < 4; i++) {
    double cross = Cross2D(poly[(i + 1) % 4] - poly[i], P - poly[i]);
    if (cross > 1e-9)
      hasPos = true;
    if (cross < -1e-9)
      hasNeg = true;
    if (hasPos && hasNeg)
      return false;  // mixed signs -> outside
  }
  return true;  // all same sign (or on edge)
}

// Calculate efficiency
std::vector<double> getEfficiency(
    TH2F *numerator,
    TH2F *denominator,
    std::vector<std::vector<TVector2>> polygons,
    bool maskHistograms = false  // The function changes the histograms (map over used polygonal regions)
) {
  std::vector<double> efficiencies;
  std::set<std::pair<int, int>> notToDelete;

  int nx = numerator->GetNbinsX();
  int ny = numerator->GetNbinsY();

  for (auto &poly : polygons) {
    poly = orderCCW(poly);  // make sure polygon is CCW
    long long numSum = 0;
    long long denSum = 0;

    for (int ix = 1; ix <= nx; ix++) {
      for (int iy = 1; iy <= ny; iy++) {
        double x = numerator->GetXaxis()->GetBinCenter(ix);
        double y = numerator->GetYaxis()->GetBinCenter(iy);
        TVector2 point(x, y);

        if (pointInConvexPolygon(poly, point)) {
          if (maskHistograms) {
            notToDelete.insert({ix, iy});
          }
          numSum += numerator->GetBinContent(ix, iy);
          denSum += denominator->GetBinContent(ix, iy);
        }
      }
    }
    if (numSum > denSum) {
      std::cout << "Numerator sum greater than denominator sum" << " " << numSum << " " << denSum << '\n';
    }
    double eff = (denSum > 0) ? (static_cast<double>(numSum) / denSum) : 0;
    efficiencies.push_back(eff);
  }

  if (maskHistograms) {
    for (int ix = 1; ix <= nx; ix++) {
      for (int iy = 1; iy <= ny; iy++) {
        if (notToDelete.count({ix, iy}) == 0) {
          numerator->SetBinContent(ix, iy, 0);
          numerator->SetBinError(ix, iy, 0);
          denominator->SetBinContent(ix, iy, 0);
          denominator->SetBinError(ix, iy, 0);
        }
      }
    }
  }

  return efficiencies;
}

class TimingEfficiencyPerRP {
public:
  TimingEfficiencyPerRP(const std::string &name,
                        TH2F *num,
                        TH2F *den,
                        const std::vector<std::vector<TVector2>> &polygons)
      : name(name) {
    if (num)
      numerator = (TH2F *)num->Clone((name + "_num").c_str());
    if (den)
      denominator = (TH2F *)den->Clone((name + "_den").c_str());

    for (auto &polygon : polygons) {
      this->polygons.push_back(orderCCW(polygon));
      if (!isConvex(this->polygons.back())) {
        Fatal("TimingEfficiencyPerRP::TimingEfficiencyPerRP", "Polygons have to be convex.");
      }
    }
  }

  std::vector<double> draw(TCanvas *c,
                           int pad = 1,
                           bool drawNumber = false,
                           bool drawPolygon = false,
                           bool maskOut = false,
                           Option_t *opt = "colz") {
    compute(maskOut);
    if (!efficiency) {  // Functions calculates efficiency, draws it on given canvas and returns the efficiency as vector
      Fatal("TimingEfficiencyPerRP::draw", "No efficiency plot!");
      return {};
    }

    c->cd(pad);
    efficiency->Draw(opt);
    efficiency->SetTitle(name.c_str());
    efficiency->SetStats(0);
    efficiency->SetMaximum(1);

    if (drawPolygon)
      drawPolygons();
    if (drawNumber)
      drawNumbers();

    return efficiencyPerDiamond;
  }

private:
  std::string name;
  TH2F *numerator = nullptr;
  TH2F *denominator = nullptr;
  std::vector<std::vector<TVector2>> polygons;
  std::vector<double> efficiencyPerDiamond;
  std::vector<TPolyLine *> polyObjects;

  TH2F *efficiency = nullptr;

  void compute(bool maskOut = false) {
    if (!numerator || !denominator) {
      Fatal("RPEfficiency::compute", "No numerator or denominator plot!");
      return;
    }
    efficiency = (TH2F *)numerator->Clone((name + "_eff").c_str());

    efficiencyPerDiamond = getEfficiency(efficiency, denominator, polygons, maskOut);

    efficiency->Divide(denominator);
    efficiency->SetStats(0);
    efficiency->SetMaximum(1);
  }

  void drawPolygons() {
    for (const auto &poly : polygons) {
      TPolyLine *pl = new TPolyLine(poly.size());
      for (int i = 0; i <= poly.size(); ++i) {
        int polygonPointIdx = i % poly.size();
        pl->SetPoint(i, poly[polygonPointIdx].X(), poly[polygonPointIdx].Y());
      }
      pl->SetLineColor(kBlack);
      pl->SetLineWidth(2);
      pl->Draw("L SAME");
      polyObjects.push_back(pl);  // keep it alive
    }
  }

  void drawNumbers() {
    for (int i = 0; i < polygons.size() && i < efficiencyPerDiamond.size(); ++i) {
      TVector2 centroid = findCentroid(polygons[i]);

      TLatex t;
      t.SetTextAlign(22);
      t.SetTextSize(0.05);
      t.DrawLatex(centroid.X(), centroid.Y(), Form("%.2f", efficiencyPerDiamond[i]));
    }
  }
};

// Vector calculations to find the diamond polygons
struct Coords {
  TVector2 topLeftPoint;
  double shiftNumber;
  double diamondAngle;
  double extraShiftNumber;

  Coords() = default;

  Coords(TVector2 topLeftPoint_, double shiftNumber_, double diamondAngle_, double extraShiftNumber_)
      : topLeftPoint(topLeftPoint_),
        shiftNumber(shiftNumber_),
        diamondAngle(diamondAngle_),
        extraShiftNumber(extraShiftNumber_) {}
};

std::vector<std::vector<TVector2>> approximateDiamondShape(Coords coords,
                                                           double sizeX = 4.5,
                                                           double sizeY = 4.5,
                                                           double shiftSize = 0.2) {  // diamondAngle in radians

  TVector2 diamondRight(sizeX, 0), diamondDown(0, -sizeY), diamondShift(0, shiftSize * coords.shiftNumber);

  if (coords.diamondAngle != 0) {
    diamondRight = diamondRight.Rotate(coords.diamondAngle);
    diamondDown = diamondDown.Rotate(coords.diamondAngle);
    diamondShift = diamondShift.Rotate(coords.diamondAngle);
  }

  coords.topLeftPoint -= diamondShift;

  return std::vector<std::vector<TVector2>>{{coords.topLeftPoint,
                                             coords.topLeftPoint + diamondRight,
                                             coords.topLeftPoint + diamondDown,
                                             coords.topLeftPoint + diamondDown + diamondRight},
                                            {coords.topLeftPoint + diamondRight,
                                             coords.topLeftPoint + diamondRight + diamondRight,
                                             coords.topLeftPoint + diamondDown + diamondRight,
                                             coords.topLeftPoint + diamondRight + diamondRight + diamondDown},
                                            {coords.topLeftPoint + 2 * (diamondRight),
                                             coords.topLeftPoint + 2 * (diamondRight) + diamondRight,
                                             coords.topLeftPoint + 2 * (diamondRight) + diamondDown,
                                             coords.topLeftPoint + 2 * (diamondRight) + diamondDown + diamondRight}};
}

using json = nlohmann::json;

std::map<std::string, Coords> readCoordsJson(std::string coordsJson, int targetRun) {
  /*
  keys: 
  45cyl
  45box
  56cyl
  56 box
  */
  std::ifstream file(coordsJson);
  if (!file) {
    Fatal("readCoordsJson", "Cannot open file %s", coordsJson.c_str());
    return {};
  }

  json data;
  try {
    file >> data;
  } catch (const std::exception &e) {
    Fatal("readCoordsJson", "%s", e.what());
    return {};
  }

  std::vector<int> runs;
  for (auto &el : data.items()) {
    try {
      runs.push_back(std::stoi(el.key()));
    } catch (...) {
      std::cout << "Warning: Skipping non-numeric key " << el.key() << "\n";
    }
  }

  if (runs.empty()) {
    Fatal("readCoordsJson", "No numeric run keys found in JSON");
    return {};
  }

  std::sort(runs.begin(), runs.end());

  int candidateCoordsNum = -1;
  int candidateShiftNum = -1;

  for (int r : runs) {
    if (r <= targetRun) {
      auto runData = data[std::to_string(r)];
      if (runData.contains("45") && runData["45"].contains("cylX")) {
        candidateCoordsNum = r;
      }
      if (runData.contains("45") && runData["45"].contains("cylShift")) {
        candidateShiftNum = r;
      }
    } else {
      break;
    }
  }

  if (candidateCoordsNum == -1) {
    Fatal("readCoordsJson", "No run with absolute coords found less than or equal to %d", targetRun);
    return {};
  }

  auto CoordsRunData = data[std::to_string(candidateCoordsNum)];

  auto ret = std::map<std::string, Coords>{{"45box", {}}, {"56box", {}}, {"45cyl", {}}, {"56cyl", {}}};

  if (candidateShiftNum != -1) {
    auto shiftRunData = data[std::to_string(candidateShiftNum)];
    if (shiftRunData.count("45") > 0) {
      if (shiftRunData["45"].count("boxShift") > 0) {
        ret["45box"].shiftNumber = shiftRunData["45"]["boxShift"].get<double>();
      } else {
        Fatal("readCoordsJson", "no boxShift in shiftData[45]");
        return {};
      }
      if (shiftRunData["45"].count("cylShift") > 0) {
        ret["45cyl"].shiftNumber = shiftRunData["45"]["cylShift"].get<double>();
      } else {
        Fatal("readCoordsJson", "no cylShift in shiftData[45]");
        return {};
      }
    } else {
      Fatal("readCoordsJson", "no 45 in shift data");
      return {};
    }

    if (shiftRunData.count("56") > 0) {
      if (shiftRunData["56"].count("boxShift") > 0) {
        ret["56box"].shiftNumber = shiftRunData["56"]["boxShift"].get<double>();
      } else {
        Fatal("readCoordsJson", "no boxShift in shiftData[56]");
        return {};
      }
      if (shiftRunData["56"].count("cylShift") > 0) {
        ret["56cyl"].shiftNumber = shiftRunData["56"]["cylShift"].get<double>();
      } else {
        Fatal("readCoordsJson", "no cylShift in shiftData[56]");
        return {};
      }
    } else {
      Fatal("readCoordsJson", "no 56 in shift data");
      return {};
    }
  }

  // Read absolute coordinates from candidateCoordsNum
  if (CoordsRunData.count("45") > 0) {
    auto run45 = CoordsRunData["45"];
    // Cylinders
    if (run45.count("cylX") > 0 && run45.count("cylY") > 0 && run45.count("cylAngle") > 0) {
      ret["45cyl"].topLeftPoint = TVector2(run45["cylX"].get<double>(), run45["cylY"].get<double>());
      ret["45cyl"].diamondAngle = run45["cylAngle"].get<double>() * TMath::DegToRad();
    } else {
      Fatal("readCoordsJson", "missing cylinder coords in 45");
      return {};
    }

    // Boxes
    if (run45.count("boxX") > 0 && run45.count("boxY") > 0 && run45.count("boxAngle") > 0) {
      ret["45box"].topLeftPoint = TVector2(run45["boxX"].get<double>(), run45["boxY"].get<double>());
      ret["45box"].diamondAngle = run45["boxAngle"].get<double>() * TMath::DegToRad();
    } else {
      Fatal("readCoordsJson", "missing box coords in 45");
      return {};
    }
  } else {
    Fatal("readCoordsJson", "no 45 section in coords data");
    return {};
  }

  if (CoordsRunData.count("56") > 0) {
    auto run56 = CoordsRunData["56"];
    // Cylinders
    if (run56.count("cylX") > 0 && run56.count("cylY") > 0 && run56.count("cylAngle") > 0) {
      ret["56cyl"].topLeftPoint = TVector2(run56["cylX"].get<double>(), run56["cylY"].get<double>());
      ret["56cyl"].diamondAngle = run56["cylAngle"].get<double>() * TMath::DegToRad();
    } else {
      Fatal("readCoordsJson", "missing cylinder coords in 56");
      return {};
    }

    // Boxes
    if (run56.count("boxX") > 0 && run56.count("boxY") > 0 && run56.count("boxAngle") > 0) {
      ret["56box"].topLeftPoint = TVector2(run56["boxX"].get<double>(), run56["boxY"].get<double>());
      ret["56box"].diamondAngle = run56["boxAngle"].get<double>() * TMath::DegToRad();
    } else {
      Fatal("readCoordsJson", "missing box coords in 56");
      return {};
    }
  } else {
    Fatal("readCoordsJson", "no 56 section in coords data");
    return {};
  }

  return ret;
}

void CalculateDiamondEfficiency(TString infile = "timingHistograms.root",
                                TString coordsJson = "coords.json",
                                TString outputJson = "defficiency.json",
                                TString runNumber = "000000") {
  gStyle->SetPalette(kBird);

  TFile *f = TFile::Open(infile);

  TH2F *den45 = (TH2F *)f->Get("diamondHistograms/deffden45");
  TH2F *den56 = (TH2F *)f->Get("diamondHistograms/deffden56");

  TH2F *num45 = (TH2F *)f->Get("diamondHistograms/deffnum45");
  TH2F *num56 = (TH2F *)f->Get("diamondHistograms/deffnum56");

  TH2F *boxnum45 = (TH2F *)f->Get("diamondHistograms/dboxeffnum45");
  TH2F *boxnum56 = (TH2F *)f->Get("diamondHistograms/dboxeffnum56");

  bool drawNumbers = true;
  bool drawPolygons = true;
  bool maskOut = false;

  double diamondSizeX = 4.6;
  double diamondSizeY = 4.65;
  double diamondShiftSize = 0.2;

  auto coords = readCoordsJson(std::string(coordsJson.Data()), std::stoi(runNumber.Data()));

  std::vector<std::string> keys = {"45cyl", "45box", "56cyl", "56box"};
  bool containsAll = std::all_of(
      keys.begin(), keys.end(), [&coords](const std::string &key) { return coords.find(key) != coords.end(); });

  if (!containsAll) {
    Fatal("CalculateDiamondEfficiency", "Some keys are missing - input JSON format not supported.");
    return;
  }

  TimingEfficiencyPerRP box45efficiency(
      "45 box time-tracks efficiency per diamond",
      boxnum45,
      den45,
      approximateDiamondShape(coords["45box"], diamondSizeX, diamondSizeY, diamondShiftSize));

  TimingEfficiencyPerRP box56efficiency(
      "56 box time-tracks efficiency per diamond",
      boxnum56,
      den56,
      approximateDiamondShape(coords["56box"], diamondSizeX, diamondSizeY, diamondShiftSize));

  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->Divide(2, 1);
  auto box45efficiencyEffV = box45efficiency.draw(c1, 1, drawNumbers, drawPolygons, maskOut);
  auto box56efficiencyEffV = box56efficiency.draw(c1, 2, drawNumbers, drawPolygons, maskOut);

  TimingEfficiencyPerRP cylindrical45efficiency(
      "45 cylindrical time-tracks efficiency per diamond",
      num45,
      den45,
      approximateDiamondShape(coords["45cyl"], diamondSizeX, diamondSizeY, diamondShiftSize));

  TimingEfficiencyPerRP cylindrical56efficiency(
      "56 cylindrical time-tracks efficiency per diamond",
      num56,
      den56,
      approximateDiamondShape(coords["56cyl"], diamondSizeX, diamondSizeY, diamondShiftSize));

  TCanvas *c2 = new TCanvas("c2", "c2");
  c2->Divide(2, 1);
  auto cylindrical45efficiencyEffV = cylindrical45efficiency.draw(c2, 1, drawNumbers, drawPolygons, maskOut);
  auto cylindrical56efficiencyEffV = cylindrical56efficiency.draw(c2, 2, drawNumbers, drawPolygons, maskOut);

  TObjArray *tx = infile.Tokenize(".");
  TString outplot = ((TObjString *)(tx->At(0)))->String();
  c1->SaveAs("BoxTrackPerDiamondEfficiency_" + outplot + ".pdf");
  c2->SaveAs("CylindricalTrackPerDiamondEfficiency_" + outplot + ".pdf");

  // save efficiency to JSON:
  json j;
  std::string runStr = std::string(runNumber.Data());
  std::string outputJsonStr = std::string(outputJson.Data());

  j[runStr]["45"]["box"] = box45efficiencyEffV;
  j[runStr]["45"]["cyl"] = cylindrical45efficiencyEffV;
  j[runStr]["56"]["box"] = box56efficiencyEffV;
  j[runStr]["56"]["cyl"] = cylindrical56efficiencyEffV;

  std::ofstream file(outputJsonStr);
  if (!file.is_open()) {
    Fatal("CalculateDiamondEfficiency", "Cannot open output file to save JSON.");
    return;
  }

  file << j.dump(4);
  std::cout << "JSON saved to " << outputJsonStr << '\n';
}
