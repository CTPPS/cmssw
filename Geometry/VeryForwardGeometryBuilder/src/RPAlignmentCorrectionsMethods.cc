/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *	Jan Kašpar (jan.kaspar@gmail.com)
 *
 ****************************************************************************/

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/typelookup.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSPixelDetId.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"

#include "Geometry/VeryForwardGeometryBuilder/interface/RPAlignmentCorrectionsMethods.h"

#include <set>

#include "TMatrixD.h"
#include "TVectorD.h"

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>

using namespace std;
using namespace xercesc;

//----------------------------------------------------------------------------------------------------

RPAlignmentCorrectionsData RPAlignmentCorrectionsMethods::GetCorrectionsDataFromFile(const string &fileName)
{
  printf(">> RPAlignmentCorrectionsMethods::LoadXMLFile(%s)\n", fileName.c_str());

  // prepend CMSSW src dir
  char *cmsswPath = getenv("CMSSW_BASE");
  size_t start = fileName.find_first_not_of("   ");
  string fn = fileName.substr(start);
  if (cmsswPath && fn[0] != '/' && fn.find("./") != 0)
    fn = string(cmsswPath) + string("/src/") + fn;

  // load DOM tree first the file
  try {
    XMLPlatformUtils::Initialize();
  }
  catch (const XMLException& toCatch)
  {
    char* message = XMLString::transcode(toCatch.getMessage());
    throw cms::Exception("RPAlignmentCorrectionsMethods") << "An XMLException caught with message: " << message << ".\n";
    XMLString::release(&message);
  }

  XercesDOMParser* parser = new XercesDOMParser();
  parser->setValidationScheme(XercesDOMParser::Val_Always);
  parser->setDoNamespaces(true);

  try {
    parser->parse(fn.c_str());
  }
  catch (...)
  {
    throw cms::Exception("RPAlignmentCorrectionsMethods") << "Cannot parse file `" << fn << "' (exception)." << endl;
  }

  if (!parser)
    throw cms::Exception("RPAlignmentCorrectionsMethods") << "Cannot parse file `" << fn << "' (parser = NULL)." << endl;

  DOMDocument* xmlDoc = parser->getDocument();

  if (!xmlDoc)
    throw cms::Exception("RPAlignmentCorrectionsMethods") << "Cannot parse file `" << fn << "' (xmlDoc = NULL)." << endl;

  DOMElement* elementRoot = xmlDoc->getDocumentElement();
  if (!elementRoot)
    throw cms::Exception("RPAlignmentCorrectionsMethods") << "File `" << fn << "' is empty." << endl;

  RPAlignmentCorrectionsData d = GetCorrectionsData(elementRoot);

  XMLPlatformUtils::Terminate();

  return d;
}

//----------------------------------------------------------------------------------------------------

RPAlignmentCorrectionsData RPAlignmentCorrectionsMethods::GetCorrectionsData(DOMNode *root)
{
  RPAlignmentCorrectionsData result;

  DOMNodeList *children = root->getChildNodes();
  for (unsigned int i = 0; i < children->getLength(); i++)
  {
    DOMNode *n = children->item(i);
    if (n->getNodeType() != DOMNode::ELEMENT_NODE)
      continue;

    // check node type
    unsigned char nodeType = 0;
    if (!strcmp(XMLString::transcode(n->getNodeName()), "det")) nodeType = 1;
    if (!strcmp(XMLString::transcode(n->getNodeName()), "rp")) nodeType = 2;

    if (!nodeType)
      throw cms::Exception("RPAlignmentCorrectionsMethods") << "Unknown node `" << XMLString::transcode(n->getNodeName()) << "'.";

    // check children
    if (n->getChildNodes()->getLength() > 0)
    {
        edm::LogProblem("RPAlignmentCorrectionsMethods") << ">> RPAlignmentCorrectionsMethods::LoadXMLFile > Warning: tag `" <<
          XMLString::transcode(n->getNodeName()) << "' has " << n->getChildNodes()->getLength() << 
          " children nodes - they will be all ignored.";
    }

    // TODO: add shr2

    // default values
    double sh_r1 = 0., sh_r2 = 0., sh_x = 0., sh_y = 0., sh_z = 0., rot_z = 0.;
    double sh_r1_e = 0., sh_r2_e =0., sh_x_e = 0., sh_y_e = 0., sh_z_e = 0., rot_z_e = 0.;
    unsigned int id = 0;
    bool idSet = false;

    // get attributes
    DOMNamedNodeMap* attr = n->getAttributes();
    for (unsigned int j = 0; j < attr->getLength(); j++)
    {    
      DOMNode *a = attr->item(j);

      //printf("\t%s\n", XMLString::transcode(a->getNodeName()));

      if (!strcmp(XMLString::transcode(a->getNodeName()), "id"))
      {
        id = atoi(XMLString::transcode(a->getNodeValue()));
        idSet = true;
      } else if (!strcmp(XMLString::transcode(a->getNodeName()), "sh_r1"))
          sh_r1 = atof(XMLString::transcode(a->getNodeValue()));
        else if (!strcmp(XMLString::transcode(a->getNodeName()), "sh_r1_e"))
          sh_r1_e = atof(XMLString::transcode(a->getNodeValue()));
          else if (!strcmp(XMLString::transcode(a->getNodeName()), "sh_r2"))
              sh_r2 = atof(XMLString::transcode(a->getNodeValue()));
            else if (!strcmp(XMLString::transcode(a->getNodeName()), "sh_r2_e"))
              sh_r2_e = atof(XMLString::transcode(a->getNodeValue()));
              else if (!strcmp(XMLString::transcode(a->getNodeName()), "sh_x"))
                sh_x = atof(XMLString::transcode(a->getNodeValue()));
                else if (!strcmp(XMLString::transcode(a->getNodeName()), "sh_x_e"))
                  sh_x_e = atof(XMLString::transcode(a->getNodeValue()));
                  else if (!strcmp(XMLString::transcode(a->getNodeName()), "sh_y"))
                    sh_y = atof(XMLString::transcode(a->getNodeValue()));
                    else if (!strcmp(XMLString::transcode(a->getNodeName()), "sh_y_e"))
                      sh_y_e = atof(XMLString::transcode(a->getNodeValue()));
                      else if (!strcmp(XMLString::transcode(a->getNodeName()), "sh_z"))
                        sh_z = atof(XMLString::transcode(a->getNodeValue()));
                        else if (!strcmp(XMLString::transcode(a->getNodeName()), "sh_z_e"))
                          sh_z_e = atof(XMLString::transcode(a->getNodeValue()));
                          else if (!strcmp(XMLString::transcode(a->getNodeName()), "rot_z"))
                            rot_z = atof(XMLString::transcode(a->getNodeValue()));
                            else if (!strcmp(XMLString::transcode(a->getNodeName()), "rot_z_e"))
                              rot_z_e = atof(XMLString::transcode(a->getNodeValue()));
                            else
                              edm::LogProblem("RPAlignmentCorrectionsMethods") << ">> RPAlignmentCorrectionsMethods::LoadXMLFile > Warning: unknown attribute `"
                                << XMLString::transcode(a->getNodeName()) << "'.";
    }

    // id must be set
    if (!idSet)
        throw cms::Exception("RPAlignmentCorrectionsMethods") << "Id not set for tag `" << XMLString::transcode(n->getNodeName()) << "'.";

    // build alignment
    RPAlignmentCorrectionData a(sh_r1*1E-3, sh_r1_e*1E-3, sh_r2*1E-3, sh_r2_e*1E-3, sh_x*1E-3, sh_x_e*1E-3, sh_y*1E-3, sh_y_e*1E-3,
      sh_z*1E-3, sh_z_e*1E-3, rot_z*1E-3, rot_z_e*1E-3);

    // add the alignment to the right list
    if (nodeType == 1)
    {
      result.AddSensorCorrection(id, a, true);
    }

    if (nodeType == 2)
    {
      result.AddRPCorrection(id, a, true);
    }
  }

  return result;
}

//----------------------------------------------------------------------------------------------------

#define WRITE(q, dig, lim) \
  if (precise) \
    fprintf(f, " " #q "=\"%.15E\"", data.q()*1E3);\
  else \
    if (fabs(data.q()*1E3) < lim && data.q() != 0) \
      fprintf(f, " " #q "=\"%+8.1E\"", data.q()*1E3);\
    else \
      fprintf(f, " " #q "=\"%+8." #dig "f\"", data.q()*1E3);

void RPAlignmentCorrectionsMethods::WriteXML(const RPAlignmentCorrectionData &data, FILE *f, bool precise, bool wrErrors, bool wrSh_r, bool wrSh_xy,
  bool wrSh_z, bool wrRot_z)
{
  if (wrSh_r)
  {
    WRITE(sh_r1, 2, 0.02);
    WRITE(sh_r2, 2, 0.02);
    if (wrErrors)
    {
      WRITE(sh_r1_e, 2, 0.02);
      WRITE(sh_r2_e, 2, 0.02);
    }
  }

  if (wrSh_xy)
  {
    WRITE(sh_x, 2, 0.02);
    WRITE(sh_y, 2, 0.02);
    if (wrErrors)
    {
      WRITE(sh_x_e, 2, 0.02);
      WRITE(sh_y_e, 2, 0.02);
    }
  }

  if (wrRot_z)
  {
    WRITE(rot_z, 3, 0.01);
    if (wrErrors)
    {
      WRITE(rot_z_e, 3, 0.01);
    }
  }

  if (wrSh_z)
  {
    WRITE(sh_z, 2, 0.1);
    if (wrErrors)
    {
      WRITE(sh_z_e, 2, 0.1);
    }
  }
}

#undef WRITE

//----------------------------------------------------------------------------------------------------

void RPAlignmentCorrectionsMethods::WriteXMLFile(const RPAlignmentCorrectionsData & data, const string &fileName,
  bool precise, bool wrErrors, bool wrSh_r, bool wrSh_xy, bool wrSh_z, bool wrRot_z)
{
  FILE *rf = fopen(fileName.c_str(), "w");
  if (!rf)
    throw cms::Exception("RPAlignmentCorrections::WriteXMLFile") << "Cannot open file `" << fileName
      << "' to save alignments." << endl;

  // TODO
  //fprintf(rf, "<!--\nShifts in um, rotations in mrad.\n\nFor more details see RPAlignmentCorrections::LoadXMLFile in\n");
  //fprintf(rf, "Alignment/RPDataFormats/src/RPAlignmentCorrectionsSequence.cc\n-->\n\n");
  fprintf(rf, "<xml DocumentType=\"AlignmentDescription\">\n");

  WriteXMLBlock(data, rf, precise, wrErrors, wrSh_r, wrSh_xy, wrSh_z, wrRot_z);

  fprintf(rf, "</xml>\n");
  fclose(rf);
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentCorrectionsMethods::WriteXMLBlock(const RPAlignmentCorrectionsData & data, FILE *rf, bool precise, bool wrErrors, bool wrSh_r,
  bool wrSh_xy, bool wrSh_z, bool wrRot_z)
{
  bool firstRP = true;
  unsigned int prevRP = 0;
  set<unsigned int> writtenRPs;

  RPAlignmentCorrectionsData::mapType sensors = data.GetSensorMap();
  RPAlignmentCorrectionsData::mapType rps = data.GetRPMap();

  for (auto it = sensors.begin(); it != sensors.end(); ++it)
  {
    CTPPSDetId sensorId(it->first);
    unsigned int rpId = sensorId.getRPId();
    unsigned int decRPId = sensorId.arm()*100 + sensorId.station()*10 + sensorId.rp();

    // start a RP block
    if (firstRP || prevRP != rpId)
    {
      if (!firstRP)
        fprintf(rf, "\n");
      firstRP = false;

      fprintf(rf, "\t<!-- RP %3u -->\n", decRPId);

      auto rit = rps.find(rpId);
      if (rit != rps.end())
      {
        fprintf(rf, "\t<rp id=\"%u\" ", rit->first);
        WriteXML( rit->second , rf, precise, wrErrors, false, wrSh_xy, wrSh_z, wrRot_z );
        fprintf(rf, "/>\n");
        writtenRPs.insert(rpId);
      }
    }
    prevRP = rpId;

    // write plane id
    unsigned int planeIdx = 1000;
    if (sensorId.subdetId() == CTPPSDetId::sdTrackingStrip) planeIdx = TotemRPDetId(it->first).plane();
    if (sensorId.subdetId() == CTPPSDetId::sdTrackingPixel) planeIdx = CTPPSPixelDetId(it->first).plane();
    if (sensorId.subdetId() == CTPPSDetId::sdTimingDiamond) planeIdx = CTPPSDiamondDetId(it->first).plane();
    fprintf(rf, "\t<!-- plane %u --> ", planeIdx);

    // write the correction
    fprintf(rf, "<det id=\"%u\"", it->first);
    WriteXML(it->second, rf, precise, wrErrors, wrSh_r, wrSh_xy, wrSh_z, wrRot_z);
    fprintf(rf, "/>\n");
  }

  // write remaining RPs
  for (auto it = rps.begin(); it != rps.end(); ++it)
  {
    set<unsigned int>::iterator wit = writtenRPs.find(it->first);
    if (wit == writtenRPs.end())
    {
      fprintf(rf, "\t<rp id=\"%u\"                                ", it->first);
      WriteXML(it->second, rf, precise, wrErrors, false, wrSh_xy, wrSh_z, wrRot_z);
      fprintf(rf, "/>\n");
    }
  }
}
