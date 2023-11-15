// =========================================================================
// File:       tledModel.cpp
// Purpose:    Class for loading and parsing an XML model file
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    July 2008
// 
// Copyright (c) 2010, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================


#define _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_DEPRECATE

#include "tledModel.h"
#include "tledMatrixFunctions.h"
#include "tledVTKMeshLoader.h"
#include "tledSubModelManager.h"
#include "tledVectorArithmetic.h"
#include "tledHelper.h"
#include "tledMeshSurface.h"
#include "tledShellMaterialLinearPlateDecorator.h"
#include "tledShellMaterialLinearThickPlateDecorator.h"
#include "tledSurface.h"
#include "tledNodeRejector.h"

#include <algorithm>
#include <sstream>
#include <iterator>

using namespace std;

XMLNode tledModel::GetMandatoryXMLNode(const std::string &nodeName, const XMLNode &parentNode) {
  XMLNode node;

  if (parentNode.nChildNode(nodeName.c_str()) > 0) {
    node = parentNode.getChildNode(nodeName.c_str());
  } else {
    tledLogErrorStream(tledHelper::FatalError() << "No node " << nodeName << " found below " << parentNode.getName());    
  }

  return node;
}

XMLNode tledModel::GetSystemParamsNode(void) const {
  return tledModel::GetMandatoryXMLNode("SystemParams", xModel);
}

void tledModel::LoadFromXML(const XMLNode &rootNode) {
  xModel = rootNode;
  if (xModel.nChildNode("SubModel")) {
    int sMInd;
       
    mp_SubModelManager = new tledSubModelManager;
    if (xModel.nChildNode("SubModelNodeMergerDistance")) {
      stringstream ss;
      float dist;

      ss << xModel.getChildNode("SubModelNodeMergerDistance").getText();
      ss >> dist;

      if (ss.fail()) {
	tledFatalError("Could not parse value of SubModelNodeMergerDistance tag.");
      }

      mp_SubModelManager->SetMinNodeDistance(dist);
    } else mp_SubModelManager->SetMinNodeDistance(-1);

    mp_SubModelManager->SetNumberOfSubModels(xModel.nChildNode("SubModel"));
    for (sMInd = 0; sMInd < xModel.nChildNode("SubModel"); sMInd++) {
      XMLNode subXMLModel;

      subXMLModel = xModel.getChildNode("SubModel", sMInd);
      mp_SubModelManager->AddSubModel(subXMLModel);

      {
	const int numElSets = mp_SubModelManager->GetNumberOfSubModelElementSets(sMInd);

	int elSetInd;

	for (elSetInd = 0; elSetInd < numElSets; elSetInd++) {
	  xModel.addChild(mp_SubModelManager->GetSubModelElementSet(sMInd, elSetInd));
	}
      }

      {
	const int numElSets = mp_SubModelManager->GetNumberOfSubModelShellElementSets(sMInd);

	int elSetInd;

	for (elSetInd = 0; elSetInd < numElSets; elSetInd++) {
	  xModel.addChild(mp_SubModelManager->GetSubModelShellElementSet(sMInd, elSetInd));
	}
      }

      {
	const int numConsts = mp_SubModelManager->GetNumberOfSubModelConstraints(sMInd);

	int constInd;

	for (constInd = 0; constInd < numConsts; constInd++) {
	  xModel.addChild(mp_SubModelManager->GetSubModelConstraint(sMInd, constInd));
	}
      }
    }       

    Mesh = new tledMesh();    
    mp_SubModelManager->WriteMesh(*Mesh);
    if (mp_SubModelManager->HasMembrane()) xModel.addChild(mp_SubModelManager->CreateShellMeshXML());
  } else {
    mp_SubModelManager = NULL;
    Mesh = new tledMesh(&xModel);      
  }
}

tledModel::tledModel(const XMLNode &rootNode) : xModel(rootNode) {
  LoadFromXML(rootNode);
}

tledModel::tledModel(const char* xmlFileName)
{
   error = 0;
   // Load XML file
   XMLResults pResults;
   xModel = XMLNode::parseFile(xmlFileName,"Model",&pResults);
   if (pResults.error == 10) {
     tledLogErrorStream(tledHelper::NonFatalError() << "File " << xmlFileName << " not found");
     error = 1;
     return;
   }

   xFName = xmlFileName;
   LoadFromXML(xModel);
}

tledModel::~tledModel(void)
{
   delete Mesh;
   if (HasSubModels()) delete mp_SubModelManager;
   xModel.deleteNodeContent();
}

bool tledModel::DoDeformableDeformableContactHandling() const {
  if (xModel.getChildNode("SystemParams").nChildNode("DoDeformableCollision") > 0) {
    istringstream iss(xModel.getChildNode("SystemParams").getChildNode("DoDeformableCollision").getText());

    return *istream_iterator<int>(iss) != 0 || this->DoSelfCollisionContactHandling();
  } else return this->DoSelfCollisionContactHandling();
}

bool tledModel::DoSelfCollisionContactHandling() const {
  if (xModel.getChildNode("SystemParams").nChildNode("DoSelfCollision") > 0) {
    istringstream iss(xModel.getChildNode("SystemParams").getChildNode("DoSelfCollision").getText());

    return *istream_iterator<int>(iss) != 0;
  } else return false;
}

std::string tledModel::GetBoundingVolumeType() const {
  XMLNode sysParms = this->GetSystemParamsNode();

  if (sysParms.nChildNode("BVType") > 0) {
    return sysParms.getChildNode("BVType").getText();
  } 
  
  return "AABB";
}

bool tledModel::DoMultiPassContactHandling() const {
  XMLNode sysParms = this->GetSystemParamsNode();

  if (sysParms.nChildNode("DoDeformableCollision") > 0) {
    XMLNode defColl = sysParms.getChildNode("DoDeformableCollision");

    return (defColl.getAttribute("Type") != NULL && std::string(defColl.getAttribute("Type")) == "MultiPass");
  } 
  if (sysParms.nChildNode("DoSelfCollision") > 0) {
    XMLNode defColl = sysParms.getChildNode("DoSelfCollision");

    return (defColl.getAttribute("Type") != NULL && std::string(defColl.getAttribute("Type")) == "MultiPass");
  } 

  return sysParms.nChildNode("DoMultiPassContacts") > 0;
}

float tledModel::GetDeformableFrictionCoefficient() const {
  if (xModel.getChildNode("SystemParams").nChildNode("DeformableFrictionCoefficient") >= 1) {
    istringstream iss(xModel.getChildNode("SystemParams").getChildNode("DeformableFrictionCoefficient").getText());

    return *istream_iterator<float>(iss);
  } else return std::numeric_limits<float>::quiet_NaN();
}

float tledModel::GetContactSafetyMargin(void) const {
  if (xModel.getChildNode("SystemParams").nChildNode("ContactSafetyMargin") >= 1) {
    istringstream iss(xModel.getChildNode("SystemParams").getChildNode("ContactSafetyMargin").getText());

    return *istream_iterator<float>(iss);
  } else return std::numeric_limits<float>::quiet_NaN();
}

float tledModel::GetRateConstraintDistance() const {
  if (xModel.getChildNode("SystemParams").nChildNode("RateConstraintDistance") >= 1) {
    istringstream iss(xModel.getChildNode("SystemParams").getChildNode("RateConstraintDistance").getText());

    return *istream_iterator<float>(iss);
  } else return std::numeric_limits<float>::quiet_NaN();
}

int tledModel::GetNumberOfRigidContactSurfaces() const {
  return xModel.nChildNode("ContactSurface");
}

XMLNode tledModel::GetRigidContactSurfaceDefinition(const int si) const {
  assert(si < this->GetNumberOfRigidContactSurfaces());

  return xModel.getChildNode("ContactSurface", si);
}

string tledModel::GetDirectory(void)
{
  if (xFName != NULL) {
    string fname(xFName);
    size_t pos = fname.rfind("/");

    if (pos < fname.length()) {
      return fname.substr(0,pos+1);
    } else return std::string("");
  } else return std::string("");
}

vector<int> tledModel::GetElNodeInds(int ElNum) const
{
   return Mesh->GetElNodeInds(ElNum);
}

vector<float> tledModel::GetNodeCds(int NodeNum) const
{
   return Mesh->GetNodeCds(NodeNum);
}

int tledModel::GetNumElSets() const
{
   return xModel.nChildNode("ElementSet");
}

vector<int> tledModel::GetElSet(int ElSetNum) const
{
   XMLNode xNode = xModel.getChildNode("ElementSet",ElSetNum);
   int sz = GetElSetSize(ElSetNum);

   if (sz == ELEMENT_SET_SIZE_ALL) {
     return tledSequenceGenerator::MakeSequence(0, GetMesh()->GetNumEls());
   } else {
     vector<int> data = GetXMLTextAsVector<int>(xNode);

     if ((int)data.size() == sz) // Elements are listed explicitly
       return data;
     else if (data.size() == 1) // Starting element number is given
       {
	 vector<int> elset(sz,0);
	 for (int i = 0; i < sz; i++)
	   elset[i] = data[0] + i;
	 return elset;
       }
     else {
       tledNonFatalError("Invalid element set listing");
       return data;
     }
   }
}

int tledModel::GetElSetSize(int ElSetNum) const
{
  istringstream iss(xModel.getChildNode("ElementSet", ElSetNum).getAttribute("Size"));
  int size;

  iss >> size;
  if (iss.fail() || size > GetMesh()->GetNumEls() || size < 0) {
    if (iss.str() == "all") return ELEMENT_SET_SIZE_ALL;
    else {
      tledLogErrorStream(tledHelper::FatalError() << "Element set size attribute has to be a numeric value 0.." << GetMesh()->GetNumEls() << " or \"all\"");
      return 0;
    }
  } else return size;
}

const char* tledModel::GetMatType(int ElSetNum) const
{
  return this->GetElSetMaterialNode(ElSetNum).getAttribute("Type");
}

double tledModel::GetTimeStep(void) const
{
  XMLNode sysNode = this->GetSystemParamsNode();
  XMLNode dtNode = tledModel::GetMandatoryXMLNode("TimeStep", sysNode);

  if (dtNode.getText() != NULL) {
    return atof(dtNode.getText());
  } else {
    tledFatalError("Missing mandatory XML element: TimeStep");

    return std::numeric_limits<double>::quiet_NaN();
  }
}

double tledModel::GetTotalTime(void) const
{
  XMLNode sysNode = this->GetSystemParamsNode();
  XMLNode tNode = tledModel::GetMandatoryXMLNode("TotalTime", sysNode);

  if (tNode.getText() != NULL) {
    return atof(tNode.getText());
  } else {
    tledFatalError("Missing mandatory XML element: TotalTime");

    return std::numeric_limits<double>::quiet_NaN();
  }
}

float tledModel::GetDensity(void) const
{
  float rho = std::numeric_limits<float>::quiet_NaN();

  if (this->GetSystemParamsNode().nChildNode("Density") > 0) {
    std::istringstream iss(this->GetSystemParamsNode().getChildNode("Density").getText());

    iss >> rho;
    if (iss.fail()) {
      tledLogErrorStream(tledHelper::Warning() << "A default mass density was set and requested, but the XML node's value of " 
			 << xModel.getChildNode("SystemParams").getChildNode("Density").getText() << " could not be parsed.");
    } 
  }

  return rho;
}

float tledModel::GetDampingCoeff(void) const
{
   return (float)atof(xModel.getChildNode("SystemParams").getChildNode("DampingCoeff").getText());
}

const char* tledModel::GetElType(void) const
{
   return Mesh->GetElType();
}

int tledModel::GetNodesPerEl(void)
{
   return Mesh->GetNodesPerEl();
}

int tledModel::GetNumElasticParams(int ElSetNum) const
{
  return atoi(this->GetElSetMaterialNode(ElSetNum).getChildNode("ElasticParams").getAttribute("NumParams"));
}

int tledModel::GetMaxNumElasticParams(void) const
{
   // Loop over element sets and find the highest occuring NumElasticParams
   int maxNumElasticParams = 0;
   int NumElSets = this->GetNumElSets();
   for (int ElSetNum = 0; ElSetNum < NumElSets; ElSetNum++)
   {
      int currNumElasticParams = this->GetNumElasticParams(ElSetNum);
      maxNumElasticParams = currNumElasticParams > maxNumElasticParams ? currNumElasticParams : maxNumElasticParams;
   }

   return maxNumElasticParams;
}

int tledModel::GetNumViscIsoTerms(int ElSetNum) const
{
   const char* MType = this->GetMatType(ElSetNum);
   if (!strcmp(MType,"NHV") || !strcmp(MType,"TIV"))
   {
     return atoi(this->GetElSetMaterialNode(ElSetNum).getChildNode("ViscoParams").getAttribute("NumIsoTerms"));
   }
   else
      return 0;
}

int tledModel::GetNumViscVolTerms(const int ElSetNum) const
{
   const char* MType = this->GetMatType(ElSetNum);
   if (!strcmp(MType,"NHV") || !strcmp(MType,"TIV"))
   {
      return atoi(this->GetElSetMaterialNode(ElSetNum).getChildNode("ViscoParams").getAttribute("NumVolTerms"));
   }
   else
      return 0;
}

int tledModel::GetMaxNumViscIsoTerms(void) const
{
   // Loop over element sets and find the highest occuring numViscIsoTerms
   int maxNumViscIsoTerms = 0;
   int NumElSets = this->GetNumElSets();
   for (int ElSetNum = 0; ElSetNum < NumElSets; ElSetNum++)
   {
      int currNumViscIsoTerms = this->GetNumViscIsoTerms(ElSetNum);
      maxNumViscIsoTerms = currNumViscIsoTerms > maxNumViscIsoTerms ? currNumViscIsoTerms : maxNumViscIsoTerms;
   }

   return maxNumViscIsoTerms;
}

int tledModel::GetMaxNumViscVolTerms(void) const
{
   // Loop over element sets and find the highest occuring numViscVolTerms
   int maxNumViscVolTerms = 0;
   int NumElSets = this->GetNumElSets();
   for (int ElSetNum = 0; ElSetNum < NumElSets; ElSetNum++)
   {
      int currNumViscVolTerms = this->GetNumViscVolTerms(ElSetNum);
      maxNumViscVolTerms = currNumViscVolTerms > maxNumViscVolTerms ? currNumViscVolTerms : maxNumViscVolTerms;
   }

   return maxNumViscVolTerms;
}

XMLNode tledModel::GetElSetMaterialNode(const int elSetNum) const {
  XMLNode xNode = xModel.getChildNode("ElementSet", elSetNum);

  if (!xNode.isEmpty()) {
     xNode = xNode.getChildNode("Material");
     if (xNode.isEmpty()) {
       tledLogErrorStream(tledHelper::FatalError() << "Could not find a material node in element set " << elSetNum);
     }     
  } else {
    tledLogErrorStream(tledHelper::FatalError() << "Invalid element set index " << elSetNum);
  }

  return xNode;
}

float tledModel::GetElSetMassDensity(const int elSetNum) const {
  XMLNode xNode = this->GetElSetMaterialNode(elSetNum);
  float rho = this->GetDensity();

  if (xNode.nChildNode("Density") > 0) {
    std::istringstream iss(xNode.getChildNode("Density").getText());

    iss >> rho;
    if (iss.fail()) {
      tledLogErrorStream(tledHelper::Warning() << "Element set " << elSetNum << "'s mass density was set and requested, but the XML node's value of " 
			 << xNode.getChildNode("Density").getText() << " could not be parsed.");
    } 
  } 

  if (std::isnan(rho)) {
    tledLogErrorStream(tledHelper::FatalError() << "Have no valid density for element set " << elSetNum);
  }

  return rho;
}

void tledModel::GetElasticParams(const int ElSetNum, float* Params) const
{
  XMLNode xNode = this->GetElSetMaterialNode(ElSetNum);

   xNode = xNode.getChildNode("ElasticParams");
   int numElasticParams = this->GetNumElasticParams(ElSetNum);
   // Get params from xModel
   vector<float> xmlMatParams = GetXMLTextAsVector<float>(xNode);
   for (int i = 0; i < numElasticParams; i++)
   {
      Params[i] = xmlMatParams[i];
   }
}

void tledModel::GetViscoParams(const int ElSetNum, float* Params) const
{
   XMLNode xNode = this->GetElSetMaterialNode(ElSetNum);
   xNode = xNode.getChildNode("ViscoParams");
   int numViscParams = 2*(this->GetNumViscIsoTerms(ElSetNum) + this->GetNumViscVolTerms(ElSetNum));
   // Get params from xModel
   if (numViscParams > 0)
   {
      vector<float> xmlMatParams = GetXMLTextAsVector<float>(xNode);
      for (int i = 0; i < numViscParams; i++)
      {
         Params[i] = xmlMatParams[i];
      }
   }
}

void tledModel::GetCombinedMatParams(const int ElSetNum, float* Params) const
{
   XMLNode xNode = this->GetElSetMaterialNode(ElSetNum);
   // Get elastic params
   int numElasticParams = this->GetNumElasticParams(ElSetNum);
   float* elasticParams = new float[numElasticParams];
   this->GetElasticParams(ElSetNum,elasticParams);
   // Get visco params (if applicable)
   int numViscIsoTerms = this->GetNumViscIsoTerms(ElSetNum);
   int numViscVolTerms = this->GetNumViscVolTerms(ElSetNum);
   int numViscParams = 2*(numViscIsoTerms + numViscVolTerms);
   float* viscParams = new float[numViscParams];
   this->GetViscoParams(ElSetNum,viscParams);
   // Combine all params
   for (int i = 0; i < numElasticParams; i++)
   {
      Params[i] = elasticParams[i];
   }
   if (numViscParams > 0)
   {
      Params[numElasticParams] = (float)this->GetTimeStep();
      Params[numElasticParams+1] = (float)numViscIsoTerms;
      Params[numElasticParams+2] = (float)numViscVolTerms;
      for (int i = 0; i < numViscParams; i++)
      {
         Params[numElasticParams+3+i] = viscParams[i];
      }
   }

   delete[] elasticParams;
   delete[] viscParams;
}

float tledModel::GetHGKappa(void) const
{
   return (float)atof(xModel.getChildNode("SystemParams").getChildNode("HGKappa").getText());
}

//const bool tledModel::GetMeshColouring(void)
//{
//   int nmshClr = xModel.getChildNode("SystemParams").nChildNode("MeshColouring");
//   if (nmshClr == 0)
//      return false;
//   else
//   {
//      int mshClr = atoi(xModel.getChildNode("SystemParams").getChildNode("MeshColouring").getText());
//      if (mshClr == 0)
//         return false;
//      else
//         return true;
//   }
//}

// Get constraint data

int tledModel::GetNumConstraints(const std::string &reqConstraintType) const {
   int nConst = xModel.nChildNode("Constraint");	// Total number of constraints
   int count = 0;
   
   for (int i = 0; i < nConst; i++)
   {
      const char* currConstraintType = xModel.getChildNode("Constraint",i).getAttribute("Type");
      if (currConstraintType == reqConstraintType)
         count++;
   }

   return count;
}

int tledModel::GetNumConstraints(void) const {
   return xModel.nChildNode("Constraint");
}

static XMLNode _GetConstraintNodeSafely(const XMLNode root, const int constraintInd) {
  XMLNode constNode;

  if (constraintInd > root.nChildNode("Constraint")) {
    tledLogErrorStream(tledHelper::FatalError() << "Request for constraint " << constraintInd << ", but there are only " << root.nChildNode("Constraint") << " constraints.");
  } else {
    constNode = root.getChildNode("Constraint", constraintInd);
  }

  return constNode;
}

static const char* _GetConstraintAttributeSafely(const XMLNode root, const char attName[], const int constraintInd) {
  char const *pc_val = _GetConstraintNodeSafely(root, constraintInd).getAttribute(attName);

  if (pc_val == NULL) {
    tledLogErrorStream(tledHelper::FatalError() << "No attribute " << attName << " set on constraint " << constraintInd << ".");
  }

  return pc_val;
}

const char* tledModel::GetConstraintType(const int constraintInd) const {
   return _GetConstraintAttributeSafely(xModel, "Type", constraintInd);
}

int tledModel::GetConstraintDOF(const int constraintInd) const  {
  XMLNode constNode = _GetConstraintNodeSafely(xModel, constraintInd);
  int dof;

  if (constNode.getAttribute("DOF") == NULL || std::string(constNode.getAttribute("DOF")) == "all") dof = CONSTRAINT_ALL_DOFS;
  else {
    dof = atoi(_GetConstraintAttributeSafely(xModel, "DOF", constraintInd));
    if (dof < 0 || dof >= 3) {
      tledFatalError("Only DOF values allowed are 0..2 or \"all\"");
    }
  }

  return dof;
}

enum loadShape tledModel::GetConstraintLoadShape(const int constraintInd) const {
  return atols(_GetConstraintAttributeSafely(xModel, "LoadShape", constraintInd));
}

/** Returns the surface facet definitions as a 1D vector of non-unique node indices */
template <const int t_numFacetVertices>
static vector<int> _FindBoundaryThroughNormal(const tledMesh &mesh, const float n[], const float tolAngle) {
  tledMeshSurface<t_numFacetVertices> surface(mesh);
  vector<int> nodeInds;  
  float fN[3];
       
  for (int fInd = 0; fInd < surface.GetNumberOfFacets(); fInd++) {
    if (tledVectorArithmetic::ComputeAngleNormalised(surface.ComputeNormalisedFacetNormal(fN, fInd), n) < tolAngle) {
      nodeInds.insert(nodeInds.end(), surface.GetFacet(fInd).NodeIndices, surface.GetFacet(fInd).NodeIndices + t_numFacetVertices);
    }
  }

  return nodeInds;
}

template <const int t_numFacetVertices>
static vector<int> _FindBoundaryThroughNormalUnique(const tledMesh &mesh, const float n[], const float tolAngle) {
  return tledHelper::MakeSortedUnique(_FindBoundaryThroughNormal<t_numFacetVertices>(mesh, n, tolAngle));
}

static std::vector<int> _RestrictBoundary(std::vector<int> nodeInds, const XMLNode xNode, const float *nodes) {
  if (xNode.nChildNode("RestrictTo") > 0) {
    tledNodeRejector *p_rejector = tledNodeRejector::CreateRejector(xNode.getChildNode("RestrictTo", 0));

    for (int r = 1; r < xNode.nChildNode("RestrictTo"); r++) {
      p_rejector->AppendRejectorToChain(tledNodeRejector::CreateRejector(xNode.getChildNode("RestrictTo", r)));
    }
    p_rejector->SetNodeIndices(nodeInds);
    p_rejector->SetNodes(nodes);
    p_rejector->RunRejection();

    delete p_rejector;
  } 

  return nodeInds;
}

static bool _LoadNormalSpec(float *p_n, float &r_tolAngle, const XMLNode xNode) {
  using namespace tledVectorArithmetic;

  XMLNode nrmlNode;
  vector<float> nVec;

  if (xNode.nChildNode("Normal") < 1) {
    tledLogErrorStream(tledHelper::NonFatalError() << "Expected a \"Normal\" element in constraint of type " << xNode.getAttribute("Type") << " based on requested SpecType, but none found.");
    return false;
  }

  nrmlNode = xNode.getChildNode("Normal");     
  nVec = GetXMLTextAsVector<float>(nrmlNode);
  if (nVec.size() != 3) {
    tledNonFatalError("Error in boundary specification: Normal vector needs 3 components!");
    return false;
  }
  
  std::copy(nVec.begin(), nVec.end(), p_n);
  if (Norm(p_n) < 1e-6f) {
    tledLogErrorStream(tledHelper::NonFatalError() << "Error in boundary specification: Normal vector has insignificant magnitude (||n|| = " << Norm(p_n) << ")!");
    return false;
  }
  ScalarDiv(p_n, Norm(p_n));

  if (nrmlNode.getAttribute("ToleranceAngle") == NULL) {
    tledLogErrorStream(tledHelper::NonFatalError() << "Error in boundary specification: tolerance angle not found. Specified with normal-tag attribute \"ToleranceAngle\", in degrees.");
    return false;
  }

  r_tolAngle = std::numeric_limits<float>::quiet_NaN();
  istringstream(nrmlNode.getAttribute("ToleranceAngle")) >> r_tolAngle;
  if (std::isnan(r_tolAngle) || r_tolAngle < 0.0f || r_tolAngle > 180.0f) {
    tledNonFatalError("Error in boundary specification: Tolerance angle needs to be in [0, 180] (degrees).");
    return false;
  }
  r_tolAngle *= float(tledPi/180.0f);

  return true;
}

vector<int> tledModel::GetConstraintInd(const int ConstraintNum) const 
{
   XMLNode xNode = xModel.getChildNode("Constraint",ConstraintNum);

   if ((xNode.getAttribute("NumNodes") != NULL || (xNode.getAttribute("SpecType") != NULL && string(xNode.getAttribute("SpecType")) == "NODES")) && !(xNode.getAttribute("SpecType") != NULL && string(xNode.getAttribute("SpecType")) == "NORMAL")) {
     int sz;
     vector<int> Ind;

     if (xNode.getAttribute("NumNodes") == NULL) {
       tledFatalError("SpecType=\"NODES\" requires a \"NumNodes\" attribute.");
     }

     if (std::string(xNode.getAttribute("NumNodes")) == "all") {
       sz = this->GetMesh()->GetNumNodes();
       Ind = std::vector<int>(1, 0);
     } else {
       sz = atoi(xNode.getAttribute("NumNodes"));
       if (sz > 0) {
	 if (xNode.nChildNode("Nodes") > 0) {
	   Ind = GetXMLTextAsVector<int>(xNode.getChildNode("Nodes"));
	 } else {
	   tledLogErrorStream(tledHelper::FatalError() << "Expected to find a \"Nodes\" element in constraint " << ConstraintNum);
	 }
       }
     }

     if ((int)Ind.size() == sz) // Nodes are listed explicitly
       return Ind;
     else if (Ind.size() == 1) // Starting element number is given
       {
	 vector<int> ind(sz,0);
	 for (int i = 0; i < sz; i++)
	   ind[i] = Ind[0] + i;
	 return ind;
       }
     else
       {
	 tledLogErrorStream(tledHelper::NonFatalError() << "Invalid node listing for Constraint #" << ConstraintNum);
	 Ind.resize(0);
	 return Ind;
       }
   } else if (xNode.getAttribute("SpecType") != NULL && string(xNode.getAttribute("SpecType")) == "NORMAL") {
     float n[3], tolAngle;

     if (!_LoadNormalSpec(n, tolAngle, xNode)) {
       return std::vector<int>(0);
     }

     if (string(GetMesh()->GetElType()) == "T4" || string(GetMesh()->GetElType()) == "T4ANP") return _RestrictBoundary(_FindBoundaryThroughNormalUnique<3>(*GetMesh(), n, tolAngle), xNode, this->GetMesh()->GetAllNodeCds());
     else return _RestrictBoundary(_FindBoundaryThroughNormalUnique<4>(*GetMesh(), n, tolAngle), xNode, this->GetMesh()->GetAllNodeCds());
   } else {
     /* if node list spec else if spec thru normal */
     tledFatalError("Don't have a valid boundary specification mode.");
   }

   return vector<int>();
}

vector<float> tledModel::GetConstraintMag(const int ConstraintNum) const {
   XMLNode xNode = xModel.getChildNode("Constraint",ConstraintNum);
   int ConstraintSize = this->GetConstraintInd(ConstraintNum).size();
   xNode = xNode.getChildNode("Magnitudes");
   const char* MagType = xNode.getAttribute("Type");
   vector<float> Mag;
   Mag.reserve(ConstraintSize);
   if (!strcmp(MagType,"UNIFORM"))
   {
      float singleMag = (float)atof(xNode.getText());
      Mag.resize(ConstraintSize,singleMag);
   }
   else if (!strcmp(MagType,"DIFFORM"))
   {
      Mag = GetXMLTextAsVector<float>(xNode);
      if ((int)Mag.size() != ConstraintSize)
      {
	tledLogErrorStream(tledHelper::NonFatalError() << "Constraint #" << ConstraintNum << ": NumNodes = " << ConstraintSize 
			   << ", but " << (int)Mag.size() << " magnitudes are listed.");
      }
   }
   else
   {
     tledLogErrorStream(tledHelper::NonFatalError() << "Unknown constraint magnitude type: " << MagType);
   }

   return Mag;
}

// vector<float>* tledModel::GetConstraintMag(const int ConstraintNum)
// {
//    XMLNode xNode = xModel.getChildNode("Constraint",ConstraintNum);
//    int ConstraintSize = atoi(xNode.getAttribute("NumNodes"));
//    loadShape LS = atols(xNode.getAttribute("LoadShape"));
//    xNode = xNode.getChildNode("Magnitudes");
//    const char* MagType = xNode.getAttribute("Type");
//    vector<float>* Mag=NULL;
//    if (!strcmp(MagType,"UNIFORM"))
//    {
//       Mag = new vector<float>;
//       Mag->reserve(ConstraintSize);
//       float singleMag = (float)atof(xNode.getText());
//       Mag->resize(ConstraintSize,singleMag);
//    }
//    else if (!strcmp(MagType,"DIFFORM"))
//    {
//       Mag = new vector<float>;
//       Mag->reserve(ConstraintSize);
//       vector<float> xmlMag = getXMLTextAsVector<float>(xNode);
//       *Mag = xmlMag;
//       if (Mag->size() != ConstraintSize)
//       {
//          cerr << "\n!!! Constraint #" << ConstraintNum << ": NumNodes = " << ConstraintSize 
//             << ", but " << Mag->size() << " magnitudes are listed." << endl;
//       }
//    }
//    else if (!strcmp(MagType,"FILE"))
//    {
//       double T = this->GetTotalTime();
//       double dt = this->GetTimeStep();
//       int NumSteps = (int)(T/dt + 0.5);	// Actually this is NumSteps + 1
//       Mag = new vector<float>[NumSteps+1];
//       for (int i = 0; i < NumSteps+1; i++)
//       {
//          Mag[i].reserve(ConstraintSize);
//       }
//       // Load all displacement vals from file
//       vector<float>* allMags = new vector<float>;
//       const char* fname = xNode.getText();
//       int retval = loadFile(fname,allMags);
//       // Check size
//       if (allMags->size() != ConstraintSize*(NumSteps+1))
//          cerr << "\n!!! Constraint #" << ConstraintNum << ": Invalid number of magnitudes in file " << fname << endl;
//       else
//       {
//          // Separate magnitudes for each time step
//          vector<float>* currMag=NULL;
//          int count = 0;
//          for (int i = 0; i < NumSteps+1; i++)
//          {
//             currMag = &Mag[i];
//             for (int j = 0; j < ConstraintSize; j++)
//                currMag->push_back((*allMags)[count++]);
//          }
//       }
//    }
//    else
//    {
//       cerr << "\n!!! Unknown constraint magnitude type: " << MagType << endl;
//    }
// 
//    return Mag;
// }

float tledModel::GetGravityMag(const int ConstraintNum) const
{
   XMLNode xNode = xModel.getChildNode("Constraint",ConstraintNum);
   const char* type = xNode.getAttribute("Type");
   if (strcmp(type,"Gravity"))
   {
     tledLogErrorStream(tledHelper::NonFatalError() << "Constraint #" << ConstraintNum << " is not of type \"Gravity\" --> cannot retrieve acceleration magnitude.");
     return 0.0f;
   }
   
   return (float)atof(xNode.getChildNode("AccelerationMagnitude").getText());
}

vector<float> tledModel::GetGravityDirec(const int ConstraintNum) const
{
   XMLNode xNode = xModel.getChildNode("Constraint",ConstraintNum);
   const char* type = xNode.getAttribute("Type");
   if (strcmp(type,"Gravity"))
   {
     tledLogErrorStream(tledHelper::NonFatalError() << "Constraint #" << ConstraintNum << " is not of type \"Gravity\" --> cannot retrieve acceleration magnitude");
     vector<float> x;
     return x;
   }
   
   return GetXMLTextAsVector<float>(xNode.getChildNode("AccelerationDirection"));
}

int tledModel::_GetSurfaceConstraintNumFaces(const int ConstraintNum, const char type[]) const {
  const XMLNode pc = xModel.getChildNode("Constraint", ConstraintNum);
  const std::string cType = pc.getAttribute("Type");

  if (cType != type) {
    tledLogErrorStream(tledHelper::FatalError() << "Constraint #" << ConstraintNum << " is not of type \"" << type << "\" --> cannot retrieve NumFaces");
  } else {
    if (pc.getAttribute("SpecType") != NULL && std::string(pc.getAttribute("SpecType")) == "NORMAL") {
      return SURFACE_CONSTRAINT_FACE_SPEC_NORMAL;
    } else {
      const char* nFSpec = pc.getAttribute("NumFaces");

      if (nFSpec == NULL) {
	tledFatalError("Need the number of faces to which the pressure is to be applied (either \"all\" or a number >= 0).");
      } else {
	if (std::string(nFSpec) == "all") return SURFACE_CONSTRAINT_FACE_SPEC_ALL;
	else return atoi(pc.getAttribute("NumFaces"));
      }
    }
  }

  return 0;
}

const char* tledModel::_GetSurfaceConstraintFaceType(const int ConstraintNum, const char type[]) const {
  const std::string cType = xModel.getChildNode("Constraint",ConstraintNum).getAttribute("Type");

  if (cType != type) {
    tledLogErrorStream(tledHelper::FatalError() << "Constraint #" << ConstraintNum << " is not of type \"" << type << "\" --> cannot retrieve FaceType");
  } else {  
    if (_GetSurfaceConstraintNumFaces(ConstraintNum, type) < 0) {
      switch (this->GetMesh()->GetNodesPerEl()) {
      case 4:
	return "Tri";

      case 8:
	return "Quad";

      default:
	tledFatalError("Unsupported solid mesh type.");
      }
    } else return xModel.getChildNode("Constraint",ConstraintNum).getAttribute("FaceType");
  }

  return NULL;       
}

template <const int t_numFacetVertices>
static std::vector<int> _ExtractSurfaceConstraintMeshSurface(const tledMesh &solidMesh) {
  tledMeshSurface<t_numFacetVertices> surface(solidMesh);
  std::vector<int> facetDefs;
  
  facetDefs.reserve(t_numFacetVertices*surface.GetNumberOfFacets());
  for (int f = 0; f < surface.GetNumberOfFacets(); f++) {
    facetDefs.insert(facetDefs.end(), surface.GetFacet(f).NodeIndices, surface.GetFacet(f).NodeIndices + t_numFacetVertices);
  }
  
  return facetDefs;
}

vector<int> tledModel::_GetSurfaceConstraintFaceNodeInds(const int ConstraintNum, const char constType[]) const {
   vector<int> faces;
   int nConstraints = xModel.nChildNode("Constraint");

   if (nConstraints > ConstraintNum) {
     if (_GetSurfaceConstraintNumFaces(ConstraintNum, constType) == SURFACE_CONSTRAINT_FACE_SPEC_ALL) {
       if (this->GetMesh()->GetNodesPerEl() == 4) {
	 faces = _ExtractSurfaceConstraintMeshSurface<3>(*this->GetMesh());
       } else if (this->GetMesh()->GetNodesPerEl() == 8) {
	 faces = _ExtractSurfaceConstraintMeshSurface<4>(*this->GetMesh());
       } else {
	 tledFatalError("Element-type is not supported by facet extraction.");
       }
     } else if (_GetSurfaceConstraintNumFaces(ConstraintNum, constType) == SURFACE_CONSTRAINT_FACE_SPEC_NORMAL) {
       const XMLNode pc = xModel.getChildNode("Constraint", ConstraintNum);

       float n[3], tolAngle;

       assert(!(pc.getAttribute("Type") == NULL || std::string(pc.getAttribute("Type")) != constType));
       assert(pc.getAttribute("SpecType") != NULL && std::string(pc.getAttribute("SpecType")) == "NORMAL");
       if (_LoadNormalSpec(n, tolAngle, pc)) {
	 if (this->GetMesh()->GetNodesPerEl() == 4) {
	   faces = _FindBoundaryThroughNormal<3>(*this->GetMesh(), n, tolAngle);
	 } else {
	   assert(this->GetMesh()->GetNodesPerEl() == 8);
	   faces = _FindBoundaryThroughNormal<4>(*this->GetMesh(), n, tolAngle);
	 }
       } else return std::vector<int>(0);
     } else {
       int nFaces = xModel.getChildNode("Constraint",ConstraintNum).nChildNode("Faces");

       if (nFaces > 0) {
	 faces = GetXMLTextAsVector<int>(xModel.getChildNode("Constraint",ConstraintNum).getChildNode("Faces"));
       }
     }
   }
   
   return faces;
}

float tledModel::GetPressureMagnitude(const int ConstraintNum) const {
   float mag = 0.0f;
   int nConstraints = xModel.nChildNode("Constraint");
   if (nConstraints > ConstraintNum)
   {
      int nMags = xModel.getChildNode("Constraint",ConstraintNum).nChildNode("Magnitude");
      if (nMags > 0) {
	mag = (float)atof(xModel.getChildNode("Constraint",ConstraintNum).getChildNode("Magnitude").getText());
      }
   }
   
   return mag;
}

std::vector<float> tledModel::GetTractionFaceTractions(const int ConstraintNum) const {
  std::vector<float> ts;

  if (xModel.nChildNode("Constraint") > ConstraintNum) {
    XMLNode xNode = xModel.getChildNode("Constraint",ConstraintNum);
     
    if (xNode.nChildNode("Magnitudes")) {
      const int numConstFacets = this->GetTractionNumFaces(ConstraintNum);
      
      xNode = xNode.getChildNode("Magnitudes");
      if (std::string(xNode.getAttribute("Type")) == "UNIFORM") {
	ts = GetXMLTextAsVector<float>(xNode);
	
	if (ts.size() != 3) {
	  tledLogErrorStream(tledHelper::FatalError() << "Need one 3-component vector in uniform traction magnitude elements. Error in constraint " << ConstraintNum);
	}
      } else if (std::string(xNode.getAttribute("Type")) == "DIFFORM") {
	ts = GetXMLTextAsVector<float>(xNode);
	if (numConstFacets < 0 || int(ts.size()) != 3*numConstFacets) {
	  tledLogErrorStream(tledHelper::FatalError() << "Need 3 values/face with traction magnitude elements. Error in constraint " << ConstraintNum);
	}
      } else {
	tledLogErrorStream(tledHelper::FatalError() << "Need a magnitude type on constraint " << ConstraintNum);
      }
    } else {
      tledLogErrorStream(tledHelper::FatalError() << "Need a magnitude spec. for constraint " << ConstraintNum);
    }
  } else {
    tledLogErrorStream(tledHelper::FatalError() << "Need invalid constraint index: " << ConstraintNum);
  }
  
  return ts;
}

void tledModel::SetTimeStep(double dt)
{
   char* str = new char[15];
   sprintf(str,"%15.14f",dt);
   xModel.getChildNode("SystemParams").getChildNode("TimeStep").updateText(str);
   delete str;
}

void tledModel::SetConstraintMag(const int ConstraintNum, const std::vector<float> &Mag)
{
   // Check that Mag.size = NumNodes
   XMLNode xNode = xModel.getChildNode("Constraint",ConstraintNum);
   int NumNodes = atoi(xNode.getAttribute("NumNodes"));
   if ((int)Mag.size() != NumNodes)
   {
     tledLogErrorStream(tledHelper::NonFatalError() << "Invalid number of magnitude vals submitted - Attribute NumNodes = " << NumNodes);
     return;
   }

   // Update Type Att
   xNode = xNode.getChildNode("Magnitudes");
   xNode.updateAttribute("DIFFORM",NULL,"Type");
   // Update Text
   stringstream ssText;
   for (int i = 0; i < (int)Mag.size(); i++)
   {
      ssText << Mag[i];
      ssText << " ";
   }
   char* chrText = new char[ssText.str().size()+1];
   strcpy(chrText, ssText.str().c_str());
   xNode.updateText(chrText);

   delete chrText;
}

void tledModel::SetConstraintInd(const int ConstraintNum, const std::vector<int> &Ind)
{
   // Update NumNodes Att
   stringstream ssAtt;
   ssAtt << Ind.size();
   char* chrAtt = new char[ssAtt.str().size()+1];
   strcpy(chrAtt,ssAtt.str().c_str());
   XMLNode xNode = xModel.getChildNode("Constraint",ConstraintNum);
   xNode.updateAttribute(chrAtt,NULL,"NumNodes");
   // Update Text
   stringstream ssText;
   for (int i = 0; i < (int)Ind.size(); i++)
   {
      ssText << Ind[i];
      ssText << " ";
   }
   char* chrText = new char[ssText.str().size()+1];
   strcpy(chrText, ssText.str().c_str());
   xNode = xNode.getChildNode("Nodes");
   xNode.updateText(chrText);

   delete chrText;
   delete chrAtt;
}

void tledModel::SetPressureMag(const int ConstraintNum, const float Mag)
{
   XMLNode xNode = xModel.getChildNode("Constraint",ConstraintNum);
   // Check if constraint is a pressure constraint
   const char* cType = xNode.getAttribute("Type");
   if (strcmp(cType,"Pressure"))
   {
     tledLogErrorStream(tledHelper::NonFatalError() << "Constraint number " << ConstraintNum << " is not a Pressure constraint");
     return;
   }
   // Update Text
   stringstream ssText;
   ssText << Mag;
   char* chrText = new char[ssText.str().size()+1];
   strcpy(chrText, ssText.str().c_str());
   xNode.getChildNode("Magnitude").updateText(chrText);
   
   delete chrText;
}

int tledModel::GetNumContactObjects(void)
{
   return xModel.nChildNode("ContactCylinder") + xModel.nChildNode("ContactUSProbe") + xModel.nChildNode("ContactPlate");
}

int tledModel::GetNumContactCyls(void)
{
   return xModel.nChildNode("ContactCylinder");
}

vector<float> tledModel::GetContactCylOrigin(int cylNum)
{
   vector<float> orig;
   int nCyls = xModel.nChildNode("ContactCylinder");
   if (nCyls >= cylNum+1)
   {
      int nOrig = xModel.getChildNode("ContactCylinder",cylNum).nChildNode("Origin");
      if (nOrig > 0)
         orig = GetXMLTextAsVector<float>(xModel.getChildNode("ContactCylinder",cylNum).getChildNode("Origin"));
   }
   return orig;
}

vector<float> tledModel::GetContactCylAxis(int cylNum)
{
   vector<float> axis;
   int nCyls = xModel.nChildNode("ContactCylinder");
   if (nCyls >= cylNum+1)
   {
      int nAxis = xModel.getChildNode("ContactCylinder",cylNum).nChildNode("Axis");
      if (nAxis > 0)
         axis = GetXMLTextAsVector<float>(xModel.getChildNode("ContactCylinder",cylNum).getChildNode("Axis"));
   }
   return axis;
}

float tledModel::GetContactCylRadius(int cylNum)
{
   int nCyls = xModel.nChildNode("ContactCylinder");
   if (nCyls >= cylNum+1)
   {
      int nRad = xModel.getChildNode("ContactCylinder",cylNum).nChildNode("Radius");
      if (nRad > 0)
	return (float)atof(xModel.getChildNode("ContactCylinder",cylNum).getChildNode("Radius").getText());
   }
   return 0.0f;
}

float tledModel::GetContactCylLength(int cylNum)
{
   int nCyls = xModel.nChildNode("ContactCylinder");
   if (nCyls >= cylNum+1)
   {
      int nLen = xModel.getChildNode("ContactCylinder",cylNum).nChildNode("Length");
      if (nLen > 0)
	return (float)atof(xModel.getChildNode("ContactCylinder",cylNum).getChildNode("Length").getText());
   }
   return 0.0f;
}

vector<float> tledModel::GetContactCylOrigDisp(int cylNum)
{
   vector<float> disp;
   int nCyls = xModel.nChildNode("ContactCylinder");
   if (nCyls >= cylNum+1)
   {
      int nDisp = xModel.getChildNode("ContactCylinder",cylNum).nChildNode("OrigDisp");
      if (nDisp > 0)
         disp = GetXMLTextAsVector<float>(xModel.getChildNode("ContactCylinder",cylNum).getChildNode("OrigDisp"));
   }
   return disp;
}

float tledModel::GetContactCylRadChange(int cylNum)
{
   int nCyls = xModel.nChildNode("ContactCylinder");
   if (nCyls >= cylNum+1)
   {
      int nRadChng = xModel.getChildNode("ContactCylinder",cylNum).nChildNode("RadChange");
      if (nRadChng > 0) {
	return (float)atof(xModel.getChildNode("ContactCylinder",cylNum).getChildNode("RadChange").getText());
      }
   }
   return 0.0f;
}

vector<int> tledModel::GetContactCylSlvs(int cylNum)
{
   vector<int> slvs;
   int nCyls = xModel.nChildNode("ContactCylinder");
   if (nCyls >= cylNum+1)
   {
      int nSlvs = xModel.getChildNode("ContactCylinder",cylNum).nChildNode("SlvNodes");
      if (nSlvs > 0)
         slvs = GetXMLTextAsVector<int>(xModel.getChildNode("ContactCylinder",cylNum).getChildNode("SlvNodes"));
   }
   return slvs;
}

int tledModel::GetNumContactPrbs(void)
{
   return xModel.nChildNode("ContactUSProbe");
}

vector<float> tledModel::GetContactPrbOrigin(int prbNum)
{
   vector<float> orig;
   int nPrbs = xModel.nChildNode("ContactUSProbe");
   if (nPrbs >= prbNum+1)
   {
      int nOrig = xModel.getChildNode("ContactUSProbe",prbNum).nChildNode("Origin");
      if (nOrig > 0)
         orig = GetXMLTextAsVector<float>(xModel.getChildNode("ContactUSProbe",prbNum).getChildNode("Origin"));
   }
   return orig;
}

vector<float> tledModel::GetContactPrbAxis(int prbNum)
{
   vector<float> axis;
   int nPrbs = xModel.nChildNode("ContactUSProbe");
   if (nPrbs >= prbNum+1)
   {
      int nAxis = xModel.getChildNode("ContactUSProbe",prbNum).nChildNode("Axis");
      if (nAxis > 0)
         axis = GetXMLTextAsVector<float>(xModel.getChildNode("ContactUSProbe",prbNum).getChildNode("Axis"));
   }
   return axis;
}

float tledModel::GetContactPrbRadius(int prbNum)
{
   int nPrbs = xModel.nChildNode("ContactUSProbe");
   if (nPrbs >= prbNum+1)
   {
      int nRad = xModel.getChildNode("ContactUSProbe",prbNum).nChildNode("Radius");
      if (nRad > 0) {
	return (float)atof(xModel.getChildNode("ContactUSProbe",prbNum).getChildNode("Radius").getText());
      }
   }
   return 0.0f;
}

float tledModel::GetContactPrbLength(int prbNum)
{
   int nPrbs = xModel.nChildNode("ContactUSProbe");
   if (nPrbs >= prbNum+1)
   {
      int nLen = xModel.getChildNode("ContactUSProbe",prbNum).nChildNode("Length");
      if (nLen > 0) {
	return (float)atof(xModel.getChildNode("ContactUSProbe",prbNum).getChildNode("Length").getText());
      }
   }
   return 0.0f;
}

vector<float> tledModel::GetContactPrbOrigDisp(int prbNum)
{
   vector<float> disp;
   int nPrbs = xModel.nChildNode("ContactUSProbe");
   if (nPrbs >= prbNum+1)
   {
      int nDisp = xModel.getChildNode("ContactUSProbe",prbNum).nChildNode("OrigDisp");
      if (nDisp > 0)
         disp = GetXMLTextAsVector<float>(xModel.getChildNode("ContactUSProbe",prbNum).getChildNode("OrigDisp"));
   }
   return disp;
}

float tledModel::GetContactPrbRadChange(int prbNum)
{
   int nPrbs = xModel.nChildNode("ContactUSProbe");
   if (nPrbs >= prbNum+1)
   {
      int nRadChng = xModel.getChildNode("ContactUSProbe",prbNum).nChildNode("RadChange");
      if (nRadChng > 0) {
	return (float)atof(xModel.getChildNode("ContactUSProbe",prbNum).getChildNode("RadChange").getText());
      }
   }
   return 0.0f;
}

vector<int> tledModel::GetContactPrbSlvs(int prbNum)
{
   vector<int> slvs;
   int nPrbs = xModel.nChildNode("ContactUSProbe");
   if (nPrbs >= prbNum+1)
   {
      int nSlvs = xModel.getChildNode("ContactUSProbe",prbNum).nChildNode("SlvNodes");
      if (nSlvs > 0)
         slvs = GetXMLTextAsVector<int>(xModel.getChildNode("ContactUSProbe",prbNum).getChildNode("SlvNodes"));
   }
   return slvs;
}

int tledModel::GetNumContactPlts(void)
{
   return xModel.nChildNode("ContactPlate");
}

vector<float> tledModel::GetContactPltCrnrA(int pltNum)
{
   vector<float> a;
   int nPlts = xModel.nChildNode("ContactPlate");
   if (nPlts >= pltNum+1)
   {
      int na = xModel.getChildNode("ContactPlate",pltNum).nChildNode("a");
      if (na > 0)
         a = GetXMLTextAsVector<float>(xModel.getChildNode("ContactPlate",pltNum).getChildNode("a"));
   }
   return a;
}

vector<float> tledModel::GetContactPltCrnrB(int pltNum)
{
   vector<float> b;
   int nPlts = xModel.nChildNode("ContactPlate");
   if (nPlts >= pltNum+1)
   {
      int nb = xModel.getChildNode("ContactPlate",pltNum).nChildNode("b");
      if (nb > 0)
         b = GetXMLTextAsVector<float>(xModel.getChildNode("ContactPlate",pltNum).getChildNode("b"));
   }
   return b;
}

vector<float> tledModel::GetContactPltCrnrC(int pltNum)
{
   vector<float> c;
   int nPlts = xModel.nChildNode("ContactPlate");
   if (nPlts >= pltNum+1)
   {
      int nc = xModel.getChildNode("ContactPlate",pltNum).nChildNode("c");
      if (nc > 0)
         c = GetXMLTextAsVector<float>(xModel.getChildNode("ContactPlate",pltNum).getChildNode("c"));
   }
   return c;
}

vector<float> tledModel::GetContactPltDisp(int pltNum)
{
   vector<float> disp;
   int nPlts = xModel.nChildNode("ContactPlate");
   if (nPlts >= pltNum+1)
   {
      int nDisp = xModel.getChildNode("ContactPlate",pltNum).nChildNode("Disp");
      if (nDisp > 0)
         disp = GetXMLTextAsVector<float>(xModel.getChildNode("ContactPlate",pltNum).getChildNode("Disp"));
   }
   return disp;
}

vector<int> tledModel::GetContactPltSlvs(int pltNum)
{
   vector<int> slvs;
   int nPlts = xModel.nChildNode("ContactPlate");
   if (nPlts >= pltNum+1)
   {
      int nSlvs = xModel.getChildNode("ContactPlate",pltNum).nChildNode("SlvNodes");
      if (nSlvs > 0)
         slvs = GetXMLTextAsVector<int>(xModel.getChildNode("ContactPlate",pltNum).getChildNode("SlvNodes"));
   }
   return slvs;
}

int tledModel::GetNumKinematicContacts(void)
{
   return xModel.nChildNode("KinematicContact");
}

vector<float> tledModel::GetKinematicContactPoints(int kinNum)
{
   vector<float> points;
   int nKin = xModel.nChildNode("KinematicContact");
   if (nKin >= kinNum+1)
   {
      int nPnts = xModel.getChildNode("KinematicContact",kinNum).nChildNode("Points");
      if (nPnts > 0)
         points = GetXMLTextAsVector<float>(xModel.getChildNode("KinematicContact",kinNum).getChildNode("Points"));
   }
   return points;
}

vector<int> tledModel::GetKinematicContactFaces(int kinNum)
{
   vector<int> faces;
   int nKin = xModel.nChildNode("KinematicContact");
   if (nKin >= kinNum+1)
   {
      int nFcs = xModel.getChildNode("KinematicContact",kinNum).nChildNode("Faces");
      if (nFcs > 0)
         faces = GetXMLTextAsVector<int>(xModel.getChildNode("KinematicContact",kinNum).getChildNode("Faces"));
   }
   return faces;
}

vector<int> tledModel::GetKinematicContactSlvs(int kinNum)
{
   vector<int> slvs;
   int nKin = xModel.nChildNode("KinematicContact");
   if (nKin >= kinNum+1)
   {
      int nSlvs = xModel.getChildNode("KinematicContact",kinNum).nChildNode("SlaveNodes");
      if (nSlvs > 0)
         slvs = GetXMLTextAsVector<int>(xModel.getChildNode("KinematicContact",kinNum).getChildNode("SlaveNodes"));
   }
   return slvs;
}

int tledModel::GetNumOutputVars(void)
{
   int nOutputElements = xModel.nChildNode("Output");
   if (nOutputElements == 0)
      return 0;
   else
      return xModel.getChildNode("Output").nChildNode("Variable");
}

int tledModel::GetOutputFreq(void)
{
   int nOutputElements = xModel.nChildNode("Output");
   if (nOutputElements == 0)
      return 0;
   else
      return atoi(xModel.getChildNode("Output").getAttribute("Freq"));
}

const char* tledModel::GetOutputVar(int varNum)
{
   int nOutputElements = xModel.nChildNode("Output");
   if (nOutputElements == 0)
      return "";
   else
      return xModel.getChildNode("Output").getChildNode("Variable",varNum).getText();
}

void tledModel::WriteModel(char* outFName)
{
   xModel.writeToFile(outFName);
}

bool tledModel::GetROM(void)
{
   int nrom = xModel.getChildNode("SystemParams").nChildNode("ROM");
   if (nrom == 0)
      return false;
   else
   {
      int rom = atoi(xModel.getChildNode("SystemParams").getChildNode("ROM").getText());
      if (rom == 0)
         return false;
      else
         return true;
   }
}

int tledModel::GetNumBasisVectors(void)
{
   int nrom = xModel.getChildNode("SystemParams").nChildNode("ROM");
   if (nrom == 0)
      return 0;
   else
      return atoi(xModel.getChildNode("ReducedBasis").getAttribute("NumVectors"));
}

float* tledModel::GetReducedBasis(void)
{
   int nrom = xModel.getChildNode("SystemParams").nChildNode("ROM");
   if (nrom == 0)
   {
      float* ret=NULL;
      return ret;
   }
   else
   {
      int numBasisVecs = GetNumBasisVectors();
      int numNodes = GetNumNodes();
      XMLNode xNode = xModel.getChildNode("ReducedBasis");
      bool isFileSet = (xNode.isAttributeSet("File") != 0);
      if (isFileSet == true)
      {
         // Get reduced basis from external file
         const char* fileName = xNode.getAttribute("File");
         return GetExternalDataAsArray<float>(fileName,numBasisVecs*numNodes*3);
      }
      else
         return GetXMLTextAsArray<float>(xNode,numBasisVecs*numNodes*3);
   }
}

/***************************************** Shell Elements *****************************************/
#include "tledMembraneMaterialLinear.h"
#include "tledMembraneMaterialNeoHookean.h"

tledSurface* tledModel::GetGenericShellMesh(void) const {
  const std::string surfaceType = GetShellMeshType();
  const std::string meshType = GetElType();

  if (surfaceType == "T3" || (surfaceType == "SURFACE" && (meshType == "T4" || meshType == "T4ANP"))) return GetShellMesh<3>();
  else if (surfaceType == "Q4" || (surfaceType == "SURFACE" && meshType == "H8")) return GetShellMesh<4>();
  else {
    tledLogErrorStream(tledHelper::FatalError() << "Shell element type " << surfaceType << " not recognised or supported.");
    return NULL;
  }
}

int tledModel::GetNumberOfShellElementSets() const {
   return xModel.nChildNode("ShellElementSet");
}

int tledModel::GetShellElementSetSize(const int ElSetNum) const {
  istringstream iss(xModel.getChildNode("ShellElementSet").getAttribute("Size"));
  int size;

  iss >> size;
  if (iss.fail() && iss.str() == "all") return -1;
  else return size;
}

vector<int> tledModel::GetShellElementSet(int ElSetNum) const {
   XMLNode xNode = xModel.getChildNode("ShellElementSet", ElSetNum);
   int sz = GetShellElementSetSize(ElSetNum);

   if (sz == -1) return tledSequenceGenerator::MakeSequence(-1, 0);
   else {
     vector<int> data = GetXMLTextAsVector<int>(xNode);

     if ((int)data.size() == sz) return data;
     else if (data.size() == 1)  return tledSequenceGenerator::MakeSequence(data[0], data[0] + sz);
     else {
       tledLogErrorStream(tledHelper::FatalError() << "Element set " << ElSetNum << " has invalid format.");
       return data;
     }
   }
}

std::vector<int> tledModel::GetClampedShellNodes(int elSetNum) const {
  std::vector<int> clampedNodeInds;

  if (xModel.nChildNode("ShellElementSet") > elSetNum) {
    const XMLNode esNode = xModel.getChildNode("ShellElementSet", elSetNum);
    
    for (int n = 0; n < esNode.nChildNode("ClampedNodes"); n++) {
      const XMLNode xNode = esNode.getChildNode("ClampedNodes", n);

      std::vector<int> xInds = GetXMLTextAsVector<int>(xNode);

      if (xInds.size() == 1) {
	std::istringstream iss(xNode.getAttribute("NumNodes"));
	int sz = *std::istream_iterator<int>(iss);

	if (iss.fail() || sz <= 0) {
	  tledLogErrorStream(tledHelper::FatalError() << "NumNodes attribute of ClampedNodes tag invalid; expected strictly positive numerical value, got " << xNode.getAttribute("NumNodes"));
	} else xInds = tledSequenceGenerator::MakeSequence(xInds[0], xInds[0] + sz);
      }

      clampedNodeInds.insert(clampedNodeInds.end(), xInds.begin(), xInds.end());
    }

    clampedNodeInds = tledHelper::MakeSortedUnique(clampedNodeInds);
  } else {
    tledLogErrorStream(tledHelper::Warning() << elSetNum << " exceeds number of available shell element sets.");
  }

  return clampedNodeInds;
}

tledShellMaterial* tledModel::GetShellMaterial(const int ElSetNum) const {
  if (xModel.nChildNode("ShellElementSet") <= ElSetNum) {
    goto __shellMaterialFail;
  } else {
    const XMLNode matNode = xModel.getChildNode("ShellElementSet", ElSetNum).getChildNode("Material");

    tledShellMaterial *p_mat;
  
    if (string("Linear") == matNode.getAttribute("Type")) {
      p_mat = new tledMembraneMaterialLinear();
    } else if (string("NeoHookean") == matNode.getAttribute("Type")) {
      p_mat = new tledMembraneMaterialNeoHookean();
    } else if (string("LinearShell") == matNode.getAttribute("Type")) {
      p_mat = new tledShellMaterialLinearPlateDecorator<tledMembraneMaterialLinear>();      
    } else if (string("LinearThickShell") == matNode.getAttribute("Type")) {
      p_mat = new tledShellMaterialLinearThickPlateDecorator<tledMembraneMaterialLinear>();      
    } else {
      tledLogErrorStream(tledHelper::FatalError() << "Material type " << matNode.getAttribute("Type") << " not supported for shells.");
    }

    if (matNode.nChildNode("Density") > 0 && matNode.nChildNode("Thickness") > 0) {
      istringstream iss;

      iss.str(matNode.getChildNode("Density").getText());
      p_mat->SetDensity(*istream_iterator<float>(iss));
      if (iss.fail()) goto  __shellMaterialFail;
      iss.clear();
      iss.str(matNode.getChildNode("Thickness").getText());
      p_mat->SetThickness(*istream_iterator<float>(iss));
      if (iss.fail()) goto  __shellMaterialFail;
      p_mat->SetParameters(GetXMLTextAsVector<float>(matNode));

      return p_mat;
    }
  }

 __shellMaterialFail:
  tledLogErrorStream(tledHelper::FatalError() << "No/Invalid material specification given for shell element set " << ElSetNum);

  return NULL;
}

bool tledModel::DoShellUseMeshSurface(void) const {  
  if (!this->HasSubModels()) {
    if (xModel.nChildNode("ShellElements") == 0) return false;
    return string(xModel.getChildNode("ShellElements").getAttribute("Type")) == "SURFACE";
  } else {
    return this->GetSubModelManager().DoShellUseMeshSurface();
  }
}

const char* tledModel::GetShellMeshType(void) const {
  if (xModel.nChildNode("ShellElements") > 0) {
    return xModel.getChildNode("ShellElements").getAttribute("Type");
  } else return NULL;
}

/***************************************** XML Import/Export of Geometrical Data *****************************************/
static const char g_GeometryExportTag[] = "RestartedSimulation";

void tledModel::SetPrecomputedGeometryFilename(const char path[]) {
  XMLResults pResults;

  m_GeomPrecompXML = XMLNode::parseFile(path, g_GeometryExportTag, &pResults);
  if (pResults.error != eXMLErrorNone) {
    tledLogErrorStream(tledHelper::FatalError() << "Parsing " << path << " @ line " << pResults.nLine << ", column " << pResults.nColumn << ":\n"
		       << XMLNode::getError(pResults.error));
  }  
}

bool tledModel::HasPrecomputedData(const char tag[]) const {
  if (m_GeomPrecompXML.isEmpty()) return false;
  else {
    assert(m_GeomPrecompXML.nChildNode(tag) == 0 || m_GeomPrecompXML.nChildNode(tag) == 1);
    return m_GeomPrecompXML.nChildNode(tag) > 0;
  }
}

XMLNode tledModel::GetPrecomputedData(const char tag[]) const {
  return m_GeomPrecompXML.getChildNode(tag);
}

void tledModel::StartXMLExport(void) {
  if (!m_GeomPrecompXML.isEmpty()) m_GeomPrecompXML.deleteNodeContent();
  m_GeomPrecompXML = XMLNode::createXMLTopNode(g_GeometryExportTag);
}

bool tledModel::WritePrecomputableDataXML(const char path[]) const {
  return m_GeomPrecompXML.writeToFile(path) == eXMLErrorNone;  
}
