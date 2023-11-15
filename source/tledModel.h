// =========================================================================
// File:       tledModel.h
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
#ifndef tledModel_H
#define tledModel_H

#include "xmlParser.h"
#include "tledMesh.h"
#include "tledConstraint.h"
#include "tledMatrixFunctions.h"

#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>

/**
 * \defgroup model Simulation Setup
 */

template <class T> std::vector<T> GetXMLTextAsVector(XMLNode xNode);
template <class T> T* GetXMLTextAsArray(XMLNode xNode, int NumElements);
template <class T> T* GetExternalDataAsArray(const char* fileName, int NumElements);

class tledSubModelManager;
template <const int t_numFacetVertices>
class tledShellMesh;
class tledShellMaterial;
class tledSurface;

/**
 * \brief In-memory representation of simulation setup
 * \ingroup model
 */
class tledModel
{
protected:
  /** Returns a mandatory child node of parentNode, throws a fatal error if it's not found. */
  static XMLNode GetMandatoryXMLNode(const std::string &nodeName, const XMLNode &parentNode);

  /** Returns the SystemParams element of the model. */
  XMLNode GetSystemParamsNode(void) const;

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledModel(void) : mp_SubModelManager(NULL), Mesh(NULL), xFName(NULL) {;}
  tledModel(const XMLNode &rootNode);
  tledModel(const char* xmlFileName);
  ~tledModel(void);
  /** @} */

  /**
   * \name XML Description
   * @{
   */
protected:
  const XMLNode GetRootNode(void) const { return xModel; }
  void ResetError(void) { error = 0; }

public:
  /** Parse error status, non-zero if parsing failed. */
  int GetError(void) {return error;}
   
  /** XML file name */
   const char* GetFileName(void) {return xFName;}
  
  /** Parent directory of XML file, default directory for output. */
  std::string GetDirectory(void);

  /** Loads the model description from a parsed XML file */
  void LoadFromXML(const XMLNode &rootNode);
  /** @} */

  /**
   * \name Mesh data
   * @{
   */
public:
   int GetNumNodes(void) const {return Mesh->GetNumNodes();}
   int GetNumEls(void) const {return Mesh->GetNumEls();}
   std::vector<int> GetElNodeInds(int ElNum) const;
   std::vector<float> GetNodeCds(int NodeNum) const;
   int* GetAllElNodeInds(void) {return Mesh->GetAllElNodeInds();}
   float* GetAllNodeCds(void) {return Mesh->GetAllNodeCds();}
   const float* GetAllNodeCds(void) const {return Mesh->GetAllNodeCds();}
   int GetNodeDOF(void) const { return Mesh->GetNumberOfDOFs(); }
   const char* GetElType(void) const;
   int GetNodesPerEl(void);
   tledMesh* GetMesh(void) {return Mesh;}
   const tledMesh* GetMesh(void) const {return Mesh;}
  /** @} */
   
  /**
   * \name Element set data
   * @{
   */
public:
  static const int ELEMENT_SET_SIZE_ALL = -1;

protected:
  XMLNode GetElSetMaterialNode(const int elSetNum) const;

public:
   int GetNumElSets(void) const;
   std::vector<int> GetElSet(const int ElSetNum) const;

  /** Returns the number of elements in the set or ELEMENT_SET_SIZE_ALL if the size attribute was set to \"all\" */ 
   int GetElSetSize(const int ElSetNum) const;

  /** 
   * Returns the mass density associated with the material of the element set with index elSetInd. 
   * Returns the global mass density if none set for this element set. Fails if neither is available.
   */
  float GetElSetMassDensity(const int elSetNume) const;

   const char* GetMatType(const int ElSetNum) const;
   int GetNumElasticParams(const int ElSetNum) const;
   int GetMaxNumElasticParams(void) const;
   int GetNumViscIsoTerms(const int ElSetNum) const;
   int GetNumViscVolTerms(const int ElSetNum) const;
   int GetMaxNumViscIsoTerms(void) const;
   int GetMaxNumViscVolTerms(void) const;
   void GetElasticParams(const int ElSetNum, float* Params) const;	// Get only elastic params
   void GetViscoParams(const int ElSetNum, float* Params) const;	// Get only visco params
   void GetCombinedMatParams(const int ElSetNum, float* Params) const;	// Get both elastic and visco params, plus time step if applicable -> for tledSolverCPU
  /** @} */
   
  /**
   * \name System data
   * @{
   */
public:
   double GetTimeStep(void) const;
   double GetTotalTime(void) const;

   /** Default material density set in SystemParams, returns NaN if not set */
   float GetDensity(void) const;
   float GetDampingCoeff(void) const;
   float GetHGKappa(void) const;
   void SetTimeStep(double dt);
//   const bool GetMeshColouring(void);
  /** @} */

  /**
   * \name Constraint data
   * @{
   */
public:
  static const int CONSTRAINT_ALL_DOFS = -1;

private:
  const char* _GetSurfaceConstraintFaceType(const int ConstraintNum, const char constType[]) const;
  std::vector<int> _GetSurfaceConstraintFaceNodeInds(const int ConstraintNum, const char constType[]) const;  
  int _GetSurfaceConstraintNumFaces(const int ConstraintNum, const char constType[]) const;

public:
   int GetNumConstraints(const std::string &reqConstraintType) const;
   int GetNumConstraints(void) const;
   const char* GetConstraintType(const int ConstraintNum) const;
   int GetConstraintDOF(const int ConstraintNum) const;
   enum loadShape GetConstraintLoadShape(const int ConstraintNum) const;
   std::vector<int> GetConstraintInd(const int ConstraintNum) const;
   std::vector<float> GetConstraintMag(const int ConstraintNum) const;
   float GetGravityMag(const int ConstraintNum) const;
   std::vector<float> GetGravityDirec(const int ConstraintNum) const;

   static const int SURFACE_CONSTRAINT_FACE_SPEC_ALL = -1;
   static const int SURFACE_CONSTRAINT_FACE_SPEC_NORMAL = -2;

   /** 
    * Returns a number >= 0 if the number is not "all". In the latter case, SURFACE_CONSTRAINT_FACE_SPEC_ALL is returned and the facets are defined by extraction of the mesh surface (works with solid meshes only) 
    * and their number has to be determined by checking the size of the vector returned by tledModel::GetPressureFaceNodeInds.<br>
    * Similarly if the SpecType type attribute is set to NORMAL, SURFACE_CONSTRAINT_FACE_SPEC_NORMAL is returned, and again the number of facets must be determined through the size of the facet vector.
    */
  int GetPressureNumFaces(const int ConstraintNum) const { return _GetSurfaceConstraintNumFaces(ConstraintNum, "Pressure"); }

  const char* GetPressureFaceType(const int ConstraintNum) const { return _GetSurfaceConstraintFaceType(ConstraintNum, "Pressure"); }
  std::vector<int> GetPressureFaceNodeInds(const int ConstraintNum) const { return _GetSurfaceConstraintFaceNodeInds(ConstraintNum, "Pressure"); }
  float GetPressureMagnitude(const int ConstraintNum) const;

  /** 
   * Returns a number >= 0 if the number is not "all". In the latter case, SURFACE_CONSTRAINT_FACE_SPEC_ALL is returned and the facets are defined by extraction of the mesh surface (works with solid meshes only) 
   * and their number has to be determined by checking the size of the vector returned by tledModel::GetTractionFaceNodeInds.<br>
   * Similarly if the SpecType type attribute is set to NORMAL, SURFACE_CONSTRAINT_FACE_SPEC_NORMAL is returned, and again the number of facets must be determined through the size of the facet vector.
   */
  int GetTractionNumFaces(const int ConstraintNum) const { return _GetSurfaceConstraintNumFaces(ConstraintNum, "Traction"); }
  
  const char* GetTractionFaceType(const int ConstraintNum) const { return _GetSurfaceConstraintFaceType(ConstraintNum, "Traction"); }
  std::vector<int> GetTractionFaceNodeInds(const int ConstraintNum) const { return _GetSurfaceConstraintFaceNodeInds(ConstraintNum, "Traction"); }
  std::vector<float> GetTractionFaceTractions(const int ConstraintNum) const;

  void SetConstraintMag(const int ConstraintNum, const std::vector<float> &Mag);
  void SetConstraintInd(const int ConstraintNum, const std::vector<int> &Ind);
  void SetPressureMag(const int ConstraintNum, const float Mag);
  /** @} */
   
  /**
   * \name Contact modelling
   * @{
   */
   int GetNumContactObjects(void);
   int GetNumContactCyls(void); // Contact cylinders
   std::vector<float> GetContactCylOrigin(int cylNum);
   std::vector<float> GetContactCylAxis(int cylNum);
   float GetContactCylRadius(int cylNum);
   float GetContactCylLength(int cylNum);
   std::vector<float> GetContactCylOrigDisp(int cylNum);
   float GetContactCylRadChange(int cylNum);
   std::vector<int> GetContactCylSlvs(int cylNum);
   int GetNumContactPrbs(void); // Contact US probes
   std::vector<float> GetContactPrbOrigin(int prbNum);
   std::vector<float> GetContactPrbAxis(int cylNum);
   float GetContactPrbRadius(int prbNum);
   float GetContactPrbLength(int prbNum);
   std::vector<float> GetContactPrbOrigDisp(int prbNum);
   float GetContactPrbRadChange(int prbNum);
   std::vector<int> GetContactPrbSlvs(int prbNum);
   int GetNumContactPlts(void); // Contact plates
   std::vector<float> GetContactPltCrnrA(int pltNum);
   std::vector<float> GetContactPltCrnrB(int pltNum);
   std::vector<float> GetContactPltCrnrC(int pltNum);
   std::vector<float> GetContactPltDisp(int pltNum);
   std::vector<int> GetContactPltSlvs(int pltNum);
   
   int GetNumKinematicContacts(void); // Kinematic contacts
   std::vector<float> GetKinematicContactPoints(int kinNum);
   std::vector<int> GetKinematicContactFaces(int kinNum);
   std::vector<int> GetKinematicContactSlvs(int kinNum);

  /**
   * \brief Friction coefficient for deformable-deformable contacts.
   *
   * Returns NaN if contacts are to be frictionless.
   */
  float GetDeformableFrictionCoefficient(void) const;

  /** User has requested handling of collisions between deformables */
  bool DoDeformableDeformableContactHandling(void) const;

  /** Number of rigid contact surfaces (moving or unmoving) */
  int GetNumberOfRigidContactSurfaces(void) const;

  /** XML root node of a rigid contact surface definition */
  XMLNode GetRigidContactSurfaceDefinition(const int surfaceIndex) const;

  /** 
   * \brief User has requested handling of collisions between deformables and self-collisions.  
   *
   * Implies DoDeformableDeformableContactHandling.
   */
  bool DoSelfCollisionContactHandling(void) const;

  /** User has requested multi-pass contact handling */
  bool DoMultiPassContactHandling(void) const;

  /** 
   * Returns a user-requested bounding volume type as a string ("AABB", "OBB").
   * Defaults to AABB.
   */
  std::string GetBoundingVolumeType(void) const;

  /** 
   * \brief User-defined contact safety margin
   * 
   * Safety margin: distance at which geometry is assumed to be in contact and margin that should be satisfied after collisions have been resolved.<br />
   * Returns NaN if the user didn't specify such a margin (default defined by tledUnstructuredContactManager).
   */
  float GetContactSafetyMargin(void) const;

  /**
   * \brief User-defined distance at which rate constraints come into force.
   *
   * Returns NaN if the user didn't specify such a margin (default defined by tledUnstructuredContactManager).
   */
  float GetRateConstraintDistance(void) const;
  /** @} */
   
  /**
   * \name Output requests
   * @{
   */
   int GetNumOutputVars(void);
   int GetOutputFreq(void);
   const char* GetOutputVar(int varNum);
  /** @} */
   
   // Write model to XML file
   void WriteModel(char* outFName);
   
  /**
   * \name Reduced order modelling
   * @{
   */
   bool GetROM(void);
   int GetNumBasisVectors(void);
   float* GetReducedBasis(void);
  /** @} */

  /**
   * \name Geometry XML Import/Export
   *
   * Allows for precomputation of material parameter-independent data in repeated simulations
   * on the same geometry.
   * @{
   */
private:  
  XMLNode m_GeomPrecompXML;

public:
  /** Setter for path to precomputed geometry data */
  void SetPrecomputedGeometryFilename(const char path[]);
  
  /** Tells if there is precomputed data for a structure corresponding to the given XML tag */
  bool HasPrecomputedData(const char tag[]) const;

  /** Returns the available precomputed data for a given XML tag, check if there's any data available with HavePrecomputedData */
  XMLNode GetPrecomputedData(const char tag[]) const;

  /** Starts the XML export process */
  void StartXMLExport(void);

  /** To be used upstream to add data to the exported XML. */
  XMLNode& GetExportXMLRootNode(void) { return m_GeomPrecompXML; }

  /** 
   * \brief Exports precomputable material-independent information to an XML file with given path. 
   *
   * Last step in the process started with tledModel::StartXMLExport.
   * \return true on successful write, false on error
   */
  bool WritePrecomputableDataXML(const char path[]) const;
  /** @} */

  /**
   * \name Sub-Models
   * @{
   */
private:
  tledSubModelManager *mp_SubModelManager;

public:
  bool HasSubModels(void) const { return mp_SubModelManager != NULL; }
  tledSubModelManager& GetSubModelManager(void) { return *mp_SubModelManager; }
  const tledSubModelManager& GetSubModelManager(void) const { return *mp_SubModelManager; }
  /** @} */

  /**
   * \name Shells/Membranes
   * @{
   */
public:
  /** Returns a shell element mesh. The client is responsible for determining the right type (template parameter), and memory clean-up. */
  template <const int t_numFacetVertices>
  tledShellMesh<t_numFacetVertices>* GetShellMesh(void) const;

  /** Returns the shell mesh through generic surface pointer, wrapper around GetShellMesh */
  tledSurface* GetGenericShellMesh(void) const;

  /** 
   * \brief Returns the type of the shell element.
   *
   * Possible values:<br>   
   * <ul>
   * <li>T3:         Triangles</li>
   * <li>Q4:         Quads</li>
   * <li>SURFACE:    T3 if the main solid mesh's type is T4 or T4ANP, Q4 if it's H8</li>
   * </ul>
   */
  const char* GetShellMeshType(void) const;

  /**
   * \brief Query if shell mesh is explicitly specified or if the entire solid mesh is to be used (ShellElement type set to SURFACE).
   */
  bool DoShellUseMeshSurface(void) const;

  int GetNumberOfShellElementSets(void) const;
  std::vector<int> GetShellElementSet(const int ElSetNum) const;
  int GetShellElementSetSize(const int ElSetNum) const;

  /** Returns a dynamically allocated shell material object for a given element set. Responsibility for freeing the memeory lies with the client. */
  tledShellMaterial* GetShellMaterial(const int ElSetNum) const;

  /** Returns the indices of the clamped nodes in a particular shell element set */
  std::vector<int> GetClampedShellNodes(int elSetNum) const;
  /** @} */
  
private:
   XMLNode xModel;
   tledMesh* Mesh;
   const char* xFName;
   int error;
};

template <const int t_numFacetVertices>
tledShellMesh<t_numFacetVertices>* tledModel::GetShellMesh() const {
  if (DoShellUseMeshSurface()) return new tledShellMesh<t_numFacetVertices>(*GetMesh());
  else {
    tledShellMesh<t_numFacetVertices> *p_surface = new tledShellMesh<t_numFacetVertices>(xModel);

    p_surface->SetNodeVector(GetAllNodeCds(), GetNumNodes());

    return p_surface;
  }
}

template <class T> std::vector<T> GetXMLTextAsVector(XMLNode xNode)
{
  std::vector<T> vec;
   std::istringstream str(xNode.getText());
   T val;

   while (!(str.eof())) {
      str >> val;
      if (!str.fail()) vec.push_back(val);
      else break;
   }

   return vec;
}

template <class T> T* GetXMLTextAsArray(XMLNode xNode, int NumElements)
{
   T* arr;
   arr = new T[NumElements];
   std::stringstream str(xNode.getText());
   int count = 0;
   while (!(str.eof()))
   {
      if (count < NumElements)
         str >> arr[count++];
      else
      {
	std::cerr << "!!!\n!!! Too many data values in XMLNode " << xNode.getName() << "\n!!!" << std::endl;
	std::abort();
      }
   }
   return arr;
}

template <class T> T* GetExternalDataAsArray(const char* fileName, int NumElements)
{
  std::ifstream file(fileName);
  if (!file)
    {
      std::cerr << "Failed to open file " << fileName << std::endl;
      return NULL;
    }
  T* data = new T[NumElements];
  int i = 0;
  while ( !(file.eof()) && (i < NumElements) )
    file >> data[i++];
  if (i < NumElements-1)
    std::cerr << "!!!\n!!! Not enough data values in file " << fileName << "\n!!!" << std::endl;
  file.close();
  return data;
}


#endif // tledModel_H
