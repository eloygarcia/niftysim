// =========================================================================
// File:       tledSubModelManager.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2011
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledSubModelManager_H
#define tledSubModelManager_H

#include "tledSolverCPU.h"
#include "tledConstraintManager.h"
#include "tledShellMesh.h"
#include "xmlParser.h"

#include <string>
#include <vector>

/**
 * \brief Helper module for storing and assembling simulations from sub-model tags
 * \ingroup model
 */
class tledSubModelManager {
  /**
   * \name Mesh-Type Queries
   * @{
   */
private:
  std::string m_MeshType;
  int m_NumElementNodes;
  std::vector<XMLNode> m_SubXMLModels;
  std::vector<tledModel> m_SubModels;

public:
  const char* GetMeshType(void) const { return m_MeshType.c_str(); }
  /** @} */
  
  /**
   * \name Sub-Mesh Management
   * @{
   */
public:
  struct SubMesh {
    int ElementStartIndex, ElementEndIndex;
    int MembraneElementStartIndex, MembraneElementEndIndex;
    std::vector<int> SubMeshToMeshNodeIndexMap;
  };

private:
  std::vector<SubMesh> m_SubMeshes;

public:
  int GetNumberOfSubMeshes(void) const { return m_SubMeshes.size(); }
  int GetSubMeshElementStartIndex(const int subMeshInd) const { return m_SubMeshes[subMeshInd].ElementStartIndex; }
  int GetSubMeshElementEndIndex(const int subMeshInd) const { return m_SubMeshes[subMeshInd].ElementEndIndex; }
  const std::vector<int>& GetSubMeshToMeshNodeIndexMap(const int subMeshInd) const { return m_SubMeshes[subMeshInd].SubMeshToMeshNodeIndexMap; }
  const std::vector<SubMesh>& GetAllSubMeshes(void) const { return m_SubMeshes; }

  int GetNumberOfSubMeshElements(const int subMeshInd) const { return this->GetSubMeshElementEndIndex(subMeshInd) - this->GetSubMeshElementStartIndex(subMeshInd); }
  int GetNumberOfSubMeshNodes(const int subMeshInd) const { return GetSubMeshToMeshNodeIndexMap(subMeshInd).size(); }

  /** Extracts a submesh as a tledMesh (mainly for export purposes) */
  tledMesh* ExportSubMesh(const int subMeshInd, const tledMesh &globalMesh) const;

  /** 
   * \brief Extracts nodal attributes from a vector defined wrt. the global mesh. 
   *
   * Assumes 3 components/node.
   */
  void ExportSubMeshNodeAttributes(float *p_dst, const float attributes[], const int subMeshInd) const;
  /** @} */

  /**
   * \name Constraints
   * @{
   */
public:
  int GetNumberOfSubModelConstraints(const int smIndex) const { return m_SubXMLModels[smIndex].nChildNode("Constraint"); }

  /** Returns the global version of a constraint defined on sub-domain. */
  XMLNode GetSubModelConstraint(const int smIndex, const int constraintIndex) const;
  /** @} */

  /**
   * \name Shell/Membrane Elements
   * @{
   */
private:
  std::vector<int> m_MembraneElements;
  std::string m_ShellElementType;
  bool m_IsSurface;

private:
  /** \internal Only call after the appropriate node-index mapping has been computed in _AddSubMesh, but before the next sub-mesh is added. */
  void _AddSubMembrane(const tledSurface &membrane, const bool isSurface);

public:
  /** Builds an XML representation for the compound membrane mesh. */
  XMLNode CreateShellMeshXML(void) const;

  /** True iff any of the sub-models contained a membrane definition */
  bool HasMembrane(void) const { return m_MembraneElements.size() > 0; }

  /** Returns an XML representation of the sub-model shell element sets, updated with global mesh indices. */
  XMLNode GetSubModelShellElementSet(const int sMIndex, const int elSetIndex) const;

  int GetNumberOfSubModelShellElementSets(const int sMIndex) const { return m_SubXMLModels[sMIndex].nChildNode("ShellElementSet"); }

  bool DoShellUseMeshSurface(void) const { return m_IsSurface; }
  /** @} */

  /**
   * \name Model, Mesh I/O
   * @{
   */
private:
  std::vector<float> m_Nodes;
  std::vector<int> m_Elements;

private:
  void _AddSubMesh(const tledMesh &mesh);

public:
  /** Always call this before AddSubModel, as otherwise memory access issues may arise. */
  void SetNumberOfSubModels(const int numSubModels);

  /** Reads the sub-model specs from XML */
  void AddSubModel(const XMLNode &subModelRoot);

  void WriteMesh(tledMesh &r_mesh);
  /** @} */

  /**
   * \name Element Sets
   * @{
   */
public:
  /** Returns an XML representation of the sub-model element sets, updated with global mesh indices. */
  XMLNode GetSubModelElementSet(const int sMIndex, const int elSetIndex) const;
  int GetNumberOfSubModelElementSets(const int sMIndex) const { return m_SubXMLModels[sMIndex].nChildNode("ElementSet"); }
  /** @} */

  /**
   * \name Node Merging
   * @{
   */
private:
  float m_MinNodeDistance;

public:
  /** 
   * \brief Distance at which a node is considered one and the same. Default: 10^-4 
   *
   *
   * Set to -1 (or other negative value), if node merging is not desired.
   */
  float GetMinNodeDistance(void) const { return m_MinNodeDistance; }

  /** \brief Setter for distance at which nodes are merged. */
  void SetMinNodeDistance(const float minNodeDist) { m_MinNodeDistance = minNodeDist; }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledSubModelManager(void) : m_IsSurface(false) {}
  ~tledSubModelManager(void);
  /** @} */
};
#endif
