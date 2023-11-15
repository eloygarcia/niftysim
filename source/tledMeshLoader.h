// =========================================================================
// File:       tledMeshLoader.h
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
#ifndef tledMeshLoader_H
#define tledMeshLoader_H

#include "tledMesh.h"
#include "tledBasicMeshFileReader.h"
#include "tledHelper.h"

#include <string>

/** 
 * \brief tledMesh Importer Interface
 * \ingroup import
 */
class tledMeshLoader : public tledBasicMeshFileReader<tledMesh> {
  /**
   * \name Mesh-Type Queries
   * @{
   */
private:
  std::string m_MeshType;
  int m_NumElementNodes, m_NumDOFs;

public:
  /** @{ */
  /** Mesh type string representation: T4, H8, etc. */
  const char* GetMeshType(void) const { return m_MeshType.c_str(); }
  void SetMeshType(const char meshType[]);
  /** @} */

  /** Nodes/element */
  int GetNumberOfElementNodes(void) const { return m_NumElementNodes; }

  int GetNumberOfDOFs(void) const { return m_NumDOFs; }
  void SetNumberOfDOFs(const int numDOFs) { m_NumDOFs = numDOFs; }
  /** @} */

  /**
   * \name Transforms
   * @{
   */
protected:
  virtual void ApplyTransforms(void);
  /** @} */

  /**
   * \name FS I/O
   * @{
   */
protected:
  /** 
   * Reads the mesh from file has to be implemented by all format-specific subclasses.<br>
   * Application of transforms, etc. happens after mesh loading, in the Read member function and does not have to be implemented in subclasses. 
  */
  virtual void ReadFile(void) = 0;

public:
  virtual void Read(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledMeshLoader(void) : m_NumDOFs(3) {}
  virtual ~tledMeshLoader(void) {}
  /** @} */
};
#endif
