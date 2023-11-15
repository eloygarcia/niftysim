// =========================================================================
// File:       tledGenericVTKMeshSource.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    February 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledGenericVTKMeshSource_H
#define tledGenericVTKMeshSource_H

#include "tledHelper.h"

#ifdef _Visualisation_
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>

#include <string>
#include <algorithm>
#include <cassert>

/**
 * \brief Converts any TLED mesh (surfacic, volumetric) into a VTK mesh 
 */
template <class TMesh, class TVTKMesh> 
class tledGenericVTKMeshSource {
  /**
   * \name Mesh Constituent Access
   * @{
   */
protected:
  const float* GetPointCoordinates(const int ptInd) const;
  
  int GetNumberOfPoints(void) const;

  int GetNumberOfElements(void) const;

  /** \return Iterator pointing to start of element definition */
  const int* GetElementBegin(const int elementIndex) const;

  /** \return Iterator pointing to end of element definition */
  const int* GetElementEnd(const int elementIndex) const;

  int GetVTKElementType(void) const;
  /** @} */

  /**
   * \name Attributes
   * @{
   */
private:
  void _AddVTKNodeAttribute(const std::string &name, const float attribs[], const int numComponents);

public:
  virtual void AddNodeVectorAttribute(const std::string &name, const float attribs[]);
  virtual void AddNodeScalarAttribute(const std::string &name, const float attribs[]);
  /** @} */

  /**
   * \name Input
   * @{
   */
private:
  const TMesh *mpc_Mesh;
  std::vector<int> m_Out2InNodeIndexMap, m_In2OutNodeIndexMap;

protected:
  const std::vector<int>& GetOut2InNodeIndexMap(void) const { return m_Out2InNodeIndexMap; }
  const std::vector<int>& GetIn2OutNodeIndexMap(void) const { return m_In2OutNodeIndexMap; }
  
public:
  void SetInput(const TMesh &mesh) { mpc_Mesh = &mesh; }
  const TMesh& GetInput(void) const { return *mpc_Mesh; }

  /** Eliminates any nodes that are not part of some element. */
  void RemoveUnreferencedNodes(void);
  /** @} */

  /**
   * \name Standard VTK-Style Algorithm Interface
   * @{
   */
private:
  vtkSmartPointer<TVTKMesh> msp_Output;

protected:
  void SetOutputObject(vtkSmartPointer<TVTKMesh> sp_out) { msp_Output = sp_out; }

public:
  vtkSmartPointer<TVTKMesh> GetOutput(void) { return msp_Output; }
  void Update(void);
  /** @} */
};

#include "tledGenericVTKMeshSource.tpp"
#else
template <class TMesh, class TVTKMesh> 
class tledGenericVTKMeshSource {};
#endif
#endif
