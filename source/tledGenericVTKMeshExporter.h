// =========================================================================
// File:       tledGenericVTKMeshExporter.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledGenericVTKMeshExporter_H
#define tledGenericVTKMeshExporter_H

#include "tledHelper.h"
#include "tledBasicMeshFileWriter.h"

#ifdef _Visualisation_
#include "tledGenericVTKMeshSource.h"

#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataWriter.h>
#include <vtkErrorCode.h>

/**
 * \brief Exporter for writing of VTK unstructured grids
 * \ingroup export
 *
 *
 * Allows for making of more sophisticated animations in Paraview.
 */
template <class TMesh, class TVTKMesh> 
class tledGenericVTKMeshExporter : public tledBasicMeshFileWriter<TMesh> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledBasicMeshFileWriter<TMesh> Superclass;
  typedef TMesh Mesh;
  typedef tledGenericVTKMeshSource<TMesh, TVTKMesh> Converter;
  /** @} */

  /**
   * \name Mesh 
   * @{
   */
private:  
  Converter m_Converter;

protected:
  Converter& GetConverter(void) { return m_Converter; }
  
public:
  virtual void SetMesh(const Mesh &mesh);
  /** @} */

  /**
   * \name Attributes
   * @{
   */
public:
  virtual void AddNodeVectorAttribute(const std::string &name, const float attribs[]);
  virtual void AddNodeScalarAttribute(const std::string &name, const float attribs[]);
  /** @} */

  /**
   * \name I/O
   */
private:
  bool m_DoNodeCompression;

public:
  /** Removes nodes that are not referenced in by any elements in the input mesh. Must be set before the mesh is inputted. */
  void SetDoNodeCompression(const bool doCompress) { m_DoNodeCompression = doCompress; }
  bool DoNodeCompression(void) const { return m_DoNodeCompression; }

  virtual bool Write(void);
  /** @} */

  /**
   * \name Constrcution, Destruction
   * @{
   */
public:
  tledGenericVTKMeshExporter(const std::string &filename) : Superclass(filename), m_DoNodeCompression(false) {} 
  tledGenericVTKMeshExporter(void) : m_DoNodeCompression(false) {}
  virtual ~tledGenericVTKMeshExporter(void) {}
  /** @} */
};

#include "tledGenericVTKMeshExporter.tpp"
#else 
/*
 * Dummy class used when VTK is not available
 */
template <class TMesh> 
class tledGenericVTKMeshExporter : public tledBasicMeshFileWriter<TMesh> {
public:
  typedef tledBasicMeshFileWriter<TMesh> Superclass;
  typedef TMesh Mesh;

public:
  virtual void SetMesh(const TMesh &mesh) {}
  virtual void AddNodeVectorAttribute(const std::string &name, const float attribs[]) {}
  virtual void AddNodeScalarAttribute(const std::string &name, const float attribs[]) {}
  virtual bool Write(void) { return false; }

public:
  tledGenericVTKMeshExporter(const std::string &filename) : Superclass(filename) {} 
  tledGenericVTKMeshExporter(void) {}
  virtual ~tledGenericVTKMeshExporter(void) {}
};
#endif
#endif
