// =========================================================================
// File:       tledBasicMeshFileWriter.h
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
#ifndef tledBasicMeshFileWriter_H
#define tledBasicMeshFileWriter_H

#include <string>

/**
 * \defgroup export Mesh Exporters
 * \ingroup fileio
 */

/**
 * \brief Basic mesh file writer API
 */
template <class TMesh>
class tledBasicMeshFileWriter {
  /**
   * \name I/O
   */
private:
  std::string m_FileName;

public:
  virtual void SetMesh(const TMesh &mesh) = 0;

  void SetFileName(const std::string &filename) { m_FileName = filename; }
  const std::string& GetFileName(void) const { return m_FileName; }

  virtual bool Write(void) = 0;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledBasicMeshFileWriter(void) {} 
  tledBasicMeshFileWriter(const std::string &filename) : m_FileName(filename) {} 
  tledBasicMeshFileWriter(void) {}
  /** @} */
};

#endif
