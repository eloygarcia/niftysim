// =========================================================================
// File:       tledMSHMeshLoader.h
// Purpose:    Loader for MSH (GMSH) files
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
#ifndef tledMSHMeshLoader_H
#define tledMSHMeshLoader_H
#include "tledMeshLoader.h"

/** 
 * \brief Reader for GMSH solid mesh files
 * \ingroup import
 */
class tledMSHMeshLoader : public tledMeshLoader {
  /**
   * \name tledMesh I/O
   * @{
   */
protected:
  virtual void ReadFile(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledMSHMeshLoader(void) {}
  /** @} */
};

#endif
