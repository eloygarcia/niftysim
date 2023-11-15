// =========================================================================
// File:       tledVTKMeshLoader.h
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
#ifndef tledVTKMeshLoader_H
#define tledVTKMeshLoader_H
#include "tledHelper.h"
#include "tledModel.h"
#include "tledMeshLoader.h"

#include <vector>
#include <string>

/**
 * \brief Allows for loading of FEM meshes from VTK files
 * \ingroup import
 */
class tledVTKMeshLoader : public tledMeshLoader {
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
  tledVTKMeshLoader(void);
  virtual ~tledVTKMeshLoader(void) {}
  /** @} */
};

#endif
