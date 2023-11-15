// =========================================================================
// File:       tledBasicMeshCreator.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledBasicMeshCreator_H
#define tledBasicMeshCreator_H

template <class TMesh>
class tledBasicMeshCreator {
  /**
   * \name Mesh 
   * @{
   */
private:
  TMesh *mp_Mesh;

public:
  /** Allocates an output buffer for the next read operation. Client is responsible for memory clean-up. */
  void AllocateOutputMesh(void) { mp_Mesh = new TMesh(); }

  /** Setter allowing the user to specify the memory location of the output mesh. Either this or AllocateOutputMesh must be called before a read. */
  void SetOutputMesh(TMesh &r_out) { mp_Mesh = &r_out; }

  const TMesh& GetOutput(void) const { return *mp_Mesh; }
  TMesh& GetOutput(void) { return *mp_Mesh; }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledBasicMeshCreator(void) {}
  /** @} */
};
#endif
