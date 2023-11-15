// =========================================================================
// File:       tledNodeRejector.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    October 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledNodeRejector_H
#define tledNodeRejector_H

#include "tledHelper.h"
#include "xmlParser.h"

#include <vector>

/**
 * \brief Node set exclusion criterion
 * \ingroup model
 *
 * Excludes nodes based on a criterion implemented in the member funcion DoReject. This class is mainly used for geometrically defining boundary conditions.
 */
class tledNodeRejector {
  /**
   * \name Input
   * @{
   */
private:
  const float *mpc_Nodes;
  std::vector<int> *mp_NodeIndices;

protected:
  const float* GetNodes(void) const { return mpc_Nodes; }

  std::vector<int>& GetNodeIndices(void) { return *mp_NodeIndices; }
  const std::vector<int>& GetNodeIndices(void) const { return *mp_NodeIndices; }

public:
  void SetNodes(const float nodes[]) { mpc_Nodes = nodes; }
  void SetNodeIndices(std::vector<int> &r_data) { mp_NodeIndices = &r_data; }

  const float* GetNode(const int nodeIndex) const { return mpc_Nodes + 3*nodeIndex; }
  /** @} */

  /**
   * \name Processing
   * @{
   */
private:
  tledNodeRejector *mp_Next;

protected:
  bool HasNext(void) const {return mp_Next != NULL; }

  tledNodeRejector& GetNext(void) { return *mp_Next; }
  const tledNodeRejector& GetNext(void) const { return *mp_Next; }

  void ProceedWithNextRejector(void);

public:
  virtual bool DoReject(const int nodeIndex) const = 0;  
  virtual void RunRejection(void);
  /** @} */
  
  /**
   * \name Construction, Init, Destruction
   * @{
   */
protected:
  tledNodeRejector(void) : mp_Next(NULL) {}

public:
  /** 
   * Appends a sub-criterion to the chain of criteria to be evaluated. 
   * The object is assumed to be residing in dynamically allocated memory and is destroyed upon destruction of the object to which this function is applied.
   */
  void AppendRejectorToChain(tledNodeRejector *p_rejector);

  virtual void InitialiseFromXMLSpec(const XMLNode rootNode) = 0;

  /**
   * Creates a rejector from an XML description
   */
  static tledNodeRejector* CreateRejector(const XMLNode rootNode);    

  virtual ~tledNodeRejector(void);
  /** @} */
};

#endif
