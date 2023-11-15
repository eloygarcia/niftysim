// =========================================================================
// File:       tledMesh.h
// Purpose:    Class for defining mesh characteristics
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    March 2007
// 
// Copyright (c) 2010, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================


#ifndef tledMesh_H
#define tledMesh_H

#include "xmlParser.h"

#include <iostream>
#include <string>
#include <vector>

/**
 * \defgroup mesh Mesh Representations
 * \brief NiftySim mesh class group.
 */


/**
 * \brief Solid-element mesh class
 * \ingroup mesh
 *
 * Supports T4 (= tetrahedral), T4ANP (= tetra w/ average nodal pressure), and H8 (= hexahedral) meshes.
 */
class tledMesh
{
public:
   tledMesh(void);
   tledMesh(XMLNode* xModel);
   ~tledMesh(void);

  void SetNumberOfElements(const int numEls, const char elType[]);
  void SetNumberOfNodes(const int numNodes, const int numNodeDOFs);

   std::vector<int> GetElNodeInds(int ElNum) const;
   std::vector<float> GetNodeCds(int NodeNum) const;
   int* GetAllElNodeInds() {return EInd;}
  const int* GetAllElNodeInds() const {return EInd;}
  float* GetAllNodeCds() {return NCds;}
  const float* GetAllNodeCds() const {return NCds;}
   int GetNumNodes() const {return NumNodes;}
   int GetNumEls() const {return NumEls;}
   int GetNodesPerEl() const {return NodesPerEl;}
   const char* GetElType() const {return EType;}
   void SetNodeCoordinates(std::vector<float> nodes);
   float GetMinElementCharacteristicLength(void);
   int GetNumberOfDOFs(void) const { return m_NumDOFs; }

  /** Element centroid */
  float* ComputeCentroid(float *p_cent, const int elementIndex) const;
   
private:
   float ComputeTetVol(float x[4][3]);
   float ComputeTriArea(float x[3][3]);
   
   const char* EType;	// Element type (T4,T4ANP,H8 currently supported)
   float* NCds;	// Node coordinates
   int* EInd;		// Element node indices
   int NumNodes;	// Number of nodes in the mesh
   int NumEls;		// Number of elements in the mesh
   int NodesPerEl; // Number of nodes per element
   int m_NumDOFs;
};

#endif // tledMesh_H
