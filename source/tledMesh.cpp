// =========================================================================
// File:       tledMesh.cpp
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


#include "tledMesh.h"
#include "tledModel.h"
#include "tledTimer.h"
#include "tledMatrixFunctions.h"
#include "tledMSHMeshLoader.h"
#include "tledVTKMeshLoader.h"
#include "tledVectorArithmetic.h"

#include <fstream>
#include <sstream>
#include <cstdlib>
#include <float.h>
#include <cmath>

using namespace std;

tledMesh::tledMesh(void)
{
  NCds = NULL;
  EInd = NULL;
}

tledMesh::tledMesh(XMLNode* xModel)
{
   XMLNode xNode;

   if (xModel->nChildNode("VTKMesh") || xModel->nChildNode("MSHMesh")) {
     tledMeshLoader *p_loader;

     if (xModel->nChildNode("VTKMesh")) {
       p_loader = new tledVTKMeshLoader();
       xNode = xModel->getChildNode("VTKMesh");
     } else {
       p_loader = new tledMSHMeshLoader();
       xNode = xModel->getChildNode("MSHMesh");
     }

     p_loader->SetMeshType(xNode.getAttribute("Type"));
     p_loader->SetFilename(xNode.getText());
     p_loader->SetOutputMesh(*this);

     if (xNode.nChildNode("Translation")) {
       float *p_trans = GetXMLTextAsArray<float>(xNode.getChildNode("Translation"), 3);

       p_loader->SetTranslation(p_trans);
       delete[] p_trans;
     }

     if (xNode.nChildNode("ScaleFactor")) {
       istringstream iss(xNode.getChildNode("ScaleFactor").getText());
       float s;

       iss >> s;
       if (iss.fail()) {
	 tledLogErrorStream(tledHelper::FatalError() << "Content of ScaleFactor tag is invalid, could not parse: " << xNode.getChildNode("ScaleFactor").getText());
       }
       p_loader->SetScaleFactor(s);
     }

     if (xNode.nChildNode("Rotation")) {
       float *p_rot;

      p_rot = GetXMLTextAsArray<float>(xNode.getChildNode("Rotation"), 6);
      p_loader->SetRotations(p_rot, p_rot[3], p_rot[4], p_rot[5]);

      delete[] p_rot;
     }

     p_loader->Read();
     EType = xNode.getAttribute("Type");

     delete p_loader;
   } else {
     // Get node data
     xNode = xModel->getChildNode("Nodes");
     NumNodes = atoi(xNode.getAttribute("NumNodes"));
     m_NumDOFs = atoi(xNode.getAttribute("DOF"));
     NCds = GetXMLTextAsArray<float>(xModel->getChildNode("Nodes"),NumNodes*3);
     // Get element data
     xNode = xModel->getChildNode("Elements");
     EType = xNode.getAttribute("Type");
     if (!strcmp(EType,"T4")) NodesPerEl = 4;
     else if (!strcmp(EType,"H8")) NodesPerEl = 8;
     else NodesPerEl = 4;	// T4ANP
     NumEls = atoi(xNode.getAttribute("NumEls"));
     if (NumEls > 0) EInd = GetXMLTextAsArray<int>(xModel->getChildNode("Elements"),NumEls*NodesPerEl);
     else EInd = NULL;
   }
}

tledMesh::~tledMesh(void)
{
  if (EInd != NULL) delete[] EInd;
  if (NCds != NULL) delete[] NCds;
}

float* tledMesh::ComputeCentroid(float *p_cent, const int elementIndex) const {
  using namespace tledVectorArithmetic;

  const int *elDef = GetAllElNodeInds() + elementIndex*GetNodesPerEl();
  
  std::copy(GetAllNodeCds() + 3*elDef[0], GetAllNodeCds() + 3*elDef[0] + 3, p_cent);
  for (int const *pc_nInd = elDef + 1; pc_nInd < elDef + GetNodesPerEl(); pc_nInd++) Add(p_cent, p_cent, GetAllNodeCds() + 3*(*pc_nInd));

  return ScalarDiv(p_cent, (float)GetNodesPerEl());
}

vector<int> tledMesh::GetElNodeInds(int ElNum) const
{
   vector<int> Inds;
   for (int i = 0; i < NodesPerEl; i++)
      Inds.push_back(EInd[NodesPerEl*ElNum + i]);
   return Inds;
}

vector<float> tledMesh::GetNodeCds(int NodeNum) const
{
   vector<float> Cds;
   for (int i = 0; i < 3; i++)
      Cds.push_back(NCds[3*NodeNum + i]);
   return Cds;
}

void tledMesh::SetNodeCoordinates(vector<float> nodes)
{
   int nnodes = nodes.size()/3;
   if (nnodes != NumNodes)
   {
      cerr << "!!! Warning: new list of node coords must be same size as existing one" << endl;
      return;
   }
   
   for (int i = 0; i < NumNodes*3; i++)
      NCds[i] = nodes[i];
}

float tledMesh::GetMinElementCharacteristicLength(void)
{
   float Lmin = FLT_MAX;
   
   vector<int> ind;
   if (NodesPerEl == 8)
   {
      // Node index pairs defining each of the 12 edges
      int a[12] = {0,1,2,3,4,5,6,7,0,1,2,3};
      int b[12] = {1,2,3,0,5,6,7,4,4,5,6,7};
      vector<float> x;
      for (int el = 0; el < NumEls; el++)
      {
         ind = GetElNodeInds(el);
         x.resize(NodesPerEl*3,0);
         for (int i = 0; i < NodesPerEl; i++)
         {
            vector<float> wkCds = GetNodeCds(ind[i]);
            for (int j = 0; j < 3; j++)
               x[3*i+j] = wkCds[j];
         }
         for (int edge = 0; edge < 12; edge++)
         {
            float L = sqrt((x[3*b[edge]] - x[3*a[edge]])*(x[3*b[edge]] - x[3*a[edge]]) +
                           (x[3*b[edge]+1] - x[3*a[edge]+1])*(x[3*b[edge]+1] - x[3*a[edge]+1]) +
                           (x[3*b[edge]+2] - x[3*a[edge]+2])*(x[3*b[edge]+2] - x[3*a[edge]+2]));
            Lmin = L < Lmin ? L : Lmin;
         }
      }
   }
   else
   {
      int A[4] = {0,0,0,1};
      int B[4] = {2,1,3,2};
      int C[4] = {1,3,2,3};
      float x[4][3];
      float tri[3][3];
      for (int el = 0; el < NumEls; el++)
      {
         ind = GetElNodeInds(el);
         for (int i = 0; i < NodesPerEl; i++)
         {
            vector<float> wkCds = GetNodeCds(ind[i]);
            for (int j = 0; j < 3; j++)
               x[i][j] = wkCds[j];
         }
         float V = ComputeTetVol(x);
         for (int face = 0; face < 4; face++)
         {
            for (int j = 0; j < 3; j++)
            {
               tri[0][j] = x[A[face]][j];
               tri[1][j] = x[B[face]][j];
               tri[2][j] = x[C[face]][j];
            }
            float area = ComputeTriArea(tri);
            float L = 3*V/area;
            Lmin = L < Lmin ? L : Lmin;
         }
      }
   }

   return Lmin;
}

float tledMesh::ComputeTetVol(float x[4][3])
{
   // Use formula from tet element classes
  float DhDr[4][3] = {{-1, -1, -1},	// Shape function natural derivatives
		      {1, 0, 0},
		      {0, 1, 0},
		      {0, 0, 1}};
   float J[3][3];
   MatMult43T43(DhDr,x,J);
   float detJ;
   MatDet33(J,&detJ);
   return fabs(detJ/6);
}

float tledMesh::ComputeTriArea(float x[3][3])
{
   float t1 = (x[0][0]*x[1][1] - x[0][1]*x[1][0] - x[0][0]*x[2][1] + x[0][1]*x[2][0] + x[1][0]*x[2][1] - x[1][1]*x[2][0])*
              (x[0][0]*x[1][1] - x[0][1]*x[1][0] - x[0][0]*x[2][1] + x[0][1]*x[2][0] + x[1][0]*x[2][1] - x[1][1]*x[2][0]);
   float t2 = (x[0][0]*x[1][2] - x[0][2]*x[1][0] - x[0][0]*x[2][2] + x[0][2]*x[2][0] + x[1][0]*x[2][2] - x[1][2]*x[2][0])*
              (x[0][0]*x[1][2] - x[0][2]*x[1][0] - x[0][0]*x[2][2] + x[0][2]*x[2][0] + x[1][0]*x[2][2] - x[1][2]*x[2][0]);
   float t3 = (x[0][1]*x[1][2] - x[0][2]*x[1][1] - x[0][1]*x[2][2] + x[0][2]*x[2][1] + x[1][1]*x[2][2] - x[1][2]*x[2][1])*
              (x[0][1]*x[1][2] - x[0][2]*x[1][1] - x[0][1]*x[2][2] + x[0][2]*x[2][1] + x[1][1]*x[2][2] - x[1][2]*x[2][1]);
   return sqrt(t1+t2+t3)/2;
}

void tledMesh::SetNumberOfElements(const int numEls, const char elType[]) {
  NodesPerEl = string(elType) == "T4" || string(elType) == "T4ANP"? 4 : 8;  
  EInd = new int[NodesPerEl*numEls];
  EType = elType;
  NumEls = numEls;
}

void tledMesh::SetNumberOfNodes(const int numNodes, const int numNodeDOFs) {
  m_NumDOFs = numNodeDOFs;
  NumNodes = numNodes;
  NCds = new float[numNodes*3];
}
