// =========================================================================
// File:       tledDeformableContactSurfaceGPU.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    November 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledDeformableContactSurfaceGPU.h"
#include "tledXMLDeformableContactSurfaceCreator.h"

#include "tledDeformableContactSurfaceGPU_kernels.cu"

tledDeformableContactSurfaceGPU* tledDeformableContactSurfaceGPU::CreateSurface(const tledMesh &mesh) {
  if (mesh.GetElType() == std::string("T4") || mesh.GetElType() == std::string("T4ANP")) {
    return new tledDeformableContactSurfaceT3GPU(mesh);
  } else if (mesh.GetElType() == std::string("H8")) {
    return new tledDeformableContactSurfaceQ4GPU(mesh);
  } else {
    tledLogErrorStream(tledHelper::FatalError() << "Mesh type " << mesh.GetElType() << " is not supported.");

    return NULL;
  }  
}

tledDeformableContactSurfaceGPU* tledDeformableContactSurfaceGPU::CreateSurface(const std::string &type) {
  if (type == "T3") {
    return new tledDeformableContactSurfaceT3GPU();
  } else if (type == "Q4") {
    return new tledDeformableContactSurfaceQ4GPU();
  } else {
    tledLogErrorStream(tledHelper::FatalError() << "Mesh type " << type << " is not supported.");

    return NULL;
  }  
}

template <class TSurface>
static TSurface* _LoadSurfaceFromXML(TSurface *p_surface, const XMLNode xmlRep) {
  tledXMLDeformableContactSurfaceCreator<TSurface> importer;

  importer.SetXMLRoot(xmlRep);
  importer.SetOutputMesh(*p_surface);
  importer.Create();

  return p_surface;
}

tledDeformableContactSurfaceGPU* tledDeformableContactSurfaceGPU::CreateSurface(const XMLNode xmlRep) {
  tledDeformableContactSurfaceGPU *p_surf = NULL;

  if (xmlRep.nChildNode("Elements") == 1) {
    const char *type = xmlRep.getChildNode("Elements").getAttribute("Type");

    if (type != NULL) {
      if (std::string("T3") == type) {
	p_surf = _LoadSurfaceFromXML(new tledDeformableContactSurfaceT3GPU(), xmlRep);
      } else if (std::string("Q4") == type) {
	p_surf = _LoadSurfaceFromXML(new tledDeformableContactSurfaceQ4GPU(), xmlRep);
      } else {
	tledLogErrorStream(tledHelper::FatalError() << "Element type \"" << type << "\" not recognised.");      
      }
    } else {
      tledFatalError("Found no type attribute on \"Elements\" node.");
    }
  } else {
    tledFatalError("Found no \"Elements\" node.");
  }

  return p_surf; 
}

void tledDeformableContactSurfaceGPU::UpdateNodePositions(float3 *dp_x, const float3 *dpc_x0, const float4 *dpc_uNext, const int *dpc_surfaceToVolumeNodeIndexMap, const int numNodes) {
  const int blockSize = 256;
  const int numBlks = tledCUDAHelpers::GetNumberOfBlocks(numNodes, blockSize);

  tledDeformableContactSurfaceGPU_kernels::UpdateNodePositionsKernel <<<numBlks, blockSize>>> (dp_x, dpc_x0, dpc_uNext, dpc_surfaceToVolumeNodeIndexMap, numNodes);
}

