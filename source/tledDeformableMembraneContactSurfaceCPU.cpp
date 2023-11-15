// =========================================================================
// File:       tledDeformableMembraneContactSurfaceCPU.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledDeformableMembraneContactSurfaceCPU.h"
#include "tledContactVolumeSurfaceExtractor.h"
#include "tledXMLDeformableMembraneContactSurfaceCreator.h"

#include <string>

tledDeformableMembraneContactSurfaceCPU* tledDeformableMembraneContactSurfaceCPU::CreateSurface(const std::string &type) {
  tledDeformableMembraneContactSurfaceCPU *p_surf = NULL;

  if (type == "T3") {
    p_surf = new tledDeformableMembraneContactSurfaceT3CPU();
  } else if (type == "Q4") {
    p_surf = new tledDeformableMembraneContactSurfaceQ4CPU();
  } else {
    tledLogErrorStream(tledHelper::FatalError() << "Type " << type << " not supported.");
  }

  return p_surf;
}

tledDeformableMembraneContactSurfaceCPU* tledDeformableMembraneContactSurfaceCPU::CreateSurface(const tledMesh &volumeMesh, const tledSurface &membraneMesh) {
  tledDeformableMembraneContactSurfaceCPU *p_surf = NULL;

  if (membraneMesh.GetNumberOfFacetVertices() == 3 && (std::string("T4") == volumeMesh.GetElType() || std::string("T4ANP") == volumeMesh.GetElType())) {
    p_surf = new tledDeformableMembraneContactSurfaceImplCPU<tledDeformableMembraneContactSurfaceImpl<3, tledDeformableMembraneContactSurfaceCPU> >(volumeMesh, membraneMesh);
  } else if (membraneMesh.GetNumberOfFacetVertices() == 4 && std::string("H8") == volumeMesh.GetElType()) {
    p_surf = new tledDeformableMembraneContactSurfaceImplCPU<tledDeformableMembraneContactSurfaceImpl<4, tledDeformableMembraneContactSurfaceCPU> >(volumeMesh, membraneMesh);
  } else {
    tledLogErrorStream(tledHelper::FatalError() << "Combination of solid mesh type " << volumeMesh.GetElType() << " + membrane mesh type " << (membraneMesh.GetNumberOfFacetVertices() == 3? "T3" : "Q4") << " not supported");
  }

  return p_surf;
}

template <class TSurface>
static TSurface* _LoadSurfaceFromXML(TSurface *p_surface, const XMLNode xmlRep) {
  tledXMLDeformableMembraneContactSurfaceCreator<TSurface> importer;

  importer.SetXMLRoot(xmlRep);
  importer.SetOutputMesh(*p_surface);
  importer.Create();

  return p_surface;
}

tledDeformableMembraneContactSurfaceCPU* tledDeformableMembraneContactSurfaceCPU::CreateSurface(const XMLNode xmlRep) {
  tledDeformableMembraneContactSurfaceCPU *p_surf = NULL;

  if (xmlRep.nChildNode("Elements") == 1) {
    const char *type = xmlRep.getChildNode("Elements").getAttribute("Type");

    if (type != NULL) {
      if (std::string("T3") == type) p_surf = _LoadSurfaceFromXML(new tledDeformableMembraneContactSurfaceT3CPU(), xmlRep);
      else if (std::string("Q4") == type) p_surf = _LoadSurfaceFromXML(new tledDeformableMembraneContactSurfaceQ4CPU(), xmlRep);
      else {
	tledLogErrorStream(tledHelper::FatalError() << "Element type \"" << type << "\" not recognised.");
      }
    } else {
      tledLogErrorStream(tledHelper::FatalError() << "Found no type attribute \"Elements\" node");
    }
  } else {
    tledLogErrorStream(tledHelper::FatalError() << "Found no \"Elements\" node");
  }

  return p_surf;
}
