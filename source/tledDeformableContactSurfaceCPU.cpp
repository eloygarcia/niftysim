// =========================================================================
// File:       tledDeformableContactSurfaceCPU.cpp
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
#include "tledXMLDeformableContactSurfaceCreator.h"
#include "tledDeformableContactSurfaceCPU.h"
#include "tledHelper.h"

#include <string>

tledDeformableContactSurfaceCPU* tledDeformableContactSurfaceCPU::CreateSurface(const std::string &type) {
  if (type == "T3") {
    return new tledDeformableContactSurfaceT3CPU();
  } else if (type == "Q4") {
    return new tledDeformableContactSurfaceQ4CPU();
  } else {
    tledLogErrorStream(tledHelper::FatalError() << "Mesh type " << type << " is not supported.");

    return NULL;
  }  
}

tledDeformableContactSurfaceCPU* tledDeformableContactSurfaceCPU::CreateSurface(const tledMesh &mesh) {
  tledDeformableContactSurfaceCPU *p_surf = NULL;

  if (std::string("T4") == mesh.GetElType() || std::string("T4ANP") == mesh.GetElType()) {
    p_surf = new tledDeformableContactSurfaceImplCPU<tledDeformableContactSurfaceImpl<3, tledDeformableContactSurfaceCPU> >(mesh);
  } else if (std::string("H8") == mesh.GetElType()) {
    p_surf = new tledDeformableContactSurfaceImplCPU<tledDeformableContactSurfaceImpl<4, tledDeformableContactSurfaceCPU> >(mesh);
  } else {
    tledLogErrorStream(tledHelper::FatalError() << "Could not determine surface type from solid-mesh type " << (mesh.GetElType() == NULL? "Unset" : mesh.GetElType()));
  }

  return p_surf;
}

template <class TSurface>
static TSurface* _LoadSurfaceFromXML(TSurface *p_surface, const XMLNode xmlRep) {
  tledXMLDeformableContactSurfaceCreator<TSurface> importer;

  importer.SetXMLRoot(xmlRep);
  importer.SetOutputMesh(*p_surface);
  importer.Create();

  return p_surface;
}

tledDeformableContactSurfaceCPU* tledDeformableContactSurfaceCPU::CreateSurface(const XMLNode xmlRep) {
  tledDeformableContactSurfaceCPU *p_surf = NULL;

  if (xmlRep.nChildNode("Elements") == 1) {
    const char *type = xmlRep.getChildNode("Elements").getAttribute("Type");

    if (type != NULL) {
      if (std::string("T3") == type) {
	p_surf = _LoadSurfaceFromXML(new tledDeformableContactSurfaceT3CPU(), xmlRep);
      } else if (std::string("Q4") == type) {
	p_surf = _LoadSurfaceFromXML(new tledDeformableContactSurfaceQ4CPU(), xmlRep);
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
