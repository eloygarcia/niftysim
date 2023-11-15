// =========================================================================
// File:       tledContactSurfaceExporter.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    May 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

void tledContactSurfaceExporter<t_numFacetVertices>::ResetContactForces() {
  if (m_ContactForceAccumulator.size() == 0) {
    m_ContactForceAccumulator.insert(m_ContactForceAccumulator.end(), mpc_Surface->GetNumberOfNodes()*3, 0.0f);
  } else std::fill(m_ContactForceAccumulator.begin(), m_ContactForceAccumulator.end(), 0.0f);
  m_NumContactForceSamples = 0;
}

void tledContactSurfaceExporter<t_numFacetVertices>::AddContactForceSample(const std::vector<float> &contactForces) {
  MatAdd(&m_ContactForceAccumulator.front(), &contactForces.front(), mpc_Surface->GetNumberOfNodes()*3, 1, &m_ContactForceAccumulator.front());
  m_NumContactForceSamples += 1;
}

void tledContactSurfaceExporter<t_numFacetVertices>::Write(const float t) {
  vtkSmartPointer<vtkPolyDataWriter> sp_writer = vtkSmartPointer<vtkPolyDataWriter>::New();
  tledGenericVTKMeshSource<Surface, vtkPolyData> converter;
  std::vector<float> scaledForces(m_ContactForceAccumulator.size());
  std::ostringstream oss;

  converter.SetInput(*mpc_Surface);  
  converter.AddNodeVectorAttribute("contact_forces", MatMultScalar(&m_ContactForceAccumulator.front(), m_ContactForceAccumulator.size(), 1, 1.0f/m_NumContactForceSamples, &scaledForces.front()));
  converter.Update();

  oss << m_FilePrefix << "_t=" << t << ".vtk";
  sp_writer->SetFileName(oss.str().c_str());
  sp_writer->SetInput(converter.GetOutput());
  sp_writer->Update();
}
