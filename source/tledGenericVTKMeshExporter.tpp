// =========================================================================
// File:       tledGenericVTKMeshExporter.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#ifdef _Visualisation_

template <class TVTKMesh>
struct _NiftySim_VTKMeshWriterAdapter {
  typedef void Writer;
};

template <>
struct _NiftySim_VTKMeshWriterAdapter<vtkUnstructuredGrid> {
  typedef vtkUnstructuredGridWriter Writer;
};

template <>
struct _NiftySim_VTKMeshWriterAdapter<vtkPolyData> {
  typedef vtkPolyDataWriter Writer;
};

template <class TMesh, class TVTKMesh> 
bool tledGenericVTKMeshExporter<TMesh, TVTKMesh>::Write() {
  typedef typename _NiftySim_VTKMeshWriterAdapter<TVTKMesh>::Writer __Writer;

  vtkSmartPointer<__Writer> sp_writer;

  sp_writer = vtkSmartPointer<__Writer>::New();
  sp_writer->SetFileName(this->GetFileName().c_str());
  tledVTK6CompatSetInput(sp_writer, m_Converter.GetOutput());
  sp_writer->Update();

  if (sp_writer->GetErrorCode() != 0) {
    tledLogErrorStream(tledHelper::Warning() << "Error writing to " + this->GetFileName() << ": " << vtkErrorCode::GetStringFromErrorCode(sp_writer->GetErrorCode()));
  }

  return sp_writer->GetErrorCode() == 0;
}

template <class TMesh, class TVTKMesh>
void tledGenericVTKMeshExporter<TMesh, TVTKMesh>::AddNodeVectorAttribute(const std::string &name, const float attribs[]) { 
  m_Converter.AddNodeVectorAttribute(name, attribs); 
}

template <class TMesh, class TVTKMesh>
void tledGenericVTKMeshExporter<TMesh, TVTKMesh>::AddNodeScalarAttribute(const std::string &name, const float attribs[]) { 
  m_Converter.AddNodeScalarAttribute(name, attribs); 
}

template <class TMesh, class TVTKMesh>
void tledGenericVTKMeshExporter<TMesh, TVTKMesh>::SetMesh(const Mesh &mesh) {  
  m_Converter.SetInput(mesh);
  if (this->DoNodeCompression()) m_Converter.RemoveUnreferencedNodes();
  m_Converter.Update();
}

#else
template <class TMesh, class TVTKMesh> void tledGenericVTKMeshExporter<TMesh, TVTKMesh>::SetMesh(const tledMesh &mesh) {}
template <class TMesh, class TVTKMesh> void tledGenericVTKMeshExporter<TMesh, TVTKMesh>::AddNodeVectorAttribute(const std::string &name, const float attribs[]) {}
template <class TMesh, class TVTKMesh> void tledGenericVTKMeshExporter<TMesh, TVTKMesh>::AddNodeScalarAttribute(const std::string &name, const float attribs[]) {}
template <class TMesh, class TVTKMesh> bool tledGenericVTKMeshExporter<TMesh, TVTKMesh>::Write() { return false; }
#endif
