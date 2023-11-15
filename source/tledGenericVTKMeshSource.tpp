// =========================================================================
// File:       tledGenericVTKMeshSource.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    February 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TMesh, class TVTKMesh> 
void tledGenericVTKMeshSource<TMesh, TVTKMesh>::Update() {
  this->SetOutputObject(vtkSmartPointer<TVTKMesh>::New());

  this->GetOutput()->SetPoints(vtkSmartPointer<vtkPoints>::New());
  assert(this->GetIn2OutNodeIndexMap().size() != 0 || this->GetOut2InNodeIndexMap().size() == 0);
  if (this->GetIn2OutNodeIndexMap().size() > 0) {
    for (std::vector<int>::const_iterator ic_p = this->GetOut2InNodeIndexMap().begin(); ic_p < this->GetOut2InNodeIndexMap().end(); ic_p++) {
      this->GetOutput()->GetPoints()->InsertNextPoint(this->GetPointCoordinates(*ic_p));
    }
  } else {
    for (int ptInd = 0; ptInd < this->GetNumberOfPoints(); ptInd++) {
      this->GetOutput()->GetPoints()->InsertNextPoint(this->GetPointCoordinates(ptInd));
    }
  }

  this->GetOutput()->Allocate();
  for (int elInd = 0; elInd < this->GetNumberOfElements(); elInd++) {
    vtkIdType vtkElNodeInds[8];

    if (this->GetIn2OutNodeIndexMap().size() > 0) {
      for (int const *pc_i = this->GetElementBegin(elInd); pc_i < this->GetElementEnd(elInd); pc_i++) {
	assert(this->GetIn2OutNodeIndexMap()[*pc_i] >= 0);
	vtkElNodeInds[pc_i-this->GetElementBegin(elInd)] = this->GetIn2OutNodeIndexMap()[*pc_i];
      }
    } else {
      std::copy(this->GetElementBegin(elInd), this->GetElementEnd(elInd), vtkElNodeInds);
    }
    this->GetOutput()->InsertNextCell(this->GetVTKElementType(), this->GetElementEnd(elInd) - this->GetElementBegin(elInd), vtkElNodeInds);
  }
  assert(this->GetOutput()->GetNumberOfCells() == this->GetNumberOfElements());
}

template <class TMesh, class TVTKMesh> 
void tledGenericVTKMeshSource<TMesh, TVTKMesh>::RemoveUnreferencedNodes() {
  m_Out2InNodeIndexMap.clear();
  m_Out2InNodeIndexMap.reserve(8*this->GetNumberOfElements());

  for (int elInd = 0; elInd < this->GetNumberOfElements(); elInd++) {
    std::copy(this->GetElementBegin(elInd), this->GetElementEnd(elInd), std::back_inserter(m_Out2InNodeIndexMap));
  }
  m_Out2InNodeIndexMap = tledHelper::MakeSortedUnique(m_Out2InNodeIndexMap);

  m_In2OutNodeIndexMap = std::vector<int>(this->GetNumberOfPoints(), -1);
  for (size_t i = 0; i < m_Out2InNodeIndexMap.size(); i++) {
    m_In2OutNodeIndexMap[m_Out2InNodeIndexMap[i]] = i;
  }
}

template <class TMesh, class TVTKMesh> 
void tledGenericVTKMeshSource<TMesh, TVTKMesh>::_AddVTKNodeAttribute(const std::string &name, const float attribs[], const int numComponents) {
  vtkSmartPointer<vtkFloatArray> sp_attribs;

  assert(this->msp_Output != NULL);
  sp_attribs = vtkSmartPointer<vtkFloatArray>::New();
  sp_attribs->SetName(name.c_str());
  sp_attribs->SetNumberOfComponents(numComponents);

  if (this->GetIn2OutNodeIndexMap().size() > 0) {
    for (std::vector<int>::const_iterator ic_p = this->GetOut2InNodeIndexMap().begin(); ic_p < this->GetOut2InNodeIndexMap().end(); ic_p++) {
      // sp_attribs->InsertNextTupleValue(attribs + numComponents*(*ic_p));
      sp_attribs->InsertNextTuple(attribs + numComponents*(*ic_p));
      // Considering https://github.com/PointCloudLibrary/pcl/issues/2060
      // sp_attribs->InsertNextTypedTuple(attribs + numComponents*(*ic_p));
    }
  } else {
    for (float const *pc_nodeAttrib = attribs; pc_nodeAttrib < attribs + numComponents*this->GetOutput()->GetNumberOfPoints(); pc_nodeAttrib += numComponents) sp_attribs->InsertNextTuple(pc_nodeAttrib);
    // sp_attribs->InsertNextTypedTuple(pc_nodeAttrib);
    //sp_attribs->InsertNextTupleValue(pc_nodeAttrib);
  }
  this->GetOutput()->GetPointData()->AddArray(sp_attribs);

  assert(this->msp_Output->GetPointData()->HasArray(name.c_str()));
}

template <class TMesh, class TVTKMesh> 
void tledGenericVTKMeshSource<TMesh, TVTKMesh>::AddNodeScalarAttribute(const std::string &name, const float attribs[]) {
  _AddVTKNodeAttribute(name, attribs, 1);
}

template <class TMesh, class TVTKMesh> 
void tledGenericVTKMeshSource<TMesh, TVTKMesh>::AddNodeVectorAttribute(const std::string &name, const float attribs[]) {
  _AddVTKNodeAttribute(name, attribs, 3);
}
