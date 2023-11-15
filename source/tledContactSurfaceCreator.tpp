template <class TSurface>
void tledContactSurfaceCreator<TSurface>::Create() {
  mp_BaseCreator->SetOutputMesh(this->GetOutput());
  mp_BaseCreator->Create();
  this->InitEdges();
}

template <class TSurface>
void tledContactSurfaceCreator<TSurface>::InitEdges() {
  tledSurfaceTopology<TSurface> topo(this->GetOutput());

  topo.ComputeEdges();
  this->GetOutput().GetAllEdges() = topo.GetEdges();
  for (int fInd = 0; fInd < this->GetOutput().GetNumberOfFacets(); fInd++) {
    std::copy(topo.GetFacetEdges()[fInd].begin(), topo.GetFacetEdges()[fInd].end(), this->GetOutput().GetFacet(fInd).EdgeIndices);
  }
}
