
template <const int t_NumElNodes>
tledMeshTopology<t_NumElNodes>::Facet::Facet(const int v0, const int v1, const int v2, const int v3) {
#if defined __GXX_EXPERIMENTAL_CXX0X__ || defined __WIN32
  static_assert(NumFacetVertices == 3 || NumFacetVertices == 4, "Only quadrangular/triangular facets are supported!");
#endif
  assert(NumFacetVertices == 4 || v3 == -1);

  (*this)[0] = v0;
  (*this)[1] = v1;
  (*this)[2] = v2;
  if (NumFacetVertices == 4) (*this)[3] = v3;
  
  if (NumFacetVertices == 3) std::sort(this->begin(), this->end());
  else {
    int sVal;

    std::rotate(this->begin(), std::min_element(this->begin(), this->end()), this->end());    
    sVal = *(this->begin() + 1);
    std::reverse(this->begin(), this->end());
    std::rotate(this->begin(), std::min_element(this->begin(), this->end()), this->end());    
    if (*(this->begin() + 1) > sVal) {
      std::reverse(this->begin(), this->end());
      std::rotate(this->begin(), std::min_element(this->begin(), this->end()), this->end());    
    }      
  }
} 

/**********************************************************************************************************/

template <const int t_NumElNodes>
void tledMeshTopology<t_NumElNodes>::ComputeNodeAdjacency() {
  const int *elNodeInds = mc_Mesh.GetAllElNodeInds();
  const int numEls = mc_Mesh.GetNumEls();
  const int numNodes = mc_Mesh.GetNumNodes();

  m_NodeEls.resize(numNodes);
  assert((int)m_NodeEls.size() == numNodes);
  for (int const *pc_elNodeInd = elNodeInds; pc_elNodeInd < elNodeInds + numEls*t_NumElNodes;) {
    for (int vidx = 0; vidx < t_NumElNodes; vidx++, pc_elNodeInd++) {
      assert((pc_elNodeInd - elNodeInds)/t_NumElNodes >= 0 && (pc_elNodeInd - elNodeInds)/t_NumElNodes < numEls);
      assert(*pc_elNodeInd >= 0 && *pc_elNodeInd < numNodes);
      m_NodeEls[*pc_elNodeInd].push_back((pc_elNodeInd - elNodeInds)/t_NumElNodes);
    }
  }
} 

template <const int t_NumElNodes>
std::vector<int> tledMeshTopology<t_NumElNodes>::GetSurfaceFacets() const {
  std::vector<int> surface;
  int facetInd;
  
  if (m_FacetEls.size() == 0) {
    std::cerr << "Need to call computeFaces before calling getSurfaceFaces!\n";
    abort();
  }

  surface.reserve(GetNumFacets());
  for (facetInd = 0; facetInd < GetNumFacets(); facetInd++) {
    const std::pair<int,int> &faceels = GetFacetElements(facetInd);

    assert(faceels.first >= 0);
    if (faceels.second == -1) surface.push_back(facetInd);
  }

  return surface;
} 

template <const int t_NumElNodes>
void tledMeshTopology<t_NumElNodes>::ComputeElementDiameters() {
  using namespace tledVectorArithmetic;

  const int *elNodeInds = mc_Mesh.GetAllElNodeInds();
  const int numEls = mc_Mesh.GetNumEls();
  const float *nodes = mc_Mesh.GetAllNodeCds();
  
  int eInd;
  
  m_MinH = std::numeric_limits<float>::max();
  m_MaxH = 0;
  for (eInd = 0; eInd < numEls; eInd++) {
    const int *currElNodes = elNodeInds + eInd*t_NumElNodes;

    int const *pc_v0Ind, *pc_v1Ind;
    float h, diff[3];

    for (pc_v0Ind = currElNodes; pc_v0Ind < currElNodes + t_NumElNodes - 1; pc_v0Ind++) for (pc_v1Ind = pc_v1Ind + 1; pc_v1Ind < currElNodes + t_NumElNodes; pc_v1Ind++) {
	h = Norm(Sub(diff, nodes + 3*(*pc_v0Ind), nodes + 3*(*pc_v1Ind)));
	m_MinH = std::min(m_MinH, h);
	m_MaxH = std::max(m_MaxH, h);
      }
  }
} 

template <>
void tledMeshTopology<4>::ComputeFacets();

template <>
void tledMeshTopology<8>::ComputeFacets();
