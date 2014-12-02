from signedind import SignedInd as SI

def polygon_vert_out_labels(w):
  Lw = len(w)
  letter_inds = dict([(W,i) for (i,W) in enumerate(w)])
  ans = []
  done_vert = [False for i in xrange(Lw)]
  for i in xrange(Lw):
    if done_vert[i]:
      continue
    done_vert[i] = True
    vL = [w[i], i]
    j=i
    while True:
      next_letter = w[(j-1)%Lw].swapcase()
      if next_letter == w[i]:
        break
      next_letter_ind = letter_inds[next_letter]
      vL.extend( [ next_letter, next_letter_ind  ] )
      done_vert[next_letter_ind] = True
      j=next_letter_ind
    ans.append(vL)
  return ans

def polygon_edge_inds(w):
  letter_inds = dict([(W,i) for (i,W) in enumerate(w)])
  gens = sorted([L for L in letter_inds if L.islower()])
  gens_to_edges = dict([(g,SI(i,1)) for (i,g) in enumerate(gens)] + \
                       [(g.swapcase(),SI(i,-1)) for (i,g) in enumerate(gens)])
  return gens_to_edges




class Vertex :
  def __init__(self, iE, iF):
    self.i_edges = iE
    self.i_faces = iF
  def __repr__(self):
    return "Vertex(" + str(self.i_edges) + "," + str(self.i_faces) + ")"
  def __str__(self):
    return repr(self)

class Edge :
  def __init__(self, v0, v1, f0, f1):
    self.source = v0
    self.dest = v1
    self.on_right = f0
    self.on_left = f1
  def __repr__(self):
    return "Edge(" + str(self.source) + "," + str(self.dest) + "," + str(self.on_right) + "," + str(self.on_left) + ")"
  def __str__(self):
    return repr(self)

class Triangle :
  def __init__(self, V=None, E=None):
    self.i_edges = E
    self.i_verts = V
  def __repr__(self):
    return "Triangle(" + str(self.i_verts) + "," + str(self.i_edges) + ")"
  def __str__(self):
    return repr(self)

class TopSurface :
  def __init__(self, method=None, w=None, g=None):
    if method == None:
      self.v = []
      self.e = []
      self.f = []
    
    elif method=='polygon':
      Lw = len(w)
      v_labels = polygon_vert_out_labels(w)
      gens_to_edges = polygon_edge_inds(w)
      num_outside_edges = max([ei.ind for ei in gens_to_edges.values()])+1
      num_verts = len(v_labels) + 1
      self.e = [Edge(None,None,None,None) for i in xrange(num_outside_edges + Lw)]
      self.f = [None for i in xrange(Lw)]
      for i in xrange(num_outside_edges, num_outside_edges + Lw):
        self.e[i].source = num_verts-1
      
      # fill in the edge data for the vertices
      # and the vertex data for the edges
      v_incident_edges = []
      for j,V in enumerate(v_labels):
        ans = []
        i=0
        while i<len(V):
          gte = gens_to_edges[V[i]]
          ans.extend( [ gte, SI(num_outside_edges + V[i+1],-1) ] )
          if gte.sign>0:
            self.e[gte.ind].source = j
          else:
            self.e[gte.ind].dest = j
          self.e[num_outside_edges+V[i+1]].dest = j
          i += 2
        v_incident_edges.append(ans)
      
      v_incident_edges.append( [SI(num_outside_edges+i,1) for i in xrange(Lw)] )
      
      self.v = [Vertex(VIE, None) for VIE in v_incident_edges]
      
      #build the faces
      #it only knows the vertices now; not where in the vertex incidence 
      #list it is
      for i in xrange(Lw):
        e = [ SI(num_outside_edges+i), gens_to_edges[w[i]], SI(num_outside_edges+((i+1)%Lw), -1) ]
        v = [num_verts-1, self.e[e[0].ind].dest, self.e[e[2].ind].dest]
        for j,ei in enumerate(e):
          if ei.sign < 0:
            self.e[ei.ind].on_right = (i, j)
          else:
            self.e[ei.ind].on_left = (i, j)
        self.f[i] = Triangle(v,e)
      
      #fill in the face data for the edges and vertices
      #and also where in the vertex each angle is
      for i,V in enumerate(self.v):
        V.i_faces = [-1 for _ in xrange(len(V.i_edges))]
        for j in xrange(len(V.i_edges)):
          ei = V.i_edges[j]
          if ei.sign > 0:
            OL = self.e[ei.ind].on_left
            V.i_faces[j] = OL
            self.f[OL[0]].i_verts[OL[1]] = (i,j)
          else:
            OR = self.e[ei.ind].on_right
            V.i_faces[j] = OR
            self.f[OR[0]].i_verts[OR[1]] = (i,j)       
    
  def __repr__(self):
    return str(self)
  
  def __str__(self):
    ans = "Vertices: \n"
    for V in self.v:
      ans += str(V) + "\n"
    ans += "\nEdges: \n"
    for E in self.e:
      ans += str(E) + "\n"
    ans += "\nTriangles: \n"
    for F in self.f:
      ans += str(F) + "\n"
    return ans

  def euler_char(self):
    return len(self.v) - len(self.e) + len(self.f)
  

























        