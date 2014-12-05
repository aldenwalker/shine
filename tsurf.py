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
    self.i_tris = iF
  def __repr__(self):
    return "Vertex(" + str(self.i_edges) + "," + str(self.i_tris) + ")"
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
  def __init__(self, v, e, t):
    self.v = v
    self.e = e
    self.t = t

  @classmethod
  def from_polygon(cls, w):
    Lw = len(w)
    v_labels = polygon_vert_out_labels(w)
    gens_to_edges = polygon_edge_inds(w)
    num_outside_edges = max([ei.ind for ei in gens_to_edges.values()])+1
    num_verts = len(v_labels) + 1
    E = [Edge(None,None,None,None) for i in xrange(num_outside_edges + Lw)]
    T = [None for i in xrange(Lw)]
    for i in xrange(num_outside_edges, num_outside_edges + Lw):
      E[i].source = num_verts-1
    
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
          E[gte.ind].source = j
        else:
          E[gte.ind].dest = j
        E[num_outside_edges+V[i+1]].dest = j
        i += 2
      v_incident_edges.append(ans)
    
    v_incident_edges.append( [SI(num_outside_edges+i,1) for i in xrange(Lw)] )
    
    V = [Vertex(VIE, None) for VIE in v_incident_edges]
    
    #build the faces
    #it only knows the vertices now; not where in the vertex incidence 
    #list it is
    for i in xrange(Lw):
      e = [ SI(num_outside_edges+i), gens_to_edges[w[i]], SI(num_outside_edges+((i+1)%Lw), -1) ]
      v = [num_verts-1, E[e[0].ind].dest, E[e[2].ind].dest]
      for j,ei in enumerate(e):
        if ei.sign < 0:
          E[ei.ind].on_right = (i, j)
        else:
          E[ei.ind].on_left = (i, j)
      T[i] = Triangle(v,e)
    
    #fill in the face data for the edges and vertices
    #and also where in the vertex each angle is
    for i,v in enumerate(V):
      v.i_tris = [-1 for _ in xrange(len(v.i_edges))]
      for j in xrange(len(v.i_edges)):
        ei = v.i_edges[j]
        if ei.sign > 0:
          OL = E[ei.ind].on_left
          v.i_tris[j] = OL
          T[OL[0]].i_verts[OL[1]] = (i,j)
        else:
          OR = E[ei.ind].on_right
          v.i_tris[j] = OR
          T[OR[0]].i_verts[OR[1]] = (i,j)   
    return cls(V,E,T)    
  
  def attach_triangle_at(self, VI, I ):
    """place a new triangle with 0th vertex at vi in position i"""
    i_verts = [(VI, I), None, None]
    i_edges = [None,None,None]
    for j in xrange(3):
      vi,i = i_verts[j]
      ei = self.v[ vi ].i_edges[i]
      e = self.e[ei.ind]
      i_edges[j] = ei
      ov = (e.dest if ei.sign>0 else e.source)
      ovi = self.v[ov].i_edges.index( -ei )
      ovi = (ovi-1)%len(self.v[ov].i_edges)
      i_verts[(j+1)%3] = (ov, ovi)
    
    for j in xrange(3):
      self.v[ i_verts[j][0] ].i_tris[ i_verts[j][1] ] = (len(self.t),j)
      if i_edges[j].sign>0:
        self.e[ i_edges[j].ind ].on_left = (len(self.t), j)
      else :
        self.e[ i_edges[j].ind ].on_right = (len(self.t), j)
    
    self.t.append( Triangle(i_verts, i_edges) )
    return
  
  def fill_in_triangles(self):
    """fill in any triangles that aren't there, as detected by looking for """
    v_to_fill = [i for i,v in enumerate(self.v) if None in v.i_tris]
    if not isinstance(self.t, list):
      self.t = []
    for vi in v_to_fill:
      v = self.v[vi]
      for i in xrange(len(v.i_tris)):
        if v.i_tris[i] == None:
          self.attach_triangle_at( vi, i )
    return
    
  
  def __repr__(self):
    return str(self)
  
  def __str__(self):
    ans = "Topological surface\n"
    ans += "Vertices: \n"
    for i,V in enumerate(self.v):
      ans += str(i) + ": " + str(V) + "\n"
    ans += "\nEdges: \n"
    for i,E in enumerate(self.e):
      ans += str(i) + ": " + str(E) + "\n"
    ans += "\nTriangles: \n"
    for i,T in enumerate(self.t):
      ans += str(i) + ": " + str(T) + "\n"
    return ans

  def euler_char(self):
    return len(self.v) - len(self.e) + len(self.t)
  
  def highest_valence_vertex(self):
    L = [(i,len(v.i_edges)) for i,v in enumerate(self.v)]
    L.sort(key=lambda x:x[1])
    return L[-1][0]

  def is_isolated_vertex(self, vi):
    """an isolated vertex has the property that none of its incident 
    edges touch it twice"""
    E = [ei.ind for ei in self.v[vi].i_edges]
    return len(E) == len(set(E))

  def find_isolated_vertex(self):
    for vi in xrange(len(self.v)):
      if self.is_isolated_vertex(vi):
        return vi
    return None























        