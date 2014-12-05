import tsurf
from signedind import SignedInd as SI

import math

class PlanarVertex:
  def __init__(self, pt, i_edges):
    self.pt = pt
    self.i_edges = i_edges
  def __repr__(self):
    return "PlanarVertex("+str(self.pt) + "," + str(self.i_edges)+")"
  def __str__(self):
    return repr(self)

class PlanarEdge:
  def __init__(self, src, dest):
    self.src = src
    self.dest = dest
  def __repr__(self):
    return "PlanarEdge" + str((self.src, self.dest))
  def __str__(self):
    return repr(self)
    
class PlanarGraph:
  def __init__(self, v, e):
    self.v = v
    self.e = e
  
  @classmethod
  def from_file(cls, graph_filename):
    f = open(graph_filename, 'r')
    lines = f.read().split('\n')
    f.close()
    num_v, num_e = map(int, lines[0])
    V = num_v*[None]
    E = num_e*[None]
    for i in xrange(num_v):
      V[i] = PlanarVertex(complex(*map(float, lines[i+1].split())), [])
    for i in xrange(num_g_edges):
      E[i] = PlanarEdge(*map(int, lines[i+num_g_verts+1].split()))
      V[E[i].src].i_edges.append(SI(i,1))
      V[E[i].dest].i_edges.append(SI(i,-1))
    
    #sort the edges around each vertex and record the angle
    for v in V:
      for i,ei in enumerate(v.i_edges):
        ovi = (E[ei.ind].dest if ei.sign>0 else E[ei.ind].src)
        vec = V[ovi].pt - v.pt
        ang = math.atan2(vec.imag, vec.real)
        v.i_edges[i] = (ei, ang)
      v.i_edges.sort(key=lambda x:x[1])
    return cls(V,E)

class EmbeddedSurface(tsurf.TopSurface):
  def __init__(self, TS, em_v, em_e, em_t):
    self.v = TS.v
    self.e = TS.e
    self.t = TS.t
    self.em_v = em_v
    self.em_e = em_e
    self.em_t = em_t
  
  @classmethod
  def from_planar_graph_file(cls, graph_filename):
    return cls.from_planar_graph(PlanarGraph.from_file(graph_filename))
  
  @classmethod
  def from_planar_graph(cls, PG):
    
    
    V = []
    






































