import tsurf
import R3
from signedind import SignedInd as SI

import math

def average_angle(a1,a2):
  """compute the average angle (halfway from a1 to a2)"""
  while a2 < a1:
    a2 += 2*math.pi
  while a2 > a1 + 2*math.pi:
    a2 -= 2*math.pi
  ans = 0.5*(a1+a2)
  while ans < 0:
    ans += 2*math.pi
  while ans > 2*math.pi:
    ans -= 2*math.pi
  return ans

class PlanarVertex:
  def __init__(self, pt, i_edges, angles):
    self.pt = pt
    self.i_edges = i_edges
    self.i_edge_angles = angles
  def __repr__(self):
    return "PlanarVertex("+str(self.pt) + "," + str(self.i_edges)+","+str(self.i_edge_angles) + ")"
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
    num_v, num_e = map(int, lines[0].split())
    V = num_v*[None]
    E = num_e*[None]
    for i in xrange(num_v):
      V[i] = PlanarVertex(complex(*map(float, lines[i+1].split())), [], [])
    for i in xrange(num_e):
      E[i] = PlanarEdge(*map(int, lines[i+num_v+1].split()))
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
      v.i_edge_angles = [x[1] for x in v.i_edges]
      v.i_edges = [x[0] for x in v.i_edges]
    return cls(V,E)
  
  def __str__(self):
    ans = "Vertices: \n"
    for i,v in enumerate(self.v):
      ans += str(i) + ": " + str(v) + "\n"
    ans += "Edges:\n"
    for i,e in enumerate(self.e):
      ans += str(i) + ": " + str(e) + "\n"
    return ans
  
  








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
    E = []
    T = []
    em_V = []
    V_from_PG = len(PG.v)*[None]
    #create vertices
    for i,pv in enumerate(PG.v):
      val = len(pv.i_edges)
      V_from_PG[i] = dict()
      
      V_from_PG[i]['top'] = len(V)
      V.append( tsurf.Vertex( 2*val*[None], 2*val*[None] ) )
      em_V.append( R3.Vector( (pv.pt.real, pv.pt.imag, 0.5)) )
      
      V_from_PG[i]['bottom'] = len(V)
      V.append( tsurf.Vertex( 2*val*[None], 2*val*[None] ) )
      em_V.append( R3.Vector( (pv.pt.real, pv.pt.imag, -0.5)) )
      
      V_from_PG[i]['around'] = val*[None]
      print "Doing around vertices for ", i, "of valence ", val
      for j in xrange(val):
        V_from_PG[i]['around'][j] = len(V)
        V.append( tsurf.Vertex( 6*[None], 6*[None] ) )
        ang = average_angle( pv.i_edge_angles[j], pv.i_edge_angles[(j+1)%val] )
        print "doing angle", ang, " average of ", pv.i_edge_angles[j], pv.i_edge_angles[(j+1)%val]
        em_V.append( R3.Vector( (pv.pt.real + 0.5*math.cos(ang),                    \
                                 pv.pt.imag + 0.5*math.sin(ang),                    \
                                 0) ))
    #create edges around the vertices
    for i,pv in enumerate(PG.v):
      val = len(pv.i_edges)
      for j in xrange(val):
        V[ V_from_PG[i]['top'] ].i_edges[2*j+1] = SI(len(E),1)
        V[ V_from_PG[i]['around'][j] ].i_edges[0] = SI(len(E),-1)
        E.append( tsurf.Edge( V_from_PG[i]['top'], V_from_PG[i]['around'][j], None, None ) )
        
        V[ V_from_PG[i]['bottom'] ].i_edges[-2*(j+1)] = SI(len(E),1)
        V[ V_from_PG[i]['around'][j] ].i_edges[3] = SI(len(E),-1)
        E.append( tsurf.Edge( V_from_PG[i]['bottom'], V_from_PG[i]['around'][j], None, None ) )
    
    #create edges from edges
    #(but not diagonal because that'll make the top difficult)
    E_from_PG = [dict() for i in xrange(len(PG.e))]
    for i,pe in enumerate(PG.e):
      src_index = PG.v[pe.src].i_edges.index( SI(i,1) )
      src_val = len(PG.v[pe.src].i_edges)
      dest_index = PG.v[pe.dest].i_edges.index( SI(i,-1) )
      dest_val = len(PG.v[pe.dest].i_edges)
      
      V[ V_from_PG[pe.src]['top'] ].i_edges[ 2*src_index ] = SI(len(E),1)
      V[ V_from_PG[pe.dest]['top'] ].i_edges[ 2*dest_index ] = SI(len(E),-1)
      E_from_PG[i]['top'] = len(E)
      E.append( tsurf.Edge( V_from_PG[pe.src]['top'], V_from_PG[pe.dest]['top'], None, None ) )
      
      V[ V_from_PG[pe.src]['bottom'] ].i_edges[ -2*(src_index+1)+1 ] = SI(len(E),1)
      V[ V_from_PG[pe.dest]['bottom'] ].i_edges[ -2*(dest_index+1)+1 ] = SI(len(E),-1)
      E_from_PG[i]['bottom'] = len(E)
      E.append( tsurf.Edge( V_from_PG[pe.src]['bottom'], V_from_PG[pe.dest]['bottom'], None, None ) )
      
      vs = V_from_PG[pe.src]['around'][src_index]
      vd = V_from_PG[pe.dest]['around'][(dest_index-1)%dest_val]
      E_from_PG[i]['left'] = len(E)
      V[ vs ].i_edges[1] = SI(len(E),1)
      V[ vd ].i_edges[4] = SI(len(E),-1)
      E.append( tsurf.Edge( vs, vd, None, None ) )
      
      vs = V_from_PG[pe.src]['around'][(src_index-1)%src_val]
      vd = V_from_PG[pe.dest]['around'][dest_index]
      E_from_PG[i]['right'] = len(E)
      V[ vs ].i_edges[4] = SI(len(E),1)
      V[ vd ].i_edges[1] = SI(len(E),-1)
      E.append( tsurf.Edge( vs, vd, None, None ) )
    
    #print "Before diagonals"
    #print V
    #print E
    
    #create the diagonal edges
    for i,pe in enumerate(PG.e):
      src_index = PG.v[pe.src].i_edges.index( SI(i,1) )
      src_val = len(PG.v[pe.src].i_edges)
      dest_index = PG.v[pe.dest].i_edges.index( SI(i,-1) )
      dest_val = len(PG.v[pe.dest].i_edges)
      
      vs = V_from_PG[pe.src]['around'][(src_index-1)%src_val]
      vd = V_from_PG[pe.dest]['top']
      vdi = V[vd].i_edges.index( SI(E_from_PG[i]['top'],-1) )
      vdi = (vdi+1)%len(V[vd].i_edges)
      V[vs].i_edges[5] = SI(len(E),1)
      V[vd].i_edges.insert(vdi, SI(len(E),-1))
      V[vd].i_tris.append(None)
      E.append( tsurf.Edge( vs, vd, None, None) )
      
      vs = V_from_PG[pe.src]['top']
      vsi = V[vs].i_edges.index( SI(E_from_PG[i]['top'],1) )
      vsi = (vsi+1)%len(V[vs].i_edges)
      vd = V_from_PG[pe.dest]['around'][ (dest_index-1)%dest_val ]
      V[vs].i_edges.insert(vsi, SI(len(E),1))
      V[vs].i_tris.append(None)
      V[vd].i_edges[5] = SI(len(E),-1)  
      E.append( tsurf.Edge( vs, vd, None, None) ) 
      
      vs = V_from_PG[pe.src]['around'][src_index]
      vd = V_from_PG[pe.dest]['bottom']
      vdi = V[vd].i_edges.index( SI(E_from_PG[i]['bottom'],-1) )
      vdi = (vdi+1)%len(V[vd].i_edges)
      V[vs].i_edges[2] = SI(len(E),1)
      V[vd].i_edges.insert(vdi, SI(len(E),-1))  
      V[vd].i_tris.append(None)
      E.append( tsurf.Edge( vs, vd, None, None) )    
      
      vs = V_from_PG[pe.src]['bottom']
      vsi = V[vs].i_edges.index( SI(E_from_PG[i]['bottom'],1) )
      vsi = (vsi+1)%len(V[vs].i_edges)
      vd = V_from_PG[pe.dest]['around'][dest_index]
      V[vs].i_edges.insert(vsi, SI(len(E),1))
      V[vs].i_tris.append(None)
      V[vd].i_edges[2] = SI(len(E),-1)
      E.append( tsurf.Edge( vs, vd, None, None) )
    
    ES = cls( tsurf.TopSurface(V,E,[]), em_V, len(E)*[None], len(T)*[None])
    #print ES
    ES.fill_in_triangles()
    
    ES.em_e = len(ES.e)*[None]
    for i,e in enumerate(ES.e):
      ES.em_e[i] = [ES.em_v[e.source], ES.em_v[e.dest]]
    ES.em_t = len(ES.t)*[None]
    for i,t in enumerate(ES.t):
      ES.em_t[i] = [ES.em_v[t.i_verts[j][0]] for j in xrange(3)]
    
    return ES
  
  def __str__(self):
    ans = "Embedded surface:\n";
    ans += "Vertices: \n"
    for i,v in enumerate(self.v):
      ans += str(i) + ": " + str(v) + "\n" + "embedded: " + str(self.em_v[i]) + "\n"
    ans += "Edges: \n"
    for i,e in enumerate(self.e):
      ans += str(i) + ": " + str(e) + "\n" + "embedded: " + str(self.em_e[i]) + "\n"
    ans += "Triangles: \n"
    for i,t in enumerate(self.t):
      ans += str(i) + ": " + str(t) + "\n" + "embedded: " + str(self.em_t[i]) + "\n"
    return ans




































