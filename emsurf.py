import tsurf
import R3
import copy
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

#############################################################################
# a planar graph class
# for technical reasons, it also can remember fundamental group
# generators of a neighborhood
#############################################################################
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
  def __init__(self, v, e, loops=None):
    self.v = v
    self.e = e
    self.loops = (dict() if loops == None else loops)
  
  @classmethod
  def from_file(cls, graph_filename):
    f = open(graph_filename, 'r')
    lines = f.read().split('\n')
    f.close()
    lines = [ell for ell in lines if len(ell)>0 and ell[0] != '#']
    num_v, num_e = map(int, lines[0].split())
    V = num_v*[None]
    E = num_e*[None]
    for i in xrange(num_v):
      V[i] = PlanarVertex(complex(*map(float, lines[i+1].split())), [], [])
    for i in xrange(num_e):
      E[i] = PlanarEdge(*map(int, lines[i+num_v+1].split()))
      V[E[i].src].i_edges.append(SI(i,1))
      V[E[i].dest].i_edges.append(SI(i,-1))
    
    #record the loops, if there are any
    loops = dict()
    for i in xrange(1+num_v+num_e, len(lines)):
      ell = lines[i].split()
      loops[ell[0]] = ell[1:]
    
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
    return cls(V,E,loops=loops)
  
  def __str__(self):
    ans = "Vertices: \n"
    for i,v in enumerate(self.v):
      ans += str(i) + ": " + str(v) + "\n"
    ans += "Edges:\n"
    for i,e in enumerate(self.e):
      ans += str(i) + ": " + str(e) + "\n"
    ans += "Loops:\n"
    for ell in self.loops:
      ans += str(ell) + ": " + str(self.loops[ell]) + "\n"
    return ans
  
  


##########################################################################
# a topological path in an embedded surface
# the map is recorded by saying how far along each edge it is
##########################################################################
class EmbeddedPath(tsurf.TopologicalPath):
  def __init__(self, edges, edge_coords):
    self.edges = edges
    self.edge_coords = edge_coords
  
  @classmethod
  def from_topological_path(cls, tp):
    return cls(tp.edges, [0.5 for e in tp.edges])

  def __repr__(self):
    return str(self)
  
  def __str__(self):
    return "EmbeddedPath(" + str(zip(self.edges, self.edge_coords)) + ")"
  
  def subdivide(self, old_TS, vertices_from_edges, edges_from_edges, edges_from_tris, tris_from_tris):
    pass
    

###########################################################################
# a topological surface, embedded in R3
###########################################################################
class EmbeddedSurface(tsurf.TopSurface):
  def __init__(self, TS, em_v, em_e, em_t, em_loops):
    self.v = TS.v
    self.e = TS.e
    self.t = TS.t
    self.loops = TS.loops
    self.em_v = em_v
    self.em_e = em_e
    self.em_t = em_t
    self.em_loops = em_loops
  
  @classmethod
  def from_planar_graph_file(cls, graph_filename):
    return cls.from_planar_graph(PlanarGraph.from_file(graph_filename))
  
  #########################################################################
  # create an embedded surface which is a neighborhood of a graph
  #########################################################################
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
      for j in xrange(val):
        V_from_PG[i]['around'][j] = len(V)
        V.append( tsurf.Vertex( 6*[None], 6*[None] ) )
        ang = average_angle( pv.i_edge_angles[j], pv.i_edge_angles[(j+1)%val] )
        em_V.append( R3.Vector( (pv.pt.real + 0.5*math.cos(ang),                    \
                                 pv.pt.imag + 0.5*math.sin(ang),                    \
                                 0) ))
    #create edges around the vertices
    E_from_PGV = [dict() for i in xrange(len(PG.v))]
    for i,pv in enumerate(PG.v):
      val = len(pv.i_edges)
      E_from_PGV[i]['around'] = val*[None]
      for j in xrange(val):
        V[ V_from_PG[i]['top'] ].i_edges[2*j+1] = SI(len(E),1)
        V[ V_from_PG[i]['around'][j] ].i_edges[0] = SI(len(E),-1)
        E_from_PGV[i]['around'][j] = SI(len(E),1)
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
      E_from_PG[i]['top_right'] = len(E)
      E.append( tsurf.Edge( vs, vd, None, None) )
      
      vs = V_from_PG[pe.src]['top']
      vsi = V[vs].i_edges.index( SI(E_from_PG[i]['top'],1) )
      vsi = (vsi+1)%len(V[vs].i_edges)
      vd = V_from_PG[pe.dest]['around'][ (dest_index-1)%dest_val ]
      V[vs].i_edges.insert(vsi, SI(len(E),1))
      V[vs].i_tris.append(None)
      V[vd].i_edges[5] = SI(len(E),-1)  
      E_from_PG[i]['top_left'] = len(E)
      E.append( tsurf.Edge( vs, vd, None, None) ) 
      
      vs = V_from_PG[pe.src]['around'][src_index]
      vd = V_from_PG[pe.dest]['bottom']
      vdi = V[vd].i_edges.index( SI(E_from_PG[i]['bottom'],-1) )
      vdi = (vdi+1)%len(V[vd].i_edges)
      V[vs].i_edges[2] = SI(len(E),1)
      V[vd].i_edges.insert(vdi, SI(len(E),-1))  
      V[vd].i_tris.append(None)
      E_from_PG[i]['bottom_left'] = len(E)
      E.append( tsurf.Edge( vs, vd, None, None) )    
      
      vs = V_from_PG[pe.src]['bottom']
      vsi = V[vs].i_edges.index( SI(E_from_PG[i]['bottom'],1) )
      vsi = (vsi+1)%len(V[vs].i_edges)
      vd = V_from_PG[pe.dest]['around'][dest_index]
      V[vs].i_edges.insert(vsi, SI(len(E),1))
      V[vs].i_tris.append(None)
      V[vd].i_edges[2] = SI(len(E),-1)
      E_from_PG[i]['bottom_right'] = len(E)
      E.append( tsurf.Edge( vs, vd, None, None) )
    
    #create the loops, if there are any
    loops = dict()
    for ell in PG.loops:
      print "Making loop: ", ell, PG.loops[ell]
      raw_loop_edges = [(SI.from_string(l) if l!='a' else 'a') for l in PG.loops[ell]]
      Lrle = len(raw_loop_edges)
      padded_loop_edges = []
      #figure out which *graph* edges need to get added to get around the vertices
      for i,ei in enumerate(raw_loop_edges):
        if ei == 'a':
          continue
        next_ei = raw_loop_edges[(i+1)%Lrle]
        if next_ei == 'a':
          next_ei = raw_loop_edges[(i+2)%Lrle]
          padded_loop_edges.extend([ei, 'a'])
        else:
          padded_loop_edges.append(ei)
        if next_ei == ei or next_ei == -ei:
          continue
        central_vi = (PG.e[ei.ind].dest if ei.sign>0 else PG.e[ei.ind].src)
        central_v = PG.v[central_vi]
        central_v_ei1 = PG.v[central_vi].i_edges.index( -ei )
        central_v_ei2 = PG.v[central_vi].i_edges.index( next_ei )
        extra_edges = (central_v.i_edges[central_v_ei1+1:central_v_ei2] \
                                          if central_v_ei1<central_v_ei2 else  \
                       central_v.i_edges[central_v_ei1+1:] + central_v.i_edges[:central_v_ei2])
        extra_edges = [SI(x.ind, s*x.sign) for x in extra_edges for s in [1,-1]]
        padded_loop_edges.extend(extra_edges)
      print "Padded loop edges: ", padded_loop_edges
      #now actually figure out which real edges we cross
      loop_edges = []
      Lple = len(padded_loop_edges)
      for i,ei in enumerate(padded_loop_edges):
        next_ei = padded_loop_edges[(i+1)%Lple]
        if next_ei == 'a':
          #add the edges around
          loop_edges.extend([ SI(E_from_PG[ei.ind][x], -1) for x in ['top_right', \
                                                                         'right', \
                                                                  'bottom_right', \
                                                                        'bottom', \
                                                                   'bottom_left', \
                                                                          'left', \
                                                                      'top_left', \
                                                                        'top']])
          next_ei = padded_loop_edges[(i+2)%Lple]
        if ei == 'a' or next_ei == ei or next_ei == -ei:
          continue
        #add the edges to get from ei to next_ei
        vi = (PG.e[ei.ind].dest if ei.sign>0 else PG.e[ei.ind].src)
        vii = PG.v[vi].i_edges.index(-ei)
        v_crossing_edge = E_from_PGV[vi]['around'][vii].ind #the sign is always positive
        if ei.sign>0:
          loop_edges.append( SI(E_from_PG[ei.ind]['top_right'],-1) )
        else:
          loop_edges.extend( [SI(E_from_PG[ei.ind][x],1) for x in ['top', 'top_left']] )
        loop_edges.append( SI(v_crossing_edge, 1) )
        if next_ei.sign>0:
          pass #don't append anything
        else:
          loop_edges.append( SI(E_from_PG[next_ei.ind]['top'],-1) )
      loops[ell] = loop_edges

    ES = cls( tsurf.TopSurface(V,E,[],loops=loops), em_V, len(E)*[None], len(T)*[None], None)
    #print ES
    ES.fill_in_triangles()
    ES.em_loops = dict([(ell, EmbeddedPath.from_topological_path(ES.loops[ell])) for ell in ES.loops])
    
    ES.em_e = len(ES.e)*[None]
    for i,e in enumerate(ES.e):
      ES.em_e[i] = [ES.em_v[e.source], ES.em_v[e.dest]]
    ES.em_t = len(ES.t)*[None]
    for i,t in enumerate(ES.t):
      ES.em_t[i] = [ES.em_v[t.i_verts[j][0]] for j in xrange(3)]
    
    return ES
  
  def subdivide(self):
    old_nv = len(self.v)
    old_ne = len(self.e)
    old_nt = len(self.t)
    old_TS, vertices_from_edges, edges_from_edges, edges_from_tris, tris_from_tris = super(EmbeddedSurface, self).subdivide()
    self.em_v = [(self.em_v[i] if i<old_nv else None) for i in xrange(len(self.v))]
    self.em_e = len(self.e)*[None]
    self.em_t = len(self.t)*[None]
    for i in xrange(old_ne):
      v0i = self.e[edges_from_edges[i][0]].source
      v1i = self.e[edges_from_edges[i][1]].dest
      new_vi = self.e[edges_from_edges[i][0]].dest
      self.em_v[new_vi] = (self.em_v[v0i] + self.em_v[v1i])/2.0
    for ei,e in enumerate(self.e):
      self.em_e[ei] = [ self.em_v[e.source], self.em_v[e.dest] ]
    for ti,t in enumerate(self.t):
      self.em_t[ti] = [ self.em_v[t.i_verts[i][0]] for i in xrange(3) ]
    for ell in self.loops:
      self.loops[ell].subdivide(old_TS, vertices_from_edges, edges_from_edges, edges_from_tris, tris_from_tris)
      self.em_loops[ell] = EmbeddedPath.from_topological_path(self.loops[ell]) #(old_TS, vertices_from_edges, edges_from_edges, edges_from_tris, tris_from_tris)
    return old_TS, vertices_from_edges, edges_from_edges, edges_from_tris, tris_from_tris
      
  def flow(self):
    for i,v in enumerate(self.v):
      av = R3.Vector([0,0,0])
      for ei in v.i_edges:
        ovi = (self.e[ei.ind].dest if ei.sign>0 else self.e[ei.ind].source)
        av = av + self.em_v[ovi]
      av = av*(1.0/len(v.i_edges))
      self.em_v[i] = self.em_v[i]*0.8 + av*0.2
    for i,e in enumerate(self.e):
      self.em_e[i] = [self.em_v[e.source], self.em_v[e.dest]]
    for i,t in enumerate(self.t):
      self.em_t[i] = [self.em_v[t.i_verts[j][0]] for j in xrange(3)]
    return
  
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




































