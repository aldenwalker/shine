import tsurf
import gsurf
import emsurf
import hyp
from signedind import SignedInd as SI

import math


def hash_complex(p):
  return ( int(1e4*p.real), int(1e4*p.imag) )

class CoveringVertex(tsurf.Vertex):
  def __init__(self, pt, covered_v, iE, iT):
    self.i_edges = iE
    self.i_tris = iT
    self.pt = pt
    self.covered_v = covered_v
  def __repr__(self):
    return "CoveringVertex(" + str(self.pt) + "," + str(self.covered_v) + "," + str(self.i_edges) + str(self.i_tris) + ")"
  def __str__(self):
    return repr(self)
 
class CoveringEdge(tsurf.Edge):
  def __init__(self, gi, covered_e, source, dest, on_left, on_right):
    self.gi, self.covered_e, self.source = gi, covered_e, source
    self.dest, self.on_left, self.on_right = dest, on_left, on_right
  def __repr__(self):
    return "CoveringEdge(" + str(self.gi) + "," + str(self.covered_e) + "," \
                           + str(self.source) + "," + str(self.dest) + "," \
                           + str(self.on_left) + "," + str(self.on_right) + ")"
  def __str__(self):
    return repr(self)
  
class CoveringTri(tsurf.Triangle):
  def __init__(self, t, covered_t, i_edges, i_verts):
    self.t, self.covered_t = t, covered_t
    self.i_edges, self.i_verts = i_edges, i_verts
  def __repr__(self):
    return "CoveringTri(" + str(self.t) + "," + str(self.covered_t) + ","  \
                          + str(self.i_edges) + str(self.i_verts) + ")"
  def __str__(self):
    return repr(self)


class LiftedSurface(gsurf.GeoSurface):
  def __init__(self, GS, ev, ee, et, v_lifts, e_lifts, t_lifts):
    self.v, self.e, self.t = GS.v, GS.e, GS.t
    self.h_tris, self.h_lengths = GS.h_tris, GS.h_lengths
    self.GS = GS
    self.em_v = ev
    self.em_e = ee
    self.em_t = et
    self.v_lifts = v_lifts
    self.e_lifts = e_lifts
    self.t_lifts = t_lifts
    self.v_hashes = dict([( hash_complex(v.pt), i) for i,v in enumerate(self.em_v)])
  
  @classmethod
  def lift_gsurf(cls, GS):
    LS = cls(GS, [], [], [], [], [], [])
    LS.relay()
    return LS
  
  def is_vertex_surrounded(self, vi):
    """are all adjacent triangles placed?"""
    return all([len(self.t_lifts[ti])>0 for ti,j in self.v[vi].i_tris])
  
  def attach_triangle_to_gi(self, ti, i, gi):
    """do the bookkeeping to attach edge index i in triangle index ti
    to the geodesic interval gi (it must be of the correct length)"""
    ei = self.t[ti].i_edges[i]
    self.em_t[ti] = self.h_tris[ti].realize_along_gi(gi, i)
    #record what we've done
    for j,ei in enumerate(self.t[ti].i_edges):
      if ei.sign > 0:
        self.em_e[ei.ind][0] = self.em_t[ti].sides[j]
        self.em_v[self.e[ei.ind].source].append( self.em_t[ti].v[j] )
      else:
        self.em_e[ei.ind][1] = self.em_t[ti].sides[j]
        self.em_v[self.e[ei.ind].dest].append( self.em_t[ti].v[j] )
    return
  
  #######################################################################
  # lifts the appropriate triangle so that it contains the lifted_ei
  # with the appropriate sign
  #######################################################################
  def lift_triangle_to_lifted_edge(self, lifted_ei):
    s = lifted_ei.sign
    lifted_e = self.em_e[lifted_ei.ind]
    lower_e = self.e[lifted_e.covered_e]
    if (lifted_e.on_left if s>0 else lifted_e.on_right) != None:
      raise ValueError("Attaching a triangle where there is one already?")
    gi = (lifted_e.gi if s>0 else lifted_e.gi.reversed())
    lower_ti,j = (lower_e.on_left if s>0 else lower_e.on_right)
    lower_t = self.t[lower_ti]
    new_t = self.h_tris[lower_ti].realize_along_gi(gi, j)
    new_t = CoveringTri(new_t, lower_ti, 3*[None], 3*[None])
    new_t.i_edges[j] = lifted_ei
    new_t.i_verts[j] = ( (lifted_e.source if s>0 else lifted_e.dest), lower_t.i_verts[j][1])
    new_t.i_verts[(j+1)%3] = ( (lifted_e.dest if s>0 else lifted_e.source), lower_t.i_verts[(j+1)%3][1] )
    ov = None
    
    #determine if the edge before exists
    pvie = self.em_v[new_t.i_verts[j][0]].i_edges
    pe = pvie[ (new_t.i_verts[j][1] + 1)%len(pvie) ] 
    jm1 = (j-1)%3
    if pe != None:
      new_t.i_edges[jm1] = -pe
      ov = (self.em_e[pe.ind].dest if pe.sign>0 else self.em_e[pe.ind].source)
    else:
      new_t.i_edges[jm1] = SI( len(self.em_e), lower_t.i_edges[jm1].sign )
      self.e_lifts[lower_t.i_edges[jm1].ind].append( len(self.em_e) )
      self.em_e.append(CoveringEdge(None,None,None,None,None,None))
    
    #determine if the next edge exists
    jp1 = (j+1)%3
    nvie = self.em_v[new_t.i_verts[jp1][0]].i_edges
    ne = nvie[ new_t.i_verts[jp1][1] ]
    if ne != None:
      new_t.i_edges[jp1] = ne
      ov = (self.em_e[ne.ind].dest if ne.sign>0 else self.em_e[ne.ind].source)
    else:
      new_t.i_edges[jp1] = SI( len(self.em_e), lower_t.i_edges[jp1].sign )
      self.e_lifts[lower_t.i_edges[jp1].ind].append( len(self.em_e) )
      self.em_e.append(CoveringEdge(None,None,None,None,None,None))
    
    #determine if the other vertex exists
    if ov == None:
      #check if it's secretly there
      ov_pt = new_t.t.sides[(j+2)%3].start
      ov = self.v_hashes.get( hash_complex(ov_pt), None )
    if ov != None:
      new_t.i_verts[(j+2)%3] = ( ov, lower_t.i_verts[(j+2)%3][1] )
    else:
      new_t.i_verts[(j+2)%3] = ( len(self.em_v), lower_t.i_verts[(j+2)%3][1] )
      lower_vi = lower_t.i_verts[(j+2)%3][0]
      self.v_lifts[lower_vi].append(len(self.em_v))
      self.em_v.append( CoveringVertex( None, None, len(self.v[lower_vi].i_edges)*[None], len(self.v[lower_vi].i_tris)*[None]) )
    
    new_ti = len(self.em_t)
    self.em_t.append(new_t)
    self.t_lifts[lower_ti].append(new_ti)
    
    #now go around and fill in all relevant info all the way around
    #this may overwrite stuff, but it *should* be the same
    IV = new_t.i_verts
    IE = new_t.i_edges
    for i in xrange(3):
      #fill in vert
      vi,j = IV[i]
      v = self.em_v[vi]
      v.pt = new_t.t.v[i]
      self.v_hashes[ hash_complex(v.pt) ] = vi
      v.covered_v = lower_t.i_verts[i][0]
      v.i_edges[j] = new_t.i_edges[i]
      v.i_edges[(j+1)%len(v.i_edges)] = -new_t.i_edges[(i-1)%3]
      v.i_tris[j] = (new_ti, i)
      #fill in edge
      ei = IE[i]
      e = self.em_e[ei.ind]
      if e.gi == None:
        e.gi = (new_t.t.sides[i] if ei.sign>0 else new_t.t.sides[i].reversed())
      e.covered_e = lower_t.i_edges[i].ind 
      e.source = (new_t.i_verts[i][0] if ei.sign>0 else new_t.i_verts[(i+1)%3][0])
      e.dest = (new_t.i_verts[(i+1)%3][0] if ei.sign>0 else new_t.i_verts[i][0])
      if ei.sign>0:
        e.on_left = (new_ti, i)
      else:
        e.on_right = (new_ti, i)

  #########################################################################
  # relift the liftedsurface by going back to the geometric surface
  # and lifting all the triangles
  #########################################################################
  def relay(self):
    num_edges = len(self.e)
    self.em_v = []
    self.v_lifts = [[] for v in self.v]
    self.em_e = []
    self.e_lifts = [[] for e in self.e]
    self.em_t = []
    self.t_lifts = [[] for t in self.t]
    
    vdata = [(i, self.is_isolated_vertex(i), len(v.i_edges)) for i,v in enumerate(self.v)]
    vdata.sort(key=lambda x:(int(x[1]), x[2]), reverse=True)
    first_vertex = vdata[0][0]
    
    #lift the first triangle so the first edge happens to go straight up
    T = self.t[0]
    IE = T.i_edges
    gi = hyp.HypGeodesicInterval.from_pt_angle_dist(1j, math.pi/2.0, self.h_lengths[IE[0].ind])
    if IE[0].sign < 0:
      giR = gi.reversed()
    new_t = CoveringTri( self.h_tris[0].realize_along_gi(gi, 0), 0, 3*[None], 3*[None] )
    new_es = [None, None, None]
    new_vs = [None, None, None]
    self.t_lifts[0].append(0)
    for i in xrange(3):
      new_es[i] = CoveringEdge( (new_t.t.sides[i] if IE[i].sign > 0 else new_t.t.sides[i].reversed()), \
                                IE[i].ind,                           \
                                (i if IE[i].sign > 0 else (i+1)%3),          \
                                ((i+1)%3 if IE[i].sign > 0 else i),          \
                                (0 if IE[i].sign > 0 else None),             \
                                (None if IE[i].sign > 0 else 0) )
      new_t.i_edges[i] = SI(i, IE[i].sign)
      self.e_lifts[IE[i].ind].append(i)
      cvi = T.i_verts[i]
      new_vs[i] = CoveringVertex( new_t.t.v[i], cvi[0],                        \
                                  len(self.v[cvi[0]].i_edges)*[None],        \
                                  len(self.v[cvi[0]].i_tris)*[None] )
      self.v_hashes[ hash_complex(new_vs[i].pt) ] = i
      new_vs[i].i_tris[cvi[1]] = (0, i)
      new_vs[i].i_edges[cvi[1]] = SI(i, IE[i].sign)
      new_vs[i].i_edges[(cvi[1]+1)%len(new_vs[i].i_edges)] = SI((i-1)%3, -IE[(i-1)%3].sign)
      new_t.i_verts[i] = (i, cvi[1])
      self.v_lifts[cvi[0]].append(i)
    self.em_v.extend(new_vs)
    self.em_e.extend(new_es)
    self.em_t.append(new_t)
    
    while True:
      #find a vertex which is (in increasing order of desirability and rarity):
      #lifted
      #partially lifted
      #appears in the cover only once so far
      #isolated topologically
      lifted = [i for i in xrange(len(self.v)) if len(self.v_lifts[i])>0]
      partially_lifted = [i for i in lifted if not self.is_vertex_surrounded(i)]
      appears_once = [i for i in partially_lifted if len(self.v_lifts[i])==1]
      isolated = [i for i in appears_once if self.is_isolated_vertex(i)]
      if len(isolated)>1:
        isolated.sort(key=lambda x:len(self.v[x].i_edges), reverse=True)
      current_v = (isolated[0] if len(isolated)>0 else                         \
                  (appears_once[0] if len(appears_once)>0 else                 \
                  (partially_lifted[0] if len(partially_lifted)>0 else None)))
      if current_v == None:
        break
      
      #go around the vertex to find a triangle which hasn't been placed
      i=0
      IT = self.v[current_v].i_tris
      IE = self.v[current_v].i_edges
      LIT = len(IT)
      while len(self.t_lifts[IT[i][0]])==0:
        i = (i+1)%LIT
      while len(self.t_lifts[IT[i][0]])>0:
        i = (i+1)%LIT
      while len(self.t_lifts[IT[i][0]])==0:
        im1 = (i-1)%LIT
        lifted_ti = self.t_lifts[IT[im1][0]][0]
        lifted_t = self.em_t[lifted_ti]
        lifted_ei = -lifted_t.i_edges[ (IT[im1][1]-1)%3 ]
        #attach this triangle
        self.lift_triangle_to_lifted_edge(lifted_ei)
        i = (i+1)%LIT
    return
  ###########################################################################
  # act by a mobius transformation
  ###########################################################################
  def act_by_mobius(self, M):
    self.v_hashes = dict()
    for i,v in enumerate(self.em_v):
      v.pt = M(v.pt)
      self.v_hashes[hash_complex(v.pt)] = i
    for e in self.em_e:
      e.gi = e.gi.act_by_mobius(M)
    for t in self.em_t:
      t.t = t.t.act_by_mobius(M)
    return None
  
  def __repr__(self):
    return str(self)
  
  def __str__(self):
    ans = "Lifted Surface:"
    ans += "Vertices: \n"
    for i,V in enumerate(self.v):
      ans += str(i) + ": " + str(V) + "\n"
      ans += "covered by: " + str(self.v_lifts[i]) + "\n"
    ans += "\nEdges: \n"
    for i,E in enumerate(self.e):
      ans += str(i) + ": " + str(E) + "\n"
      ans += "covered by: " + str(self.e_lifts[i]) + "\n"
    ans += "\nTriangles: \n"
    for i,T in enumerate(self.t):
      ans += str(i) + ": " + str(T) + "\n"
      ans += "covered by: " + str(self.t_lifts[i]) + "\n"
    ans += "Covering vertices:\n"
    for i,CV in enumerate(self.em_v):
      ans += str(i) + ": " + str(self.em_v[i]) + "\n"
    ans += "Covering edges:\n"
    for i,CE in enumerate(self.em_e):
      ans += str(i) + ": " + str(self.em_e[i]) + "\n"
    ans += "Covering triangles:\n"
    for i,CT in enumerate(self.em_t):
      ans += str(i) + ": " + str(self.em_t[i]) + "\n"
    return ans
    
    