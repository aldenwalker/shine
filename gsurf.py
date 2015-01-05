import math
import random
import scipy.optimize
import copy

import hyp
import mobius
import tsurf
import emsurf
from signedind import SignedInd as SI

##############################################################################
# Geometric (hyperbolic) surface
##############################################################################
class GeometricSurface(tsurf.TopSurface):
  def __init__(self, TS, lens):
    self.v = copy.deepcopy(TS.v)
    self.e = copy.deepcopy(TS.e)
    self.t = copy.deepcopy(TS.t)
    self.loops = dict()
    self.h_lengths = lens
    self.h_tris = [hyp.HypTri([lens[ei.ind] for ei in T.i_edges]) for T in self.t]
  
  ##########################################################################
  # find a hyperbolic structure on a topological surface
  ##########################################################################
  @classmethod
  def geometrize_tsurf(cls, TS, verbose=0):
    if TS.euler_char() >= 0:
      print "Only hyperbolic supported"
      return
    num_edges = len(TS.e)
    num_verts = len(TS.v)
    num_tris = len(TS.t)
    
    #optimization inputs
    #edge lengths
    #    def obj_fun(x):
    #      return sum([(x[i]-desired_edge_lengths[i])**2 for i in xrange(len(x))])
    #    def obj_fun_grad(x):
    #      return [2*(x[i]-desired_edge_lengths[i]) for i in xrange(len(x))]
    #    def f(x):
    #      return [obj_fun(x), obj_fun_grad(x)]
    
    #angle deviations
    def obj_fun(x):
      ans = 0
      v_vals = [float(len(v.i_edges)) for v in TS.v]
      for i,t in enumerate(TS.t):
        Tlens = [x[ei.ind] for ei in t.i_edges]
        Tangles = hyp.tri_angles(Tlens)
        for j in xrange(3):
          ans += (Tangles[j] - (6*math.pi/v_vals[t.i_verts[j][0]]))**2
      return ans
    def obj_fun_grad(x):
      ans = [0 for _ in xrange(len(x))]
      v_targets = [6*math.pi/float(len(v.i_edges)) for v in TS.v]
      for i,t in enumerate(TS.t):
        Tlens = [x[ei.ind] for ei in t.i_edges]
        Tangles = hyp.tri_angles(Tlens)
        eis = t.i_edges
        for j in xrange(3):
          for k in xrange(3):
            ans[eis[k].ind] += 2*(Tangles[j] - v_targets[t.i_verts[j][0]])*hyp.tri_angle_deriv(Tlens, j, k) 
      return ans
    def f(x):
      return [obj_fun(x), obj_fun_grad(x)]



    cons = [dict() for _ in xrange(num_verts+3*num_tris )]
    #constraints from vertices (angles = 0)
    for i,v in enumerate(TS.v):
      def this_func(x,v=v):
        ans = -2*math.pi
        for ti,j in v.i_tris:
          ans += hyp.tri_angle( [x[ei.ind] for ei in TS.t[ti].i_edges], j )
        return ans
      def this_jac(x,v=v):
        ans = [0 for _ in xrange(len(x))]
        incident_faces = len(v.i_tris)
        for ti,j in v.i_tris:   #ti is the triangle index, j is the index of this angle in the triangle
          T = TS.t[ti]
          Tlens = [x[ei.ind] for ei in T.i_edges]
          for k,ei in enumerate(T.i_edges):  #k is the side index in the triangle, ei is the edge index
            ans[ei.ind] += hyp.tri_angle_deriv( Tlens, j, k )
        return ans
      cons[i]['type'] = 'eq'
      cons[i]['fun'] = this_func
      cons[i]['jac'] = this_jac
    #constraints from faces (triangle inequalities)
    for j in xrange(num_tris):
      Te = TS.t[j].i_edges
      for k in xrange(3):
        def this_func(x,Te=Te,k=k):
          return -x[Te[k].ind] + x[Te[(k+1)%3].ind] + x[Te[(k+2)%3].ind] - 1e-5
        def this_jac(x,Te=Te,k=k):
          ans = [0 for _ in xrange(num_edges)]
          ans[Te[k].ind], ans[Te[(k+1)%3].ind], ans[Te[(k+2)%3].ind] = -1, 1, 1
          return ans
        cons[num_verts + 3*j + k]['type'] = 'ineq'
        cons[num_verts + 3*j + k]['fun'] = this_func
        cons[num_verts + 3*j + k]['jac'] = this_jac
    
    x0 = [1 for k in xrange(num_edges)]
    bounds = [(0,None) for i in xrange(num_edges)]
    res = scipy.optimize.minimize(f,                              \
                                  x0,           \
                                  jac=True,                       \
                                  bounds=bounds, \
                                  constraints=cons,               \
                                  method='SLSQP',                 \
                                  options={'disp':True,'iprint':verbose, 'maxiter':200})
    if not res.success:
      print "Failed to find structure"
      raise ValueError("Couldn't find structure")
      return None
    if verbose>0:
      print "Found structure"
    gs = cls(TS, res.x)
    return gs
  
  ########################################################################
  # subdivide a geometric surface
  ########################################################################
  def subdivide(self):
    old_nv = len(self.v)
    old_ne = len(self.e)
    old_nt = len(self.t)
    old_TS, vertices_from_edges, edges_from_edges, edges_from_tris, tris_from_tris = super(GeometricSurface, self).subdivide()
    
    self.old_h_lengths = self.h_lengths
    self.h_lengths = len(self.e)*[None]
    self.old_h_tris = self.h_tris
    self.h_tris = len(self.t)*[None]
    
    #divide each of the previous edges in two
    for i in xrange(old_ne):
      e1i, e2i = edges_from_edges[i]
      self.h_lengths[e1i] = self.old_h_lengths[i]/2.0
      self.h_lengths[e2i] = self.old_h_lengths[i]/2.0
    
    #the remaining edges will have their lengths fixed
    for i in xrange(old_nt):
      for j in xrange(3):
        ell1 = self.old_h_tris[i].lengths[(j-1)%3] / 2.0
        angle = self.old_h_tris[i].angles[j]
        ell2 = self.old_h_tris[i].lengths[j] / 2.0
        new_edge_i = edges_from_tris[i][j]
        self.h_lengths[new_edge_i] = hyp.tri_opposite_length(ell1, angle, ell2)
    
    if None in self.h_lengths:
      print "Haven't found all the edges?"
      raise ValueError("yo moma")
    
    #now the triangles get their lengths from the edges
    for i in xrange(len(self.t)):
      self.h_tris[i] = hyp.HypTri( [ self.h_lengths[ self.t[i].i_edges[j].ind ] for j in xrange(3) ] )

    return old_TS, vertices_from_edges, edges_from_edges, edges_from_tris, tris_from_tris
  
  ########################################################################
  # given three edges and coordinates (in [0,1]) along them, compute
  # the angle at the middle edge (this is the external angle, so 0 means
  # that the edge is straight)
  ########################################################################
  def angle_at_edge(self, e1, e1_t, e2, e2_t, e3, e3_t):
    #find out the exiting angle from the first triangle
    t1i, t1_e1 = (self.e[e1.ind].on_left if e1.sign>0 else self.e[e1.ind].on_right)
    t1 = self.h_tris[t1i]
    t1_e1_t = (e1_t if e1.sign>0 else 1-e1_t)
    t1_e2 = (self.e[e2.ind].on_right if e2.sign>0 else self.e[e2.ind].on_left)[1]
    t1_e2_t = (1-e2_t if e2.sign>0 else e2_t)
    
    a1 = t1.inside_final_angle(t1_e1, t1_e1_t, t1_e2, t1_e2_t)
    
    t2i, t2_e2 = (self.e[e2.ind].on_left if e2.sign>0 else self.e[e2.ind].on_right)
    t2 = self.h_tris[t2i]
    t2_e2_t = (e2_t if e2.sign>0 else 1-e2_t)
    t2_e3 = (self.e[e3.ind].on_right if e3.sign>0 else self.e[e3.ind].on_left)[1]
    t2_e3_t = (1-e3_t if e3.sign>0 else e3_t)
    
    a2 = t2.inside_final_angle(t2_e3, t2_e3_t, t2_e2, t2_e2_t) #note backwards
    
    return a2 - a1
    
  ########################################################################
  # return the polygon (group of triangles) which surround the vertex
  # opposite the given edge
  ########################################################################
  def polygon_with_entrance_edge(self, ei):
    t0, ei_in_t0 = (self.e[ei.ind].on_left if ei.sign>0 else self.e[ei.ind].on_right)
    vi, t0i_in_v = self.t[t0].i_verts[(ei_in_t0+2)%3]
    tris = self.v[vi].i_tris[t0i_in_v:] + self.v[vi].i_tris[:t0i_in_v]
    num_sides = len(tris)
    eedges = num_sides*[None]
    angles = num_sides*[0]
    sides = num_sides*[None]
    for i,(ti,tii) in enumerate(tris):
      eedges[i] = -self.t[ti].i_edges[(tii+1)%3] 
      sides[i] = self.h_lengths[eedges[i].ind]
      angles[i] += self.h_tris[ti].angles[(tii+1)%3] 
      angles[(i+1)%num_sides] += self.h_tris[ti].angles[(tii+2)%3]
    return GeometricSurfacePolygon(self, vi, tris, eedges, sides, angles)
      
  ########################################################################
  # find the maximal angle deviation along an embedded path
  ########################################################################
  def maximal_angle_deviation(self, EP):
    ma = 0
    lep = len(EP.edges)
    for i in xrange(lep):
      angle = self.angle_at_edge( EP.edges[(i-1)%lep], EP.edge_coords[(i-1)%lep], \
                                  EP.edges[i], EP.edge_coords[i], \
                                  EP.edges[(i+1)%lep], EP.edge_coords[(i+1)%lep] )
      if abs(angle) > ma:
        ma = abs(angle)
    return ma
  
  ########################################################################
  # given an embedded path, find the geodesic geometric path
  # if it's just a topological path, it'll convert it to a default path first
  ########################################################################
  def geodesicify(self, TP, tol=None):
  
    #this tolerance says how big an angle is considered "straight"
    if tol==None:
      tol = 1e-3
      
    #make a default embedded path, if necessary
    new_TP = copy.deepcopy(TP)
    if isinstance(TP, emsurf.EmbeddedPath):
      #print "Already an embedded path; no need to default it"
      EP = new_TP
    else:
      #print "Default pathing topological path"
      EP = emsurf.EmbeddedPath.from_topological_path(new_TP)
    
    #for each edge, do the following:
    # - simplify the loop (remove edge followed by opposite edge)
    # - find a location at which there is a nonzero angle (if there 
    #   is no such location, then the loop is a geodesic)
    # - our edge to focus on will be the one *before* the kink
    # - find the vertex opposite from the edge
    # - build the hyperbolic polygon we get by tiling around this vertex
    # - find out where the loop enters and exits this polygon
    # - replace the loop inside the polygon with the geodesic between the 
    #   entry and exit points
    
    #it's slightly complicated because the number of edges can change!
    
    #print "Geodesicifying embedded path: ", EP
    
    EP.simplify()
    
    #print "Simplified: ", EP
    
    current_EP_index = 0
    
    while True:
      lep = len(EP.edges)
      
      #print "Current EP", EP
      #print "EP index:", current_EP_index
      
      # #find an edge with the kink
#       print "Looking for edge with a kink"
#       kinked_edge_i = None
#       for i in xrange(lep):
#         angle = self.angle_at_edge( EP.edges[(i-1)%lep], EP.edge_coords[(i-1)%lep], \
#                                     EP.edges[i], EP.edge_coords[i], \
#                                     EP.edges[(i+1)%lep], EP.edge_coords[(i+1)%lep] )
#         print "Found angle", angle, "at edge", i, "i.e.", EP.edges[i]
#         if abs(angle) > tol:
#           kinked_edge_i = i
#           break
#       if kinked_edge_i == None:
#         break
#       
#       #step back to the previous edge
#       poly_entrance_edge_i = (kinked_edge_i-1)%len(EP.edges)
      
      poly_entrance_edge_i = current_EP_index
      
      #build the polygon
      HP = self.polygon_with_entrance_edge(EP.edges[poly_entrance_edge_i])
      
      #print "Built the polygon\n", HP
      
      #find the exit edge
      poly_exit_edge_i = poly_entrance_edge_i
      while EP.edges[poly_exit_edge_i] not in HP.exit_edges:
        poly_exit_edge_i = (poly_exit_edge_i + 1 if poly_exit_edge_i<len(EP.edges)-1 else 0)
      
      #print "Found the entrance and exit edge indices", poly_entrance_edge_i, " ", poly_exit_edge_i
      #print "Which is edges", EP.edges[poly_entrance_edge_i], "and", EP.edges[poly_exit_edge_i]
      
      #find the list of edges in the straightened geodesic
      replacement_edges, replacement_edge_coords = HP.geodesic_path(EP.edges[poly_entrance_edge_i], \
                                                                    EP.edge_coords[poly_entrance_edge_i], \
                                                                    EP.edges[poly_exit_edge_i], \
                                                                    EP.edge_coords[poly_exit_edge_i])
      
      #print "Found replacement edges: ", replacement_edges
      #print "And coords: ", replacement_edge_coords
      
      #delete the interior edges and replace them; note there 
      #are two cases depending on whether we overlap the end of the list
      if poly_entrance_edge_i < poly_exit_edge_i: #we don't overlap the end
        peei = poly_entrance_edge_i
        del EP.edges[peei+1:poly_exit_edge_i]
        del EP.edge_coords[peei+1:poly_exit_edge_i]
        EP.edges[peei+1:peei+1] = replacement_edges
        EP.edge_coords[peei+1:peei+1] = replacement_edge_coords
      else: #we do overlap
        del EP.edges[poly_entrance_edge_i+1:]
        del EP.edges[:poly_exit_edge_i]
        del EP.edge_coords[poly_entrance_edge_i+1:]
        del EP.edge_coords[:poly_exit_edge_i]
        EP.edges.extend(replacement_edges)
        EP.edge_coords.extend(replacement_edge_coords)
      
      current_EP_index = poly_entrance_edge_i+1
      if current_EP_index >= len(EP.edges):
        current_EP_index = 0
        #check to see if we are done
        mad = self.maximal_angle_deviation(EP)
        print "Maximum angle deviation:", mad
        if mad < tol:
          break
        
    #########
    return EP
  
  #########################################################################
  # Print a geometric surface
  #########################################################################
  def __repr__(self):
    return str(self)
  
  def __str__(self):
    ans = "Geometric surface:\n"
    ans += "Vertices: \n"
    for i,V in enumerate(self.v):
      ans += str(i) + ": " + str(V) + "\n"
    ans += "\nEdges: \n"
    for i,E in enumerate(self.e):
      ans += str(i) + ": " + str(E) + " Length: " + str(self.h_lengths[i]) + "\n"
    ans += "\nTriangles: \n"
    for i,T in enumerate(self.t):
      ans += str(i) + ": " + str(T) + "\n  " + str(self.h_tris[i]) + "\n"
    return ans




############################################################################
# this is a polygon in a surface (a bunch of triangles), and it 
# builds the model polygon
############################################################################
class GeometricSurfacePolygon:
  def __init__(self, GS, vi, tris, eedges, sides, angles):
    self.GS = GS
    self.central_vert = vi
    self.tris = tris
    self.exit_edges = eedges
    self.HP = hyp.HypPolygon(sides, angles)
  
  ########################################################################
  # return the edges and edge coords of the geodesic which goes from edge
  # e1 at e1t to edge e2 at e2t.  Note the edges are edge indices in the 
  # main surface, so we have to search for them
  ########################################################################
  def geodesic_path(self, e1, e1t, e2, e2t):
    nv = len(self.tris)
    enter_index = self.exit_edges.index(-e1)
    exit_index = self.exit_edges.index(e2)
    enter_t = (e1t if e1.sign>0 else 1-e1t)
    exit_t = (e2t if e2.sign<0 else 1-e2t)
    enter_angle = self.HP.enter_angle(enter_index, enter_t, exit_index, exit_t)
    
    #print "Found the entrance angle: ", enter_angle
    
    #now walk through the polygon
    added_edges = []
    added_edge_coords = []
    final_tri = self.tris[exit_index][0] #if we get this, we are done
    current_edge = e1
    current_edge_t = e1t
    current_angle = enter_angle
    current_tri = self.tris[enter_index][0]
    while current_tri != final_tri:
      #print "edge, t, angle, tri =", current_edge, current_edge_t, current_angle, current_tri
    
      ind_in_tri = (self.GS.e[current_edge.ind].on_left if current_edge.sign>0 else \
                    self.GS.e[current_edge.ind].on_right)[1]
      t_in_tri = (current_edge_t if current_edge.sign>0 else 1-current_edge_t)
      #print "ind in tri, t in tri:", ind_in_tri, t_in_tri
      exit_i_in_tri, exit_t, exit_angle = self.GS.h_tris[current_tri].exit_t_and_angle(ind_in_tri, t_in_tri, current_angle)
      #print "Exit: ", exit_i_in_tri, exit_t, exit_angle
      next_edge = -self.GS.t[current_tri].i_edges[exit_i_in_tri]
      added_edges.append(next_edge)
      added_edge_coords.append( (1-exit_t if next_edge.sign>0 else exit_t) )
      ###
      current_edge = next_edge
      current_edge_t = (1-exit_t if next_edge.sign>0 else exit_t)
      current_angle = exit_angle
      current_tri = (self.GS.e[next_edge.ind].on_left if next_edge.sign>0 else \
                     self.GS.e[next_edge.ind].on_right)[0]
      if current_edge in self.exit_edges:
        print "We've found an exit edge without hitting the edge we were supposed to!"
        return 
    return added_edges, added_edge_coords

  ###########################################################################
  # print it out
  ###########################################################################
  def __repr__(self):
    return str(self)
  def __str__(self):
    ans = ''
    ans += "Polygon with " + str(len(self.tris)) + " sides and triangles: " + str(self.tris) + '\n'
    ans += 'The exit edges are ' + str(self.exit_edges) + "\n"
    ans += 'Around vertex ' + str(self.central_vert) + "\n"
    ans += 'The hyperbolic polygon is: \n' + str(self.HP)
    return ans

#############################################################################
# Realized vertices, edges, triangles in hyperbolic space which cover a 
# geometric surface
#############################################################################
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


############################################################################
# A cover of a geometric surface
############################################################################
class LiftedSurface(GeometricSurface):
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
  
  ##########################################################################
  # lift a geometric surface
  ##########################################################################
  @classmethod
  def lift_gsurf(cls, GS):
    LS = cls(GS, [], [], [], [], [], [])
    LS.relay()
    return LS
  
  ##########################################################################
  # for a vertex index downstairs, are there lifts of every adjacent triangle?
  ##########################################################################
  def is_vertex_surrounded(self, vi):
    """are all adjacent triangles placed?"""
    return all([len(self.t_lifts[ti])>0 for ti,j in self.v[vi].i_tris])
  
  #########################################################################
  # for a signed lifted edge index, lift the appropriate triangle so that 
  #the given (signed index) edge now has a triangle attached
  #########################################################################
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

  ##########################################################################
  # go back to the underlying geometric surface and lift it again
  ##########################################################################
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
  
  ##########################################################################
  # hit everything by a mobius transformation
  ##########################################################################
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
  
  ##########################################################################
  # given a topological path in the surface, find a 
  # geodesic representative of the homotopy class
  ##########################################################################
  def geodesicify(self, TP):
    #find an acceptable lift of the first edge
    lifted_edges = []
    for pli in self.e_lifts[TP.edges[0].ind]:
      pl = self.em_e[pli]
      if (pl.on_left if TP.edges[0].sign>0 else pl.on_right) != None:
        lifted_edges.append( SI(pli, TP.edges[0].sign) )
        break
    #now go through, lifting as necessary
    #when we look at the ith index, we lift the *next* edge
    #this includes the last edge, which is all we care about
    ne = len(TP.edges)
    for i in xrange(ne):
      e1 = TP.edges[i]
      lifted_e1 = lifted_edges[i]
      e2 = TP.edges[(i+1)%ne]
      tri_ind = (self.em_e[lifted_e1.ind].on_left if e1.sign>0 else self.em_e[lifted_e1.ind].on_right)
      if tri_ind == None:
        print "Something is wrong -- should be a triangle here"
      tri = self.em_t[tri_ind[0]]
      next_edge_ind_in_tri = (self.e[e2.ind].on_right if e2.sign>0 else self.e[e2.ind].on_left)[1]
      next_edge = -tri.i_edges[next_edge_ind_in_tri]
      if next_edge.sign != e2.sign:
        print "Something is wrong -- edges should have the opposite sign"
      next_tri_ind = (self.em_e[next_edge.ind].on_left if e2.sign>0 else self.em_e[next_edge.ind].on_right)
      if next_tri_ind == None:
        self.lift_triangle_to_lifted_edge(next_edge)
      lifted_edges.append(next_edge)
    print "Should take", lifted_edges[0], "to", lifted_edges[-1]
    print "i.e.", self.em_e[lifted_edges[0].ind], "to", self.em_e[lifted_edges[-1].ind]
    GI1 = self.em_e[lifted_edges[0].ind].gi
    GI2 = self.em_e[lifted_edges[-1].ind].gi
    M = mobius.MobiusTrans.unit_tangent_action(GI1.start, GI1.initial_angle, GI2.start, GI2.initial_angle)
    print "Got mobius: ", M
    GA = M.geodesic_axis()
    print "Got geodesic axis: ", GA
    return TP
  
  ##########################################################################
  # print out a lifted surface
  ##########################################################################
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
    











