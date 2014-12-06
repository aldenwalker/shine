import math
import random
import scipy.optimize

import hyp
import tsurf

class GeoSurface(tsurf.TopSurface):
  def __init__(self, TS, lens):
    self.v = TS.v
    self.e = TS.e
    self.t = TS.t
    self.h_lengths = lens
    self.h_tris = [hyp.HypTri([lens[ei.ind] for ei in T.i_edges]) for T in self.t]
  
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
    