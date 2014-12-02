import math
import scipy.optimize

import models
import tsurf

class GeoSurface(tsurf.TopSurface):
  def __init__(self, TS, lens):
    self.v = TS.v
    self.e = TS.e
    self.f = TS.f
    self.h_lengths = lens
    self.h_tris = [models.HypTri([lens[ei.ind] for ei in F.i_edges]) for F in self.f]
  
  @classmethod
  def geometrize_tsurf(cls, TS, edge_hints=None):
    num_edges = len(TS.e)
    if edge_hints == None:
      desired_edge_lengths = [1 for i in xrange(num_edges)]
    else:
      desired_edge_lengths = [eh for eh in edge_hints]
    
    #optimization inputs
    def obj_fun(x):
      return sum([(x[i]-desired_edge_lengths[i])**2 for i in xrange(len(x))])
    def obj_fun_grad(x):
      return [2*(x[i]-desired_edge_lengths[i]) for i in xrange(len(x))]
    def f(x):
      return [obj_fun(x), obj_fun_grad(x)]

    cons = [dict() for v in TS.v]
    for i,v in enumerate(TS.v):
      def this_func(x):
        ans = -2*math.pi
        for fi,j in v.i_faces:
          ans += models.hyp_tri_angle( [x[ei.ind] for ei in TS.f[fi].i_edges], j )
        return ans
      def this_jac(x):
        ans = [0 for _ in xrange(len(x))]
        inc_edges = len(v.i_edges)
        for j in xrange(inc_edges):
          f_before = v.i_faces[(j-1)%len(v.i_faces)]
          f_before_lengths = [x[ei.ind] for ei in TS.f[f_before[0]]]
          f_after = v.i_faces[j]
          f_after_lengths = [x[ei.ind] for ei in TS.f[f_after[0]]]
          ans[v.i_edges[j]] += models.hyp_tri_angle_deriv( f_before_lengths, f_before[1], (f_before[1]-1)%3 )
          ans[v.i_edges[j]] += models.hyp_tri_angle_deriv( f_after_lengths, f_after[1], (f_after[1]+1)%3 )
        return ans
      cons[i]['type'] = 'eq'
      cons[i]['fun'] = this_func
      cons[i]['jac'] = this_jac
    
    res = scipy.optimize.minimize(f,                              \
                                  desired_edge_lengths,           \
                                  jac=True,                       \
                                  bounds=[(0,None) for i in xrange(num_edges)], \
                                  constraints=cons,               \
                                  method='SLSQP',                 \
                                  )#disp=True )
    if not res.success:
      print "Failed to find structure"
      return None
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
    for i,F in enumerate(self.f):
      ans += str(i) + ": " + str(F) + " HypTri: " + str(self.h_tris[i]) + "\n"
    return ans
    