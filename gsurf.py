import scipy.optimization

import tsurf

PI = 3.14159265358979323846


class Structure:
  def __init__(self, T):
    self.T = T
    self.e = [None for i in xrange(len(T.e))]
    self.f = [[None for i in xrange(len(F.i_verts))] for F in T.f]
    self.has_structure = False
  
  def excess_angle_hyp(self, x):
    """Returns the vector of excess angles over all vertices, assuming that 
       the edge lengths are x"""
    nv = len(T.v)
    ans = [-2*PI for i in xrange(nv)]
    for F in T.f:
      edge_lengths = [x[ei.ind] for ei in F.i_edges] 
      A = hyperbolic_angles(edge_lengths)
      for i,a in enumerate(A):
        ans[ F.i_verts[i][0] ] += a
    return ans
  
  def excess_angle_hyp_jac(self, x):
    
  
  def find_structure(self, requested_edge_ratios=None):
    """Attemps to find a geometric structure on the triangulation"""
    chi = T.euler_char()
    if chi >= 0:
      print "Surface must be hyperbolic"
      return
    


class GeoSurface(tsurf.TopSurface):
  