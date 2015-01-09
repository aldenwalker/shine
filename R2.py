import math
import copy
import numpy as np

class Vector:
  def __init__(self, L):
    if isinstance(L, int):
      self.x = L*[0]
    else:
      self.x = [i for i in L]
    
  
  def __sub__(self, other):
    return Vector([self.x[i]-other.x[i] for i in xrange(len(self.x))])
  
  def __add__(self, other):
    return Vector([self.x[i]+other.x[i] for i in xrange(len(self.x))])
  
  def __mul__(self, a):
    return Vector([a*self.x[i] for i in xrange(2)])  
  
  def __div__(self, a):
    return Vector([self.x[i]/a for i in xrange(2)])
  
  def __getitem__(self, i):
    return self.x[i]
  
  def __repr__(self):
    return "Vector(" + str(self.x) + ")"
  
  def __str__(self):
    return str(self.x)
  
  def __len__(self):
    return len(self.x)

  def norm(self):
    return math.sqrt(sum([v**2 for v in self.x]))
  
  def scaled_to_len(self, ell):
    n = self.norm()
    return Vector( [(ell/n)*v for v in self.x] )
  
  def dot(self, other):
    return sum([self.x[i]*other.x[i] for i in xrange(len(self.x))])
  
  def cross(self, other):
    return self.x[0]*other.x[1]-self.x[1]*other.x[0]
  
  def copy(self):
    return Vector([i for i in self.x])
  
  def angle(self, other):
    try:
      return math.acos( self.dot(other) / (self.norm()*other.norm()) )
    except ValueError:
      if abs( -1 - (self.dot(other) / (self.norm()*other.norm())) ) < 1e-10:
        return math.pi
      else:
        raise ValueError("math domain error")


def along_segment(s, t):
  return s[1]*t + s[0]*(1-t)

def intersect_segments(s1, s2):
  """return (t,b) where t is the fraction along s1 that it intersects
  s2, and b is a boolean that says whether s2, s1 direction is a 
  positive basis"""
  #print "Intersecting ", s1, s2
  s1x1, s1y1 = s1[0]
  s1x2, s1y2 = s1[1]
  s2x1, s2y1 = s2[0]
  s2x2, s2y2 = s2[1]
  t1x, t1y = [s1x2-s1x1, s1y2-s1y1]
  t2x, t2y = [s2x2-s2x1, s2y2-s2y1]
  RHSx, RHSy = [s2x1-s1x1, s2y1-s1y1]
  a,b,c,d = t1x, -t2x, t1y, -t2y
  det = a*d-b*c
  if det == 0:
    #print "det 0"
    return None
  Ia, Ib, Ic, Id = d/det, -b/det, -c/det, a/det
  ansx = Ia*RHSx + Ib*RHSy
  ansy = Ic*RHSx + Id*RHSy
  if 0.0-1e-10 <= ansx and ansx <= 1.0+1e-10 and 0.0-1e-10 <= ansy and ansy <= 1.0+1e-10:
    return (ansx, det > 0)
  else:
    return None

###########################################################################
# fit a single cubic bezier curve to data.  Each point should be of the 
# form (float t, R2.vector v), where the curve wants to be at v at time t
# If fix_endpoints==True, the first and last points must have t=0,1
###########################################################################
def cubic_bezier_fit(pts, fix_endpoints=True):
  #if we're fixing the endpoints, we need to subtract off 
  #those constant quantities from the data points
  pts = list(pts)
  lpts = len(pts)
  if lpts == 3:
    return pts[0][1], pts[1][1], pts[1][1], pts[2][1]
  if fix_endpoints:
    p1 = pts[0][1]
    p2 = pts[-1][1]
    for i in xrange(1,lpts-1):
      t,v = pts[i]
      x_amount = (1-t)**3*p1[0] + t**3*p2[0]
      y_amount = (1-t)**3*p1[1] + t**3*p2[1]
      pts[i] = (t, Vector([v[0]-x_amount, v[1]-y_amount]))
    M = np.matrix([[3*(1-t)**2*t, 3*(1-t)*t**2] for t,v in pts[1:-1]])
    RHS = np.matrix([[v[0], v[1]] for t,v in pts[1:-1]])
    print "Points: ", pts
    print pts
    print M
    print RHS
    MT = M.getT()
    ans = (MT*M).getI()*MT
    print ans.shape, RHS.shape
    ans = ans*RHS
    ans = ans.tolist()
    cp1 = Vector(ans[0])
    cp2 = Vector(ans[1])
  else:
    M = [[(1-t)**3, 3*(1-t)**2*t, 3*(1-t)*t**2, t**3] for t,v in pts]
    RHS = [[v[0], v[1]] for t,v in pts]
    MT = M.getT()
    ans = (MT*M).getI()*MT*RHS
    ans = ans.tolist()
    p1, cp1, cp2, p2 = [Vector(a) for a in ans]
  return p1, cp1, cp2, p2

############################################################################
# fit a succession of cubic bezier curves to the sequence of points
# note the points are traversed at unit speed (equal time is given
# to each gap).  The endpoints are matched exactly
#
# the return result is of the form [[p11, p12, p13, p14],[p21,p22,p23,p24]...]
# where it is guaranteed that pn4=p(n+1)1 and that p(n+1)2 is the 
# reflection of pn3 through pn4
#
# it will make sure num_breakpoints is large enough so that at least 4 points 
# are in every interval
############################################################################
def fit_polybezier(pts, num_breakpoints=None):
  n = len(pts)
  if n<4:
    if n==1:
      return [[pts[0], pts[0], pts[0], pts[0]]]
    elif n==2:
      return [[pts[0], (pts[0]+pts[1])*0.5, (pts[0]+pts[1])*0.5, pts[1]]]
    else: # n==3
      return [[pts[0], pts[1], pts[1], pts[2]]]
  if num_breakpoints == None or 3*num_breakpoints > n-4:
    num_breakpoints = ((n-4) / 3)+1
  M_cols = 2*num_breakpoints + 2
  #each gap is defined to have length 1
  #each cubic bezier therefore covers this much distance:
  distance_per_curve = float(n-1)/float(num_breakpoints+1)
  #we'll have to provide a scaled t values for each point
  M = [M_cols * [0] for i in xrange(n-2)]
  RHS = [ [0,0] for i in xrange(n-2)]
  #print "n:",n
  #print "num_breakpoints:",num_breakpoints
  #print "distance_per_curve:",distance_per_curve
  for i in xrange(1,n-1):
    this_point_dist = float(i)
    breakpoint_index = int(this_point_dist / distance_per_curve)
    if breakpoint_index > num_breakpoints:
      breakpoint_index -= 1
    curve_start_dist = breakpoint_index*distance_per_curve
    t = (this_point_dist - curve_start_dist) / distance_per_curve
    row_ind = i-1
    #print "Doing point index ", i
    #print "t value:", t
    #print "breakpoint index:", breakpoint_index
    if breakpoint_index == 0:
      M[row_ind][0] = 3*(1-t)**2*t
      M[row_ind][1] = 3*(1-t)*t**2
      RHS[row_ind][0] = pts[i][0] - (1-t)**3*pts[0][0]
      RHS[row_ind][1] = pts[i][1] - (1-t)**3*pts[0][1]
      if num_breakpoints > 0:
        M[row_ind][2] = t**3
      else:
        RHS[row_ind][0] -= t**3*pts[-1][0]
        RHS[row_ind][1] -= t**3*pts[-1][1]
    else:
      variable_index_offset = 1 + (breakpoint_index-1)*2
      vio = variable_index_offset
      M[row_ind][vio] = -3*(1-t)**2*t
      M[row_ind][vio+1] = (1-t)**3 + 6*(1-t)**2*t
      M[row_ind][vio+2] = 3*(1-t)*t**2
      RHS[row_ind][0] = pts[i][0]
      RHS[row_ind][1] = pts[i][1]
      if breakpoint_index < num_breakpoints:
        M[row_ind][vio+3] = t**3
      else:
        RHS[row_ind][0] -= t**3*pts[-1][0]
        RHS[row_ind][1] -= t**3*pts[-1][1]
  M = np.matrix(M)
  RHS = np.matrix(RHS)
  MT = M.getT()
  #print M
  #print RHS
  #print MT*M
  ans = (MT*M).getI()*MT*RHS
  ans = ans.tolist()
  points = [[pts[0], Vector(ans[0])] ]
  for i in xrange(num_breakpoints):
    cp2 = Vector(ans[1+2*i])
    p2 = Vector(ans[1+2*i+1])
    points[-1].extend( [cp2, p2] )
    reflection = p2 + (p2-cp2)
    points.append( [ p2, reflection ] )
  points[-1].extend( [Vector(ans[-1]), pts[-1]] )
  return points
    
########################################################################
# fit a polybezier, but figure out how many control points to use 
# and how much to subdivide
######################################################################### 
def bezier_approximation(pts):
  if len(pts) < 3:
    return fit_polybezier(pts)
  new_pts = list(pts)
  #clean the points
  i=0
  while i<len(new_pts)-1:
    if (new_pts[(i+1)%len(new_pts)]-new_pts[i]).norm() < 1e-8:
      del new_pts[i+1]
    else:
      i += 1
  pts = new_pts  
  n = len(pts)
  if n < 3:
    return fit_polybezier(pts)
  total_angle = sum([(pts[i]-pts[i-1]).angle(pts[i+1]-pts[i]) for i in xrange(1,n-1)])
  #each pi of angle should give about one breakpoint
  num_breakpoints = int(1.5*(total_angle / math.pi)+1)
  #each breakpoint should have about 10 points
  n_needed = 8*num_breakpoints + 8
  cut_number = int(n_needed / float(n))+1
  new_pts = []
  for i in xrange(n-1):
    for j in xrange(cut_number):
      t = j/float(cut_number)
      new_pts.append( pts[i]*(1-t) + pts[i+1]*t )
  new_pts.append( pts[-1] )
  #print "total angle:", total_angle
  #print "num breakpoints:", num_breakpoints
  #print "n_needed:", n_needed
  #print "n now:", n
  #print "cut number:", cut_number
  return fit_polybezier(new_pts, num_breakpoints=num_breakpoints)
  
  
 
  
def bezier_approximation_cubic(pts, fix_endpoints=True):
  n = len(pts)
  if n < 3:
    return None
  denom = float(n-1)
  augmented_pts = [(i/denom, v) for i,v in enumerate(pts)]
  return cubic_bezier_fit(augmented_pts, fix_endpoints)












  