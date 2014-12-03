import math
import mobius

def same_float(x,y, tol=1e-10):
  return abs(x-y) < tol

def rough_acos(x,cutoff=1e-5):
  try:
    return math.acos(x)
  except ValueError:
    if x > 1.0 and x-1.0 < cutoff:
      return 0.0
    elif x < -1.0 and -1.0-x < cutoff:
      return math.pi/2.0
    else:
      raise ValueError('math domain error')

def dist(v1,v2):
  x = 1 + (((v2.real-v1.real)**2 + (v2.imag-v1.imag)**2)/(2*v1.imag*v2.imag))
  return math.acosh(x)

def point_at_vertical_distance(x, d):
  """return the point obtained by going vertically distance d from x"""
  a = 1
  b = -2*x.imag*math.cosh(d)
  c = x.imag**2
  height = (-b + math.sqrt(b**2-4*a*c))/(2.0*a)
  return complex(x.real, height)

def point_angle_dist(p, A, d):
  """return the point obtained by going distance d from point p at angle A"""
  vp = point_at_vertical_distance(p, d)
  if same_float(A, math.pi/2.0):
    return vp
  M = mobius.MobiusTrans.unit_tangent_action(p, math.pi/2, p, A)
  return M(vp)

def tri_angles(L):
  ans = [0, 0, 0]
  for i in xrange(3):
    top = math.cosh(L[i])*math.cosh(L[(i+2)%3]) - math.cosh( L[(i+1)%3] )
    bottom = math.sinh(L[i])*math.sinh(L[(i+2)%3])
    ans[i] = math.acos(top/bottom)
  return ans

def tri_angle(L, ind):
  """return the angle at position ind"""
  top = math.cosh(L[ind])*math.cosh(L[(ind+2)%3]) - math.cosh( L[(ind+1)%3] )
  bottom = math.sinh(L[ind])*math.sinh(L[(ind+2)%3])
  r = top/bottom
  return rough_acos(r)

def tri_angle_deriv(L, i, j):
  """compute the derivative of angle i with respect to side length j"""
  a = L[i]
  b = L[(i+2)%3]
  c = L[(i+1)%3]
  cscha = 1/math.sinh(a)
  cschb = 1/math.sinh(b)
  cothb = math.cosh(b)/math.sinh(b)
  cotha = math.cosh(a)/math.sinh(a)
  coshc = math.cosh(c)
  den = math.sqrt(1- (cotha*cothb - coshc*cscha*cschb)**2)
  if i==j:           #deriv of angle with respect to a
    num = cscha*(cothb*cscha - coshc*cotha*cschb)
  elif j==(i+2)%3:    #deriv of angle with respect to b
    num = cschb*(cotha*cschb - coshc*cothb*cscha)
  else:              #deriv of angle with respect to c
    num = cscha*cschb*math.sinh(c)
  return num/den

class HypGeodesicInterval:
  def __init__(self, v1, v2):
    self.start = v1
    self.end = v2
    self.length = dist(v1,v2)
    if same_float(v1.real, v2.real):
      self.vertical = True
      self.circ_center = v1.real
      self.circ_radius = 'inf'
      self.circ_angle1 = v1.imag
      self.circ_angle2 = v2.imag
      if v2.imag > v1.imag:
        self.initial_angle = math.pi/2.0
        self.final_angle = math.pi/2.0
      else:
        self.initial_angle = -math.pi/2.0
        self.final_angle = -math.pi/2.0
    else:
      av = (v1+v2)/2.0
      vec = v2 - v1
      perp = ( complex(-vec.imag, vec.real) if vec.real > 0 else complex(vec.imag, -vec.real) )
      t = av.imag / perp.imag
      self.vertical = False
      self.circ_center = (av - t*perp).real
      self.circ_radius = abs(v2-self.circ_center)
      self.circ_angle1 = math.atan2(v1.imag, v1.real - self.circ_center)
      self.circ_angle2 = math.atan2(v2.imag, v2.real - self.circ_center)
      if self.circ_angle2 > self.circ_angle1:
        self.initial_angle = self.circ_angle1 + math.pi/2.0
        self.final_angle = self.circ_angle2 + math.pi/2.0
      else:
        self.initial_angle = self.circ_angle1 - math.pi/2.0
        self.final_angle = self.circ_angle2 - math.pi/2.0   
  
  @classmethod
  def from_pt_angle_dist(cls, v, A, d):
    return cls(v, point_angle_dist(v, A, d))
  
  def reversed(self):
    return HypGeodesicInterval(self.end, self.start)
  
  def Euclidean_length(self):
    if self.vertical:
      return abs(self.end - self.start)
    else:
      return abs(self.circ_angle2-self.circ_angle1)*self.circ_radius
  
  def same_but_reversed(self, other):
    return same_float(self.start, other.end) and same_float(self.end, other.start)
  
  def act_by_mobius(self, M):
    return HypGeodesicInterval(M(self.start), M(self.end))
  
  def __repr__(self):
    return str(self)
  
  def __str__(self):
    return "GI(" + str(self.start) + "," +  str(self.end) + ")"

class HypTri:
  def __init__(self, lengths):
    self.lengths = lengths
    self.angles = tri_angles(self.lengths)
  
  def __repr__(self):
    return "HypTri(" + str(self.lengths) + ")"
  
  def __str__(self):
    ans = "HypTri(lengths: " + str(self.lengths) + " angles: " + str(self.angles) + ")"
    return ans
  
  def realize_along_gi(self, gi, side_ind):
    """return an embedded hyperbolic triangle in which side side_ind
    lies along the geodesic interval gi"""
    if not same_float(gi.length, self.lengths[side_ind]):
      raise ValueError("Can't realize it because the edge isn't the right size")
    #print "Realizing triangle ", self
    #print "With edge", side_ind
    #print "Along gi:", gi, " of length", gi.length
    p1,A1,d1 = gi.start, gi.initial_angle+self.angles[side_ind], self.lengths[(side_ind-1)%3]
    p2,A2,d2 = gi.end, (math.pi - self.angles[(side_ind+1)%3]) + gi.final_angle, self.lengths[(side_ind+1)%3]
    pp1 = point_angle_dist(p1, A1, d1)
    pp2 = point_angle_dist(p2, A2, d2)
    #print "First point ", p1, A1, d1, " -> ", pp1
    #print "Second point ", p2, A2, d2, " -> ", pp2
    if (not same_float(pp1.real, pp2.real)) or (not same_float(pp1.imag, pp2.imag)):
      raise ValueError("Doesn't seem to be the same going from both sides?")
    vs = [p1, p2, pp1]
    #rotate right so it aligns with the topological triangle
    vs = vs[-side_ind:] + vs[:-side_ind]
    return EmHypTri.from_vertices(vs)

class EmHypTri(HypTri):
  def __init__(self, GI):
    self.sides = [gi for gi in GI]
    self.v = [gi.start for gi in GI]
    self.lengths = [gi.length for gi in GI]
    self.angles = tri_angles(self.lengths)
  
  @classmethod
  def from_vertices(cls, em_V):
    GIs = [HypGeodesicInterval(em_V[i], em_V[(i+1)%3]) for i in xrange(3)]
    return cls(GIs)
  
  def __repr__(self):
    return str(self)
  
  def __str__(self):
    return "EmHypTri("+ str(self.v) +  ")"
    
  