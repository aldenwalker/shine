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

def tri_opposite_length(ell1, angle, ell2):
  """returns the length of the side of the triangle opposite angle"""
  return math.acosh(  math.cosh(ell1)*math.cosh(ell2) -                  \
                      math.sinh(ell1)*math.sinh(ell2)*math.cos(angle) )

def tri_angle_side_angle(a1, s, a2):
  """if the sides go [s,X,Y] with angles a1, a2, on either side of s, this 
  function returns X,Y"""
  #find the remaining angle
  S = math.acos(math.sin(a1)*math.sin(a2)*math.cosh(s) - math.cos(a1)*math.cos(a2))
  #use the law of sines to get the other sides
  sin_ratio = math.sinh(s)/math.sin(S)
  X = math.asinh(sin_ratio*math.sin(a1))
  Y = math.asinh(sin_ratio*math.sin(a2))
  return X,Y
  


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


class HypGeodesic:
  def __init__(self, center, radius):
    self.center = center
    self.radius = radius
  def __repr__(self):
    return 'HypGeodesic(' + str(self.center) + ',' + str(self.radius) + ')'
  def __str__(self):
    return repr(self)
  @classmethod
  def from_real_endpoints(cls, p1, p2):
    c = (p1 + p2)/2.0
    r = abs(p2-p1) / 2.0
    return cls(c,r)


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
  
  @classmethod
  def orthogonal_to_geodesics(cls, G1, G2):
    c1, r1 = G1
    c2, r2 = G2
    #print "Find geodesic orthogonal to ", (c1, r1), (c2, r2)
    if r1 == 'inf':
      if r2 == 'inf':
        raise ValueError("No orthogonal geodesic segment between limiting geodesics")
      d1 = abs(c2-c1)
      d2 = math.sqrt(d1**2 - r2**2)
      far_angle = math.asin(d2/d1)
      top_small_angle = math.pi/2.0 - far_angle
      bottom_length = r2*math.sin(top_small_angle)
      height = r2*math.cos(top_small_angle)
      geo_seg_rad = math.sqrt(bottom_length**2 + height**2)
      dest = complex(r2 + (-bottom_length if c2 > c1 else bottom_length), height)
      return HypGeodesicInterval( complex(c1, geo_seg_rad), dest)
    if r2 == 'inf':
      return cls.orthogonal_to_geodesics(G2, G1).reversed()
    #they're both circles
    d = abs(c2-c1)
    h = math.sqrt( (d-r1-r2)*(d-r1+r2)*(d+r1-r2)*(d+r1+r2) ) / (2*d)
    a = math.sqrt(h**2 + r1**2)
    inner_angle1 = math.acos(r1/a)
    height1 = r1*math.sin(inner_angle1)
    base1 = r1*math.cos(inner_angle1)
    #print "Stuff1: ", d, h, a, inner_angle1, height1, base1
    if (c1<c2 and c1+r1>c2) or (c2>c1 and c2-r2>c1+r1) or (c1<c2 and c2-r2<c1):
      dir1 = 'right'
    else:
      dir1 = 'left'
    source = complex( (c1 + base1 if dir1=='right' else c1-base1), height1)
        
    b = math.sqrt(h**2 + r2**2)
    inner_angle2 = math.acos(r2/b)
    height2 = r2*math.sin(inner_angle2)
    base2 = r2*math.cos(inner_angle2)
    #print "Stuff2: ", h, b, inner_angle2, height2, base2
    if (c1<c2 and c1+r1>c2) or (c2<c1 and c2+r2<c1-r1) or (c2<c1 and c2+r2>c1):
      dir2 = 'right'
    else:
      dir2 = 'left'
    dest = complex( (c2 + base2 if dir2=='right' else c2-base2), height2)
    return HypGeodesicInterval(source, dest)
    
    
  
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
  
  def pt_along_euclidean(self, t):
    """t in [0,1]; this gives a parameterization"""
    if self.vertical:
      return self.start + t*(self.end - self.start)
    a = self.circ_angle1 + t*(self.circ_angle2 - self.circ_angle1)
    cc = self.circ_center
    cr = self.circ_radius
    return complex( cc + cr*math.cos(a), cr*math.sin(a) )
  
  def pt_along(self, t):
    """gives the point a fraction t of the way along (in hyperbolic metric)"""
    #get the mobius transporting the interval to vertical
    M = mobius.MobiusTrans.unit_tangent_to_I_vert(self.start, self.initial_angle)
    #p1 = 0+j #M(self.start)
    #p2 = M(self.end)
    total_dist = self.length
    d = t*total_dist
    #note the start is at I, with imaginary part 1
    #the formula for a point distance d up from (0,y1) is:
    #y1 Cosh[d] + Sqrt[-y1^2 + y1^2 Cosh[d]^2]
    p = complex(0, math.cosh(d) + math.sqrt(math.cosh(d)**2 - 1))
    return M.inverse()(p)
  
  def __repr__(self):
    return str(self)
  
  def __str__(self):
    return "GI(" + str(self.start) + "," +  str(self.end) + ")"
  
  def mathematica_string(self, opt='Black'):
    if self.vertical:
      return 'Graphics[{' + opt + ',Line[{{' + str(self.circ_center.real) + ',' + str(self.circ_angle1) + '},{' + str(self.circ_center) + ',' + str(self.circ_angle2) + '}}]}]'
    else:
      return 'Graphics[{' + opt + ',Circle[{' + str(self.circ_center) + ',0},' + str(self.circ_radius) + ',{' + str(self.circ_angle1) + ',' + str(self.circ_angle2) + '}]}]'

#########################################################################
# a hyperbolic triangle
#########################################################################
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
    if not same_float(gi.length, self.lengths[side_ind],tol=1e-7):
      raise ValueError("Can't realize it because the edge" + str(gi.length) + "isn't the right size" + str(self.lengths[side_ind]))
    #print "Realizing triangle ", self
    #print "With edge", side_ind
    #print "Along gi:", gi, " of length", gi.length
    p1,A1,d1 = gi.start, gi.initial_angle+self.angles[side_ind], self.lengths[(side_ind-1)%3]
    p2,A2,d2 = gi.end, (math.pi - self.angles[(side_ind+1)%3]) + gi.final_angle, self.lengths[(side_ind+1)%3]
    pp1 = point_angle_dist(p1, A1, d1)
    pp2 = point_angle_dist(p2, A2, d2)
    #print "First point ", p1, A1, d1, " -> ", pp1
    #print "Second point ", p2, A2, d2, " -> ", pp2
    if (not same_float(pp1.real, pp2.real,tol=1e-7)) or (not same_float(pp1.imag, pp2.imag,tol=1e-7)):
      raise ValueError("Doesn't seem to be the same going from both sides?")
    vs = [p1, p2, pp1]
    #rotate right so it aligns with the topological triangle
    vs = vs[-side_ind:] + vs[:-side_ind]
    return EmHypTri.from_vertices(vs)
  
  #########################################################################
  # returns the angle between (the geodesic joining e1t along edge e1i
  # to e2t along edge e2i) and edge e2i.  it returns the angle we get to the
  # *left* of the incoming geodesic
  #########################################################################
  def inside_final_angle(self, e1i, e1t, e2i, e2t):
    if e2i == (e1i+1)%3:
      x = (1-e1t)*self.lengths[e1i]
      a = self.angles[e2i]
      y = e2t*self.lengths[e2i]
      arc_len = tri_opposite_length(x,a,y)
      wrong_angle = tri_angle([arc_len, x, y], 0)
      angle = math.pi - wrong_angle
    else:
      x = (1-e2t)*self.lengths[e2i]
      a = self.angles[e1i]
      y = e1t*self.lengths[e1i]
      arc_len = tri_opposite_length(x,a,y)
      angle = tri_angle([x,y,arc_len],0)
    return angle
  
  #########################################################################
  # given a side on which we enter, and an angle, figure out which edge it 
  # hits, where on the edge, and with what angle it leaves (angle to the left)
  #########################################################################
  def exit_t_and_angle(self, i, t, angle_in):
    ip1 = (i+1)%3
    ip2 = (i+2)%3
    central_cut_len = tri_opposite_length( (1-t)*self.lengths[i], \
                                           self.angles[ip1],        \
                                           self.lengths[ip1] )
    central_left_angle = tri_angle( [ (1-t)*self.lengths[i],        \
                                      self.lengths[ip1],            \
                                      central_cut_len ], 0 )
    if angle_in < central_left_angle:
      out_i = ip1
      outside_len, inside_len = tri_angle_side_angle( angle_in,                \
                                                      (1-t)*self.lengths[i], \
                                                      self.angles[ip1] )
      out_t = outside_len / self.lengths[ip1]
      out_angle = math.pi - tri_angle( [inside_len, (1-t)*self.lengths[i], outside_len], 0)
    else:
      out_i = ip2
      alpha = math.pi - angle_in
      inside_len, outside_len = tri_angle_side_angle( self.angles[i],         \
                                                      t*self.lengths[i],      \
                                                      alpha )
      out_t = 1 - (outside_len / self.lengths[ip2])
      out_angle = tri_angle( [outside_len, t*self.lengths[i], inside_len], 0)
    return out_i, out_t, out_angle
      
  
#########################################################################
# a hyperbolic polygon
#########################################################################
class HypPolygon:
  def __init__(self, sides, angles):
    self.nsides = len(sides)
    self.sides = sides
    self.angles = angles
  
  ######################################################################
  # return the angle (to the right) of the geodesic connecting the 
  # point a fraction enter_t of the way along edge enter_index to the point
  # a fraction exit_t along the edge exit_index
  ######################################################################
  def enter_angle(self, enter_index, enter_t, exit_index, exit_t):
    #draw edges to each of the intermediate vertices
    current_intermediate_vertex = (enter_index+1)%self.nsides
    current_extra_edge_len = self.sides[enter_index]*(1-enter_t)
    current_angle_cut_off = 0
    cumulative_enter_angle = 0
    while current_intermediate_vertex != exit_index:
      remaining_vertex_angle = self.angles[current_intermediate_vertex] - current_angle_cut_off
      next_extra_edge_len = tri_opposite_length(current_extra_edge_len,       \
                                                remaining_vertex_angle,       \
                                                self.sides[current_intermediate_vertex] )
      current_tri_sides = [current_extra_edge_len,                          \
                           self.sides[current_intermediate_vertex],         \
                           next_extra_edge_len]
      ##
      current_intermediate_vertex = (current_intermediate_vertex+1)%self.nsides
      current_extra_edge_len = next_extra_edge_len
      current_angle_cut_off = tri_angle(current_tri_sides, 2)
      cumulative_enter_angle += tri_angle(current_tri_sides, 0)
    #now do the final triangle, 
    remaining_vertex_angle = self.angles[current_intermediate_vertex] - current_angle_cut_off
    next_extra_edge_len = tri_opposite_length( current_extra_edge_len,          \
                                              remaining_vertex_angle,         \
                                              self.sides[exit_index]*exit_t )
    cumulative_enter_angle += tri_angle( [current_extra_edge_len,            \
                                          self.sides[exit_index]*exit_t,     \
                                          next_extra_edge_len], 0 )
    return cumulative_enter_angle
    
  ########################################################################
  # print
  ########################################################################
  def __str__(self):
    ans = 'Hyperbolic polygon with ' + str(self.nsides) + ' sides\n'
    ans += 'Sides: ' + str(self.sides) + '\n'
    ans += 'Angles: ' + str(self.angles)
    return ans

  
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
  
  def gi_between_points(self, e1i, e1t, e2i, e2t):
    p1 = self.sides[e1i].pt_along(e1t)
    p2 = self.sides[e2i].pt_along(e2t)
    return HypGeodesicInterval(p1, p2)
  
  def act_by_mobius(self, M):
    return EmHypTri([gi.act_by_mobius(M) for gi in self.sides])
  
  def __repr__(self):
    return str(self)
  
  def __str__(self):
    return "EmHypTri("+ str(self.v) +  ")"
    
  