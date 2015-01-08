import math
import R2

def triangle_normal(L):
  v1, v2, v3 = L
  return (v2-v1).cross(v3-v1)

def along_segment(s, t):
  return s[0]*(1-t) + s[1]*t

def segment_intersect_plane_t_value(s, pt, n):
  d = s[1] - s[0]
  dn = d.dot(n)
  if abs(dn) < 1e-8:
    return None
  return (pt.dot(n) - s[0].dot(n)) / dn

def subsegment_from_t_values(s, T):
  return [along_segment(s, T[0]), along_segment(s, T[1])]

def remove_duplicate_floats(L):
  """removes duplicate floats in a list -- it assumes they are sorted"""
  i = 0
  while i < len(L)-1:
    if L[i+1]-L[i] < 1e-10:
      del L[i+1]
    else:
      i += 1
  return None

class Vector:
  def __init__(self, L):
    if isinstance(L, int):
      self.x = L*[0]
    elif isinstance(L, list) or isinstance(L, tuple):
      self.x = [i for i in L]
  
  def __neg__(self):
    return Vector([-v for v in self.x])
  
  def __sub__(self, other):
    return Vector([self.x[i]-other.x[i] for i in xrange(len(self.x))])
  
  def __add__(self, other):
    return Vector([self.x[i]+other.x[i] for i in xrange(len(self.x))])
  
  def __mul__(self, a):
    return Vector([a*self.x[i] for i in xrange(3)])  
  
  def __div__(self, a):
    return Vector([self.x[i]/a for i in xrange(3)])
  
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
    x1, x2, x3 = self.x
    o1, o2, o3 = other.x
    return Vector([ x2*o3-x3*o2, x3*o1-x1*o3, x1*o2-x2*o1 ])
  
  def copy(self):
    return Vector([i for i in self.x])
  

class Matrix:
  def __init__(self, L):
    self.M = L
  def __call__(self, x):
    return Vector( [ sum([self.M[i][j]*x[j] for j in xrange(len(x))]) for i in xrange(len(self.M))] )
  def __mul__(self, other):
    return Matrix( [ [ sum([self.M[i][k]*other.M[k][j] for k in xrange(len(self.M[0]))]) for j in xrange(len(other.M[0]))] for i in xrange(len(self.M))])




class PolygonalPath:
  def __init__(self, L):
    self.L = list(L)
  def __str__(self):
    return repr(self)
  def __repr__(self):
    return str(self.L)



class ProjectionViewer:
  def __init__(self, eye, towards, lights):
    self.lights = lights
    self.eye = eye
    self.unscaled_towards = towards.copy()
    self.towards = towards.scaled_to_len(2.0)  # this is the plane normal
    self.plane_origin = self.eye + self.towards
    self.right = self.towards.cross( Vector([0,0,1]) )
    self.right = self.right.scaled_to_len(1.0)
    self.up = self.right.cross(self.towards)
    self.up = self.up.scaled_to_len(1.0)
    self.viewer_grid = None
  
  def zoom(self, factor):
    self.eye = self.eye + self.unscaled_towards*factor
    self.unscaled_towards = self.unscaled_towards*(1.0-factor)
    self.plane_origin = self.eye + self.towards
  
  def is_segment_hidden(self, T, v1, v2):
    return self.is_point_hidden(T, v1) and self.is_point_hidden(T, v2)
  
  def is_point_hidden(self, T, p):
    #print "\n\nChecking if", p, "is hidden"
    for t in T:
      if self.triangle_hides_point(t, p):
        #print "Yes; it's hidden"
        return True
    #print "No, it's not hidden"
    return False
    
  def triangle_hides_point(self, t, p):
    hyp_pt = [t[0]]
    hyp_n = [-(t[1]-t[0]).cross(t[2]-t[0])]
    hyp_pt.extend( [t[i] for i in xrange(3)] )
    hyp_n.extend( [(self.eye-t[i]).cross(t[(i+1)%3]-t[i]) for i in xrange(3)] )
    inclusions = [ (p-hyp_pt[i]).dot(hyp_n[i]) > 1e-10 for i in xrange(4)]
    return all(inclusions)
    
  def faces_eye(self, pt, v):
    return v.dot(self.eye-pt) > 0
  
  def project_point(self, pt):
    tp = pt - self.eye
    t = 1.0/self.towards.dot(tp)
    pt_in_plane = self.eye + tp*t
    X = (pt_in_plane - self.plane_origin).dot(self.right)
    Y = (pt_in_plane - self.plane_origin).dot(self.up)
    return (X,Y)
  
  def project_point_flat(self, pt):
    v = pt - self.eye
    X = v.dot(self.right)
    Y = v.dot(self.up)
    return (X,Y)
  
  def project_triangle(self, T):
    points = [self.project_point(p) for p in T]
    tnorm = (T[1]-T[0]).cross( T[2]-T[0] )
    tnorm = tnorm.scaled_to_len(1.0)
    to_lights = [(ell-T[0]).scaled_to_len(1.0) for ell in self.lights]
    amount_from_lights = [tnorm.dot(x) for x in to_lights]
    amount = max(amount_from_lights)
    return (points, amount)
  
  def sort_triangles(self, T):
    sorted_T = sorted(T, key=lambda x: ((x[0]+x[1]+x[2])/3.0 -self.eye).norm(), reverse=True)
    return sorted_T
  
  def project_triangles(self, T):
    """project a list of triangles, *not* sorting them"""
    #sorted_T = sorted(T, key=lambda x: ((x[0]+x[1]+x[2])/3.0 -self.eye).norm(), reverse=True)
    return [self.project_triangle(t) for t in T]
  
  def segments_cross_t_value(self, s0, s1):
    ps0 = map(self.project_point, s0)
    ps1 = map(self.project_point, s1)
    t = R2.intersect_segments(ps0, ps1)
    if t == None:
      return None
    else:
      return t[0]
  
  def visible_subsegments_one_triangle_t_values(self, segment, t):
    hyp_pt = [t[0]]
    hyp_n = [-(t[1]-t[0]).cross(t[2]-t[0])]
    hyp_pt.extend( [t[i] for i in xrange(3)] )
    hyp_n.extend( [(self.eye-t[i]).cross(t[(i+1)%3]-t[i]) for i in xrange(3)] )
    inclusion_dots = [ [ (s-hyp_pt[i]).dot(hyp_n[i]) for i in xrange(4)] for s in segment ]
    inclusion_status = [ [(0 if abs(d)<1e-10 else (1 if d>0 else -1)) for d in sd] for sd in inclusion_dots]
    #print "Intersecting ", segment, "with", t
    t_values = []
    for i in xrange(4):
      is0 = inclusion_status[0][i]
      is1 = inclusion_status[1][i]
      if is0<=0 and is1<=0:
        return [[0.0,1.0]]
      if is0*is1 < 0:
        pierce_t = segment_intersect_plane_t_value(segment, hyp_pt[i], hyp_n[i])
        if pierce_t == None:
          print "Inclusions: ", inclusion_dots
          print "Inclusion_status: ", inclusion_status
          print "Index ", i
        pierce_point = along_segment(segment, pierce_t)
        #print pierce_t
        #print pierce_point
        #print [ (pierce_point-hyp_pt[j]).dot(hyp_n[j])  for j in xrange(4) if j!=i]
        if any([ (pierce_point-hyp_pt[j]).dot(hyp_n[j]) < 0 for j in xrange(4) if j != i]):
          continue
        t_values.append(pierce_t)
        continue
    t_values.sort()
    #print t_values
    remove_duplicate_floats(t_values)
    all_in_0 = all([x>=0 for x in inclusion_status[0]])
    if all_in_0:
      ltv = len(t_values)
      if ltv == 0:
        return None
      elif ltv > 1:
        raise ValueError("Wrong number of t values?")
      else:
        return [[t_values[0], 1.0]]
    all_in_1 = all([x>=0 for x in inclusion_status[1]])
    if all_in_1:
      if len(t_values) != 1:
        print inclusion_dots
        print inclusion_status
        print t_values
        raise ValueError("Wrong number of t values?")
      return [[0.0,t_values[0]]]
    ltv = len(t_values)
    if ltv == 0 or ltv == 1:
      return [[0.0,1.0]]
    if ltv == 2:
      return [[0.0,t_values[0]],[t_values[1],1.0]]
    print inclusion_dots
    print inclusion_status
    print t_values
    raise ValueError("Wrong number of t values?")

  
  def visible_subsegments(self, segment, T):
    """returns a list of 3d segments which remain after cutting with all 
    the triangles in T"""
    segs = [segment]
    #print "Cutting ", segs
    for t in T:
      #print "With ", t
      new_segments = []
      for s in segs:
        s_cut = self.visible_subsegments_one_triangle_t_values(s, t)
        #print "Got ", s_cut
        if s_cut == None:
          continue
        new_segments.extend( [ [along_segment(s, sc[0]), along_segment(s, sc[1])] for sc in s_cut] )
      segs = new_segments
      #print "Current segs:", segs
      if len(segs) == 0:
        return None
    #print "Returning ", segs
    return segs

  def visible_subsegments_t_values(self, segment, T):
    VSS = self.visible_subsegments(segment, T)
    if VSS == None:
      return VSS
    ans = []
    v = segment[1] - segment[0]
    vdv = v.dot(v)
    for vss in VSS:
      vssv1 = vss[0]-segment[0]
      t1 = v.dot(vssv1)/vdv
      vssv2 = vss[1]-segment[0]
      t2 = v.dot(vssv2)/vdv
      ans.append([t1, t2])
    return ans
    
  
  def viewer_grid_init_triangles(self, T):
    T_projected = [map(self.project_point, t) for t in T]
    box_ll = list(T_projected[0][0])
    box_ur = list(T_projected[0][0])
    max_segment_height = 0
    max_segment_width = 0
    for tp in T_projected:
      for i in xrange(3):
        slx = tp[(i+1)%3][0] - tp[i][0]
        sly = tp[(i+1)%3][1] - tp[i][1]
        if slx > max_segment_width:
          max_segment_width = slx
        if sly > max_segment_height:
          max_segment_height = sly
        tv = tp[i]
        if tv[0] < box_ll[0]:
          box_ll[0] = tv[0]
        elif tv[0] > box_ur[0]:
          box_ur[0] = tv[0]
        if tv[1] < box_ll[1]:
          box_ll[1] = tv[1]
        elif tv[1] > box_ur[1]:
          box_ur[1] = tv[1]
    self.viewer_grid_ll = box_ll
    self.viewer_grid_height = box_ur[1] - box_ll[1]
    self.viewer_grid_ur = box_ur
    self.viewer_grid_width = box_ur[0] - box_ll[0]
    self.viewer_grid_num_horiz_boxes = int( self.viewer_grid_width/max_segment_width )
    self.viewer_grid_num_vert_boxes = int( self.viewer_grid_height/max_segment_height )
    self.viewer_grid_box_width = self.viewer_grid_width / self.viewer_grid_num_horiz_boxes
    self.viewer_grid_box_height = self.viewer_grid_height / self.viewer_grid_num_vert_boxes
    self.viewer_grid_grid = [ [ [set(),set()] for y in xrange(self.viewer_grid_num_vert_boxes)] \
                                        for x in xrange(self.viewer_grid_num_horiz_boxes) ]
    #print "made grid:", self.viewer_grid_num_horiz_boxes, self.viewer_grid_num_vert_boxes, self.viewer_grid_box_width, self.viewer_grid_box_height, self.viewer_grid_ll, self.viewer_grid_ur
    for i,t in enumerate(T_projected):
      self.viewer_grid_add_projected_triangle(t, i)
    #print "Grid:", self.viewer_grid_grid
  
  def viewer_grid_near_segment(self, s):
    grid_indices = self.viewer_grid_projected_segment_indices( map(self.project_point, s) )
    nearby_triangles = set()
    nearby_segments = set()
    for i,j in grid_indices:
      nearby_triangles.update(self.viewer_grid_grid[i][j][0])
      nearby_segments.update(self.viewer_grid_grid[i][j][1])
    return nearby_triangles, nearby_segments
  
  def viewer_grid_near_point(self, pt):
    grid_index = self.viewer_grid_projected_point_indices( self.project_point(pt) )
    return self.viewer_grid_grid[grid_index[0]][grid_index[1]]
  
  def viewer_grid_add_segment(self, s, ind):
    grid_indices = self.viewer_grid_projected_segment_indices( map(self.project_point, s) )
    for i,j in grid_indices:
      self.viewer_grid_grid[i][j][1].add(ind)
  
  def viewer_grid_add_projected_triangle(self, t, ind):
    grid_indices = self.viewer_grid_projected_triangle_indices( t )
    #print "Adding triangle -- got grid indices: ", grid_indices
    for i,j in grid_indices:
      self.viewer_grid_grid[i][j][0].add(ind)
  
  def viewer_grid_projected_point_indices(self, pt):
    return int( (pt[0] - self.viewer_grid_ll[0] - 1e-12) / self.viewer_grid_box_width  ), \
           int( (pt[1] - self.viewer_grid_ll[1] - 1e-12) / self.viewer_grid_box_height )
  
  def viewer_grid_projected_segment_indices(self, s):
    I1 = self.viewer_grid_projected_point_indices(s[0])
    I2 = self.viewer_grid_projected_point_indices(s[1])
    ans = [ (i,j) for i in xrange(min(I1[0], I2[0]), max(I1[0],I2[0])+1) \
                  for j in xrange(min(I1[1], I2[1]), max(I1[1],I2[1])+1) ]
    return ans
  
  def viewer_grid_projected_triangle_indices(self, t):
    I1 = self.viewer_grid_projected_point_indices(t[0])
    I2 = self.viewer_grid_projected_point_indices(t[1])
    I3 = self.viewer_grid_projected_point_indices(t[2])
    return [ (i,j) for i in xrange(min(I1[0],I2[0],I3[0]), max(I1[0],I2[0],I3[0])+1) \
                   for j in xrange(min(I1[1],I2[1],I3[1]), max(I1[1],I2[1],I3[1])+1) ]
  
    
    
























    
    
    