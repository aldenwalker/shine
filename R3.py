import math
import R2

def triangle_normal(L):
  v1, v2, v3 = L
  return (v2-v1).cross(v3-v1)

def along_segment(s, t):
  return s[0]*(1-t) + s[1]*t

class Vector:
  def __init__(self, L):
    if isinstance(L, int):
      self.x = L*[0]
    elif isinstance(L, list) or isinstance(L, tuple):
      self.x = [i for i in L]
  
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
    #zoom until we hit the plane containing the triangle
    #print "\nChecking if ", t, "covers", p
    n = (t[1]-t[0]).cross( t[2]-t[0] )
    x = t[0]
    e = self.eye
    d = p - e
    ddn = d.dot(n)
    if abs(ddn) < 1e-6:
      return False
    alpha = (x-e).dot(n) / ddn
    #print "Got ", n, x, e, d, ddn, alpha
    if alpha >= 1.0-1e-6:
      #if the point is closer than the plane, then it can't be hidden
      #print "Point is closer; can't be hidden"
      return False
    p_in_plane = e + d*alpha
    #check if the point lies on the correct side of the three triangle sides
    for i in xrange(3):
      v = p_in_plane - t[i]
      c = (t[(i+1)%3] - t[i]).cross(v)
      #print "Checking against hyperplane", i, "got", v, c, c.dot(n)
      if c.dot(n) < -1e-4:
        return False
    return True
    
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
    
  def project_triangles(self, T):
    """project a list of triangles, sorting them"""
    sorted_T = sorted(T, key=lambda x: ((x[0]+x[1]+x[2])/3.0 -self.eye).norm(), reverse=True)
    return [self.project_triangle(t) for t in sorted_T]
  
  def visible_subsegments(self, segment, T):
    """returns a list of 3d segments which remains after cutting with all 
    the triangles in T"""
    s = [segment]
    print "Cutting ", s
    for t in T:
      print "With ", t
      n = (t[1]-t[0]).cross(t[2]-t[0])
      #if (self.eye - t[0]).dot(n) <= 0:
      #  continue
      t_projected = [R2.Vector(self.project_point(v)) for v in t]
      new_segments = []
      for seg in s:
        print "Cutting subsegment", seg
        if (seg[0]-t[0]).dot(n) >= 0 and (seg[1]-t[0]).dot(n) >= 0:
          new_segments.append(seg)
          continue
        seg_projected = [R2.Vector(self.project_point(v)) for v in seg]
        seg_projected_cut = R2.cut_segment_with_triangle_t_values(seg_projected, t_projected)
        print "Got ", seg_projected_cut
        if seg_projected_cut == None:
          continue
        new_segments.extend( [ [along_segment(seg, spc[0]), along_segment(seg, spc[1])] for spc in seg_projected_cut] )
      s = new_segments
      print "Current s:", s
      if len(s) == 0:
        return None
    return s

  
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
    self.viewer_grid_grid = [ [ [[],[]] for y in xrange(self.viewer_grid_num_vert_boxes)] \
                                        for x in xrange(self.viewer_grid_num_horiz_boxes) ]
    print "made grid:", self.viewer_grid_num_horiz_boxes, self.viewer_grid_num_vert_boxes, self.viewer_grid_box_width, self.viewer_grid_box_height, self.viewer_grid_ll, self.viewer_grid_ur
    for i,t in enumerate(T_projected):
      self.viewer_grid_add_projected_triangle(t, i)
  
  def viewer_grid_near_segment(self, s):
    grid_indices = self.viewer_grid_projected_segment_indices( map(self.project_point, s) )
    nearby_triangles = []
    nearby_segments = []
    for i,j in grid_indices:
      nearby_triangles.extend(self.viewer_grid_grid[i][j][0])
      nearby_segments.extend(self.viewer_grid_grid[i][j][1])
    return nearby_triangles, nearby_segments
  
  def viewer_grid_add_segment(self, s, ind):
    grid_indices = self.viewer_grid_projected_segment_indices( map(self.project_point, s) )
    for i,j in grid_indices:
      self.viewer_grid_grid[i][j][1].append(ind)
  
  def viewer_grid_add_projected_triangle(self, t, ind):
    grid_indices = self.viewer_grid_projected_triangle_indices( t )
    print "Adding triangle -- got grid indices: ", grid_indices
    for i,j in grid_indices:
      self.viewer_grid_grid[i][j][0].append(ind)
  
  def viewer_grid_projected_point_indices(self, pt):
    return int( (pt[0] - self.viewer_grid_ll[0] - 1e-12) / self.viewer_grid_box_width  ), \
           int( (pt[1] - self.viewer_grid_ll[1] - 1e-12) / self.viewer_grid_box_height )
  
  def viewer_grid_projected_segment_indices(self, s):
    I1 = self.viewer_grid_projected_point_indices(s[0])
    I2 = self.viewer_grid_projected_point_indices(s[1])
    return [ (i,j) for i in xrange(min(I1[0], I2[0]), max(I1[0],I2[0])+1) \
                   for j in xrange(min(I1[1], I2[1]), max(I1[1],I2[1])+1) ]
  
  def viewer_grid_projected_triangle_indices(self, t):
    I1 = self.viewer_grid_projected_point_indices(t[0])
    I2 = self.viewer_grid_projected_point_indices(t[1])
    I3 = self.viewer_grid_projected_point_indices(t[1])
    return [ (i,j) for i in xrange(min(I1[0],I2[0],I3[0]), max(I1[0],I2[0],I3[0])+1) \
                   for j in xrange(min(I1[1],I2[1],I3[1]), max(I1[1],I2[1],I3[1])+1) ]
  
    
    
























    
    
    