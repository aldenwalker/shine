import math

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
  
  def __getitem__(self, i):
    return self.x[i]
  
  def __repr__(self):
    return "Vector(" + str(self.x) + ")"
  
  def __str__(self):
    return str(self.x)

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

class ProjectionViewer:
  def __init__(self, eye, towards, lights):
    self.lights = lights
    self.eye = eye
    self.towards = towards.scaled_to_len(1.0)  # this is the plane normal
    self.plane_origin = self.eye + self.towards
    self.right = self.towards.cross( Vector([0,0,1]) )
    self.right = self.right.scaled_to_len(1.0)
    self.up = self.right.cross(self.towards)
    self.up = self.up.scaled_to_len(1.0)
  
  def project_point(self, pt):
    tp = pt - self.eye
    t = 1.0/self.towards.dot(tp)
    pt_in_plane = self.eye + tp*t
    X = (pt_in_plane - self.plane_origin).dot(self.right)
    Y = (pt_in_plane - self.plane_origin).dot(self.up)
    return (X,Y)
  
  def project_triangle(self, T):
    points = [self.project_point(p) for p in T]
    tnorm = (T[2]-T[0]).cross( T[1]-T[0] )
    tnorm = tnorm.scaled_to_len(1.0)
    to_lights = [(ell-T[0]).scaled_to_len(1.0) for ell in self.lights]
    amount_from_lights = [tnorm.dot(x) for x in to_lights]
    amount = max(amount_from_lights)
    return (points, amount)
    




























    
    
    