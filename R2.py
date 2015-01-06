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

def cut_segment_with_triangle(s, t):
  diffs = [t[(i+1)%3]-t[i] for i in xrange(3)]
  ws = [ [diffs[i].cross(s[j]-t[i]) > 1e-8 for i in xrange(3)] for j in xrange(2)]
  if (not ws[0][0] and not ws[1][0]) or (not ws[0][1] and not ws[1][1]) or (not ws[0][2] and not ws[1][2]):
    return [s]
  if all(ws[0]):
    for i in xrange(3):
      if not ws[1][i]:
        I = intersect_segments(s, [t[i], t[(i+1)%3]])
        if I != None:
          return [[along_segment(s, I[0]), s[1]]]
    #print "all in"
    return None
  if all(ws[1]):
    for i in xrange(3):
      if not ws[0][i]:
        I = intersect_segments(s, [t[i], t[(i+1)%3]])
        if I != None:
          return [[s[0], along_segment(s, I[0])]]
    #print "all in"
    return None
  I = [intersect_segments(s, [t[i], t[(i+1)%3]]) for i in xrange(3)]
  I = [i for i in I if i != None]
  if len(I) == 0:
    return [s]
  I.sort()
  #print I
  return [[s[0], along_segment(s, I[0][0])],[along_segment(s, I[1][0]), s[1]]]
  
#the same except it returns fractions
def cut_segment_with_triangle_t_values(s, t):
  diffs = [t[(i+1)%3]-t[i] for i in xrange(3)]
  ws = [ [diffs[i].cross(s[j]-t[i]) > 1e-8 for i in xrange(3)] for j in xrange(2)]
  if (not ws[0][0] and not ws[1][0]) or (not ws[0][1] and not ws[1][1]) or (not ws[0][2] and not ws[1][2]):
    return [[0.0,1.0]]
  if all(ws[0]):
    for i in xrange(3):
      if not ws[1][i]:
        I = intersect_segments(s, [t[i], t[(i+1)%3]])
        if I != None:
          if I < 1-1e-10:
            return [[I[0], 1.0]]
          else:
            return None
    #print "all in"
    return None
  if all(ws[1]):
    for i in xrange(3):
      if not ws[0][i]:
        I = intersect_segments(s, [t[i], t[(i+1)%3]])
        if I != None:
          if I[0] > 1e-10:
            return [[0.0, I[0]]]
          else:
            return None
    #print "all in"
    return None
  I = [intersect_segments(s, [t[i], t[(i+1)%3]]) for i in xrange(3)]
  I = [i for i in I if i != None]
  if len(I) < 2:
    return [[0.0,1.0]]
  I.sort()
  #print I
  if I[1][0] - I[0][0] < 1e-10:
    return None
  return [ [0.0, I[0][0]], [I[1][0], 1.0] ]



######################################################################
# given a point in R2 and a list of triangles, return the list of 
# segments which remain visible after all the triangles have covered it
######################################################################
def visible_subsegment(segment, T):
  current_segments = [segment]
  #print "Cutting ", segment, "with", T
  for t in T:
    new_segments = []
    #print "Cutting with ", t
    #print "Segments: ", current_segments
    for s in current_segments:
      s_cut = cut_segment_with_triangle(s, t)
      #print "Cut ", s, "to", s_cut
      if s_cut == None:
        continue
      new_segments.extend(s_cut)
    current_segments = new_segments
    if len(current_segments) == 0:
      return None
  return current_segments
  
  
  