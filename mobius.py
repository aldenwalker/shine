import math

class MobiusTrans:
  def __init__(self, a,b,c,d):
    self.a, self.b, self.c, self.d = a,b,c,d
  
  def inverse(self):
    return MobiusTrans(self.d, -self.b, -self.c, self.a)
  
  def compose(self, other):
    return MobiusTrans(self.a*other.a+self.b*other.c, \
                       self.a*other.b+self.b*other.d, \
                       self.c*other.a+self.d*other.c, \
                       self.c*other.b+self.d*other.d)
  
  @classmethod
  def points_to_01inf(cls, z1, z2, z3):
    if z1 == 'inf':
      m = cls(0, z3-z2, -1, z3)
    elif z2 == 'inf':
      m = cls(1, -z1, 1, z3)
    elif z3 == 'inf':
      m = cls(-1, z1, 0, z1-z2)
    else:
      m = cls(z2-z3, -z1*(z2-z3), z2-z1, z3*z1-z2)
    sdet = math.sqrt(abs(m.a*m.d - m.b*m.c))
    m.a, m.b, m.c, m.d = m.a/sdet, m.b/sdet, m.c/sdet, m.d/sdet
    return m
  
  @classmethod
  def points_to_points(cls, [z1, z2, z3], [w1, w2, w3]):
    return cls.points_to_01inf(w1, w2, w3).inverse().compose( cls.points_to_01inf(z1, z2, z3) )
  
  def __call__(self, z):
    if z == 'inf':
      if self.c == 0:
        return 'inf'
      else:
        return self.a/self.c
    elif z == -self.d/self.c:
      return 'inf'
    else:
      return (self.a*z + self.b)/(self.c*z + self.d)
      