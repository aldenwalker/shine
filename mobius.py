import math
import hyp

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
  
  def __repr__(self):
    return "MobiusTrans("+str(self.a)+","+str(self.b)+","+str(self.c)+","+str(self.d)+")"
  def __str__(self):
    return repr(self)
  
  @classmethod
  def points_to_01inf(cls, z1, z2, z3):
    if z1 == 'inf':
      m = cls(0, z3-z2, -1, z3)
    elif z2 == 'inf':
      m = cls(1, -z1, 1, z3)
    elif z3 == 'inf':
      m = cls(-1, z1, 0, z1-z2)
    else:
      m = cls(z2-z3, -z1*(z2-z3), z2-z1, -z3*(z2-z1))
    sdet = math.sqrt(abs(m.a*m.d - m.b*m.c))
    m.a, m.b, m.c, m.d = m.a/sdet, m.b/sdet, m.c/sdet, m.d/sdet
    return m
  
  @classmethod
  def points_to_points(cls, Z, W):
    return cls.points_to_01inf(*W).inverse().compose( cls.points_to_01inf(*Z) )
  
  @classmethod
  def unit_tangent_to_I_vert(cls, p, A):
    #move to the I axis
    Mpar = cls(1, -p.real, 0, 1)
    #scale to i
    Mhyp = cls(1/math.sqrt(p.imag), 0, 0, math.sqrt(p.imag))
    #rotate to vertical
    theta = (A-math.pi/2.0)/2.0
    Mell = cls(math.cos(theta), -math.sin(theta), math.sin(theta), math.cos(theta))
    return Mell.compose( Mhyp.compose( Mpar ) )
  
  @classmethod
  def unit_tangent_action(cls, p1, A1, p2, A2):
    M1 = cls.unit_tangent_to_I_vert(p1, A1)
    M2 = cls.unit_tangent_to_I_vert(p2, A2)
    return M2.inverse().compose( M1 )
  
  def tr(self):
    return self.a + self.d
  
  def geodesic_axis(self):
    if abs(self.tr() <= 2):
      return None
    if hyp.same_float(0, self.c):
      return hyp.HypGeodesic( self.b/(self.d-self.a), 'inf' )
    ad = self.a-self.d
    bc = self.b*self.c
    p1 = ( ad-math.sqrt( ad**2 + 4*bc ) ) / 2*self.c
    p2 = ( ad+math.sqrt( ad**2 + 4*bc ) ) / 2*self.c
    if abs( 1.0 / (self.c*p1 + self.d)**2 ) > 1:
      return hyp.HypGeodesic.from_real_endpoints( p1, p2 )
    else:
      return hyp.HypGeodesic.from_real_endpoints( p2, p1 )
  
  def __call__(self, z):
    if z == 'inf':
      if self.c == 0:
        return 'inf'
      else:
        return self.a/self.c
    elif self.c != 0 and z == -self.d/self.c:
      return 'inf'
    else:
      return (self.a*z + self.b)/(self.c*z + self.d)
      