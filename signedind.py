class SignedInd:
  def __init__(self, ind, sign=1):
    self.ind = ind
    self.sign = sign
  @classmethod
  def from_string(cls, s):
    if len(s) == 0:
      return cls(0,1)
    if s[0].isdigit() or s[0] == '+':
      sign = 1
    else:
      sign = -1
    ind = int(s[(0 if s[0].isdigit() else 1):])
    return cls(ind, sign) 
  def __hash__(self):
    return hash(self.sign*(self.ind+1))
  def __repr__(self):
    return str(self) #"SignedInd(" + str(self.ind) + "," + str(self.sign) + ")"
  def __str__(self):
    return ("-" if self.sign < 0 else "") + str(self.ind)
  def __neg__(self):
    return SignedInd(self.ind, -self.sign)
  def __eq__(self, other):
    return isinstance(other, SignedInd) and other != None and self.ind == other.ind and self.sign == other.sign
  def __ne__(self, other):
    return not (self == other)
  