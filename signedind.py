class SignedInd:
  def __init__(self, ind, sign=1):
    self.ind = ind
    self.sign = sign
  def __repr__(self):
    return str(self) #"SignedInd(" + str(self.ind) + "," + str(self.sign) + ")"
  def __str__(self):
    return ("-" if self.sign < 0 else "") + str(self.ind)
  def __neg__(self):
    return SignedInd(self.ind, -self.sign)
  def __eq__(self, other):
    return other != None and self.ind == other.ind and self.sign == other.sign
  