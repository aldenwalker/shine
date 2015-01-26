##########################################################################
# create a complement of the Cantor set to depth n (with 1+2^n boundaries)
##########################################################################

import hyp
import gsurf
import math
from signedind import SignedInd as SI


def xor(a,b):
  return a!=b


def t_along_circle(C,t):
  c,r = C
  return (c + r*math.cos(2.0*math.pi*t), r*math.sin(2.0*math.pi*t))

def circle_between_R_and_point(x,p):
  a,b = p
  if abs(a-x) < 0.001:
    return None
  f = abs(a-x)
  d = (b*b-f*f)/(2*f)
  if x < a:
    return (a+d, a+d-x, math.atan2(b, -d), math.pi)
  else:
    return (a-d, x-a+d, 0, math.atan2(b,d))
    

def circle_between_points(p1, p2):
  x1,y1 = p1
  x2,y2 = p2
  if abs(x1-x2) < 0.001:
    return None
  if x2 < x1:
    return circle_between_points(p2, p1)
  Q1 = x2-x1
  Q2 = (y2**2-y1**2)/(x2-x1)
  a = (Q1+Q2)/2.0
  b = (Q1-Q2)/2.0
  b = (-b if b<0 else b)
  c = x1 + a
  r = math.sqrt(a**2 + y1**2)
  if x1 + a < x2:
    a1 = math.atan2(y2, b)
    a2 = math.pi-math.atan2(y1,a)
  else:
    a1 = math.pi-math.atan2(y2,b)
    a2 = math.pi-math.atan2(y1,a)
  #print "Circle angles:", a1, a2, "Diff: ", (a2-a1)%(2*math.pi)
  if (a2-a1)%(2*math.pi) > math.pi+0.001:
    #print "Swapped"
    temp = a1
    a1 = a2
    a2 = temp+2*math.pi
  return (c,r,a1,a2)
    


###########################################################################
# return a string which graphics a mathematica drawing of the path between 
# the points 
###########################################################################
def mathematica_path(pts, col, start_up):
  ans = []
  for i in xrange(len(pts)-1):
    p1 = pts[i]
    p2 = pts[i+1]
    p1R = not isinstance(p1, tuple)
    p2R = not isinstance(p2, tuple)
    if p1R and p2R:
      c = (p1+p2)/2.0
      r = abs(p1-p2)/2.0
      ans.append( 'Circle[{' + str(c) + ',0}, ' + str(r) + (',{0,Pi}' if xor(not start_up, i%2==0) else ',{Pi,2*Pi}') + ']')
    elif p1R or p2R:
      if p1R:
        C = circle_between_points((p1,0), p2)
      else:
        C = circle_between_points((p2,0), p1)
      if C == None:
        ans.append('Line[{{' + str(p1) + ',0},{' + str(p2[0]) + ',' + str(p2[1]) + '}}]')
      else:
        #we don't need to negate because the points must already have taken this into account
        c,r,a1,a2 = C
        ans.append('Circle[{' + str(c) + ',0},' + str(r) + ',{' + str(a1) + ',' + str(a2) + '}]')
    else:
      C = circle_between_points(p1, p2)
      if C == None:
        ans.append('Line[{{' + str(p1[0]) + ',' + str(p1[1]) + '},{' + str(p2[0]) + ',' + str(p2[1]) + '}}]')
      else:
        c,r,a1,a2 = C
        if i%2==1:
          a1 = -a1
          a2 = -a2
        ans.append('Circle[{' + str(c) + ',0},' + str(r) + ',{' + str(a1) + ',' + str(a2) + '}]')
  return 'Graphics[{' + col + ',' + ','.join(ans) + '}]'
      
###########################################################################
# create 2^n circles which go over the standard Cantor set intervals
# they are returned as a (center, radius) list
###########################################################################
def cantor_set_circles(n, removal_fraction=1/3.0):
  return subdivide_circle(0.5, 0.5, n, removal_fraction)
  
def subdivide_circle(center, radius, depth, removal_fraction=1/3.0):
  if depth == 0:
    return [(center, radius)]
  new_r = radius*(1.0-removal_fraction)/2.0
  new_c1 = center - radius + new_r
  new_c2 = center + radius - new_r
  return subdivide_circle(new_c1, new_r, depth-1, removal_fraction) + \
         subdivide_circle(new_c2, new_r, depth-1, removal_fraction)

##########################################################################
# create uniform circles
##########################################################################
def uniform_circles(n, radius_fraction):
  center_gap = 1.0/float(n+1.0)
  rad = radius_fraction * center_gap
  ans = [((i+1)*center_gap, rad) for i in xrange(n)]
  return ans

##########################################################################
# braid 
##########################################################################
class DiskComplementHomeo(object):
  def __init__(self, i,j):
    self.i = i
    self.j = j
    self.min = min(i,j)
    self.max = max(i,j)
    self.rightward = (self.i < self.j)
  
  def __repr__(self):
    return 'DiskComplementHomeo(' + str(self.i) + ',' + str(self.j) + ')'
  
  def __str__(self):
    return repr(self)
  
  def cuff_destination(self, x):
    return (self.j if x==self.i else (self.i if x==self.j else x))
  
  def cuff_to_cuff_image_seams(self, c1, c2, ontop):
    if not self.rightward:
      ontop = not ontop
    start_inside = self.min < c1 and c1 < self.max
    start_outside = c1 < self.min or self.max < c1
    end_inside = self.min < c2 and c2 < self.max
    end_outside = c2 < self.min or self.max < c2
    if (start_inside and end_inside) or (start_outside and end_outside):
      return [c2]
    if start_inside:
      insertions = ([self.max, self.max+1] if ontop else [self.min+1, self.min])
      if end_outside:
        return insertions + [c2]
      else:
        return [insertions[0], self.cuff_destination(c2)]
    elif start_outside:
      insertions = ([self.max+1, self.max] if ontop else [self.min, self.min+1])
      if end_inside:
        return insertions + [c2]
      else:
        return [insertions[0], self.cuff_destination(c2)]
    else:
      if end_inside:
        return [(self.max if ontop else self.min+1), c2]
      elif end_outside:
        return [(self.max+1 if ontop else self.min), c2]
      else:
        return [self.cuff_destination(c2)]
  
  def cuff_to_seam_image_seams(self, c, s, ontop):
    start_inside = self.min < c and c < self.max
    start_outside = c < self.min or self.max < c
    end_inside = self.min < s and s <= self.max
    if start_inside or start_outside:
      return self.seam_to_seam_image_seams(c,s,ontop)
    if not self.rightward:
      ontop = not ontop
    if end_inside:
      return [(self.max if ontop else self.min+1), s]
    else:
      return [(self.max+1 if ontop else self.min), s]
  
  def seam_to_cuff_image_seams(self, s, c, ontop):
    start_inside = self.min < s and s <= self.max 
    end_inside = self.min < c and c < self.max
    end_outside = c < self.min or self.max < c
    if end_inside or end_outside:
      return self.seam_to_seam_image_seams(s,c,ontop)
    if not self.rightward:
      ontop = not ontop
    if start_inside:
      return [(self.max if ontop else self.min+1), self.cuff_destination(c)]
    else:
      return [(self.max+1 if ontop else self.min), self.cuff_destination(c)]
  
  def seam_to_seam_image_seams(self, s1, s2, ontop):
    if not self.rightward:
      ontop = not ontop
    start_inside = (self.min < s1 and s1 <= self.max)
    end_inside = (self.min < s2 and s2 <= self.max)
    if start_inside == end_inside:
      return [s2]
    if ontop:
      if start_inside:
        return [self.max, self.max+1, s2]
      else:
        return [self.max+1, self.max, s2]
    else:
      if start_inside:
        return [self.min+1, self.min, s2]
      else:
        return [self.min, self.min+1, s2]
    
  def __call__(self, P):
    #print "Acting on ", P
    if P.start in [self.i, self.j]:
      new_start = self.cuff_destination(P.start)
      new_start_up = not P.start_up
    else:
      new_start = P.start
      new_start_up = P.start_up
    new_seams = []
    if len(P.seams) == 0:
      new_seams.extend( self.cuff_to_cuff_image_seams(P.start, P.end, P.start_up) )
    else:
      new_seams.extend( self.cuff_to_seam_image_seams(P.start, P.seams[0], P.start_up) )
      for i in xrange(0,len(P.seams)-1):
        ontop = xor((i%2!=0), not P.start_up)
        s = P.seams[i]
        d = P.seams[i+1]
        new_seams.extend( self.seam_to_seam_image_seams(s,d,ontop) )
      ontop = xor(len(P.seams)%2==0, not P.start_up)
      new_seams.extend( self.seam_to_cuff_image_seams(P.seams[-1], P.end, ontop) ) 
    new_p = DiskComplementBraid(new_start, new_seams[:-1], new_seams[-1], start_up=new_start_up)
    #print "Got:", new_p
    new_p.simplify(rotate_ends=True)
    #print "Simplified:", new_p
    return new_p


############################################################################
# a path in a disk complement
############################################################################
class DiskComplementBraid(object):
  def __init__(self, start, seams, end, start_up=True):
    self.start = start
    self.seams = seams
    self.end = end
    self.start_up = (start==-1 or start_up)
    self.start_t = 0.5
    self.seam_ts = [0.5 for s in seams]
    self.end_t = 0.5
  
  def __str__(self):
    return repr(self)
  
  def __repr__(self):
    return 'DiskComplementBraid(' + str(self.start) + ',' + str(self.seams) + ',' + str(self.end) + ',start_up:' + str(self.start_up) + ')'
  
  def simplify(self,rotate_ends=False):
    i=0
    while i<len(self.seams)-1:
      if self.seams[i] == self.seams[i+1]:
        del self.seams[i+1]
        del self.seams[i]
      else:
        i += 1
    if rotate_ends:
      while len(self.seams)>0 and (self.seams[-1] == self.end or self.seams[-1] == self.end+1):
        del self.seams[-1]
      if self.start != -1:
        while len(self.seams)>0 and (self.seams[0] == self.start or self.seams[0] == self.start+1):
          del self.seams[0]
          self.start_up = not self.start_up
  
  ######################################################################
  # sort several paths
  ######################################################################
  @staticmethod
  def sort(L):
    def circle_cmp(a,b,c):
      return (a<b and b<c) or (c<a and a<b) or (b<c and c<a)
    def dir_cmp(v1,v2):
      i1, j1, dir1 = v1
      p11 = (1 if dir1 else -1)
      i2, j2, dir2 = v2
      p12 = (1 if dir2 else -1)
      if L[i1].seams[j1] != L[i2].seams[j2]:
        return None
      steps = 0
      while True:
        #print "Current:", i1, j1, i2, j2, L[i1].seams[j1], L[i2].seams[j2]
        done_1 = False
        done_2 = False
        if j1 == len(L[i1].seams)-1 and dir1:
          next_1 = 2*L[i1].end+1
          done_1 = True
        elif j1 == 0 and not dir1:
          next_1 = 2*L[i1].start+1
          done_1 = True
        else:
          next_1 = 2*L[i1].seams[j1+p11]
        if j2 == len(L[i2].seams)-1 and dir2:
          next_2 = 2*L[i2].end+1
          done_2 = True
        elif j2 == 0 and not dir2:
          next_2 = 2*L[i2].start+1
          done_2 = True
        else:
          next_2 = 2*L[i2].seams[j2+p12]
        if done_1 and done_2 and next_1 == next_2:
          return None, steps
        if next_1 != next_2:
          break
        j1 += p11
        j2 += p12
        steps += 1
      prev_loc = 2*L[i1].seams[j1]
      #print "Comparing ", next_1, prev_loc, next_2, "got", circle_cmp(next_1, prev_loc, next_2)
      return ( circle_cmp(next_1, prev_loc, next_2), steps )
    def cmp_places(v1, v2):
      i1,j1,dir1 = v1
      i2,j2,dir2 = v2
      F,Fsteps = dir_cmp(v1, v2)
      B,Bsteps = dir_cmp((i1,j1, not dir1), (i2, j2, not dir2))
      #print "Comparing ", v1, v2, " got ", F,Fsteps, B,Bsteps
      if F == None:
        Bswapped = (Bsteps%2==1) != B
        return (-1 if Bswapped else 1)
      elif B == None:
        Fswapped = (Fsteps%2==1) != F
        return (-1 if Fswapped else 1)
      Fswapped = (Fsteps%2==1) != F #xor
      Bswapped = (Bsteps%2==1) != B #xor
      #print "Swapping to: ", Fswapped, Bswapped
      F = Fswapped
      B = Bswapped
      if F and B:
        return -1
      elif (not F) and (not B):
        return 1
      elif Fsteps <= Bsteps:
        return (-1 if F else 1)
      else:
        return (-1 if B else 1)
    #True means going up in the direction of travel
    decorated = [(i, j, xor( (j%2!=0), not L[i].start_up)  ) for i in xrange(len(L)) for j in xrange(len(L[i].seams))]
    #print "Decorated: ", decorated
    groups = dict()
    for (i,j,dir) in decorated:
      if L[i].seams[j] in groups:
        groups[L[i].seams[j]].append( (i,j,dir) )
      else:
        groups[L[i].seams[j]] = [ (i,j,dir) ]
    #print "Groups: ", groups
    for si in groups:
      groups[si].sort(cmp=cmp_places)
      #print "Sorted ", si, groups[si]
      for ii,(i,j,dir) in enumerate(groups[si]):
        L[i].seam_ts[j] = (float(ii)+1.0)/(len(groups[si])+1.0)
      
    
    
    
#########################################################################
# a disk complement
#########################################################################
class DiskComplement:
  def __init__(self, kind='uniform', **kwargs):
    self.kind = kind
    if kind == 'cantor':
      self.depth = kwargs.get('depth', 4)
      self.removal_fraction = kwargs.get('removal_fraction', 1/3.0)
      self.C = cantor_set_circles(self.depth, self.removal_fraction)
      self.big_circle = (0.5, 0.6)
    else:
      self.n = kwargs.get('n', 4)
      self.radius_fraction = kwargs.get('radius_fraction', 0.15)
      self.C = uniform_circles(self.n, self.radius_fraction)
      self.big_circle = (0.5, 0.6)
    self.braids = []
    
  def __repr__(self):
    ans = 'DiskComplement(kind=' + str(self.kind)
    if self.kind == 'uniform':
      ans += ', n=' + str(self.n) + ', radius_fraction=' + str(self.radius_fraction) + ')'
    else:
      ans += ', depth=' + str(self.depth) + ', removal_fraction=' + str(self.removal_fraction) + ')'
    return ans
  
  def __str__(self):
    ans = repr(self)
    if len(self.braids) > 0:
      ans += '\nBraids (' + str(len(self.braids)) + '):\n'
      for i,b in enumerate(self.braids):
        ans += str(i) + ': ' + str(b) + '\n'
    return ans
  
  def add_braid(self, *args):
    for b in args:
      if isinstance(b,list):
        self.braids.append(DiskComplementBraid(b[0], b[1:-1], b[-1]))
      else:
        self.braids.append(b)
    self.braids[-1].simplify(rotate_ends=True)
  
  def apply_homeo(self, *args):
    for h in reversed(args):
      self.braids = [h(b) for b in self.braids]
  
  def mathematica_picture(self, paths=None, big_circle=False):
    
    if paths == None:
      paths = self.braids
    else:
      paths = paths.extend(self.braids)
    
    DiskComplementBraid.sort(paths)
    
    draw_big_circle = self.big_circle
    cuff_circles = ['Circle[{'+str(c[0])+',0},'+str(c[1])+']' for c in self.C]
    outside_circle = 'Circle[{'+str(draw_big_circle[0])+',0},'+str(draw_big_circle[1])+']'
    if big_circle:
      complement_picture = 'Graphics[{' + ','.join(cuff_circles+[outside_circle]) + '}]'
    else:
      complement_picture = 'Graphics[{' + ','.join(cuff_circles) + '}]'
    path_pictures = []
    colors = ['Blue','Red','Green', 'Yellow','Cyan','Magenta']
    for pi,p in enumerate(paths):
      if p.start == -1 or p.start == len(self.C):
        start_point = t_along_circle(draw_big_circle, p.start_t/2.0)
      else:
        start_point = t_along_circle(self.C[p.start], (1 if p.start_up else -1)*p.start_t/2.0)
      mid_points = []
      for i in xrange(len(p.seams)):
        si = p.seams[i]
        if si==0:
          p1 = draw_big_circle[0]-draw_big_circle[1]
          p2 = self.C[0][0]-self.C[0][1]
          mid_points.append( (1-p.seam_ts[i])*p1 + p.seam_ts[i]*p2 )
        elif si == len(self.C):
          p1 = self.C[-1][0]+self.C[-1][1]
          p2 = draw_big_circle[0]+draw_big_circle[1]
          mid_points.append( (1-p.seam_ts[i])*p1 + p.seam_ts[i]*p2 )
        else:
          p1 = self.C[si-1][0]+self.C[si-1][1]
          p2 = self.C[si][0]-self.C[si][1]
          mid_points.append( (1-p.seam_ts[i])*p1 + p.seam_ts[i]*p2 )
      if p.end == -1 or p.end == len(self.C):
        end_point = t_along_circle(draw_big_circle, p.end_t/2.0)
      else:
        end_point = t_along_circle(self.C[p.end], (1 if xor(not p.start_up, len(p.seams)%2==0) else -1)*p.end_t/2.0)
      all_points = [start_point] + mid_points + [end_point]
      this_path_picture = mathematica_path(all_points, colors[pi%len(colors)], p.start_up)
      path_pictures.append(this_path_picture)
    return 'Show[' + complement_picture + ',' + ','.join(path_pictures) + ']'
      
      
      





#########################################################################
# a cantor set complement with a hyperbolic structure
#########################################################################
class CantorSetComplement:
  def __init__(self, n, removal_fraction=1/3.0):
    self.C = cantor_set_circles(n, removal_fraction)
    self.big_circle = (0.5, 5.0)
    #print self.C
    bottom_seams = [hyp.HypGeodesicInterval.orthogonal_to_geodesics(self.C[i], self.C[i+1]) for i in xrange(len(self.C)-1)]
    #print bottom_seams
    left_seam = hyp.HypGeodesicInterval.orthogonal_to_geodesics(self.big_circle, self.C[0])
    right_seam = hyp.HypGeodesicInterval.orthogonal_to_geodesics(self.C[-1], self.big_circle)
    #print left_seam, right_seam
    self.seams = [left_seam] + bottom_seams + [right_seam]
    self.cuffs = [hyp.HypGeodesicInterval(self.seams[i].end, self.seams[i+1].start) for i in xrange(len(self.seams)-1)]
    top_cuff = hyp.HypGeodesicInterval(self.seams[-1].end, self.seams[0].start)
    self.cuffs.append(top_cuff)
    #print "Cuff lengths: ", [x.length for x in self.cuffs]
    #print "Seam lengths: ", [x.length for x in self.seams]
    #print ','.join([x.mathematica_string() for x in self.seams] + [x.mathematica_string('Red') for x in self.cuffs])
    
    #the polygon's 0th side is the 0th seam
    sides = []
    for i in xrange(len(self.seams)):
      sides.append(self.seams[i].length)
      sides.append(self.cuffs[i].length)
    self.polygon = hyp.HypPolygon( sides, [math.pi/2 for i in xrange(2*len(self.seams))] )
  
  def geodesicify_path(self, p):
    pass
  
  def mathematica_picture(self, paths):
    draw_big_circle = (0.5,0.55)
    cuff_circles = ['Circle[{'+str(c[0])+',0},'+str(c[1])+']' for c in self.C]
    outside_circle = 'Circle[{'+str(draw_big_circle[0])+',0},'+str(draw_big_circle[1])+']'
    complement_picture = 'Graphics[{' + ','.join(cuff_circles+[outside_circle]) + '}]'
    path_pictures = []
    colors = ['Blue','Red','Green', 'Yellow','Cyan','Magenta']
    for pi,p in enumerate(paths):
      start_point = t_along_circle((0.5,0.55), p.start_t/2.0)
      mid_points = []
      path_drawing = []
      for i in xrange(len(p.seams)):
        si = p.seams[i]
        if si==0:
          p1 = draw_big_circle[0]-draw_big_circle[1]
          p2 = self.C[0][0]-self.C[0][1]
          mid_points.append( (1-p.seam_ts[i])*p1 + p.seam_ts[i]*p2 )
        elif si == len(self.seams)-1:
          p1 = self.C[-1][0]+self.C[-1][1]
          p2 = draw_big_circle[0]+draw_big_circle[1]
          mid_points.append( (1-p.seam_ts[i])*p1 + p.seam_ts[i]*p2 )
        else:
          p1 = self.C[si-1][0]+self.C[si-1][1]
          p2 = self.C[si][0]-self.C[si][1]
          mid_points.append( (1-p.seam_ts[i])*p1 + p.seam_ts[i]*p2 )
      end_point = t_along_circle(self.C[p.end], p.start_t/2.0)
      if len(p.seams)%2 == 1:
        end_point = (end_point[0], -end_point[1])
      if len(p.seams) == 0:
        C = circle_between_R_and_point(end_point[0], start_point)
        if C == None:
          path_drawing.append( 'Line[{{'+str(start_point[0])+','+str(start_point[1])+'},{'+str(end_point[0])+','+str(end_point[1])+'}}]')
        else:
          c,r,a1,a2 = C
          path_drawing.append( 'Circle[{'+str(c)+',0},'+str(r)+',{' + str(a1) + ','+str(a2) + '}]' )
      else:
        C = circle_between_R_and_point(mid_points[0], start_point)
        if C==None:
          path_drawing.append('Line[{{'+str(start_point[0]) + ',' + str(start_point[1])+'},{'+str(mid_points[0])+',0}}]')
        else:
          c,r,a1,a2 = C
          path_drawing.append( 'Circle[{'+str(c)+',0},'+str(r)+',{' + str(a1) + ','+str(a2) + '}]' )
        for i in xrange(len(mid_points)-1):
          c = (mid_points[i]+mid_points[i+1])/2.0
          r = (abs(mid_points[i+1]-mid_points[i]))/2.0
          path_drawing.append('Circle[{'+str(c)+',0},'+str(r)+',' + ('{Pi,2*Pi}]' if i%2==0 else '{0,Pi}]'))
        c = (mid_points[-1]+end_point[0])/2.0
        r = (abs(end_point[0]-mid_points[-1]))/2.0
        #path_drawing.append('Line[{{'+str(mid_points[-1])+',0},{'+str(end_point[0])+','+str(end_point[1])+'}}]')
        path_drawing.append('Circle[{'+str(c)+',0},'+str(r)+',' + ('{Pi,2*Pi}]' if len(p.seams)%2==1 else '{0,Pi}]'))
      path_pictures.append( 'Graphics[{' + colors[pi%len(colors)] + ',' + ','.join(path_drawing) + '}]' )
    return 'Show[' + complement_picture + ',' + ','.join(path_pictures) + ']'
    
  def geoify(self, p):
    while True:
      max_t_diff = 0
      for i in xrange(len(p.seams)):
        if i==0:
          prev_geo = self.cuffs[p.start]
          prev_pt = prev_geo.pt_along_euclidean(p.start_t)
        else:
          prev_geo = self.seams[p.seams[i-1]]
          prev_pt = prev_geo.pt_along_euclidean(p.seam_ts[i-1])
        if i==len(p.seams)-1:
          next_geo = self.cuffs[p.end]
          next_pt = next_geo.pt_along_euclidean(p.end_t)
        else:
          next_geo = self.seams[p.seams[i+1]]
          next_pt = next_geo.pt_along_euclidean(p.seam_ts[i+1])
        geo = self.seams[p.seams[i]]
        current_t = p.seam_ts[i]
        print "Seam ", i, " current t: ", current_t
        print "Geo: ", geo
        print "prev pt: ", prev_pt, " next_pt: ", next_pt
        reflected_next_pt = hyp.reflect_in_geodesic( geo.circ_center, geo.circ_radius, next_pt )
        print "Reflected: ", reflected_next_pt
        crossing_geo = hyp.HypGeodesicInterval(prev_pt, reflected_next_pt)
        print "Crossing geo: ", crossing_geo
        next_t = geo.euclidean_intersection_t(crossing_geo)
        print "Next t: ", next_t
        t_diff = abs(next_t-current_t)
        print "t_diff: ", t_diff
        if t_diff > max_t_diff:
          max_t_diff = t_diff
        p.seam_ts[i] = next_t
      break
    
    















