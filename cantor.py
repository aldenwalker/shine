##########################################################################
# create a complement of the Cantor set to depth n (with 1+2^n boundaries)
##########################################################################

import hyp
import gsurf
import math
from signedind import SignedInd as SI



####
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



class CantorSetPathBraid(object):
  def __init__(self, i,j):
    self.i = i
    self.j = j
    self.min = min(i,j)
    self.max = max(i,j)
    self.rightward = (self.i < self.j)
      
    
  def new_crossing_seams(self, x1, x2, ontop, end_on_cuff=False):
    if not self.rightward:
      ontop = not ontop
    start_inside = (self.min < x1 and x1 <= self.max)
    if end_on_cuff:
      if x2 != self.i and x2 != self.j:
        pass #it acts like any other one
      else:
        if ontop:
          if start_inside:
            return ([self.max, self.max] if x2 == self.min else [self.max, self.min])
          else:
            return ([self.max+1, self.max] if x2 == self.min else [self.max+1, self.min])
        else:
          if start_inside:
            return ([self.min+1, self.min] if x2 == self.max else [self.min+1, self.max])
          else:
            return ([self.min, self.min] if x2 == self.max else [self.min, self.max])
    end_inside = (self.min < x2 and x2 <= self.max)
    if start_inside == end_inside:
      return [x2]
    if ontop:
      if start_inside:
        return [self.max, self.max+1, x2]
      else:
        return [self.max+1, self.max, x2]
    else:
      if start_inside:
        return [self.min+1, self.min, x2]
      else:
        return [self.min, self.min+1, x2]
    
  def __call__(self, P):
    print "Acting on ", P
    new_seams = []
    new_start = P.start
    for i in xrange(len(P.seams)):
      ontop = (i%2==0)
      s = (P.seams[i-1] if i>0 else P.start)
      d = P.seams[i]
      print "Acting on seams", s,d
      new_seams.extend(self.new_crossing_seams(s,d, ontop))
      print "New seams total: ", new_seams
    if len(P.seams)>0:
      new_seams.extend(self.new_crossing_seams(P.seams[-1],P.end, len(P.seams)%2==0, end_on_cuff=True))
    else:
      new_seams.extend(self.new_crossing_seams(P.start,P.end, True, end_on_cuff=True))
    new_p = CantorSetPath(new_start, new_seams[:-1], new_seams[-1])
    new_p.simplify()
    return new_p



class CantorSetPath(object):
  def __init__(self, start, seams, end):
    self.start = start
    self.seams = seams
    self.end = end
    self.geometrized = False
    self.start_t = 0.5
    self.seam_ts = [0.5 for s in seams]
    self.end_t = 0.5
  
  def __str__(self):
    return repr(self)
  
  def __repr__(self):
    return 'CantorSetPath(' + str(self.start) + ',' + str(self.seams) + ',' + str(self.end) + ')'
  
  def simplify(self):
    i=0
    while i<len(self.seams)-1:
      if self.seams[i] == self.seams[i+1]:
        del self.seams[i+1]
        del self.seams[i]
      else:
        i += 1
  
  def sort(self):
    def circle_cmp(a,b,c):
      return (a<b and b<c) or (c<a and a<b) or (b<c and c<a)
    def dir_cmp(v1,v2):
      i1, dir1 = v1
      p11 = (1 if dir1 else -1)
      i2, dir2 = v2
      p12 = (1 if dir2 else -1)
      if self.seams[i1]!=self.seams[i2]:
        return None
      steps = 0
      while True:
        print "Current:", i1, i2, self.seams[i1], self.seams[i2]
        if i1 == len(self.seams)-1 and dir1:
          next_1 = 2*self.end+1
        elif i1 ==0 and not dir1:
          next_1 = -1
        else:
          next_1 = 2*self.seams[i1+p11]
        if i2 == len(self.seams)-1 and dir2:
          next_2 = 2*self.end+1
        elif i2 ==0 and not dir2:
          next_2 = -1
        else:
          next_2 = 2*self.seams[i2+p12]
        if next_1 != next_2:
          break
        i1 += p11
        i2 += p12
        steps += 1
      ontop = ((i1%2 != 0) == dir1)
      print "Comparing ", next_1, 2*self.seams[i1], next_2, "got", circle_cmp(next_1, 2*self.seams[i1], next_2)
      return ( circle_cmp(next_1, 2*self.seams[i1], next_2), steps )
    def cmp_places(v1, v2):
      i1,dir1 = v1
      i2,dir2 = v2
      F,Fsteps = dir_cmp(v1, v2)
      B,Bsteps = dir_cmp((i1,not dir1), (i2,not dir2))
      print "Comparing ", v1, v2, " got ", F,Fsteps, B,Bsteps
      Fswapped = (Fsteps%2==1) != F #xor
      Bswapped = (Bsteps%2==1) != B #xor
      print "Swapping to: ", Fswapped, Bswapped
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
    decorated = [(i, i%2!=0) for i in xrange(len(self.seams))]
    print "Decorated: ", decorated
    groups = dict()
    for (i,dir) in decorated:
      if self.seams[i] in groups:
        groups[self.seams[i]].append( (i,dir) )
      else:
        groups[self.seams[i]] = [ (i,dir) ]
    print "Groups: ", groups
    for si in groups:
      groups[si].sort(cmp=cmp_places)
      print "Sorted ", si, groups[si]
      for ii,(i,dir) in enumerate(groups[si]):
        self.seam_ts[i] = (float(ii)+1.0)/(len(groups[si])+1.0)
    
      
    
    
    

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
    
    















