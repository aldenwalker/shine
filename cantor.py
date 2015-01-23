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




class CantorSetPath:
  def __init__(self, start, seams, end):
    self.start = start
    self.seams = seams
    self.end = end
    self.geometrized = False
    self.start_t = 0.5
    self.seam_ts = [0.5 for s in seams]
    self.end_t = 0.5
    

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
    
    




def cantor_set_complement(n, removal_fraction=1/3.0):
  C = cantor_set_circles(n, removal_fraction)
  big_circle = (0.5, 5.0)
  print C
  bottom_seams = [hyp.HypGeodesicInterval.orthogonal_to_geodesics(C[i], C[i+1]) for i in xrange(len(C)-1)]
  print bottom_seams
  left_seam = hyp.HypGeodesicInterval.orthogonal_to_geodesics(big_circle, C[0])
  right_seam = hyp.HypGeodesicInterval.orthogonal_to_geodesics(C[-1], big_circle)
  print left_seam, right_seam
  seams = [left_seam] + bottom_seams + [right_seam]
  cuffs = [hyp.HypGeodesicInterval(seams[i].end, seams[i+1].start) for i in xrange(len(seams)-1)]
  top_cuff = hyp.HypGeodesicInterval(seams[-1].end, seams[0].start)
  print "Cuff lengths: ", [x.length for x in cuffs+[top_cuff]]
  print "Seam lengths: ", [x.length for x in seams]
  print ','.join([x.mathematica_string() for x in seams] + [x.mathematica_string('Red') for x in cuffs + [top_cuff]])
  #there's one vertex for every seam and cuff
  #starts at the left seam beginning, and goes around
  #make an edge for the each seam
  V = []
  E = []
  h_lens = []
  cuff_top_edges = []
  cuff_bottom_edges = []
  
  for i in xrange(len(seams)):
    E.append( tsurf.Edge(len(V),len(V)+1, None, None) )
    h_lens.append( seams[i].length )
    V.append( tsurf.Vertex([],[]) )
    V.append( tsurf.Vertex([],[]) )
  for i in xrange(len(cuffs)):
    cuff_top_edges.append( len(E) )
    E.append( tsurf.Edge(2*i+1, 2*i+2, tsurf.BD, None) )
    h_lens.append( cuffs[i].length )
    cuff_bottom_edges.append( len(E) )
    E.append( tsurf.Edge(2*i+2, 2*i+1, tsurf.BD, None) )
    h_lens.append(cuffs[i].length )
  cuff_top_edges.append( len(E) )
  E.append( tsurf.Edge(len(V)-1, 0, tsurf.BD, None) )
  h_lens.append( top_cuff.length )
  cuff_bottom_edges.append( len(E) )
  E.append( tsurf.Edge(0, len(V)-1, tsurf.BD, None) )
  h_lens.append( top_cuff.length )
  
    
  [tsurf.Vertex([],[]) for i in xrange(2*len(seams))]



















