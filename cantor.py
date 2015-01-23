##########################################################################
# create a complement of the Cantor set to depth n (with 1+2^n boundaries)
##########################################################################

import hyp
import gsurf



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
  return ','.join([x.mathematica_string() for x in seams] + [x.mathematica_string('Red') for x in cuffs + [top_cuff]])


