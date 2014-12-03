import tsurf
import gsurf
import models


class LaidOutSurface(gsurf.GeoSurface):
  def __init__(self, GS):
    self.v, self.e, self.t = GS.v, GS.e, GS.t
    self.h_tris, self.h_lengths = GS.h_tris, GS.h_lengths
    self.relay()
  
  def attach_triangle_to_gi(self, ti, i, gi):
    """do the bookkeeping to attach edge index i in triangle index ti
    to the geodesic interval gi (it must be of the correct length)"""
    ei = self.t[ti].i_edges[i]
    #G = (gi if ei.sign>0 else gi.reversed())
    self.em_t[ti] = self.h_tris[ti].realize_along_gi(gi, i)
    #record what we've done
    for j,ei in enumerate(self.t[ti].i_edges):
      if ei.sign > 0:
        self.em_e[ei.ind][0] = self.em_t[ti].sides[j]
        self.em_v[self.e[ei.ind].source].append( self.em_t[ti].v[j] )
      else:
        self.em_e[ei.ind][1] = self.em_t[ti].sides[j].reversed()
        self.em_v[self.e[ei.ind].dest].append( self.em_t[ti].v[j] )
    return
        
        
  def relay(self):
    num_edges = len(self.e)
    self.em_v = [[] for _ in xrange(len(self.v))]
    #this gives the geodesic intervals on the left, right, respectively
    self.em_e = [[None,None] for _ in xrange(len(self.e))] 
    self.em_t = [None for _ in xrange(len(self.t))]
    #place edge 0 going straight up from 0
    p1 = complex(0,1)
    p2 = models.hyp_point_at_vertical_distance(complex(0,1), self.h_lengths[0])
    gi = models.HypGeodesicInterval(p1,p2).reversed()
    ti,j = self.e[0].on_right
    self.attach_triangle_to_gi(ti, j, gi)
    print "Attached triangle ", ti, " to ", gi
    while True:
      #find the edge which is (1) half glued and (2) the median Euclidean size
      half_glued_edges = []
      for i in xrange(num_edges):
        if None in self.em_e[i] and self.em_e[i] != [None,None]:
          if self.em_e[i][0] == None:
            half_glued_edges.append( (i, self.em_e[i][1].Euclidean_length()) )
          else:
            half_glued_edges.append( (i, self.em_e[i][0].Euclidean_length()) )
      if len(half_glued_edges)==0:
        break
      half_glued_edges.sort(key=lambda x: x[1])
      edge_to_glue = half_glued_edges[len(half_glued_edges)/2][0]
      side_to_glue = (0 if self.em_e[edge_to_glue][0] == None else 1) #0 means glue left
      ti,j = (self.e[edge_to_glue].on_right if side_to_glue == 1 else self.e[edge_to_glue].on_left)
      gi = self.em_e[edge_to_glue][1-side_to_glue].reversed()
      print "Attaching triangle", ti, j, "to", gi, "of length", gi.length
      print "Triangle: ", self.h_tris[ti]
      self.attach_triangle_to_gi(ti, j, gi)
    return
      
  def __repr__(self):
    return str(self)
  
  def __str__(self):
    ans = "Laid Out Surface:"
    ans += "Vertices: \n"
    for i,V in enumerate(self.v):
      ans += str(i) + ": " + str(V) + " Realized: " + str(self.em_v[i]) + "\n"
    ans += "\nEdges: \n"
    for i,E in enumerate(self.e):
      ans += str(i) + ": " + str(E) + " Length: " + str(self.h_lengths[i]) + "\n"
      ans += "  Realized: " + str(self.em_e[i][0]) + "\n"
    ans += "\nTriangles: \n"
    for i,T in enumerate(self.t):
      ans += str(i) + ": " + str(T) + "\n"
      ans += "    HypTri: " + str(self.h_tris[i]) + "\n"
      ans += "    Realized: " + str(self.em_t[i]) + "\n"
    return ans
    
    