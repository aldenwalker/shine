import tsurf
import gsurf
import models


class LaidOutSurface(gsurf.GeoSurface):
  def __init__(self, GS):
    self.v, self.e, self.f = GS.v, GS.e, GS.f
    self.h_tris, self.h_lengths = GS.h_tris, GS.h_lengths
    self.relay()
  
  def relay(self):
    num_edges = len(self.e)
    self.em_v = [[] for _ in xrange(len(self.v))]
    self.em_e = [[None,None] for _ in xrange(len(self.e))]
    self.em_f = [None for _ in xrange(len(self.f))]
    #place edge 0 going straight up from 0
    p1 = complex(0,1)
    p2 = models.hyp_point_at_vertical_distance(complex(0,1), self.h_lengths[0])
    self.em_e[0] = models.HypGeodesicInterval(p1,p2)
    #place the triangle on the right
    fi,j = self.e[0].on_right
    self.em_f[fi] = self.h_tris[fi].realize_along_gi(self.em_e[0].reversed(), j)
    #record what we've done
    for ei,j in enumerate(self.f[fi].i_edges):
      if ei.sign > 0:
        self.em_e[ei.ind][0] = self.em_f[fi].em_sides[j]
        self.em_v[self.e[ei.ind].source].append( self.em_f[fi].em_v[j] )
      else:
        self.em_e[ei.ind][1] = self.em_f[fi].sides[j].reversed()
        self.em_v[self.e[ei.ind].dest].append( self.em_f[fi].em_v[j] )
    
    while True:
      #find the edge which is (1) half glued and (2) the median Euclidean size
      half_glued_edges = []
      for i in xrange(num_edges):
        if None in self.em_e[i] and self.em_e[i] != [None,None]:
          if self.em_e[i][0] == None:
            half_glued_edges.append( (i, self.em_e[1].Euclidean_length()) )
          else:
            half_glued_edges.append( (i, self.em_e[0].Euclidean_length()) )
      half_glued_edges.sort(key=lambda x: x[1])
      edge_to_glue = half_glued_edges[len(half_glued_edges)/2][0]
      side_to_glue = (0 if self.em_e[edge_to_glue][0] == None else 1) #0 means glue left
      
      
      
    
    
    