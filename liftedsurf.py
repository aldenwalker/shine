import gsurf
import hyp

import math

class LiftedSurface(gsurf.GeoSurface):
  def __init__(self, GS, ev, ee, et, e_single_lifts):
    self.v, self.e, self.t = GS.v, GS.e, GS.t
    self.h_tris, self.h_lengths = GS.h_tris, GS.h_lengths
    self.em_v = ev
    self.em_e = ee
    self.em_t = et
    self.e_single_lift = e_single_lifts
  
  @classmethod
  def lift_gsurf(cls, GS, method='vertices'):
    LS = cls(GS, None, None, None, None)
    LS.relay(method=method)
    return LS
  
  def clean_vertices(self):
    for i in xrange(len(self.em_v)):
      if self.em_v[i] == None:
        continue
      new_points = []
      for p in self.em_v[i]:
        if all([not hyp.same_float(p, x) for x in new_points]):
          new_points.append(p)
      self.em_v[i] = new_points
  
  def is_vertex_surrounded(self, vi):
    """are all adjacent triangles placed?"""
    return all([self.em_t[ti] != None for ti,j in self.v[vi].i_tris])
  
  def attach_triangle_to_gi(self, ti, i, gi):
    """do the bookkeeping to attach edge index i in triangle index ti
    to the geodesic interval gi (it must be of the correct length)"""
    ei = self.t[ti].i_edges[i]
    self.em_t[ti] = self.h_tris[ti].realize_along_gi(gi, i)
    #record what we've done
    for j,ei in enumerate(self.t[ti].i_edges):
      if ei.sign > 0:
        self.em_e[ei.ind][0] = self.em_t[ti].sides[j]
        self.em_v[self.e[ei.ind].source].append( self.em_t[ti].v[j] )
      else:
        self.em_e[ei.ind][1] = self.em_t[ti].sides[j]
        self.em_v[self.e[ei.ind].dest].append( self.em_t[ti].v[j] )
    return
        
        
  def relay(self, method='vertices'):
    num_edges = len(self.e)
    self.em_v = [[] for _ in xrange(len(self.v))]
    #this gives the geodesic intervals on the left, right, respectively
    self.em_e = [[None,None] for _ in xrange(len(self.e))] 
    self.em_t = [None for _ in xrange(len(self.t))]
    self.e_single_lifts = [None for _ in xrange(len(self.e))]
    
    if method=='edges':
      #place edge 0 going straight up from 0
      p1 = complex(0,1)
      p2 = hyp.point_at_vertical_distance(complex(0,1), self.h_lengths[0])
      gi = hyp.HypGeodesicInterval(p1,p2).reversed()
      ti,j = self.e[0].on_right
      self.attach_triangle_to_gi(ti, j, gi)
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
        #print "Edge to glue: ", edge_to_glue
        #print "Side to glude: ", side_to_glue
        #print "Attaching triangle", ti, j, "to", gi, "of length", gi.length
        #print "Triangle: ", self.h_tris[ti]
        self.attach_triangle_to_gi(ti, j, gi)
    
    elif method=='vertices':
      first_vertex = self.find_isolated_vertex()
      if first_vertex == None:
        first_vertex = self.highest_valence_vertex()
      gi = hyp.HypGeodesicInterval.from_pt_angle_dist(1j, math.pi/2.0, self.h_lengths[self.v[first_vertex].i_edges[0].ind])
      for i,ei in enumerate(self.v[first_vertex].i_edges):
        ti,j = self.v[first_vertex].i_tris[i]
        if self.em_t[ti] != None:
          break
        self.attach_triangle_to_gi(ti,j,gi)
        gi = self.em_t[ti].sides[(j-1)%3].reversed()
        #print "Placed triangle", ti,j, " around first vertex", first_vertex
      while True:
        self.clean_vertices()
        #vertex preference (1) isolated and appears once (2) appears once
        current_v = None
        partially_done = [i for i in xrange(len(self.v)) if not self.is_vertex_surrounded(i)]
        appears_once = [i for i in partially_done if self.em_v[i] != None and len(self.em_v[i])==1]
        isolated = [i for i in appears_once if self.is_isolated_vertex(i)]
        #print "Found partially done, appears once, and isolated:"
        #print partially_done, appears_once, isolated
        if len(isolated)>0:
          current_v = iso_appear_once[0]
        elif len(appears_once)>0:
          current_v = iso_appear_once[0]
        elif len(partially_done)>0:
          current_v = partially_done[0]
        else:
          break
        #find the first unplaced triangle, and go from there
        IT = self.v[current_v].i_tris
        LIT = len(IT)
        i=0
        while self.em_t[IT[i][0]]==None:
          i = (i+1)%LIT
        while self.em_t[IT[i][0]]!=None:
          i = (i+1)%LIT
        while self.em_t[ IT[i][0] ] == None:
          ti,j = IT[i]
          prev_ti, prev_j = IT[(i-1)%LIT]
          gi = self.em_t[prev_ti].sides[(prev_j-1)%3].reversed()
          self.attach_triangle_to_gi(ti, j, gi)
          #print "Placed triangle", ti,j, " around vertex", current_v
          i = (i+1)%LIT
      #end of vertices method
    self.clean_vertices()
    for ei in xrange(len(self.e)):
      self.e_single_lifts[ei] = hyp.same_float(self.em_e[ei][0].start, self.em_e[ei][1].end,tol=1e-6)
    return
      
  def __repr__(self):
    return str(self)
  
  def __str__(self):
    ans = "Lifted Surface:"
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
    
    