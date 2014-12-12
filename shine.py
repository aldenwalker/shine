import hyp
import tsurf
import gsurf
import liftedsurf
import emsurf
import mobius
import R3

import math
import copy
import random
import Tkinter as tk

def arc_disjoint_from_box(cc, cr, ca1, ca2, br, bh ) :
  """is the circular arc at center cc, radius cr, between angles ca1, ca2
  disjoint from the box with bottom center at 0, horiz radius br, and 
  height bh?"""
  if cr == 'inf':
    return not (abs(cc) < br and min(ca1,ca2) < bh)
  left_e = cc - cr
  right_e = cc + cr
  left_in = abs(left_e) < br
  right_in = abs(right_e) < br
  if left_in and right_in:
    return False
  Ma = max(ca1, ca2)
  ma = min(ca1, ca2)
  if left_in:
    x,y = cc + cr*math.cos(Ma), cr*math.sin(Ma)
    return not (x < br and y < bh)
  elif right_in:
    x,y = cc + cr*math.cos(ma), cr*math.sin(ma)
    return not (-br < x and y < bh)
  else:
    return (math.sqrt( (br-cc)**2 + bh**2 ) < cr and math.sqrt( (-br-cc)**2 + bh**2 ) < cr)
    
def rand_bright_color():
  patt = random.choice( [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1)] )
  x = random.uniform( *( (0.3,0.6) if patt[0]==0 else (0.9,1)) )
  y = random.uniform( *( (0.3,0.6) if patt[1]==0 else (0.9,1)) )
  z = random.uniform( *( (0.3,0.6) if patt[2]==0 else (0.9,1)) )
  c = (int(255*x), int(255*y), int(255*z) )
  return '#%02x%02x%02x' % c
  

#############################################################################
#  Hyperbolic surface visualizer
#############################################################################

class HypSurfaceVisualizer:
  def __init__(self, master, LS):
    self.master = master
    self.master.geometry('500x500+100+100')
    
    self.canvas = tk.Canvas(self.master, borderwidth=0)
    self.canvas.bind('<Configure>', self.canvas_resize)
    self.canvas.bind('<Button-1>', self.canvas_click)
    
    self.button_quit = tk.Button(self.master, text = 'Quit', command=self.quit)
    self.button_rotate_left = tk.Button(self.master, text="RotL", command=lambda : self.rotate('left'))
    self.button_rotate_right = tk.Button(self.master, text='RotR', command=lambda : self.rotate('right'))
    self.draw_do_propagate = tk.IntVar()
    self.check_propagate = tk.Checkbutton(self.master, text='Propagate', variable=self.draw_do_propagate, command=self.redraw)
    self.draw_do_propagate.set(0)
    
    self.master.rowconfigure(3, weight=1)
    self.master.columnconfigure(0, weight=1)
    
    self.canvas.grid(column=0, row=0, rowspan=4, sticky=tk.W+tk.E+tk.N+tk.S)
    self.button_quit.grid(column=1, row=0, sticky=tk.N)
    self.button_rotate_left.grid(column=1, row=1)
    self.button_rotate_right.grid(column=2, row=1)
    self.check_propagate.grid(column=1, row=2)
    
    #remember the surface 
    self.LS = LS
    
    #set up the drawing
    #a complex number is scaled by self.draw_scale and then the y coordinate 
    #is taken from the bottom of the screen
    #and the x coord is taken from the middle
    self.draw_scale = 50.0
    self.draw_width = None           #All will be set by Configure
    self.draw_height = None        
    self.draw_middle = None         
    self.draw_complex_width = None  
    self.draw_complex_height = None 
    self.draw_extended_LS = None
    #everything is acted upon by self.draw_trans before drawing
    self.draw_trans = mobius.MobiusTrans(1,0,0,1)
    #the drawing the currently empty (but it will be filled by the canvas Configure)
    self.drawing_items = []
    
    self.draw_colors = [rand_bright_color() for _ in xrange(100)]
  
  
  ######################### drawing functions
  
  def draw_complex_to_canvas(self, z):
    z *= self.draw_scale
    y = self.draw_height - z.imag
    x = self.draw_middle + z.real
    return (x,y)
  
  def draw_canvas_to_complex(self, x,y):
    return (1.0/self.draw_scale)*complex(x-self.draw_middle, self.draw_height - y)
  
  def canvas_resize(self, event):
    self.canvas.config(background='#FFFFFF')
    self.canvas.config(height=event.height, width=event.width)
    self.draw_width = event.width
    self.draw_height = event.height
    self.draw_middle = self.draw_width/2.0
    self.draw_complex_width = self.draw_width / self.draw_scale
    self.draw_complex_height = self.draw_height / self.draw_scale
    self.redraw()
        
  def draw_point(self, p, col, r=2, do_trans=False):
    t_p = (p if not do_trans else self.draw_trans(p))
    pp = self.draw_complex_to_canvas(t_p)
    di = self.canvas.create_oval(pp[0]-r,pp[1]-r,pp[0]+r,pp[1]+r,fill=col)
    self.drawing_items.append(di)
  
  def draw_check_pt_in_window(self, p):
    br = self.draw_complex_width/2.0
    bh = self.draw_complex_height
    return -br < p.real and p.real < br and 0 < p.imag and p.imag < bh
  
  def draw_geodesic_segment(self, gi, thickness=1, do_trans=False):
    trans_gi = (gi if not do_trans else gi.act_by_mobius(self.draw_trans))
    #print "Drawing geodesic segment: ", gi
    #print "After trans: ", trans_gi
    if trans_gi.vertical:
      s = self.draw_complex_to_canvas(trans_gi.start)
      e = self.draw_complex_to_canvas(trans_gi.end)
      di = self.canvas.create_line(s[0], s[1], e[0], e[1], width=2)
      self.drawing_items.append(di)
      return
    ma = min(trans_gi.circ_angle1, trans_gi.circ_angle2)
    ma *= 180.0/math.pi
    Ma = max(trans_gi.circ_angle1, trans_gi.circ_angle2)
    Ma *= 180.0/math.pi
    c = trans_gi.circ_center
    r = trans_gi.circ_radius
    #scale to the canvas
    cc = self.draw_complex_to_canvas(c)
    rr = self.draw_scale * r
    bbox = (cc[0]-rr, cc[1]-rr, cc[0]+rr, cc[1]+rr)
    bbox = [int(b) for b in bbox]
    di = self.canvas.create_arc(bbox[0], bbox[1], bbox[2], bbox[3], style=tk.ARC, start=ma, extent=(Ma-ma), width=thickness)
    #print "Drew arc at ", bbox, " with start, extent", ma, Ma-ma
    self.drawing_items.append(di)
  
  def draw_triangle(self, t, col):
      #make a polygon with several points per side
      pts = []
      for i in xrange(3):
        gi = t.sides[i]
        if gi.vertical:
          num_joints = 1
        else:
          total_angle = abs(gi.circ_angle2 - gi.circ_angle1)
          num_joints = int((total_angle/(2*math.pi))*20 + 1)
          EL = gi.Euclidean_length()
          num_joints = int(num_joints*EL+1)
        for j in xrange(num_joints):
          pts.append( gi.pt_along( float(j)/float(num_joints) ) )
      cpts = [self.draw_complex_to_canvas(p) for p in pts]
      cpts = [x for P in cpts for x in P]
      di = self.canvas.create_polygon(*cpts, fill=col)
      self.drawing_items.append(di)
  
  def propagate_surface(self):
    v_stack = [(i,v) for i,v in enumerate(self.extended_LS.em_v) if None in v.i_tris]
    while len(v_stack)>0:
      #print "Stack: ", v_stack
      (i,v) = v_stack.pop()
      if (not self.draw_check_pt_in_window(v.pt)) or         \
         (v.pt.imag < 0.25)                         or         \
         (None not in v.i_tris):
        continue
      j = 0
      IT = v.i_tris
      LIT = len(IT)
      while IT[j] == None:
        j = (j+1)%LIT
      while IT[j] != None:
        j = (j+1)%LIT
      while IT[j] == None:
        ei = v.i_edges[j] #note this is the edge of the previous triangle
        old_n_vs = len(self.extended_LS.em_v)
        self.extended_LS.lift_triangle_to_lifted_edge(ei)
        if len(self.extended_LS.em_v) > old_n_vs:
          v_stack.append( (old_n_vs, self.extended_LS.em_v[old_n_vs]) )
        j = (j+1)%LIT
    return
  
  def redraw(self):
    for di in self.drawing_items:
      self.canvas.delete(di)
    self.drawing_items = []
    
    self.extended_LS = copy.deepcopy(self.LS)
    self.extended_LS.act_by_mobius(self.draw_trans)
    
    if self.draw_do_propagate.get()==1:
      self.propagate_surface()
    
    for ti,t in enumerate(self.LS.t):
      #for i,lti in enumerate(self.extended_LS.t_lifts[ti]):
      #  self.draw_triangle(self.extended_LS.em_t[lti].t, ('#FFA0A0' if i==0 else '#FFF0F0'))
      lti = self.extended_LS.t_lifts[ti][0]
      self.draw_triangle(self.extended_LS.em_t[lti].t, '#FFA0A0')
    
    for i,e in enumerate(self.extended_LS.em_e):
      if len(self.LS.e_lifts[e.covered_e])==1:
        self.draw_geodesic_segment(e.gi, thickness=1)
      else:
        self.draw_geodesic_segment(e.gi, thickness=2)
    
    for i,v in enumerate(self.extended_LS.em_v):
      col = self.draw_colors[v.covered_v%len(self.draw_colors)]
      self.draw_point(v.pt, col, r=3)
    
  def disjoint_from_drawing(self, gi, do_trans=True):
    """True if the geodesic interval is disjoint from the drawing 
    (note this applies transformations, etc) (also true if it would be 
    really small)"""
    t_gi = (gi if not do_trans else gi.act_by_mobius(self.draw_trans))
    if t_gi.Euclidean_length() < 1:
      return True
    return arc_disjoint_from_box(t_gi.circ_center, t_gi.circ_radius,        \
                                 t_gi.circ_angle1, t_gi.circ_angle2,       \
                                 self.draw_complex_width/2.0,               \
                                 self.draw_complex_height)    
  
  ####################### signals
  
  def canvas_click(self, event):
    p = (event.x, event.y)
    #print "Clicked on ", p
    can_p = (self.canvas.canvasx(p[0]), self.canvas.canvasy(p[1]))
    #print "Canvas: ", can_p
    com_p = self.draw_canvas_to_complex(*can_p)
    can_center = (self.draw_width/2.0, self.draw_height/2.0)
    com_center = self.draw_canvas_to_complex(*can_center)
    #print can_center, com_center
    #print "Moving point", com_p, "to", com_center
    M = mobius.MobiusTrans.unit_tangent_action(com_p, 0, com_center, 0)
    self.draw_trans = M.compose(self.draw_trans)
    self.redraw()
  
  def rotate(self, dir):
    theta = (math.pi/40.0 if dir=='left' else -math.pi/40.0)
    M = mobius.MobiusTrans.unit_tangent_action(1j, 0, 1j, theta)
    self.draw_trans = self.draw_trans.compose(M)
    self.redraw()
  
  def quit(self):
    self.master.destroy()





############################################################################
# Embedded surface visualizer
#############################################################################

class EmSurfaceVisualizer:
  def __init__(self, master, ES):
    self.ES = ES
    
    self.master = master
    self.master.geometry('500x500+100+100')
    
    self.canvas = tk.Canvas(self.master, borderwidth=0)
    self.canvas.bind('<Configure>', self.canvas_resize)
    
    #self.button_quit = tk.Button(self.master, text = 'Quit', command=self.quit)
    self.button_rotate_left = tk.Button(self.master, text="Rot left", command=lambda : self.rotate('left'))
    self.button_rotate_right = tk.Button(self.master, text='Rot right', command=lambda : self.rotate('right'))
    self.button_rotate_forward = tk.Button(self.master, text='Rot forw', command=lambda : self.rotate('forward'))
    self.button_rotate_back = tk.Button(self.master, text='Rot back', command=lambda : self.rotate('back'))
    self.button_subdivide = tk.Button(self.master, text='Subdivide', command=self.subdivide)
    self.button_flow = tk.Button(self.master, text='Flow', command=self.flow)
    self.draw_do_mesh = tk.IntVar()
    self.check_do_mesh = tk.Checkbutton(self.master, text='Draw mesh', variable=self.draw_do_mesh, command=self.redraw)
    self.draw_do_mesh.set(1)
    
    self.master.rowconfigure(4, weight=1)
    self.master.columnconfigure(0, weight=1)
    
    self.canvas.grid(column=0, row=0, rowspan=5, sticky=tk.W+tk.E+tk.N+tk.S)
    #self.button_quit.grid(column=1, row=0, sticky=tk.N)
    self.button_rotate_left.grid(column=1, row=0)
    self.button_rotate_right.grid(column=2, row=0)
    self.button_rotate_forward.grid(column=1, row=1)
    self.button_rotate_back.grid(column=2, row=1)
    self.button_subdivide.grid(column=1, row=2)
    self.button_flow.grid(column=2, row=2)
    self.check_do_mesh.grid(column=1, row=3)
    
    self.draw_viewer = R3.ProjectionViewer( R3.Vector([2,2,1]),              \
                                            R3.Vector([-2,-2,-1]),              \
                                            [R3.Vector([0,0,3]), R3.Vector([1,1,-3])] )
    self.draw_transformation = R3.Matrix([[1,0,0],[0,1,0],[0,0,1]])
    self.draw_center = None
    self.draw_scale = 100.0
    self.draw_plane_hr = None
    self.draw_plane_wr = None
    self.draw_width = None
    self.draw_height = None
    self.drawing_items = []
    
  def draw_plane_to_canvas(self, pt):
    scaled = (self.draw_scale*pt[0], self.draw_scale*pt[1])
    return (self.draw_center[0] + scaled[0], self.draw_center[1] - scaled[1])
    
  def quit(self):
    self.master.destroy()
  
  def canvas_resize(self, event):
    self.canvas.config(background='#FFFFFF')
    self.canvas.config(height=event.height, width=event.width)
    self.draw_width = event.width
    self.draw_height = event.height
    self.draw_plane_wr = event.width/(2.0*self.draw_scale)
    self.draw_plane_hr = event.height/(2.0*self.draw_scale)
    self.draw_center = (self.draw_width/2, self.draw_height/2)
    self.redraw()
  
  def redraw(self):
    for di in self.drawing_items:
      self.canvas.delete(di)
    self.drawing_items = []
    acted_on_T = [[self.draw_transformation(x) for x in t] for t in self.ES.em_t]
    pT = self.draw_viewer.project_triangles(acted_on_T)
    outline = ('black' if self.draw_do_mesh.get()==1 else '')
    for pt,am in pT:
      canvas_coords = [self.draw_plane_to_canvas(x) for x in pt]
      flat_coord_list = [x for p in canvas_coords for x in p]
      grayscale = int(128*(am+1))
      rgb = '#%02x%02x%02x' % (grayscale, grayscale, grayscale)
      #print "Drawing triangle: ", canvas_coords
      #print "Amount: ", rgb
      di = self.canvas.create_polygon(*flat_coord_list, fill=rgb, outline=outline)
      self.drawing_items.append(di)
  
  def subdivide(self):
    self.ES.subdivide()
    self.redraw()
    
  def flow(self):
    self.ES.flow()
    self.redraw()
  
  def rotate(self, dir):
    if dir == 'left' or dir == 'left':
      ang = (math.pi/6 if dir=='left' else -math.pi/6)
      self.draw_transformation = R3.Matrix([[math.cos(ang), -math.sin(ang), 0], \
                                            [math.sin(ang), math.cos(ang), 0],  \
                                            [0,0,1]])*self.draw_transformation
    else:
      ang = (-math.pi/6 if dir=='forward' else math.pi/6)
      self.draw_transformation = R3.Matrix([[math.cos(ang), 0, -math.sin(ang)], \
                                            [0,1,0],  \
                                            [math.sin(ang), 0, math.cos(ang)]])*self.draw_transformation
    self.redraw()

def visualize_em_surface(ES):
  root = tk.Tk()
  vs = EmSurfaceVisualizer(root, ES)
  root.mainloop()



def visualize_hyp_surface(LS):
  root = tk.Tk()
  vs = HypSurfaceVisualizer(root, LS)
  root.mainloop()



def test():
  T = tsurf.TopSurface(method='polygon', w='abABcdCD')
  G = gsurf.GeoSurface.geometrize_tsurf(T)
  L = liftedsurf.LiftedSurface.lift_gsurf(G)
  visualize_surface(L)
  




