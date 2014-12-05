import hyp
import tsurf
import gsurf
import liftedsurf
import emsurf
import mobius

import math
import copy
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
    


class SurfaceVisualizer:
  def __init__(self, master, LS):
    self.master = master
    self.master.geometry('500x500+100+100')
    
    self.canvas = tk.Canvas(self.master, borderwidth=0)
    self.canvas.bind('<Configure>', self.canvas_resize)
    self.canvas.bind('<Button-1>', self.canvas_click)
    
    self.button_quit = tk.Button(self.master, text = 'Quit', command=self.quit)
    self.button_rotate_left = tk.Button(self.master, text="RotL", command=lambda : self.rotate('left'))
    self.button_rotate_right = tk.Button(self.master, text='RotR', command=lambda : self.rotate('right'))
    
    self.master.rowconfigure(2, weight=1)
    self.master.columnconfigure(0, weight=1)
    
    self.canvas.grid(column=0, row=0, rowspan=3, sticky=tk.W+tk.E+tk.N+tk.S)
    self.button_quit.grid(column=1, row=0, sticky=tk.N)
    self.button_rotate_left.grid(column=1, row=1)
    self.button_rotate_right.grid(column=2, row=1)
    
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
    self.draw_colors = ['red','green','blue','cyan','yellow','magenta']
  
  
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
        (v.pt.imag < 0.1)                         or         \
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
    
    self.propagate_surface()
    
    for ti,t in enumerate(self.extended_LS.t):
      lti = self.extended_LS.t_lifts[ti][0]
      lt = self.extended_LS.em_t[lti]
      col = '#FFA0A0'
      self.draw_triangle(lt.t, col)
    
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



def visualize_surface(LS):
  root = tk.Tk()
  vs = SurfaceVisualizer(root, LS)
  root.mainloop()



def test():
  T = tsurf.TopSurface(method='polygon', w='abABcdCD')
  G = gsurf.GeoSurface.geometrize_tsurf(T)
  L = liftedsurf.LiftedSurface.lift_gsurf(G)
  visualize_surface(L)
  




