import hyp
import tsurf
import gsurf
import liftedsurf
import mobius

import math
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
    x,y = cc + cr*math.cos(Ma), cc + cr*math.sin(Ma)
    return not (x < br and y < bh)
  elif right_in:
    x,y = cc + cr*math.cos(ma), cc + cr*math.sin(ma)
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
    self.draw_width = None           #will be set by Configure
    self.draw_height = None          #will be set by Configure
    self.draw_middle = None          #will be set by Configure
    self.draw_complex_width = None   #will be set by Configure
    self.draw_complex_height = None  #will be set by Configure
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
        
  def draw_point(self, p, col, do_trans=True):
    t_p = (p if not do_trans else self.draw_trans(p))
    pp = self.draw_complex_to_canvas(t_p)
    di = self.canvas.create_oval(pp[0]-2,pp[1]-2,pp[0]+2,pp[1]+2,fill=col)
    self.drawing_items.append(di)
  
  def draw_geodesic_segment(self, gi, thickness=1, do_trans=True):
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
  
  def redraw(self):
    for di in self.drawing_items:
      self.canvas.delete(di)
    self.drawing_items = []
    
    self.drawn_v = [None for v in self.v]
    self.drawn_e = []
    self.drawn_t = []
    
    for vi in xrange(len(self.v)):
      for lifted_v in self.em_v[vi]:
        
    for ti in xrange(len(self.t)):
      
    
    
    for ei, (gi_left, gi_right) in enumerate(self.LS.em_e):
      if self.LS.e_single_lifts[ei]:
        self.draw_geodesic_segment(gi_left, thickness=1)
      else:
        self.e_lifts.append( (ei, gi_left.act_by_mobius(self.draw_trans), [1,0]) )
        self.draw_geodesic_segment(gi_left, thickness=2)
        self.e_lifts.append( (ei, gi_right.act_by_mobius(self.draw_trans), [0,1]) )
        self.draw_geodesic_segment(gi_right, thickness=2)
    for vi,V in enumerate(self.LS.em_v):
      c = self.draw_colors[vi%len(self.draw_colors)]
      for p in V:
        self.draw_point(p,c)
    
    self.propagate_lifts()
    
    for ei, gi in self.e_lifts:
      self.draw_geodesic_segment(gi, thickness=1, do_trans=False)
    
  
  def propagate_lifts(self):
    """using the current drawing stuff, propagate the triangles so that they 
    will cover the entire region; this assume self.e_lifts contains the 
    current boundary"""
    self.t_lifts = []
    self.v_lifts = []
    putative_lifts = [x for x in self.e_lifts]
    self.e_lifts = []
    e_lift_points = [[] for x in self.LS.e]
    while len(putative_lifts)>0:
      ei, gi, [on_L, on_R] = putative_lifts.pop()
      if self.disjoint_from_drawing(gi, do_trans=False):
        continue
      if any([hyp.same_float(gi.start, elsp[0], tol=1e-4) and \
              hyp.same_float(gi.end, elsp[1], tol=1e-4) for elsp in e_lift_points[ei]]):
        continue
      if on_L!=0:
        gLr = gi.reversed()
        ti,j = self.LS.e[ei].on_right
        print "Adding triangle ", self.LS.h_tris[ti]
        print "To edge", ti,j
        print "Length", gLr.length
        print "Elength", gLr.Euclidean_length()
        print "All start points: ", e_lift_points[ei]
        new_tri = self.LS.h_tris[ti].realize_along_gi(gLr, j)
        self.e_lifts.append( (ei, gi) )
        e_lift_points[ei].append( (gi.start, gi.end) )
      else:
        ti,j = self.LS.e[ei].on_left
        new_tri = self.LS.h_tris[ti].realize_along_gi(gi, j)
        self.e_lifts.append( (ei, gi )  )   
        e_lift_points[ei].append( (gi.start, gi.end) ) 
      self.t_lifts.append( (ti, new_tri) )
      jp1, jp2 = (j+1)%3, (j+2)%3
      pl1 = self.LS.t[ti].i_edges[jp1]
      pl1 = ((pl1.ind, new_tri.sides[jp1], [1,0]) if pl1.sign>0 else \
             (pl1.ind, new_tri.sides[jp1].reversed(), [0,1]))
      pl2 = self.LS.t[ti].i_edges[jp2]
      pl2 = ((pl2.ind, new_tri.sides[jp2], [1,0]) if pl2.sign>0 else \
             (pl2.ind, new_tri.sides[jp2].reversed(), [0,1]))
      putative_lifts.extend([pl1,pl2])
    
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

