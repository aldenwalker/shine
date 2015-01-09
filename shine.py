import hyp
import tsurf
import gsurf
import emsurf
import mobius
import R2
import R3
from signedind import SignedInd as SI

import math
import copy
import random
import collections

import Tkinter as tk
import tkFileDialog
import tkMessageBox
import tkColorChooser

def clamp(x):
  """clamp values between 0 and 255"""
  if x < 0:
    return 0
  if x > 255:
    return 255
  return x

def whiten_color(col):
  """make a color string more white"""
  c1 = int(col[1:3], 16)
  c2 = int(col[3:5], 16)
  c3 = int(col[5:7], 16)
  c1 = c1 + int((255-c1)/1.5)
  c2 = c2 + int((255-c2)/1.5)
  c3 = c3 + int((255-c3)/1.5)
  return '#%02x%02x%02x' % (c1, c2, c3)

def get_rgb_01(s):
  """get rgb 01 values from a string"""
  c1 = int(s[1:3], 16)
  c2 = int(s[3:5], 16)
  c3 = int(s[5:7], 16)
  return [float(x)/256 for x in [c1,c2,c3]]

def remove_duplicate_floats(L):
  """removes duplicate floats in a list -- it assumes they are sorted"""
  i = 0
  while i < len(L)-1:
    if L[i+1]-L[i] < 1e-10:
      del L[i+1]
    else:
      i += 1
  return None
  

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
          pts.append( gi.pt_along_euclidean( float(j)/float(num_joints) ) )
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
    
    self.draw_viewer = R3.ProjectionViewer( R3.Vector([2,-2,1]),              \
                                            R3.Vector([-2,2,-1]),              \
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
    for ell in self.ES.loops:
      this_loop = self.ES.loops[ell]
      for i,ei in enumerate(this_loop.edges):
        #print ei
        v11 = self.ES.e[ei.ind].source
        v12 = self.ES.e[ei.ind].dest
        v1 = self.ES.em_v[v11]*this_loop.edge_coords[i] + self.ES.em_v[v12]*(1-this_loop.edge_coords[i])
        v1_a = self.draw_transformation(v1)
        v1_p = self.draw_viewer.project_point_flat(v1_a)
        v1_pc = self.draw_plane_to_canvas(v1_p)
        ip1 = (i+1)%len(this_loop.edges)
        eip1 = this_loop.edges[ip1]
        #print eip1, eip1.ind, type(eip1.ind)
        v21 = self.ES.e[eip1.ind].source
        v22 = self.ES.e[eip1.ind].dest
        v2 = self.ES.em_v[v21]*this_loop.edge_coords[ip1] + self.ES.em_v[v22]*(1-this_loop.edge_coords[ip1])
        v2_a = self.draw_transformation(v2)
        v2_p = self.draw_viewer.project_point_flat(v2_a)
        v2_pc = self.draw_plane_to_canvas(v2_p)
        #print "Drawing line ", v1, " to ", v2, "transformed: ", v1_pc, " to ", v2_pc
        di = [self.canvas.create_line(v1_pc[0], v1_pc[1], v2_pc[0], v2_pc[1], width=3), \
              self.canvas.create_oval(v1_pc[0]-2, v1_pc[1]-2, v1_pc[0]+2, v1_pc[1]+2, fill='#FF0000'), \
              self.canvas.create_oval(v2_pc[0]-2, v2_pc[1]-2, v2_pc[0]+2, v2_pc[1]+2, fill='#FF0000')]
        self.drawing_items.extend(di)
        
  
  def subdivide(self):
    self.ES.subdivide()
    self.redraw()
    
  def flow(self):
    self.ES.flow()
    self.redraw()
  
  def rotate(self, dir):
    if dir == 'left' or dir == 'right':
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





#############################################################################
# Combination visualizer
#############################################################################
class Shine:
  def __init__(self, parent):
    
    #set up the data variables
    self.ES = None
    self.GS = None
    self.LS = None
    
    #set up the window
    self.parent = parent
    self.parent.geometry('500x500+100+100')
    self.parent.title('Shine')
    
    self.emsurf_frame = tk.Frame(self.parent)#, bg='#FF0000')
    self.emsurf_displayer = ShineEmSurfDisplay(self.emsurf_frame, self)
    
    self.loop_frame = tk.Frame(self.parent)#, bg='#0000FF')
    self.loop_displayer = ShineLoopDisplay(self.loop_frame, self)
    
    self.do_show_lift = tk.IntVar()
    self.do_show_lift.set(0)
    self.liftedsurf_displayer = None
    
    self.parent.rowconfigure(0, weight=1)
    self.parent.columnconfigure(0, weight=1)
    
    self.emsurf_frame.grid(column=0, row=0, sticky=tk.W+tk.E+tk.N+tk.S)
    self.loop_frame.grid(column=1, row=0, sticky=tk.W+tk.E+tk.N+tk.S)
    
    #create the menu bar
    self.menubar = tk.Menu(parent)
    self.filemenu = tk.Menu(self.menubar, tearoff=0)
    self.filemenu.add_command(label='Open', command=self.open_file)
    self.filemenu.add_command(label='Save session')
    self.filemenu.add_command(label='Export', command=self.export)
    self.filemenu.add_command(label='Quit', command=self.parent.destroy)
    self.menubar.add_cascade(label='File', menu=self.filemenu)
    self.viewmenu = tk.Menu(self.menubar, tearoff=0)
    self.viewmenu.add_checkbutton(label='Mesh', variable=self.emsurf_displayer.draw_do_mesh, command=self.emsurf_displayer.canvas_redraw)
    self.viewmenu.add_checkbutton(label='Shading', variable=self.emsurf_displayer.draw_do_shading, command=self.emsurf_displayer.canvas_redraw)
    self.viewmenu.add_checkbutton(label='Lifted surface', variable=self.do_show_lift, command=self.swap_show_lift) 
    self.menubar.add_cascade(label='View', menu=self.viewmenu)
    self.actionmenu = tk.Menu(self.menubar, tearoff=0)
    self.actionmenu.add_command(label='Subdivide', command=self.subdivide)
    self.actionmenu.add_command(label='Flow', command=self.flow)
    self.actionmenu.add_command(label='Add loop', command=self.loop_displayer.add_loop)
    self.menubar.add_cascade(label='Actions', menu=self.actionmenu)
    
    self.parent.config(menu=self.menubar)
  
  def open_file(self):
    filename = tkFileDialog.askopenfilename(parent=self.parent)
    if len(filename) < 5:
      return
    if filename[-3:] == 'pgr':
      self.ES = emsurf.EmbeddedSurface.from_planar_graph_file(filename)
      self.GS = gsurf.GeometricSurface.geometrize_tsurf(self.ES)
      self.LS = gsurf.LiftedSurface.lift_gsurf(self.GS)
    else:
      tkMessageBox.showerror(title='Error', message='File does not have a .pgr extension', parent=self.parent)
      return
    
    self.reset()
  
  def swap_show_lift(self):
    if self.do_show_lift.get() == 1:
      self.LS = gsurf.LiftedSurface.lift_gsurf(self.GS)
      self.liftedsurf_displayer = ShineHypSurfaceDisplay(self)
      self.liftedsurf_displayer.window.protocol("WM_DELETE_WINDOW", self.kill_lift)
    else:
      self.kill_lift()
  
  def kill_lift(self):
    self.liftedsurf_displayer.window.destroy()
    self.liftedsurf_displayer = None
    self.do_show_lift.set(0)
    
  def reset(self):
    self.emsurf_displayer.reset()
    self.loop_displayer.reset()
    if self.liftedsurf_displayer != None:
      self.liftedsurf_displayer.reset()
  
  def subdivide(self):
    if self.ES == None:
      return
    sub_data = self.ES.subdivide()
    self.GS.subdivide()
    for ell in self.loop_displayer.loops:
      ell.subdivide(*sub_data)
      
    self.emsurf_displayer.canvas_redraw()
    if self.do_show_lift.get() == 1:
      self.LS = gsurf.LiftedSurface.lift_gsurf(self.GS)
      self.liftedsurf_displayer.canvas_redraw()
  
  def flow(self):
    if self.ES == None:
      return
    self.ES.flow()
    self.emsurf_displayer.canvas_redraw()
  
  def export(self):
    if self.ES == None:
      return
    dialog = tk.Toplevel(master=self.parent)
    dialog.title('Export')
    parent_location = (self.parent.winfo_rootx(), self.parent.winfo_rooty())
    dialog.geometry('+%d+%d' % (parent_location[0], parent_location[1]))
    dialog.focus_set()
    dialog.grab_set()
    
    export_as = tk.StringVar()
    export_as.set('eps')
    W_export_as_frame = tk.Frame(dialog)
    W_export_as_label = tk.Label(W_export_as_frame, text='Export as:')
    W_export_svg = tk.Radiobutton(W_export_as_frame, text='svg', variable=export_as, value='svg')
    W_export_eps = tk.Radiobutton(W_export_as_frame, text='eps', variable=export_as, value='eps')
    W_export_as_label.grid(row=0,column=0)
    W_export_eps.grid(row=0,column=1)
    W_export_svg.grid(row=0,column=2)
    
    
    W_surface_option_frame = tk.LabelFrame(dialog, text='Surface options')
    surface_outline = tk.IntVar()
    surface_outline.set(0)
    
    def change_surface_style():
      W_surface_outline_smooth.config(state=(tk.DISABLED if surface_outline.get() == 0 else tk.NORMAL))
      W_surface_mesh_shading.config(state=(tk.DISABLED if surface_outline.get() == 1 else tk.NORMAL))
      W_surface_mesh_mesh.config(state=(tk.DISABLED if surface_outline.get() == 1 else tk.NORMAL))
    
    W_surface_mesh = tk.Radiobutton(W_surface_option_frame, text='Mesh', variable=surface_outline, value=0, command=change_surface_style)
    surface_mesh_shading = tk.IntVar()
    surface_mesh_shading.set(1)
    W_surface_mesh_shading = tk.Checkbutton(W_surface_option_frame, text='Shading', variable=surface_mesh_shading)
    surface_mesh_mesh = tk.IntVar()
    surface_mesh_mesh.set(1)
    W_surface_mesh_mesh = tk.Checkbutton(W_surface_option_frame, text='Mesh', variable=surface_mesh_mesh)
    surface_mesh_outline = tk.IntVar()
    surface_mesh_outline.set(1)
    W_surface_mesh_outline = tk.Checkbutton(W_surface_option_frame, text='Thick outline', variable=surface_mesh_outline)
    
    W_surface_outline = tk.Radiobutton(W_surface_option_frame, text='Outline', variable=surface_outline, value=1, command=change_surface_style)
    surface_outline_smooth = tk.IntVar()
    surface_outline_smooth.set(1)
    W_surface_outline_smooth = tk.Checkbutton(W_surface_option_frame, text='Smooth', variable=surface_outline_smooth, state=tk.DISABLED)
    
    W_surface_mesh.grid(row=0, column=0, columnspan=2,sticky=tk.W)
    W_surface_mesh_shading.grid(row=1, column=1, sticky=tk.W, padx=10)
    W_surface_mesh_mesh.grid(row=2, column=1,sticky=tk.W, padx=10)
    W_surface_mesh_outline.grid(row=3, column=1,sticky=tk.W, padx=10)
    W_surface_outline.grid(row=4, column=0, columnspan=2,sticky=tk.W)
    W_surface_outline_smooth.grid(row=5, column=1,sticky=tk.W, padx=10)
    
    W_loop_option_frame = tk.LabelFrame(dialog, text='Loop options')
    loop_smooth = tk.IntVar()
    loop_smooth.set(0)
    W_loop_smooth = tk.Checkbutton(W_loop_option_frame, text='Smooth', variable=loop_smooth)
    
    W_loop_smooth.grid(row=0, column=0)
    
    dialog.rowconfigure(2, weight=1)
    dialog.columnconfigure(0, weight=1)
    dialog.columnconfigure(1, weight=1)
    
    W_export_as_frame.grid(row=0, column=0, sticky=tk.W)
    W_surface_option_frame.grid(row=1, column=0, sticky=tk.N+tk.E+tk.S+tk.W)
    W_loop_option_frame.grid(row=1, column=1, sticky=tk.N+tk.E+tk.S+tk.W)
    
    OK = tk.IntVar()
    OK.set(0)
    def ok():
      OK.set(1)
      dialog.destroy()
    
    W_OK = tk.Button(dialog, text='OK', command=ok)
    W_cancel = tk.Button(dialog, text='Cancel', command=dialog.destroy)
    
    W_OK.grid(row=2, column=0)
    W_cancel.grid(row=2, column=1)
    
    dialog.wait_window(dialog)
    
    if OK.get()==0:
      return
    
    filename = tkFileDialog.asksaveasfilename(parent=self.parent)
    
    f = (self.emsurf_displayer.write_svg if export_as.get() == 'svg' else self.emsurf_displayer.write_eps)
    
    f(filename, surface_outline=(surface_outline.get()==1), \
                                 surface_mesh_mesh=(surface_mesh_mesh.get()==1), \
                                 surface_mesh_shading=(surface_mesh_shading.get()==1), \
                                 surface_mesh_outline=(surface_mesh_outline.get()==1), \
                                 surface_outline_smooth=(surface_outline_smooth.get()==1), \
                                 loop_smooth=(loop_smooth.get()==1))
    
    

############################################################################
# subvisualizer based on the embedded surface visualizer
############################################################################
class ShineEmSurfDisplay:
  def __init__(self, tk_parent, shine_parent):
    self.tk_parent = tk_parent
    self.shine_parent = shine_parent
    
    # set up the window
    self.tk_parent.rowconfigure(5, weight=1)
    self.tk_parent.columnconfigure(0, weight=1)
    
    self.canvas = tk.Canvas(self.tk_parent, borderwidth=0)
    self.canvas.bind('<Configure>', self.canvas_resize)
    self.canvas.grid(column=0, row=0, rowspan=6, columnspan=4, sticky=tk.W+tk.E+tk.N+tk.S)
    
    self.draw_do_mesh = tk.IntVar()
    self.draw_do_shading = tk.IntVar()
    self.draw_do_mesh.set(1)
    self.draw_do_shading.set(1)
    
    self.rotate_vert_ccw_button = tk.Button(self.tk_parent, text='>', command=lambda : self.rotate('vert_ccw'))
    self.rotate_vert_cw_button = tk.Button(self.tk_parent, text='<', command=lambda : self.rotate('vert_cw'))
    self.rotate_horiz_ccw_button = tk.Button(self.tk_parent, text='v', command=lambda : self.rotate('horiz_ccw'))
    self.rotate_horiz_cw_button = tk.Button(self.tk_parent, text='^', command=lambda : self.rotate('horiz_cw'))
    self.zoom_in_button = tk.Button(self.tk_parent, text='+', command=lambda : self.zoom('in'))
    self.zoom_out_button = tk.Button(self.tk_parent, text='-', command=lambda : self.zoom('out'))
    self.subdivide_button = tk.Button(self.tk_parent, text='S', command=self.shine_parent.subdivide)
    self.flow_button = tk.Button(self.tk_parent, text='F', command=self.shine_parent.flow)
    self.draw_do_mesh_check = tk.Checkbutton(self.tk_parent, text='Mesh', variable=self.draw_do_mesh, command=self.canvas_redraw)
    self.draw_do_shading_check = tk.Checkbutton(self.tk_parent, text='Shading', variable=self.draw_do_shading, command=self.canvas_redraw)
    self.rotate_vert_ccw_button.grid(row=1, column=3)
    self.rotate_vert_cw_button.grid(row=1, column=1)
    self.rotate_horiz_ccw_button.grid(row=2, column=2)
    self.rotate_horiz_cw_button.grid(row=0, column=2)
    self.zoom_in_button.grid(row=0, column=1)
    self.zoom_out_button.grid(row=0, column=3)
    self.subdivide_button.grid(row=2, column=1)
    self.flow_button.grid(row=2, column=3)
    #self.draw_do_mesh_check.grid(row=3, column=1, columnspan=2, sticky=tk.W)
    #self.draw_do_shading_check.grid(row=4, column=1, columnspan=2, sticky=tk.W)
    
    
    self.drawing_visible_triangles = []
    self.drawing_boundary = []
    self.drawing_items = []
    self.reset()
    
  def reset(self):
    self.draw_viewer = R3.ProjectionViewer( R3.Vector([0,-6,2]),              \
                                            R3.Vector([0,6,-2]),              \
                                            [R3.Vector([0,0,3]), R3.Vector([1,1,-3])] )
    self.draw_transformation = R3.Matrix([[1,0,0],[0,1,0],[0,0,1]])
    self.canvas_width = int(self.canvas.config()['width'][-1])
    self.canvas_height = int(self.canvas.config()['height'][-1])
    self.draw_canvas_center = (self.canvas_width/2, self.canvas_height/2)
    self.draw_plane_to_canvas_scale = 700.0
    self.draw_do_mesh.set(1)
    self.draw_do_shading.set(1)
    self.canvas_redraw()
    
  def canvas_resize(self, event):
    self.canvas.config(background='#FFFFFF')
    self.canvas.config(height=event.height, width=event.width)
    self.canvas_width = event.width
    self.canvas_height = event.height
    self.draw_canvas_center = (self.canvas_width/2, self.canvas_height/2)
    self.canvas_redraw()
  
#   def swap_mesh(self):
#     self.draw_do_mesh.set(1-self.draw_do_mesh.get())
#     self.canvas_redraw()
#     
#   def swap_shading(self):
#     self.draw_do_shading.set(1-self.draw_do_shading.get())
#     self.canvas_redraw()
    
  def draw_plane_to_canvas(self, pt):
    scaled = (self.draw_plane_to_canvas_scale*pt[0], self.draw_plane_to_canvas_scale*pt[1])
    return (self.draw_canvas_center[0] + scaled[0], self.draw_canvas_center[1] - scaled[1])
  
  #########################################################################
  # this function recomputes the triangles and boundary
  #########################################################################
  def recompute_drawing(self):
    #act on everything
    T_acted_on = [[self.draw_transformation(x) for x in t] for t in self.shine_parent.ES.em_t]
    T_normals = [R3.triangle_normal(t) for t in T_acted_on]
    T_visible = [self.draw_viewer.faces_eye(t[0], T_normals[i]) for i,t in enumerate(T_acted_on)]
    T_visible_only = [T_acted_on[i] for i in xrange(len(T_visible)) if T_visible[i]]
    T_visible_only = self.draw_viewer.sort_triangles(T_visible_only)
    self.drawing_visible_triangles = T_visible_only
    self.draw_viewer.viewer_grid_init_triangles(T_visible_only)
    
    #find the edges which lie on the visible boundary
    edges_on_boundary = dict()
    visible_segments = []
    for ei in xrange(len(self.shine_parent.ES.e)):
      left_i, left_i_in_t = self.shine_parent.ES.e[ei].on_left
      right_i, right_i_in_t = self.shine_parent.ES.e[ei].on_right
      left_faces_eye = T_visible[left_i]
      right_faces_eye = T_visible[right_i]
      if left_faces_eye == right_faces_eye:
        continue
      #now find the visible part
      segment = [ T_acted_on[left_i][left_i_in_t], T_acted_on[left_i][(left_i_in_t+1)%3] ]
      Ti_near_segment, unused_nearby_segments = self.draw_viewer.viewer_grid_near_segment(segment)
      visible_subsegments = self.draw_viewer.visible_subsegments_t_values(segment, [T_visible_only[i] for i in Ti_near_segment])
      if visible_subsegments == None:
        continue
      if len(visible_subsegments) > 1:
        print "We have more than 1 visible subsegment?"
        print visible_subsegments
        raise ValueError("help")
      edges_on_boundary[ei] = len(visible_segments)
      #self.draw_viewer.viewer_grid_add_segment(visible_subsegments[0], len(visible_segments))
      visible_segments.append( ( visible_subsegments[0], R3.subsegment_from_t_values(segment, visible_subsegments[0]) ) )
    #join up the boundary edges into connected groups
    self.drawing_boundary = []
    while True:
      #get a starting edge
      #print "Current edges_on_boundary:", edges_on_boundary
      try:
        ei, si = edges_on_boundary.popitem()
        seg_t, seg = visible_segments[si]
      except KeyError:
        break
      #the stack records; do we try to go forwards or backwards along the 
      #edge (True=forward), and do we add the next one on to the end (True) or beginning
      stack = []
      if abs(seg_t[0]) < 1e-10:
        stack.append( (ei, False, False) )
      if abs(seg_t[1]-1.0) < 1e-10:
        stack.append( (ei, True, True) )
      current_path = collections.deque(seg)
      #print "Starting on edge", ei
      while len(stack)>0:
        #print "Current stack: ", stack
        #print "Current path:", current_path
        ei, direc, append_to_end = stack.pop()
        vi = (self.shine_parent.ES.e[ei].dest if direc else self.shine_parent.ES.e[ei].source)
        this_ei_in_v = (SI(ei,-1) if direc else SI(ei,1))
        v = self.shine_parent.ES.v[vi]
        putative_edge = None
        for di in v.i_edges:
          if di.ind in edges_on_boundary and di != this_ei_in_v:
            putative_edge = di
            break
        if putative_edge == None:
          #print "we've come to the end of this part"
          continue
        #print "Next putative edge", putative_edge
        putative_seg_t, putative_seg = visible_segments[edges_on_boundary[putative_edge.ind]]
        del edges_on_boundary[putative_edge.ind]
        #print "Putative seg_t and seg:", putative_seg_t, putative_seg
        if putative_edge.sign > 0:
          if putative_seg_t[0] > 1e-10:
            #print "Not appending"
            continue
          if abs(putative_seg_t[1]-1)<1e-10:
            stack.append( (putative_edge.ind, (putative_edge.sign==1), append_to_end) )
          if append_to_end:
            current_path.append( putative_seg[1] )
          else:
            current_path.appendleft( putative_seg[1] )
        else:
          if putative_seg_t[1] < 1-1e-10:
            #print "Not appending"
            continue
          if abs(putative_seg_t[0]) < 1e-10:
            stack.append( (putative_edge.ind, (putative_edge.sign==1), append_to_end) )
          if append_to_end:
            current_path.append( putative_seg[0] )
          else:
            current_path.appendleft( putative_seg[0] )  
      #we're done, so the current list is a good polygonal path
      self.drawing_boundary.append(R3.PolygonalPath(current_path))
    #print "Drawing boundary:\n"
    #for db in self.drawing_boundary:
    #  print db
      
      
    #recompute all the loop data
    for ell in self.shine_parent.loop_displayer.loops:
      if ell.show.get() == 0:
        continue
      ell.recompute_drawing()

  
  
  #######################################################################
  # recompute the drawing and redraw the canvas
  #######################################################################
  def canvas_redraw(self):
    #erase the canvas
    for di in self.drawing_items:
      self.canvas.delete(di)
    self.drawing_items = []
    
    #if there's no surface, show nothing
    if self.shine_parent.ES == None:
      self.drawing_items.append( self.canvas.create_text(*self.draw_canvas_center, text='(No surface; open a new surface with the file menu)') )
      return
    
    #this recomputes the visible items and the grid
    self.recompute_drawing()
    
    #now draw them to the screen
    # draw all the triangles
    triangle_outline = ('black' if self.draw_do_mesh.get()==1 else '')
    triangle_shading = (self.draw_do_shading.get()==1)
    if triangle_outline != '' or triangle_shading:
      for t in self.drawing_visible_triangles:
        pt, am = self.draw_viewer.project_triangle(t)
        canvas_coords = [self.draw_plane_to_canvas(x) for x in pt]
        flat_coord_list = [x for p in canvas_coords for x in p]
        grayscale = int(128*(am+1))
        rgb = ('#%02x%02x%02x' % (grayscale, grayscale, grayscale) if triangle_shading else '')
        #print "Drawing triangle: ", canvas_coords
        #print "Amount: ", rgb
        di = self.canvas.create_polygon(*flat_coord_list, fill=rgb, outline=triangle_outline)
        self.drawing_items.append(di)
    
    #Draw the hidden loop segments
    for ell in self.shine_parent.loop_displayer.loops:
      if ell.show.get() == 0:
        continue
      for kind, pp in ell.PPs:
        if kind == 'visible':
          continue
        col = whiten_color(ell.color)
        projected_v = [self.draw_plane_to_canvas(self.draw_viewer.project_point(v)) for v in pp.L]
        coords = [x for v in projected_v for x in v]
        di = self.canvas.create_line(*coords, width=3, fill=col)
        self.drawing_items.append(di)
    
    #draw the boundary
    for be in self.drawing_boundary:
      col = '#000000' #rand_bright_color()
      projected_v = [self.draw_plane_to_canvas(self.draw_viewer.project_point(v)) for v in  be.L]
      coords = [x for v in projected_v for x in v]
      di = self.canvas.create_line(*coords, width=3, fill=col) 
      self.drawing_items.append(di)
    
    #draw the visible loop segments
    for ell in self.shine_parent.loop_displayer.loops:
      if ell.show.get() == 0:
        continue
      for kind, pp in ell.PPs:
        if kind == 'hidden':
          continue
        col = ell.color
        projected_v = [self.draw_plane_to_canvas(self.draw_viewer.project_point(v)) for v in pp.L]
        coords = [x for v in projected_v for x in v]
        di = self.canvas.create_line(*coords, width=3, fill=col)
        self.drawing_items.append(di)
    
    #draw the grid
    # for i in xrange(self.draw_viewer.viewer_grid_num_horiz_boxes):
#       for j in xrange(self.draw_viewer.viewer_grid_num_vert_boxes):
#         pts = [ (self.draw_viewer.viewer_grid_ll[0] + ia*self.draw_viewer.viewer_grid_box_width, \
#                  self.draw_viewer.viewer_grid_ll[1] + ja*self.draw_viewer.viewer_grid_box_height) 
#                  for (ia,ja) in [(i,j), (i+1, j), (i+1,j+1), (i,j+1)]]
#         cpts = [self.draw_plane_to_canvas(x) for x in pts]
#         di = self.canvas.create_polygon(*cpts, outline='#0000FF', fill='')
#         self.drawing_items.append(di)
    
    
    
  
  ########################################################################
  # rotate the surface
  ########################################################################
  def rotate(self, dir):
    if dir == 'vert_ccw' or dir == 'vert_cw':
      ang = (math.pi/12 if dir=='vert_ccw' else -math.pi/12)
      self.draw_transformation = R3.Matrix([[math.cos(ang), -math.sin(ang), 0], \
                                            [math.sin(ang), math.cos(ang), 0],  \
                                            [0,0,1]])*self.draw_transformation
    else:
      ang = (math.pi/12 if dir=='horiz_ccw' else -math.pi/12)
      self.draw_transformation = R3.Matrix([[1,0,0],\
                                            [0, math.cos(ang), -math.sin(ang)], \
                                            [0, math.sin(ang), math.cos(ang)]])*self.draw_transformation
    self.canvas_redraw()
  
  #######################################################################
  # this computes the sum of the squares of the difference between the 
  # euclidean and hyeprbolic lengths of the edges
  #######################################################################
  def triangle_deviation(self, i):
    embedded_lengths = [ (self.shine_parent.ES.em_t[i][(j+1)%3] - self.shine_parent.ES.em_t[i][j]).norm() for j in xrange(3) ]
    hyp_lengths = [ self.shine_parent.GS.h_tris[i].lengths[j] for j in xrange(3) ]
    scaling_factor = hyp_lengths[0] / embedded_lengths[0]
    scaled_lens = [scaling_factor*x for x in embedded_lengths]
    deviation = sum([ (hyp_lengths[j]-scaled_lens[j])**2 for j in xrange(3) ])
    return deviation
  
  def zoom(self, dir):
    #self.draw_viewer.zoom( (0.1 if dir=='in' else 1.1) )
    self.draw_plane_to_canvas_scale = self.draw_plane_to_canvas_scale * (0.9 if dir=='out' else 1.1)
    self.canvas_redraw()
  
  ##########################################################################
  # write out an svg 
  ##########################################################################
  def write_svg(self, filename, surface_outline=True, \
                                surface_mesh_mesh=False, \
                                surface_mesh_shading=True, \
                                surface_mesh_outline=True, \
                                surface_outline_smooth=True, \
                                loop_smooth=True):
    f = open(filename, 'w')
    f.write('<?xml version="1.0"?>\n')
    f.write('<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN"\n')
    f.write('"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">\n')
    f.write('<svg width="' + str(self.canvas_width) + '" height="' + str(self.canvas_height) + '">\n')
    # draw all the triangles
    if not surface_outline:
      for t in self.drawing_visible_triangles:
        pt, am = self.draw_viewer.project_triangle(t)
        canvas_coords = [self.draw_plane_to_canvas(x) for x in pt]
        SVG_d = 'M ' + str(canvas_coords[0][0]) + ' ' + str(canvas_coords[0][1])
        for i in xrange(1,3):
          SVG_d += ' L ' + str(canvas_coords[i][0]) + ' ' + str(canvas_coords[i][1])
        SVG_d += ' z'
        SVG_command = '<path d="' + SVG_d + '"'
        if surface_mesh_mesh:
          SVG_command += ' stroke="black" stroke-width="1" stroke-linejoin="round"'
        else:
          SVG_command += ' stroke="none"'
        if surface_mesh_shading:
          grayscale = int(128*(am+1))
          rgb = '#%02x%02x%02x' % (grayscale, grayscale, grayscale)
          SVG_command += ' fill="' + rgb + '"'
        else:
          SVG_command += ' fill="none"'
        SVG_command += '/>\n'
        f.write(SVG_command)
    
    #Draw the hidden loop segments
    for ell in self.shine_parent.loop_displayer.loops:
      if ell.show.get() == 0:
        continue
      for kind, pp in ell.PPs:
        if kind == 'visible':
          continue
        col = whiten_color(ell.color)
        projected_v = [self.draw_plane_to_canvas(self.draw_viewer.project_point(v)) for v in pp.L]
        if loop_smooth:
          bp = R2.bezier_approximation([R2.Vector(v) for v in projected_v])
          #print "Doing hidden part; got ", len(bp), "breakpoints"
          p1, cp1, cp2, p2 = bp[0]
          SVG_d = 'M ' + str(p1[0]) + ' ' + str(p1[1]) + \
                  ' C ' + str(cp1[0]) + ' ' + str(cp1[1]) + ' ' +\
                          str(cp2[0]) + ' ' + str(cp2[1]) + ' ' +\
                          str(p2[0]) + ' ' + str(p2[1])
          for i in xrange(1,len(bp)):
            p1, cp1, cp2, p2 = bp[i]
            SVG_d += ' S ' + ' '.join(map(str, [cp2[0], cp2[1], p2[0], p2[1]]))
        else:
          SVG_d = 'M ' + str(projected_v[0][0]) + ' ' + str(projected_v[0][1])
          for i in xrange(1,len(projected_v)):
            SVG_d += ' L ' + str(projected_v[i][0]) + ' ' + str(projected_v[i][1])
        SVG_command = '<path d="' + SVG_d + '" stroke="' + col + '" stroke-width="3" stroke-linejoin="round" fill="none"/>\n'
        f.write(SVG_command)
    
    #draw the boundary
    if surface_outline or surface_mesh_outline:
      for be in self.drawing_boundary:
        #print "Drawing ", be
        col = '#000000' #rand_bright_color()
        projected_v = [self.draw_plane_to_canvas(self.draw_viewer.project_point(v)) for v in  be.L]
        if surface_outline and surface_outline_smooth:
          bp = R2.bezier_approximation([R2.Vector(v) for v in projected_v])
          #print "Doing boundary part; got ", len(bp), "breakpoints"
          p1, cp1, cp2, p2 = bp[0]
          SVG_d = 'M ' + str(p1[0]) + ' ' + str(p1[1]) + \
                  ' C ' + str(cp1[0]) + ' ' + str(cp1[1]) + ' ' +\
                          str(cp2[0]) + ' ' + str(cp2[1]) + ' ' +\
                          str(p2[0]) + ' ' + str(p2[1])
          for i in xrange(1,len(bp)):
            p1, cp1, cp2, p2 = bp[i]
            SVG_d += ' S ' + ' '.join(map(str, [cp2[0], cp2[1], p2[0], p2[1]]))
        else:
          SVG_d = 'M ' + str(projected_v[0][0]) + ' ' + str(projected_v[0][1])
          for i in xrange(1,len(projected_v)):
            SVG_d += ' L ' + str(projected_v[i][0]) + ' ' + str(projected_v[i][1])
        SVG_command = '<path d="' + SVG_d + '" stroke="' + col + '" stroke-width="3" stroke-linejoin="round" fill="none"/>\n'
        f.write(SVG_command)
    
    #draw the visible loop segments
    for ell in self.shine_parent.loop_displayer.loops:
      if ell.show.get() == 0:
        continue
      for kind, pp in ell.PPs:
        if kind == 'hidden':
          continue
        col = ell.color
        projected_v = [self.draw_plane_to_canvas(self.draw_viewer.project_point(v)) for v in pp.L]
        if loop_smooth:
          bp = R2.bezier_approximation([R2.Vector(v) for v in projected_v])
          #print "Doing visible loop; got ", len(bp), "breakpoints"
          p1, cp1, cp2, p2 = bp[0]
          SVG_d = 'M ' + str(p1[0]) + ' ' + str(p1[1]) + \
                  ' C ' + str(cp1[0]) + ' ' + str(cp1[1]) + ' ' +\
                          str(cp2[0]) + ' ' + str(cp2[1]) + ' ' +\
                          str(p2[0]) + ' ' + str(p2[1])
          for i in xrange(1,len(bp)):
            p1, cp1, cp2, p2 = bp[i]
            SVG_d += ' S ' + ' '.join(map(str, [cp2[0], cp2[1], p2[0], p2[1]]))
        else:
          SVG_d = 'M ' + str(projected_v[0][0]) + ' ' + str(projected_v[0][1])
          for i in xrange(1,len(projected_v)):
            SVG_d += ' L ' + str(projected_v[i][0]) + ' ' + str(projected_v[i][1])
        SVG_command = '<path d="' + SVG_d + '" stroke="' + col + '" stroke-width="3" stroke-linejoin="round" fill="none"/>\n'
        f.write(SVG_command)
    
    f.write('</svg>\n')
    f.close()
  
  
  #########################################################################
  # write out an eps
  #########################################################################
  def write_eps(self, filename, surface_outline=True, \
                                surface_mesh_mesh=False, \
                                surface_mesh_shading=True, \
                                surface_mesh_outline=True, \
                                surface_outline_smooth=True, \
                                loop_smooth=True):
    f = open(filename, 'w')
    f.write('%!PS-Adobe-2.0 EPSF-2.0\n')
    f.write('%%BoundingBox: ' + str(0) + ' ' +str(0) + ' ' + str(self.canvas_width) + ' ' + str(self.canvas_height) + '\n')
    f.write('1 setlinejoin\n')
    
    h = self.canvas_height
    
    # draw all the triangles
    if (not surface_outline) and (surface_mesh_mesh or surface_mesh_shading):
      for t in self.drawing_visible_triangles:
        pt, am = self.draw_viewer.project_triangle(t)
        canvas_coords = [self.draw_plane_to_canvas(x) for x in pt]
        EPS_com = str(canvas_coords[0][0]) + ' ' + str(h-canvas_coords[0][1]) + ' moveto\n'
        for i in xrange(1,3):
          EPS_com += str(canvas_coords[i][0]) + ' ' + str(h-canvas_coords[i][1]) + ' lineto\n'
        EPS_com += 'closepath\n'
        if surface_mesh_mesh and surface_mesh_shading:
          EPS_com += 'gsave\n'
        if surface_mesh_mesh:
          EPS_com += '0 0 0 setrgbcolor\n'
          EPS_com += 'stroke\n'
        if surface_mesh_shading:
          if surface_mesh_mesh:
            EPS_com += 'grestore\n'
          grayscale = 0.5*(am+1)
          EPS_com += str(grayscale) + ' ' + str(grayscale) + ' ' + str(grayscale) + ' setrgbcolor\n'
          EPS_com += 'fill\n'
        f.write(EPS_com)
    
    #Draw the hidden loop segments
    for ell in self.shine_parent.loop_displayer.loops:
      if ell.show.get() == 0:
        continue
      col = get_rgb_01(whiten_color(ell.color))
      f.write(str(col[0]) + ' ' + str(col[1]) + ' ' + str(col[2]) + ' setrgbcolor\n')
      for kind, pp in ell.PPs:
        if kind == 'visible':
          continue
        projected_v = [self.draw_plane_to_canvas(self.draw_viewer.project_point(v)) for v in pp.L]
        if loop_smooth:
          bp = R2.bezier_approximation([R2.Vector(v) for v in projected_v])
          p1, cp1, cp2, p2 = bp[0]
          EPS_com = str(p1[0]) + ' ' + str(h-p1[1]) + ' moveto\n'
          EPS_com += ' '.join(map(str, [cp1[0], h-cp1[1], cp2[0], h-cp2[1], p2[0], h-p2[1]])) + ' curveto\n'
          for i in xrange(1,len(bp)):
            p1, cp1, cp2, p2 = bp[i]
            EPS_com +=  ' '.join(map(str, [cp1[0], h-cp1[1], cp2[0], h-cp2[1], p2[0], h-p2[1]])) + ' curveto\n'
        else:
          EPS_com = str(projected_v[0][0]) + ' ' + str(h-projected_v[0][1]) + ' moveto\n'
          for i in xrange(1,len(projected_v)):
            EPS_com += str(projected_v[i][0]) + ' ' + str(h-projected_v[i][1]) + ' lineto\n'
        EPS_com += 'stroke\n'
        f.write(EPS_com)
    
    #draw the boundary
    if surface_outline or surface_mesh_outline:
      for be in self.drawing_boundary:
        #print "Drawing ", be
        f.write('0 0 0 setrgbcolor\n')
        projected_v = [self.draw_plane_to_canvas(self.draw_viewer.project_point(v)) for v in  be.L]
        if surface_outline and surface_outline_smooth:
          bp = R2.bezier_approximation([R2.Vector(v) for v in projected_v])
          p1, cp1, cp2, p2 = bp[0]
          EPS_com = str(p1[0]) + ' ' + str(h-p1[1]) + ' moveto\n'
          EPS_com += ' '.join(map(str, [cp1[0], h-cp1[1], cp2[0], h-cp2[1], p2[0], h-p2[1]])) + ' curveto\n'
          for i in xrange(1,len(bp)):
            p1, cp1, cp2, p2 = bp[i]
            EPS_com +=  ' '.join(map(str, [cp1[0], h-cp1[1], cp2[0], h-cp2[1], p2[0], h-p2[1]])) + ' curveto\n'
        else:
          EPS_com = str(projected_v[0][0]) + ' ' + str(h-projected_v[0][1]) + ' moveto\n'
          for i in xrange(1,len(projected_v)):
            EPS_com += str(projected_v[i][0]) + ' ' + str(h-projected_v[i][1]) + ' lineto\n'
        EPS_com += 'stroke\n'
        f.write(EPS_com)
    
    #draw the visible loop segments
    for ell in self.shine_parent.loop_displayer.loops:
      if ell.show.get() == 0:
        continue
      col = get_rgb_01(ell.color)
      f.write(str(col[0]) + ' ' + str(col[1]) + ' ' + str(col[2]) + ' setrgbcolor\n')
      for kind, pp in ell.PPs:
        if kind == 'hidden':
          continue
        col = ell.color
        projected_v = [self.draw_plane_to_canvas(self.draw_viewer.project_point(v)) for v in pp.L]
        if loop_smooth:
          bp = R2.bezier_approximation([R2.Vector(v) for v in projected_v])
          p1, cp1, cp2, p2 = bp[0]
          EPS_com = str(p1[0]) + ' ' + str(h-p1[1]) + ' moveto\n'
          EPS_com += ' '.join(map(str, [cp1[0], h-cp1[1], cp2[0], h-cp2[1], p2[0], h-p2[1]])) + ' curveto\n'
          for i in xrange(1,len(bp)):
            p1, cp1, cp2, p2 = bp[i]
            EPS_com +=  ' '.join(map(str, [cp1[0], h-cp1[1], cp2[0], h-cp2[1], p2[0], h-p2[1]])) + ' curveto\n'
        else:
          EPS_com = str(projected_v[0][0]) + ' ' + str(h-projected_v[0][1]) + ' moveto\n'
          for i in xrange(1,len(projected_v)):
            EPS_com += str(projected_v[i][0]) + ' ' + str(h-projected_v[i][1]) + ' lineto\n'
        EPS_com += 'stroke\n'
        f.write(EPS_com)
    
    f.write('\n')
    f.close()
  
  
  
  
###########################################################################
# the list of loops
###########################################################################
class ShineLoopDisplay:
  def __init__(self, tk_parent, shine_parent):
    self.tk_parent = tk_parent
    self.shine_parent = shine_parent
    
    self.title = tk.Label(self.tk_parent, text='Loops:')
    self.add_loop_button = tk.Button(self.tk_parent, text="+", command=self.add_loop)

    self.title.grid(column=0, row=0, sticky=tk.W)
    self.add_loop_button.grid(column=0, row=1, sticky=tk.W)
    
    self.loops = []
    self.frames = []
  
  def add_loop(self):
    if self.shine_parent.ES == None:
      return
    word, EP = self.add_loop_dialog()
    #print EP
    if EP == None:
      return
    self.add_loop_button.grid_remove()
    self.frames.append( tk.Frame(self.tk_parent, bd=2, relief=tk.RIDGE))#, bg = '#00FF00' ) )
    self.frames[-1].grid(column=0, row=len(self.frames))
    self.add_loop_button.grid(column=0, row=len(self.frames)+1, sticky=tk.W)
    self.loops.append( ShineLoop(self.frames[-1], self, EP, word=word) )
    self.shine_parent.emsurf_displayer.canvas_redraw()
    if self.shine_parent.liftedsurf_displayer != None:
      self.shine_parent.liftedsurf_displayer.canvas_redraw()
      
  
  def delete_loop(self, loop_to_delete):
    for i in xrange(len(self.loops)):
      if self.loops[i] == loop_to_delete:
        self.frames[i].destroy()
        self.add_loop_button.grid_remove()
        for j in xrange(i+1, len(self.frames)):
          self.frames[j].grid_remove()
          self.frames[j].grid(column=0, row=i+j)
        del self.loops[i]
        del self.frames[i]
        self.add_loop_button.grid(column=0, row=len(self.frames)+1, sticky=tk.W)
        break
    self.shine_parent.emsurf_displayer.canvas_redraw()
  
  def reset(self):
    for i in xrange(len(self.loops)-1, 0, -1):
      self.delete_loop(self.loops[i])

  def add_loop_dialog(self):
    dialog = tk.Toplevel(master=self.tk_parent)
    dialog.title('Add loop')
    parent_location = (self.tk_parent.winfo_rootx(), self.tk_parent.winfo_rooty())
    dialog.geometry('+%d+%d' % (parent_location[0], parent_location[1]))
    dialog.focus_set()
    dialog.grab_set()
    
    W_word_label = tk.Label(dialog, text='Create a loop as a product:')
    W_known_loops = tk.Label(dialog, text='Known loops: ' + str([ell for ell in self.shine_parent.ES.loops]) )
    W_word_input = tk.Entry(dialog)
    word_input = [None] #making it a list makes it visible to the function set_word
    def set_word():
      word_input[0] = W_word_input.get()
      dialog.destroy()
    W_word_go = tk.Button(dialog, text='Add from word', command=set_word )
    W_cancel = tk.Button(dialog, text='Cancel', command=dialog.destroy)
    
    W_word_label.grid(column=0, row=0, sticky=tk.W)
    W_known_loops.grid(column=0, row=1, sticky=tk.W)
    W_word_input.grid(column=0, row=2, sticky=tk.W)
    W_word_go.grid(column=0, row=3, sticky=tk.W)
    W_cancel.grid(column=1, row=3)
    
    dialog.wait_window(dialog)
    
    if word_input[0] != None:
      word_input = word_input[0]
      EP = self.shine_parent.ES.loop_from_word(word_input)
      return (word_input, EP)
    
    return None, None

############################################################################
# a single loop displayer
############################################################################
class ShineLoop:
  def __init__(self, tk_parent, shine_parent, EP, word=None):
    self.tk_parent = tk_parent
    self.shine_parent = shine_parent
    self.shine_main = shine_parent.shine_parent
    
    self.EP = EP
    self.PPs = None
    self.word = word
    self.color = '#000000'
    
    self.W_label = tk.Label(self.tk_parent, text='Word: ' + ('(unknown)' if self.word==None else self.word))
    self.show = tk.IntVar()
    self.show.set(1)
    self.W_show = tk.Checkbutton(self.tk_parent, text='Show', variable=self.show, command=self.shine_main.emsurf_displayer.canvas_redraw)
    self.W_color = tk.Button(self.tk_parent, command=lambda :self.set_color(tkColorChooser.askcolor()[1]), bg=self.color )
    self.W_geodesicify = tk.Button(self.tk_parent, text='Geoify', command=self.geodesicify)
    self.W_delete = tk.Button(self.tk_parent, text='X', command=lambda : self.shine_parent.delete_loop(self))
    
    self.W_label.grid(row=0, column=0, columnspan=3, sticky=tk.W)
    self.W_show.grid(row=1, column=0, columnspan=3, sticky=tk.W)
    self.W_color.grid(row=2, column=0, sticky=tk.W)
    self.W_geodesicify.grid(row=2, column=1)
    self.W_delete.grid(row=2, column=2)
  
  def set_color(self, c):
    self.color = c
    self.W_color.config(bg=self.color)
    self.shine_main.emsurf_displayer.canvas_redraw()
    if self.shine_main.liftedsurf_displayer != None:
      self.shine_main.liftedsurf_displayer.canvas_redraw()
      
  
  def subdivide(self, old_TS, vertices_from_edges, edges_from_edges, edges_from_tris, tris_from_tris):
    self.EP.subdivide(old_TS, vertices_from_edges, edges_from_edges, edges_from_tris, tris_from_tris)
  
  def geodesicify(self):
    #self.EP = self.shine_main.LS.geodesicify(self.EP)
    self.EP = self.shine_main.GS.geodesicify(self.EP)
    self.shine_main.emsurf_displayer.canvas_redraw()
    if self.shine_main.liftedsurf_displayer != None:
      self.shine_main.liftedsurf_displayer.canvas_redraw()
    
  ########################################################################
  # recompute the polygonal paths; this assumes that the shine_main 
  # triangles and draw viewer are current
  ########################################################################
  def recompute_drawing(self):
    self.PPs = []
    T = self.shine_main.emsurf_displayer.drawing_visible_triangles
    
    V = [self.shine_main.ES.along_edge(self.EP.edges[i].ind, self.EP.edge_coords[i]) for i in xrange(len(self.EP.edges))]
    V = [self.shine_main.emsurf_displayer.draw_transformation(v) for v in V]
    while (V[0]-V[-1]).norm() < 1e-10:
      del V[-1]
    i=0
    while i<len(V)-1:
      if (V[i+1]-V[i]).norm() < 1e-10:
        del V[i+1]
      else:
        i += 1
    lv = len(V)
    visible_V = lv*[None]
    for i,v in enumerate(V):
      Ti_near_point, unused_segments = self.shine_main.emsurf_displayer.draw_viewer.viewer_grid_near_point(v)
      visible_V[i] = (not self.shine_main.emsurf_displayer.draw_viewer.is_point_hidden([T[j] for j in Ti_near_point], v))
    #go until we first hit a gap between a hidden and non-hidden vertex
    first_vert = 0
    while first_vert < lv and (visible_V[first_vert] or not visible_V[(first_vert+1)%lv]):
      first_vert += 1
    if first_vert == lv:
      #there's no such gap
      if visible_V[0]:
        self.PPs.append( ['visible', R3.PolygonalPath(V)] )
      else:
        self.PPs.append( ['hidden', R3.PolygonalPath(V)] )
    fvp1 = (first_vert+1)%lv
    seg = [V[first_vert], V[fvp1]]
    Ti_near_segment, unused_segments = self.shine_main.emsurf_displayer.draw_viewer.viewer_grid_near_segment( seg )
    visible_subsegments = self.shine_main.emsurf_displayer.draw_viewer.visible_subsegments_t_values(seg, [T[i] for i in Ti_near_segment])
    if visible_subsegments == None:
      current_path = [R3.along_segment(seg, 1.0)]
    else:
      if len(visible_subsegments) != 1:
        raise ValueError("Wrong number of visible subsegments?")
      current_path = R3.subsegment_from_t_values(seg, visible_subsegments[0])   
    currently_visible = True
    current_vert = fvp1
    while True:
      cvp1 = (current_vert+1)%lv
      if visible_V[cvp1] != currently_visible:
        seg = [V[current_vert], V[cvp1]]
        Ti_near_segment, unused_segments = self.shine_main.emsurf_displayer.draw_viewer.viewer_grid_near_segment( seg )
        visible_subsegments = self.shine_main.emsurf_displayer.draw_viewer.visible_subsegments_t_values(seg, [T[i] for i in Ti_near_segment])
        if visible_subsegments == None:
          t_cut_point = (0.0 if currently_visible else 1.0)
        else:
          if len(visible_subsegments) != 1:
            raise ValueError("Wrong number of visible subsegments?")
          t_cut_point = (visible_subsegments[0][1] if currently_visible else visible_subsegments[0][0])
        current_path.append( R3.along_segment(seg, t_cut_point) )
        self.PPs.append( [ ('visible' if currently_visible else 'hidden'), R3.PolygonalPath(current_path) ] )
        if cvp1 == fvp1:
          break
        current_path = R3.subsegment_from_t_values(seg, [t_cut_point,1.0])
        currently_visible = visible_V[cvp1]
      else:
        current_path.append( V[cvp1] )
      current_vert = cvp1
    #print "Loop paths:"
    #for pp in self.PPs:
    #  print pp
          
    
############################################################################
# a hyperbolic lifted surface displayer
############################################################################
class ShineHypSurfaceDisplay:
  def __init__(self, shine_parent):
    self.shine_parent = shine_parent
    self.tk_parent = shine_parent.parent
    self.window = tk.Toplevel(self.tk_parent)
    self.window.title('Hyperbolic structure')
    self.window.geometry('500x500+%d+%d' % (self.tk_parent.winfo_rootx()+self.tk_parent.winfo_width(), self.tk_parent.winfo_rooty()) )
    
    self.canvas = tk.Canvas(self.window, borderwidth=0)
    self.canvas.bind('<Configure>', self.canvas_resize)
    self.canvas.bind('<Button-1>', self.canvas_click)
    
    self.button_rotate_left = tk.Button(self.window, text="<", command=lambda : self.rotate('left'))
    self.button_rotate_right = tk.Button(self.window, text='>', command=lambda : self.rotate('right'))
    self.draw_do_propagate = tk.IntVar()
    self.check_propagate = tk.Checkbutton(self.window, text='Propagate', variable=self.draw_do_propagate, command=self.canvas_redraw)
    self.draw_do_propagate.set(0)
    
    self.window.rowconfigure(2, weight=1)
    self.window.columnconfigure(0, weight=1)
    
    self.canvas.grid(column=0, row=0, rowspan=3, columnspan=3, sticky=tk.W+tk.E+tk.N+tk.S)
    self.button_rotate_left.grid(column=1, row=0, sticky=tk.W+tk.E)
    self.button_rotate_right.grid(column=2, row=0,sticky=tk.W+tk.E)
    self.check_propagate.grid(column=1, row=1, columnspan=2, sticky=tk.W)
    
    self.draw_transformation = mobius.MobiusTrans(1,0,0,1)
    self.draw_scale = 50.0
    self.drawing_items = []
    self.draw_colors = [rand_bright_color() for _ in xrange(100)]
  
  ###########################################################################
  # resize the canvas and figure out how big things are
  ###########################################################################
  def canvas_resize(self, event):
    self.canvas.config(background='#FFFFFF')
    self.canvas.config(height=event.height, width=event.width)
    self.draw_canvas_width = event.width
    self.draw_canvas_height = event.height
    self.draw_canvas_center = (self.draw_canvas_width/2, self.draw_canvas_height/2)
    self.draw_canvas_middle = self.draw_canvas_width/2
    self.draw_complex_width = self.draw_canvas_width/self.draw_scale
    self.draw_complex_height = self.draw_canvas_height/self.draw_scale
    self.canvas_redraw()
  
  ###########################################################################
  # reset
  ###########################################################################
  def reset(self):
    #erase the canvas
    for di in self.drawing_items:
      self.canvas.delete(di)
    self.drawing_items = []
    self.draw_transformation = mobius.MobiusTrans(1,0,0,1)
    self.canvas_redraw()
  ###########################################################################
  # draw the surface
  ###########################################################################
  def canvas_redraw(self):
    #erase the canvas
    for di in self.drawing_items:
      self.canvas.delete(di)
    self.drawing_items = []
    
    #if there's no surface, show nothing
    if self.shine_parent.LS == None:
      self.drawing_items.append( self.canvas.create_text(*self.draw_canvas_center, text='(No surface; open a new surface with the file menu)') )
      return
    
    #if there is a surface, show it
    for di in self.drawing_items:
      self.canvas.delete(di)
    self.drawing_items = []
    
    self.extended_LS = copy.deepcopy(self.shine_parent.LS)
    
    if self.draw_do_propagate.get()==1:
      self.propagate_surface()
    
    #draw all the first lifts of the triangles
    for ti in xrange(len(self.extended_LS.t)):
      lti = self.extended_LS.t_lifts[ti][0]
      self.draw_triangle(self.extended_LS.em_t[lti].t, '#FFA0A0')
    
    #draw *all* the edges, including the new ones
    #if an edge only has one lift in the original LS, then
    #it should be thin
    for i,e in enumerate(self.extended_LS.em_e):
      if len(self.shine_parent.LS.e_lifts[e.covered_e])==1:
        self.draw_geodesic_segment(e.gi, thickness=1)
      else:
        self.draw_geodesic_segment(e.gi, thickness=2)
    
    #draw all the lifts of edges
    for ell in self.shine_parent.loop_displayer.loops:
      if ell.show.get() == 0:
        continue
      lep = len(ell.EP.edges)
      ES = self.shine_parent.ES
      for i in xrange(lep):
        ip1 = (i+1)%lep
        e1 = ell.EP.edges[i]
        e1t = (ell.EP.edge_coords[i] if e1.sign>0 else 1-ell.EP.edge_coords[i])
        e2 = ell.EP.edges[ip1]
        e2t = (1-ell.EP.edge_coords[ip1] if e2.sign>0 else ell.EP.edge_coords[ip1])
        ti, e1_in_t = (ES.e[e1.ind].on_left if e1.sign>0 else ES.e[e1.ind].on_right)
        e2_in_t = (ES.e[e2.ind].on_right if e2.sign>0 else ES.e[e2.ind].on_left)[1]
        GIs = []
        for lifted_ti in self.extended_LS.t_lifts[ti]:
          lifted_t = self.extended_LS.em_t[lifted_ti]
          GIs.append( lifted_t.t.gi_between_points(e1_in_t, e1t, e2_in_t, e2t) )
        for gi in GIs:
          self.draw_geodesic_segment(gi, thickness=3, color=ell.color)
        
    
    #draw all the vertices, giving the same color to lifts
    #of the same vertex
    for i,v in enumerate(self.extended_LS.em_v):
      col = self.draw_colors[v.covered_v%len(self.draw_colors)]
      self.draw_point(v.pt, col, r=3)
  ###########################################################################
  # convert a complex point to a pixel and vice versa
  ###########################################################################
  def draw_complex_to_canvas(self, z):
    z *= self.draw_scale
    y = self.draw_canvas_height - z.imag
    x = self.draw_canvas_middle + z.real
    return (x,y)
  
  def draw_canvas_to_complex(self, x,y):
    return (1.0/self.draw_scale)*complex(x-self.draw_canvas_middle, self.draw_canvas_height - y)
  
  ###########################################################################
  # draw a point (really, an oval)
  ###########################################################################
  def draw_point(self, p, col, r=2):
    t_p = self.draw_transformation(p)
    pp = self.draw_complex_to_canvas(t_p)
    di = self.canvas.create_oval(pp[0]-r,pp[1]-r,pp[0]+r,pp[1]+r,fill=col)
    self.drawing_items.append(di)
  
  ###########################################################################
  # check if a complex point lies in the canvas
  ###########################################################################
  def check_pt_in_window(self, p):
    tp = self.draw_transformation(p)
    br = self.draw_complex_width/2.0
    bh = self.draw_complex_height
    return -br < tp.real and tp.real < br and 0 < tp.imag and tp.imag < bh
  
  ##########################################################################
  # check if a geodesic segment is disjoint from the drawing, or very small
  ##########################################################################
  def check_gi_in_window(self, gi, do_trans=True):
    t_gi = gi.act_by_mobius(self.draw_trans)
    if t_gi.Euclidean_length() < 1:
      return True
    return arc_disjoint_from_box(t_gi.circ_center, t_gi.circ_radius,        \
                                 t_gi.circ_angle1, t_gi.circ_angle2,       \
                                 self.draw_complex_width/2.0,               \
                                 self.draw_complex_height)   
  ##########################################################################
  # draw a geodesic segment (arc of a circle)
  ##########################################################################
  def draw_geodesic_segment(self, gi, thickness=1, color='#000000'):
    trans_gi = gi.act_by_mobius(self.draw_transformation)
    #print "Drawing geodesic segment: ", gi
    #print "After trans: ", trans_gi
    if trans_gi.vertical:
      s = self.draw_complex_to_canvas(trans_gi.start)
      e = self.draw_complex_to_canvas(trans_gi.end)
      di = self.canvas.create_line(s[0], s[1], e[0], e[1], width=2, fill=color)
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
    di = self.canvas.create_arc(bbox[0], bbox[1], bbox[2], bbox[3], style=tk.ARC, start=ma, extent=(Ma-ma), width=thickness, outline=color)
    #print "Drew arc at ", bbox, " with start, extent", ma, Ma-ma
    self.drawing_items.append(di)
  
  #########################################################################
  # draw a triangle (and fill it, using a polygon)
  #########################################################################
  def draw_triangle(self, t, col):
      #make a polygon with several points per side
      pts = []
      for i in xrange(3):
        gi = t.sides[i].act_by_mobius(self.draw_transformation)
        if gi.vertical:
          num_joints = 1
        else:
          total_angle = abs(gi.circ_angle2 - gi.circ_angle1)
          num_joints = int((total_angle/(2*math.pi))*20 + 1)
          EL = gi.Euclidean_length()
          num_joints = int(num_joints*EL+1)
        for j in xrange(num_joints):
          pts.append( gi.pt_along_euclidean( float(j)/float(num_joints) ) )
      cpts = [self.draw_complex_to_canvas(p) for p in pts]
      cpts = [x for P in cpts for x in P]
      di = self.canvas.create_polygon(*cpts, fill=col)
      self.drawing_items.append(di)
  
  ########################################################################
  # propogate the surface (i.e. lift more triangles)
  ########################################################################
  def propagate_surface(self):
    v_stack = [(i,v) for i,v in enumerate(self.extended_LS.em_v) if None in v.i_tris]
    while len(v_stack)>0:
      #print "Stack: ", v_stack
      (i,v) = v_stack.pop()
      if (not self.check_pt_in_window(v.pt))           or         \
         (self.draw_transformation(v.pt).imag < 0.25)  or         \
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
  
  ####################### signals
  
  def canvas_click(self, event):
    if self.shine_parent.LS == None:
      return
    p = (event.x, event.y)
    #print "Clicked on ", p
    can_p = (self.canvas.canvasx(p[0]), self.canvas.canvasy(p[1]))
    #print "Canvas: ", can_p
    com_p = self.draw_canvas_to_complex(*can_p)
    can_center = (self.draw_canvas_width/2.0, self.draw_canvas_height/2.0)
    com_center = self.draw_canvas_to_complex(*can_center)
    #print can_center, com_center
    #print "Moving point", com_p, "to", com_center
    M = mobius.MobiusTrans.unit_tangent_action(com_p, 0, com_center, 0)
    self.draw_transformation = M.compose(self.draw_transformation)
    self.canvas_redraw()
  
  def rotate(self, dir):
    if self.shine_parent.LS == None:
      return
    theta = (math.pi/40.0 if dir=='left' else -math.pi/40.0)
    M = mobius.MobiusTrans.unit_tangent_action(1j, 0, 1j, theta)
    self.draw_transformation = self.draw_transformation.compose(M)
    self.canvas_redraw()











#############################################################################
#############################################################################
def visualize_em_surface(ES):
  root = tk.Tk()
  vs = EmSurfaceVisualizer(root, ES)
  root.mainloop()



def visualize_hyp_surface(LS):
  root = tk.Tk()
  vs = HypSurfaceVisualizer(root, LS)
  root.mainloop()

def run_shine():
  root = tk.Tk()
  s = Shine(root)
  root.mainloop()


def test():
  T = tsurf.TopSurface(method='polygon', w='abABcdCD')
  G = gsurf.GeometricSurface.geometrize_tsurf(T)
  L = liftedsurf.LiftedSurface.lift_gsurf(G)
  visualize_surface(L)
  
if __name__ == '__main__':
  run_shine()
  



