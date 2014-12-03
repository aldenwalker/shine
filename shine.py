import tsurf
import gsurf
import liftedsurf
import mobius

import math
import Tkinter as tk

class SurfaceVisualizer:
  def __init__(self, master, LS):
    self.master = master
    self.master.geometry('500x500+100+100')
    
    self.canvas = tk.Canvas(self.master, borderwidth=0)
    self.canvas.bind('<Configure>', self.canvas_resize)
    self.canvas.bind('<Button-1>', self.canvas_click)
    
    self.button1 = tk.Button(self.master, text = 'Quit', command=self.quit)
    
    self.master.rowconfigure(0, weight=1)
    self.master.columnconfigure(0, weight=1)
    
    self.button1.grid(column=1, row=0, sticky=tk.N)
    self.canvas.grid(column=0, row=0, sticky=tk.W+tk.E+tk.N+tk.S)
    
    #remember the surface 
    self.LS = LS
    
    #set up the drawing
    #a complex number is scaled by self.draw_scale and then the y coordinate 
    #is taken from the bottom of the screen
    #and the x coord is taken from the middle
    self.draw_scale = 50.0
    self.draw_width = None    #will be set by Configure
    self.draw_height = None   #will be set by Configure
    self.draw_middle = None   #will be set by Configure
    #everything is acted upon by self.draw_trans before drawing
    self.draw_trans = mobius.MobiusTrans(1,0,0,1)
    #the drawing the currently empty (but it will be filled by the canvas Configure)
    self.drawing_items = []
    self.draw_colors = ['red','green','blue','cyan','yellow','magenta']
  
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
    self.redraw_surface()
  
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
    self.redraw_surface()
  
  def redraw_surface(self):
    for di in self.drawing_items:
      self.canvas.delete(di)
    self.drawing_items = []
    for gi_left, gi_right in self.LS.em_e:
      self.draw_geodesic_segment(gi_left)
      if not gi_left.same_but_reversed(gi_right):
        self.draw_geodesic_segment(gi_right)
    for vi,V in enumerate(self.LS.em_v):
      c = self.draw_colors[vi%len(self.draw_colors)]
      for p in V:
        self.draw_point(p,c)
        
  def draw_point(self, p, col):
    t_p = self.draw_trans(p)
    pp = self.draw_complex_to_canvas(t_p)
    di = self.canvas.create_oval(pp[0]-2,pp[1]-2,pp[0]+2,pp[1]+2,fill=col)
    self.drawing_items.append(di)
  
  def draw_geodesic_segment(self, gi):
    trans_gi = gi.act_by_mobius(self.draw_trans)
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
    di = self.canvas.create_arc(bbox[0], bbox[1], bbox[2], bbox[3], style=tk.ARC, start=ma, extent=(Ma-ma), width=2)
    #print "Drew arc at ", bbox, " with start, extent", ma, Ma-ma
    self.drawing_items.append(di)
  
  def quit(self):
    self.master.destroy()



def visualize_surface(LS):
  root = tk.Tk()
  vs = SurfaceVisualizer(root, LS)
  root.mainloop()

