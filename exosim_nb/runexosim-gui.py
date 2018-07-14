import numpy           as     np
import pylab           as     pl
import sys, time, os, glob
import exosim


import wx
import wx.aui
import matplotlib.pyplot as pl
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as Canvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2Wx as Toolbar


class MyFrame(wx.Frame):
    def __init__(self, parent, title):
        wx.Frame.__init__(self, parent, wx.ID_ANY, title=title, size=(600,600))
        #self.control = wx.TextCtrl(self, style=wx.TE_MULTILINE)
        self.CreateStatusBar()
        self.Show(True)
        filemenu = wx.Menu()
        menuItem = filemenu.Append(wx.ID_ABOUT, "&About","ExoSim GUI")
        self.Bind(wx.EVT_MENU, self.OnAbout, menuItem)
        menuItem = filemenu.AppendSeparator()
        menuItem = filemenu.Append(wx.ID_EXIT, "&Quit","Quit ExoSim GUI")
        self.Bind(wx.EVT_MENU, self.OnExit, menuItem)
        
        menuBar = wx.MenuBar()
        menuBar.Append(filemenu, "&File")
        self.SetMenuBar(menuBar)
    
    def OnAbout(self, event):
        dlg = wx.MessageDialog(self, "ExoSim GUI\n (c) 2015 The ExoSim Authors", "About this")
        dlg.ShowModal()
        dlg.Destroy()
    
    def OnExit(self, event):
        #self.Close(True)
        self.Destroy()
        
class Plot(wx.Panel):
    def __init__(self, parent, id = -1, dpi = None, ncols=1, nrows=1, sharex=False, **kwargs):
        wx.Panel.__init__(self, parent, id=id, **kwargs)
        
        self.figure, self.axarr = pl.subplots(ncols=ncols, nrows=nrows, dpi=dpi, figsize=(2,2), sharex=sharex)
        self.canvas = Canvas(self, -1, self.figure)
        self.toolbar = Toolbar(self.canvas)
        self.toolbar.Realize()

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas,1,wx.EXPAND)
        sizer.Add(self.toolbar, 0 , wx.LEFT | wx.EXPAND)
        self.SetSizer(sizer)

class PlotNotebook(wx.Panel):
    def __init__(self, parent, id = -1):
        wx.Panel.__init__(self, parent, id=id)
        self.nb = wx.aui.AuiNotebook(self)
        sizer = wx.BoxSizer()
        sizer.Add(self.nb, 1, wx.EXPAND)
        self.SetSizer(sizer)

    def add(self,name="plot", nrows=1, ncols=1, sharex=False):
       page = Plot(self.nb, nrows=nrows, ncols=ncols, sharex=sharex)
       self.nb.AddPage(page,name)
       return page.axarr, page.figure

if __name__ == "__main__":

  # Run simulation
  exosim.exolib.exosim_msg('Reading options from file ... \n')
  opt = exosim.Options(filename='exosim_twinkle.xml').opt #, default_path = exosim.__path__[0]).opt
  star, planet, zodi, channel = exosim.run_exosim(opt)
  ####
  
  # Initialize GUI
  app = wx.App(False)
  #frame = wx.Frame(None,-1,'Plotter')
  frame = MyFrame(None, 'Plotter')
  plotter = PlotNotebook(frame)
  

  ax1, fig1 = plotter.add('Common', nrows=1, ncols=2)
  ax1[0].plot(star.sed.wl, star.sed.sed)
  ax1[0].set_xlabel(star.sed.wl.units.dimensionality.latex)
  ax1[0].set_ylabel(star.sed.sed.units.dimensionality.latex)
  ax1[1].plot(planet.cr.wl, planet.cr.sed)
  ax1[1].set_xlabel(planet.cr.wl.units.dimensionality.latex)
  ax1[1].set_ylabel('contrast ratio')
  
  
  channel_ax = {}
  for key in channel.keys():
    
    channel_ax[key], fig = plotter.add(key, nrows=2, ncols=2)
    
    # ZODI
    #channel_ax[key][0].plot(channel[key].zodi.wl, channel[key].zodi.sed, 'c')
    
    # STAR
    i0 = np.argmax(channel[key].fp.sum(axis=1))
    channel_ax[key][0,0].plot(channel[key].star.wl, channel[key].star.sed, 'r', 
			      label = 'Star')
    channel_ax[key][0,0].plot(channel[key].wl_solution[1::channel[key].opt.osf()],
			      channel[key].fp[i0,channel[key].offs::channel[key].osf], 
			      'g', label = 'Star on FP')
    
    
    # EMISSION
    channel_ax[key][0,0].plot(channel[key].emission.wl, channel[key].emission.sed,
			      'k', label = 'emission')
    channel_ax[key][0,0].plot(channel[key].zodi.wl, channel[key].zodi.sed, 
			      'b', label = 'zodi')
    channel_ax[key][0,0].legend()
    channel_ax[key][0,0].set_title('signals')
    channel_ax[key][0,0].set_ylabel(r'e$^-$ s$^{-1}$/pix')

    # PLANET CR
    channel_ax[key][0,1].plot(planet.cr.wl, planet.cr.sed, 'c')
    channel_ax[key][0,1].plot(channel[key].planet.wl, channel[key].planet.sed, 'r')
    channel_ax[key][0,1].set_title('planet')
    

    # TRANSMISSION
    channel_ax[key][1,1].plot(channel[key].transmission.wl, channel[key].transmission.sed, 'k')
    channel_ax[key][1,1].set_title('instrument transmission')

    channel_ax[key][1,0].imshow(np.array(channel[key].fp), 
	  extent=[channel[key].wl_solution.min(),
	  channel[key].wl_solution.max(),0,1])
    channel_ax[key][1,0].set_title('focal plane')
    
    for ax_r in channel_ax[key]:
      for ax in ax_r:
	ax.set_xbound(channel[key].wl_solution.min(),channel[key].wl_solution.max())
	ax.set_xlabel(channel[key].wl_solution.units.dimensionality.latex)
	
  timeline_ax = {}
  
  for key in channel.keys():
    timeline_ax[key], fig = plotter.add("Timeline - "+key, nrows=3, ncols=1, sharex=False)
    #fig.subplots_adjust(hspace=0)
    wl_solution = channel[key].wl_solution[channel[key].offs::channel[key].osf]
    tl = channel[key].tl 
    noise = channel[key].noise 
    cds_so = tl[..., opt.timeline.multiaccum()-1::opt.timeline.multiaccum()] - \
	  tl[..., 0::opt.timeline.multiaccum()]
    cds = cds_so+noise[..., opt.timeline.multiaccum()-1::opt.timeline.multiaccum()] - \
	  noise[..., 0::opt.timeline.multiaccum()]
    wl_check_points = np.percentile(wl_solution, [20, 50, 80])*wl_solution.units
    for iwl, wl in enumerate(wl_check_points):
      idx = np.argmin( np.abs(wl - wl_solution) )
      timeline_ax[key][iwl].plot(cds.sum(axis=0)[idx, ...])
      timeline_ax[key][iwl].plot(cds_so.sum(axis=0)[idx, ...])
      timeline_ax[key][iwl].set_ylabel(r'counts [e$^-$@{:0.1f}]'.format(wl.item()))
    
    
  frame.Show()
  app.MainLoop()
  #app.Destroy()
  