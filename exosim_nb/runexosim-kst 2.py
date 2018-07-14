import numpy           as     np
import pykstplot           as     plt
import pykst as kst
import sys, time, os, glob
import exosim


if __name__ == "__main__":

  # Run simulation
  exosim.exolib.exosim_msg('Reading options from file ... \n')
  opt = exosim.Options(filename='exosim_ariel.xml').opt #, default_path = exosim.__path__[0]).opt
  star, planet, zodi, channel = exosim.run_exosim(opt)
  ####

  client = kst.Client('ExoSim')
  tab = {}
  client.set_tab_text('Common')
  tab['Common'] = {'p1': {}}
  tab['Common']['p1']['V1']  = client.new_editable_vector(star.sed.wl, name=r'\lambda')
  tab['Common']['p1']['V2']  = client.new_editable_vector(star.sed.sed, name='Star SED')
  tab['Common']['p1']['c1']  = client.new_curve(tab['Common']['p1']['V1'], tab['Common']['p1']['V2'])
  tab['Common']['p1']['plt'] = client.new_plot()
  tab['Common']['p1']['plt'].add(tab['Common']['p1']['c1'])
  tab['Common']['p1']['plt'].set_bottom_label('\lambda \[{\mu}m\]')
  
  tab['Common'] = {'p2': {}}
  tab['Common']['p2']['V1']  = client.new_editable_vector(planet.cr.wl, name=r'\lambda')
  tab['Common']['p2']['V2']  = client.new_editable_vector(planet.cr.sed, name='Planet CR')
  tab['Common']['p2']['c1']  = client.new_curve(tab['Common']['p2']['V1'], tab['Common']['p2']['V2'])
  tab['Common']['p2']['plt'] = client.new_plot()
  tab['Common']['p2']['plt'].add(tab['Common']['p2']['c1'])
  tab['Common']['p2']['plt'].set_bottom_label('\lambda \[{\mu}m\]')
 
  for key in channel.keys():
    tab = client.new_tab()
    client.set_tab_text(key)
    # signals
    #  STAR
    p = client.new_plot()
    i0 = np.argmax(channel[key].fp.sum(axis=1))
    vx = client.new_editable_vector(channel[key].star.wl, name=r'\lambda')
    vy = client.new_editable_vector(channel[key].star.sed, name='Star')
    c  = client.new_curve(vx, vy)
    p.add(c)
    
    vx = client.new_editable_vector(channel[key].wl_solution[1::channel[key].opt.osf()], name=r'\lambda')
    vy = client.new_editable_vector(channel[key].fp[i0, channel[key].offs::channel[key].osf], name='Star on FP')
    c  = client.new_curve(vx, vy)
    p.add(c)
    #  EMISSION
    vx = client.new_editable_vector(channel[key].emission.wl, name=r'\lambda')
    vy = client.new_editable_vector(channel[key].emission.wl, name='Emission')
    c  = client.new_curve(vx, vy)
    p.add(c)
    #  ZODI
    vx = client.new_editable_vector(channel[key].zodi.wl, name=r'\lambda')
    vy = client.new_editable_vector(channel[key].zodi.sed, name='Emission')
    c  = client.new_curve(vx, vy)
    p.add(c)

    #  PLANET CR
    p = client.new_plot()
    vx = client.new_editable_vector(planet.cr.wl, name=r'\lambda')
    vy = client.new_editable_vector(planet.cr.sed, name='Planet/cr')
    c  = client.new_curve(vx, vy)
    p.add(c)
    vx = client.new_editable_vector(channel[key].planet.wl, name=r'\lambda')
    vy = client.new_editable_vector(channel[key].planet.sed, name='Planet/cr \[FP\]')
    c  = client.new_curve(vx, vy)
    p.add(c)
    
    #  PLANET CR
    p = client.new_plot()
    vx = client.new_editable_vector(channel[key].transmission.wl, name=r'\lambda')
    vy = client.new_editable_vector(channel[key].transmission.sed, name='Transmission')
    c  = client.new_curve(vx, vy)
    p.add(c)
  
  for key in channel.keys():
    wl_solution = channel[key].wl_solution[channel[key].offs::channel[key].osf]
    tl = channel[key].tl 
    noise = channel[key].noise 
    cds_so = tl[..., opt.timeline.multiaccum()-1::opt.timeline.multiaccum()] - \
	  tl[..., 0::opt.timeline.multiaccum()]
    cds = cds_so+noise[..., opt.timeline.multiaccum()-1::opt.timeline.multiaccum()] - \
	  noise[..., 0::opt.timeline.multiaccum()]
    wl_check_points = np.percentile(wl_solution, [20, 40, 60, 80])*wl_solution.units

    tab = client.new_tab()
    client.set_tab_text('Timeline - '+key)
    for iwl, wl in enumerate(wl_check_points):
      idx = np.argmin( np.abs(wl - wl_solution) )
      p = client.new_plot()
      time_grid = np.arange(tl.shape[0]).astype(np.float)
      vx = client.new_editable_vector(time_grid, name='Time')
      
      vy = client.new_editable_vector(cds.sum(axis=0)[idx, ...], name='Light curve')
      c  = client.new_curve(vx, vy)
      p.add(c)  
      vy = client.new_editable_vector(cds_so.sum(axis=0)[idx, ...], name='Light curve')
      c  = client.new_curve(vx, vy)
      p.add(c)  
    