from pyspecdata import *
from align_slice import align_and_slice
fl = figlist_var()
apply_correction = True
for filename,nodename,thisdir,vdrange in [
        #('200212','calibrate_clock_6'),
        #('200917_calibrate_clock_1','signal','test_equip',slice(None,-1)),
        ('201105_calibrate_clock_1','signal','test_equip',slice(None,-1)),
        #('200922_calibrate_clock_2','signal','test_equip',slice(None)),
        ]:
    fl.basename = filename
    s = find_file(filename, expno=nodename, exp_type=thisdir)
    # I don't understand how to get the carrier frequency
    # carrier = s.get_prop('acq')['carrierFreq_MHz']
    # I just set this to zero, so I can include the code below
    carrier = 0
    s.rename('t','t2').set_units('t2','s')
    s.reorder('vd') # based on field variation, it does seem the repeats were
    # done in the outer loop, but I still like to present them as though from an
    # inner loop
    try:
        centerpoint = abs(s).mean('nScans').mean('vd').argmax('t2').item()
    except:
        centerpoint = abs(s).mean('vd').argmax('t2').item()
    s.setaxis('t2', lambda x: x-centerpoint)
    try:
        s.mean('nScans')
    except:
        print("Did not find nScans axis")
    s.setaxis('vd',r_[0:len(s.getaxis('vd'))])
    fl.next('image raw -- time domain')
    fl.image(s, interpolation='bicubic')
    s.ft('t2', shift=True)
    fl.next('image raw -- f domain')
    fl.image(s, interpolation='bicubic')
    s = s['t2':(-0.03e3,0.03e3)]
    fl.next('image raw -- f domain, zoom')
    fl.image(s, interpolation='bicubic')
    if apply_correction:
        #clock_correction = 1.692
        #clock_correction = -1.089
        clock_correction = -1.21
        s *= exp(-1j*s.fromaxis('vd')*clock_correction)
    fl.next('image raw -- frequency domain')
    fl.image(s)
    s = align_and_slice(s, convwidth=0, fl=fl)
    fl.next('image shifted and sliced')
    fl.image(s)
    fl.next('mean and repeat align and slice')
    s = align_and_slice(s, convwidth=0)
    fl.image(s)
    fl.next('phase error vs vd')
    s_forcalc = s.C.sum('t2')
    fl.plot(s_forcalc.angle, 'o')
    fl.next('phase error, unwrapped vs vd - sig avg')
    s_forcalc = s_forcalc['vd',1:]/s_forcalc['vd',:-1]
    s_angle = s_forcalc.angle.name('signal phase').set_units('rad').run_nopop(cumsum,'vd')
    fl.plot(s_angle,'o')
    # begin fit to return clock correction
    c,thisfit = s_angle['vd',vdrange].polyfit('vd')
    fl.plot(thisfit, alpha=0.5)
    zeroth_order = c[0]
    clock_correction = c[1]
    # if I have the carrier:
    text(0.5,0.5, '%#0.4g rad/s, %#0.4g Hz, %#0.8g Hz carrier = %#0.8g $\\times$ 75 MHz'%(clock_correction,clock_correction/2/pi,
        carrier*1e6,carrier/75),
            transform=gca().transAxes)
    s *= exp(-1j*s.fromaxis('vd')*clock_correction)
    s *= exp(-1j*zeroth_order)
    fl.next('corrected')
    fl.image(s)
fl.show()
