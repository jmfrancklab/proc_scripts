#!/usr/bin/env python
# coding: utf-8

# In[5]:


from pyspecdata import *
from utility import dBm2power
fl = figlist_var()
for date,id_string in [
        ('200115','CPMG_DNP_1')
        ]:
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(
                exp_type = 'test_equip'))
    #{{{ pulling acq params
    s.setaxis('power',r_[
        0,dBm2power(array(s.get_prop('meter_powers'))+20)]
        ).set_units('power','W')
    SW_kHz = s.get_prop('acq_params')['SW_kHz']
    nPoints = s.get_prop('acq_params')['nPoints']
    nEchoes = s.get_prop('acq_params')['nEchoes']
    nPhaseSteps = s.get_prop('acq_params')['nPhaseSteps']
    nScans = s.get_prop('acq_params')['nScans']
    p90_s = s.get_prop('acq_params')['p90_us']*1e-6
    deadtime_s = s.get_prop('acq_params')['deadtime_us']*1e-6
    deblank_s = s.get_prop('acq_params')['deblank_us']*1e-6
    marker_s = s.get_prop('acq_params')['marker_us']*1e-6
    tau1_s = s.get_prop('acq_params')['tau1_us']*1e-6
    pad_start_s = s.get_prop('acq_params')['pad_start_us']*1e-6
    pad_end_s = s.get_prop('acq_params')['pad_end_us']*1e-6
    #}}}
    orig_t = s.getaxis('t')
    acq_time_s = orig_t[nPoints]
    s.set_units('t','s')
    twice_tau = deblank_s + 2*p90_s + deadtime_s + pad_start_s + acq_time_s + pad_end_s + marker_s
    t2_axis = linspace(0,acq_time_s,nPoints)
    tE_axis = r_[1:nEchoes+1]*twice_tau
    s.setaxis('t',None)
    #s.setaxis('nScans',r_[0:nScans])
    s.chunk('t',['ph1','tE','t2'],[nPhaseSteps,nEchoes,-1])
    s.setaxis('ph1',r_[0.,2.]/4)
    s.setaxis('tE',tE_axis)
    s.setaxis('t2',t2_axis)
    fl.next(id_string+'raw data - chunking')
    fl.image(s)
    s.ft('t2', shift=True)
    fl.next(id_string+'raw data - chunking ft')
    fl.image(s)
    s.ft(['ph1'])
    s = s['ph1',1].C
    fl.next(id_string+' image plot coherence-- ft ')
    fl.image(s)
    s.ift('t2')
    fl.next(id_string+' image plot coherence ')
    fl.image(s)
    #s.mean('nScans',return_error=False)
    s.reorder('t2',first=True)
    t2_max = zeros_like(s.getaxis('power'))
    for x in range(len(s.getaxis('power'))):
        t2_max[x] = abs(s['power',x]['tE',0]).argmax('t2',raw_index=True).data
    s.setaxis('t2',lambda t: t -s.getaxis('t2')[int(t2_max.mean())])
    s.rename('tE','nEchoes').setaxis('nEchoes',r_[1:nEchoes+1])
    s.reorder('nEchoes',first=True)
    s.ft('t2')
    # as of right now, this plot doesn't show anything meaningful
    fl.next('check center')
    fl.image(s)
    f_axis = s.fromaxis('t2')
    def costfun(p):
        zeroorder_rad,firstorder = p
        phshift = exp(-1j*2*pi*f_axis*(firstorder*1e-6))
        phshift *= exp(-1j*2*pi*zeroorder_rad)
        corr_test = phshift * s
        return (abs(corr_test.data.imag)**2)[:].sum()
    iteration = 0
    def print_fun(x, f, accepted):
        global iteration
        iteration += 1
        #print (iteration, x, f, int(accepted))
        return
    sol = basinhopping(costfun, r_[0.,0.],
            minimizer_kwargs={"method":'L-BFGS-B'},
            callback=print_fun,
            stepsize=100.,
            niter=100,
            T=1000.
            )
    zeroorder_rad, firstorder = sol.x
    phshift = exp(-1j*2*pi*f_axis*(firstorder*1e-6))
    phshift *= exp(-1j*2*pi*zeroorder_rad)
    s *= phshift
    #print("RELATIVE PHASE SHIFT WAS {:0.1f}\us and {:0.1f}$^\circ$".format(
    #        firstorder,angle(zeroorder_rad)/pi*180))
    #if s['nEchoes':0]['t2':0]['power',-4].item().real > 0:
        #print s['nEchoes':0]['t2':0]['power',-4].item().real
        #print "Sign correction"
        #s *= -1
        #print s['nEchoes':0]['t2':0]['power',-4].item().real
    #print ndshape(s)
    s.reorder('power',first=True)
    fl.next('after phased - real ft')
    fl.image(s.real)
    fl.next('after phased - imag ft')
    fl.image(s.imag)


# In[6]:


d = s['power',:-3]['t2':(-100,150)].C.sum('t2')
data = d['power',-1].C
data = data.data
tau = tE_axis
T2_len = 50
T2 = linspace(log10(3e-3),log10(10.0),T2_len)
tau_2d = reshape(tau,(shape(tau)[0],1))
T2_2d = reshape(10**T2,(1,shape(T2)[0]))
kernel = exp(-tau_2d/T2_2d)
x,rnorm = nnls_regularized(kernel,data,l=0.)
fl.next('plot data and fit')
plot(data,'.')
plot(kernel.dot(x))


# In[7]:


def L_curve(l,r_norm,x_norm,**kwargs):
    """plot L-curve using
    
    Parameters
    ==========
    l: double
        lambda values
    r_norm: double
        norm of the residual
    x_norm: double
        norm of solution vector"""
    plot(log10(r_norm),log10(x_norm),'o',**kwargs)
    for j,this_l in enumerate(l):
        annotate('%5g'%this_l, (log10(r_norm[j]),log10(x_norm[j])),
                ha='left',va='bottom',rotation=45)
    ylabel('xnorm')
    xlabel('residual')


# In[ ]:


l = sqrt(logspace(-10,1,25))

def nonvec_lcurve(A,l):
    x_norm = empty_like(l)
    r_norm = empty_like(l)
    for j,this_l in enumerate(l):
        x,r_norm[j] = nnls_regularized(A,data,l=this_l)
        x_norm[j] = linalg.norm(x)
    return x,x_norm,r_norm
x,x_norm,r_norm = nonvec_lcurve(kernel,l)

fl.next('L-curve')
L_curve(l,r_norm,x_norm,markersize=10,alpha=0.5,label='manual loop')


# In[ ]:


x,rnorm = nnls_regularized(kernel,data,l=0.13)


# In[ ]:


this_data = ndshape([shape(d.getaxis('power'))[0],T2_len],['power','T2']).alloc()
this_data.setaxis('power',d.getaxis('power'))
this_data.setaxis('T2',T2)
for y in range(shape(d.getaxis('power'))[0]):
    temp = d['power',y].C
    if temp.data[0].real < 0:
        temp *= -1
    x,rnorm = nnls_regularized(kernel,temp.data,l=0.13)
    this_data['power',y]['T2',:] = x[:]

print(this_data.getaxis('power'))
this_data.rename('T2','log(T2)')

fl.next('T2 distribution vs power')
fl.image(abs(this_data))
clim(0,max(this_data['power',0].data).real)

fl.show();quit()
