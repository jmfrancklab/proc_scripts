#{Imports & Initialization
from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import postproc_dict
from sympy import symbols, Symbol, latex,limit,init_printing
from numpy import *
import time
import matplotlib.pyplot as plt
plt.rcParams.update({
    "figure.facecolor":  (1.0, 1.0, 1.0, 0.0),  # clear
    "axes.facecolor":    (1.0, 1.0, 1.0, 0.0),  # clear ## 1,1,1,0.9 == 90% transparent white
    "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
})
logger = init_logging("info")
t2 = symbols('t2')
fl = fl_mod()
measured_vs_actual = 22. # how many dB down the split + measured power is from
                         # the forward power
#}
#{Input parameters & Load Data
date = time.strftime('%y%m%d') 
file_location = 'odnp'
postproc = 'spincore_ODNP_v2' # 'spincore_ODNP_v1'
plotname = 'Q183R1a pR Hoffmeister Series Anions'
names = ['KCl','KI','$KH_{2}PO_{4}$']
plot_all = True 
save_figs = True
curves = nddata(zeros((len(names),18),dtype='complex128'),['sample','power'])
enhancements = nddata(zeros((len(names),18),dtype='complex128'),['sample','power'])
    # { Load T1_values
#import pandas as pd
#path = 'C:/Users/saman/Research/Sam_Notebooks/ODNP_proc'
#T1_df = pd.read_csv('%s/210707_Q183R1a_pR_DDM_T1s.csv'%path)
#print(T1_df);quit()
    # }
for (idx,filename,nodename,f_range,C,T1_0,T1_vals,ppt) in [
#        (0,'210714_A174R1a_pR_DDM_ODNP','enhancement1',
#            (-250,100),240e-6,1.49,None,1.5154),
#        (1,'210715_A174R1a_pR_KI_ODNP','enhancement',
#            (-200,200),236.9e-6,1.45,None,1.5154),
#        (2,'210715_A174R1a_pR_KH2PO4_ODNP','enhancement',
#            (-200,200),162.1e-6,2.05,None,1.5154),
#        (3,'210714_A174R1a_pR_DHPC_ODNP','enhancement',
#            (-250,75),238.5e-6,1.47,None,1.5154),
    
        (0,'210707_Q183R1a_pR_DDM_ODNP','enhancement',
            (-225,75),207.4e-6,1.00,None,1.5154), # have T1_values for this data
        (1,'210708_Q183R1a_pR_KI_ODNP','enhancement',
            (-200,100),113.1e-6,1.12,None,1.5154),
        (2,'210708_Q183R1a_pR_KH2PO4_ODNP','enhancement',
            (-225,75),115.2e-6,1.9,None,1.5154),
#        (3,'210707_Q183R1a_pR_DHPC_ODNP','enhancement',
#            (-225,75),103.2e-6,2.15,None,1.5154)
        ]:
    titl = '%s %s'%(' '.join(filename.split('_')[1:-1]),filename.split('_')[0])
    s = find_file(filename,exp_type=file_location,expno=nodename,
            postproc=postproc,lookup=postproc_dict,fl=fl)
#}
    fl.next('raw data')
    s = s['t2':(-1e3,1e3)]
    fl.image(s)
    ph0 = s['power',-4].sum('t2')
    ph0 /= abs(ph0)
    s /= ph0
    fl.next('phased')
    fl.image(s)
    s = s['t2':f_range]
    fl.next('sliced to integration bounds')
    fl.image(s)
    if 'ph2' in s.dimlabels:
        s = s['ph1',1]['ph2',0]
    else:
        s = s['ph1',1]
    fl.next('sliced to integration bounds\nselect coh. pathway\nreal')
    fl.image(s.real)
    s.integrate('t2')
    s /= max(s.data.real)
#{getting power axis
    s.setaxis('power',r_[-9999,array(s.get_prop('meter_powers'))])
    print("here are the dBm",s.getaxis('power'))
    s.setaxis('power', lambda x:
            1e-3*10**((x+measured_vs_actual)/10.))
    print("here are the powers",s.getaxis('power'))
    s.set_units('power','W')
#}
#{Integral
    fl.next('simple integral')
    fl.plot(s['power',:-3],'ko', human_units=False) # human_units = False because the range for the red points tries to force mW, which is incompatible with W
    fl.plot(s['power',-3:],'ro', human_units=False)
    plt.ylabel('integrated $^1$H NMR signal')
#}
#{Enhancement
    s /= s['power',0]
    fl.next('enhancement curve')
    fl.plot(s['power',:-3],'ko', human_units=False)
    fl.plot(s['power',-3:],'ro', human_units=False)
    plt.xlabel('power (W)')
    plt.ylabel('enhancement')
    if plot_all:
        enhancements['sample',idx] = s
    if C is not None:
        epsilon = 1 - s.C
        ks = epsilon/C
        if T1_0 is None:
            fl.next('%s\n$k_{\sigma}s(p)T_{1}(p)$'%titl)
        else:
            if T1_vals is None: # if T1_0 is not None and T1_vals is None
                ks /= T1_0
                fl.next('%s\n$k_{\sigma}s(p)T_{1}(p)$ $\div$ $T_{1}(0)$'%titl)
            else:
                m,b = np.polyfit(T1_vals.getaxis('power'),T1_vals.data,1)
                T1_p = lambda p: (m*p)+b
                T1_p = nddata(T1_p(ks.getaxis('power')),ks.getaxis('power')).labels('power',ks.getaxis('power'))
                ks /= T1_p
                fl.next('%s\n$k_{\sigma}s(p)$'%titl)
        if ppt is not None:
            ks *= ppt*1e-3 # gets you in units of k_sigma!!!
        fl.plot(ks['power',:-3],'ko', human_units=False)
        fl.plot(ks['power',-3:],'ro', human_units=False)
        plt.xlabel('power (W)')
        if T1_0 is None:
            plt.ylabel('$k_{\sigma}s(p)T_{1}(p)$')
        else:
            if T1_vals is None:
                plt.ylabel('$k_{\sigma}s(p)T_{1}(p)$ $\div$ $T_{1}(0)$')
            else:
                plt.ylabel('$k_{\sigma}s(p)$')
        if plot_all:
            curves['sample',idx] = ks
#}
###############################################################################
#    if C is not None:
#        epsilon = 1 - s
#        ks = epsilon/C
#        if T1_0 is not None:
#            ks /= T1_0
#            fl.next('$k_{\sigma}s(p)T_{1}(p)$ $\div$ $ T_{1}(0)$')
#        elif T1_vals is not None:
#            m,b = np.polyfit(T1_vals.getaxis('power'),T1_vals.data,1)
#            T1_p = lambda p: (m*p)+b
#            T1_p = nddata(T1_p(ks.getaxis('power')),ks.getaxis('power')).labels('power',ks.getaxis('power'))
#            ks /= T1_p
#            fl.next('$k_{\sigma}s(p)$')
#        if ppt is not None:
#            ks *= ppt*1e-3 # gets you in units of k_sigma!!!
#            fl.plot(ks['power',:-3],'ko', human_units=False)
#            fl.plot(ks['power',-3:],'ro', human_units=False)
#            if plot_all:
#                curves['sample',idx] = ks

fl.show()
#{ Plotting all results together
if save_figs:
    import os
    os.chdir('C:/Users/saman/Research/Sam_Notebooks/ODNP_proc/')
if plot_all:
    figure(figsize=(7,5))
    title('%s %s\nenhancements'%(filename.split('_')[0],plotname))
    enhancements.labels('power',s.getaxis('power'))
    colors = ['dark magenta','brick red','dark orange','gold','greenish','dark cyan','grape','wine']   
    for i,name in enumerate(names):
        plot(enhancements['sample',i]['power',:-3],color='xkcd:%s'%colors[i],marker='o',ls='',human_units=False,
                label='%s'%name,alpha=0.5)
        plot(enhancements['sample',i]['power',-3:],color='xkcd:%s'%colors[i],marker='x',ls='',human_units=False,
                label='%s back'%name,markersize=7)
    plt.ylabel('E')
    plt.xlabel('power (W)')
    plt.legend()#bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.)
    if save_figs:
        plt.savefig('%s_%s_enhancement.png'%(date,plotname),transparent=True,overwrite=True,bbox_inches='tight')
    plt.show()
    plt.close()
    
if plot_all and C is not None:
    figure(figsize=(7,5))
    title('%s %s\n$k_{\sigma}s(p)$'%(filename.split('_')[0],plotname))
    curves.labels('power',s.getaxis('power'))
    colors = ['dark magenta','brick red','dark orange','gold','greenish','dark cyan','grape','wine']
    for i,name in enumerate(names):
        plot(curves['sample',i]['power',:-3],color='xkcd:%s'%colors[i],marker='o',ls='',human_units=False,
                label='%s'%name,alpha=0.5)
        plot(curves['sample',i]['power',-3:],color='xkcd:%s'%colors[i],marker='x',ls='',human_units=False,
                label='%s back'%name,markersize=7)
    plt.ylabel('$k_{\sigma}s(p)$ $(s^{-1})$')
    plt.xlabel('power (W)')
    plt.legend()#bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.)
    if save_figs:
        plt.savefig('%s_%s_ksgima_sp.png'%(date,plotname),transparent=True,overwrite=True,bbox_inches='tight')
    plt.show()
    plt.close()
#}
