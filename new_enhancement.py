#{Imports & Initialization
from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import postproc_dict
from sympy import symbols, Symbol, latex,limit,init_printing
from numpy import *
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
file_location = 'odnp'
postproc = 'spincore_ODNP_v2' # 'spincore_ODNP_v1'
plot_ks = True
plot_enhancements = True
plotname = 'Q183R1a pR (210708) vs TEMPOL (210714)'
names = ['TEMPOL','pR in DDM','pR in DHPC']#['DDM KCl','DHPC KCl','DDM KI', 'DDM $KH_{2}PO_{4}$']
names = ['TEMPOL','DDM KCl','DHPC KCl','DDM KI', 'DDM $KH_{2}PO_{4}$']
curves = nddata(zeros((len(names),18),dtype='complex128'),['sample','power'])
enhancements = nddata(zeros((len(names),18),dtype='complex128'),['sample','power'])

for (idx,filename,nodename,outname,f_range,C,T1_0,ppt,export_csv) in [
        (0,'210714_150uM_TEMPOL_SMB_ODNP','enhancement_real',
            '210714_150uM_TEMPOL_SMB_enhancement',
            (-200,200),150e-6,3.56,1.5163,False),
#        (1,'210714_A174R1a_pR_DDM_ODNP','enhancement1',
#            '210714_A174R1a_pR_DDM_enhancement',
#            (-250,100),240e-6,1.49,1.5154,False),
#        (2,'210714_A174R1a_pR_DHPC_ODNP','enhancement',
#            '210714_A174R1a_pR_DHPC_enhancement',
#            (-250,75),238.5e-6,1.47,1.5154,False),
    
        (1,'210707_Q183R1a_pR_DDM_ODNP','enhancement',
            '210707_Q183R1a_pR_DDM_enhancement',
            (-225,75),207.4e-6,1.89,1.5154,True),
        (2,'210707_Q183R1a_pR_DHPC_ODNP','enhancement',
            '210707_Q183R1a_pR_DHPC_enhancement',
            (-225,75),103.2e-6,2.15,1.5154,True),
        (3,'210708_Q183R1a_pR_KI_ODNP','enhancement',
            '210708_Q183R1a_pR_KI_enhancement',
            (-200,100),113.1e-6,1.12,1.5154,True),
        (4,'210708_Q183R1a_pR_KH2PO4_ODNP','enhancement',
            '210708_Q183R1a_pR_KH2PO4_enhancement',
            (-225,75),115.2e-6,1.9,1.5154,True)
        ]:
    s = find_file(filename,exp_type=file_location,expno=nodename,
            postproc=postproc,lookup=postproc_dict,fl=fl)
    outname = filename.split('_')
    outname = '_'.join(outname[:-1])
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
    #s *= (207.4/113.1)
    s /= s['power',0]
    fl.next('enhancement curve')
    fl.plot(s['power',:-3],'ko', human_units=False)
    fl.plot(s['power',-3:],'ro', human_units=False)
    if plot_enhancements:
        enhancements['sample',idx] = s
    if C is not None:
        epsilon = 1 - s
        ks = epsilon/C
        if T1_0 is not None:
            ks /= T1_0
            if ppt is not None:
                ks *= ppt*1e-3 # gets you in units of k_sigma!!!
                fl.next('$k_{\sigma}s(p)$ using $T_{1}(0)$')
                fl.plot(ks['power',:-3],'ko', human_units=False)
                fl.plot(ks['power',-3:],'ro', human_units=False)
                if plot_ks:
                    curves['sample',idx] = ks
#}
# { Exporting data from the enhancement curve as csv
    if export_csv:
        import pandas as pd
        reenhancementdf = pd.DataFrame(data=s.real.data,columns=['Re[E(p)]'])
        imenhancementdf = pd.DataFrame(data=s.imag.data,columns=['Im[E(p)]'])
        enhancementdf = pd.DataFrame(data=s.data,columns=['enhancement'])
        powerdf = pd.DataFrame(data=s.getaxis('power'),columns=['power'])
        enhancement = powerdf.join(reenhancementdf)
        enhancement2 = enhancement.join(imenhancementdf)
        enhancement3 = enhancement2.join(enhancementdf)
        enhancement3.to_csv('%s_enhancement.csv'%outname,index=False)
#}
#    fl.show()
#{ Plotting all results together
if plot_enhancements:
    figure(figsize=(7,5))
    title('%s enhancements'%plotname)
    enhancements.labels('power',s.getaxis('power'))
    colors = ['dark magenta','brick red','dark orange','gold','greenish','dark cyan','grape','wine']   
    for i,name in enumerate(names):
        plot(enhancements['sample',i]['power',:-3],color='xkcd:%s'%colors[i],marker='o',ls='',human_units=False,
                label='%s'%name,alpha=0.5)
        plot(enhancements['sample',i]['power',-3:],color='xkcd:%s'%colors[i],marker='x',ls='',human_units=False,
                label='%s back'%name,markersize=7)
    plt.ylabel('E')
    plt.xlabel('power (W)')
    plt.legend(bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.)
    plt.savefig('%s_enhancement.png'%('_'.join(plotname.split(' '))),transparent=True,overwrite=True,bbox_inches='tight')
    plt.close()

if plot_ks:
    figure(figsize=(7,5))
    title('%s $k_{\sigma}s(p)$'%plotname)
    curves.labels('power',s.getaxis('power'))
    colors = ['dark magenta','brick red','dark orange','gold','greenish','dark cyan','grape','wine']
    for i,name in enumerate(names):
        plot(curves['sample',i]['power',:-3],color='xkcd:%s'%colors[i],marker='o',ls='',human_units=False,
                label='%s'%name,alpha=0.5)
        plot(curves['sample',i]['power',-3:],color='xkcd:%s'%colors[i],marker='x',ls='',human_units=False,
                label='%s back'%name,markersize=7)
    plt.ylabel('$k_{\sigma}s(p)$ $(s^{-1})$')
    plt.xlabel('power (W)')
    plt.legend(bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.)
    plt.savefig('%s_ksgima_sp.png'%('_'.join(plotname.split(' '))),transparent=True,overwrite=True,bbox_inches='tight')
    plt.close()
plt.show()
#}
