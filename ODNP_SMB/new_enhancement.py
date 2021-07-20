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
plotname = 'Q183R1a pR Hoffmeister Series Anions'# 'Q183R1a pR DDM' 
names = ['KCl','KI','$KH_{2}PO_{4}$']#['TEMPOL','DDM KCl','DHPC KCl','DDM KI', 'DDM $KH_{2}PO_{4}$']
plot_all_ks = True 
plot_all_enhancements = True
save_figs = True
if save_figs:
    import os
    os.chdir('C:/Users/saman/Research/Sam_Notebooks/ODNP_proc/')
curves = nddata(zeros((len(names),18),dtype='complex128'),['sample','power'])
enhancements = nddata(zeros((len(names),18),dtype='complex128'),['sample','power'])
    # { Load T1_values
#import pandas as pd
#path = 'C:/Users/saman/Research/Sam_Notebooks/ODNP_proc'
#T1_df = pd.read_csv('%s/210707_Q183R1a_pR_DDM_T1s.csv'%path)
#print(T1_df);quit()
    # }
for (idx,filename,nodename,outname,f_range,C,T1_0,T1_vals,ppt,export_csv) in [
#        (0,'210714_150uM_TEMPOL_SMB_ODNP','enhancement_real',
#            '210714_150uM_TEMPOL_SMB_enhancement',
#            (-200,200),150e-6,3.56,None,1.5163,False),

#        (0,'210714_A174R1a_pR_DDM_ODNP','enhancement1',
#            '210714_A174R1a_pR_DDM_enhancement',
#            (-250,100),240e-6,1.49,None,1.5154,False),
#        (1,'210715_A174R1a_pR_KI_ODNP','enhancement',
#            '210715_A174R1a_pR_KI_enhancement',
#            (-200,200),236.9e-6,1.45,None,1.5154,False),
#        (2,'210715_A174R1a_pR_KH2PO4_ODNP','enhancement',
#            '210715_A174R1a_pR_KH2PO4_enhancement',
#            (-200,200),162.1e-6,2.05,None,1.5154,False),
#        (3,'210714_A174R1a_pR_DHPC_ODNP','enhancement',
#            '210714_A174R1a_pR_DHPC_enhancement',
#            (-250,75),238.5e-6,1.47,None,1.5154,False),
    
        (0,'210707_Q183R1a_pR_DDM_ODNP','enhancement',
            '210707_Q183R1a_pR_DDM_enhancement',
            (-225,75),207.4e-6,None,None,1.5154,True), # have T1_values for this data 
        (1,'210708_Q183R1a_pR_KI_ODNP','enhancement',
            '210708_Q183R1a_pR_KI_enhancement',
            (-200,100),113.1e-6,1.12,None,1.5154,True),
        (2,'210708_Q183R1a_pR_KH2PO4_ODNP','enhancement',
            '210708_Q183R1a_pR_KH2PO4_enhancement',
            (-225,75),115.2e-6,1.9,None,1.5154,True)
#        (3,'210707_Q183R1a_pR_DHPC_ODNP','enhancement',
#            '210707_Q183R1a_pR_DHPC_enhancement',
#            (-225,75),103.2e-6,2.15,None,1.5154,True),

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
    if plot_all_enhancements:
        enhancements['sample',idx] = s
    if C is not None:
        epsilon = 1 - s
        ks = epsilon/C
        if T1_0 is not None:
            ks /= T1_0
            fl.next('$k_{\sigma}s(p)T_{1}(p)$ $\div$ $ T_{1}(0)$')
        elif T1_vals is not None:
            m,b = np.polyfit(T1_vals.getaxis('power'),T1_vals.data,1)
            T1_p = lambda p: (m*p)+b
            T1_p = nddata(T1_p(ks.getaxis('power')),ks.getaxis('power')).labels('power',ks.getaxis('power'))
            ks /= T1_p
            fl.next('$k_{\sigma}s(p)$')
        if ppt is not None:
            ks *= ppt*1e-3 # gets you in units of k_sigma!!!
            fl.plot(ks['power',:-3],'ko', human_units=False)
            fl.plot(ks['power',-3:],'ro', human_units=False)
            if plot_all_ks:
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
fl.show()
#{ Plotting all results together
if plot_all_enhancements:
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
    plt.legend()#bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.)
    if save_figs:
        plt.savefig('%s_enhancement.png'%('_'.join(plotname.split(' '))),transparent=True,overwrite=True,bbox_inches='tight')
    plt.show()
    plt.close()
    
if plot_all_ks:
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
    plt.legend()#bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.)
    if save_figs:
        plt.savefig('%s_ksgima_sp.png'%('_'.join(plotname.split(' '))),transparent=True,overwrite=True,bbox_inches='tight')
    plt.show()
    plt.close()
#}
