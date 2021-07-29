from pyspecdata import *
import os
os.chdir('C:/Users/saman/proc_scripts/ODNP_SMB')
from pyspecProcScripts import *
from pyspecProcScripts import postproc_dict
from pyspecProcScripts.third_level.process_IR import process_IR
import time
import matplotlib.pyplot as plt
plt.rcParams.update({
    "figure.facecolor":  (1.0, 1.0, 1.0, 0.0),  # clear
    "axes.facecolor":    (1.0, 1.0, 1.0, 0.25),  # clear ## 1,1,1,0.9 == 90% transparent white
    "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
})
fl = fl_mod()
# {{{ input parameters
filename = '210715_A174R1a_pR_ODNP' #'210707_Q183R1a_pR_DDM_ODNP' #'210714_150uM_TEMPOL_SMB_ODNP'
powers = np.r_[0,0,0,0]#np.r_[0,2,2.5,3.16,3.98] #W -> [0,33,34,35,36] dBm
exp_type = 'odnp'
postproc = 'spincore_IR_v1'
export_csv = False

date = time.strftime('%y%m%d')
outname = '%s_%s'%(date,'_'.join(filename.split('_')[1:-1]))
titl = ' '.join(filename.split('_')[:-1])
T1_vals = nddata(np.zeros((len(powers))),'power').labels('power',powers)
T1w = 2.6 #s
coherence_pathway = {'ph1':0,'ph2':1}

for idx,filename,nodename,f_range,t_range,rep,clock_correction,IR,ILT in [
    (0,'210729_T177R1a_pR_KI_ODNP','FIR_36dBm',(-250,250),(None,83e-3),4,True,False,False),
#    (0,'210729_T177R1a_pR_DHPC_ODNP','FIR_36dBm',(-250,250),(None,83e-3),4,True,False,False),
#    (0,'210728_T177R1a_pR_DDM_ODNP','FIR_36dBm',(-250,250),(None,83e-3),4,False,False,False),
#    (0,'210714_A174R1a_pR_DDM_ODNP','FIR_0dBm',(-200,0),(None,83e-3),4,False,False,False),
#    (1,'210715_A174R1a_pR_KI_ODNP','FIR_0dBm',(-125,275),(None,83e-3),4,False,False,False),
#    (2,'210715_A174R1a_pR_KH2PO4_ODNP','FIR_0dBm',(-200,150),(None,83e-3),4,False,False,False),
#    (3,'210714_A174R1a_pR_DHPC_ODNP','FIR_0dBm',(-200,75),(None,83e-3),4,False,False,False),
# {210704_150uM_TEMPOL_SMB_ODNP
#    (0,filename,'FIR_0dBm',(-300,300),(None,83e-3),6.5,False,False,False),
#    (1,filename,'FIR_33dBm',(-300,400),(None,83e-3),6.5,False,False,False),#(-75,175)
#    (2,filename,'FIR_34dBm',(-300,400),(None,83e-3),6.5,False,False,False),#(-75,175)
#    (3,filename,'FIR_35dBm',(-300,300),(None,83e-3),6.5,False,False,False),#(-75,175)
#    (4,filename,'FIR_36dBm',(-300,300),(None,83e-3),6.5,False,False,False),#(-100,125),(None,75e-3)
# }
# {210707_Q183R1a_pR_DDM_ODNP
#    (0,filename,'FIR_noPower_real_newvd2',(-200,25),(None,60e-3),3,
#        False,False,False),
#    (1,filename,'FIR_33dBm_real_newvd2',(-175,150),(None,83e-3),3,
#        False,False,False),
#    (2,filename,'FIR_34dBm_real_newvd2',(-200,150),(None,65e-3),3,
#        False,False,False),
#    (3,filename,'FIR_35dBm_real_newvd2',(-200,150),(None,65e-3),3,
#        False,False,False),
#    (4,filename,'FIR_36dBm_real_newvd2',(-200,175),(None,50e-3),3,
#        False,False,False),
# }
        ]:
#}}}
    fl.basename = filename
    s = find_file(filename,exp_type=exp_type,expno=nodename,
            postproc=postproc,lookup=postproc_dict,fl=fl)
    myslice = s['t2':f_range]
    mysgn = select_pathway(myslice,coherence_pathway).real.sum('t2').run(np.sign)
#    mysgn = nddata(np.ones((len(myslice.getaxis('vd')))),'vd').labels('vd',myslice.getaxis('vd'))
    mysgn = nddata(r_[-1.,-1.,-1.,-1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'vd').labels('vd',myslice.getaxis('vd'))

#    if len(s.getaxis('vd')) == 12:
#        mysgn = nddata(r_[-1.,-1.,-1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'vd').labels('vd',s.getaxis('vd'))
#    if len(s.getaxis('vd')) == 16:
#        mysgn = nddata(r_[-1.,-1.,1.,-1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'vd').labels('vd',s.getaxis('vd'))
#    mysgn = nddata(np.ones((16)),'vd').labels('vd',s.getaxis('vd'))
    fl.next('\nshowing coherence channel zoom')
    fl.image(s.C['ph1',0]['ph2',1]['t2':(-750,750)]*mysgn)
    T1 = process_IR(s,rd=rep,f_range=f_range,t_range=t_range,clock_correction=clock_correction,IR=IR,
            flip=True,sgn=mysgn,fl=fl,hermitian_phasing=True,best_shift=0.004167) 
#    T1_vals['power',idx] = T1
    fl.show()
#print(T1_vals*1e6)
print(T1*1e6)
#{ Fit T1 curve
#m,b = np.polyfit(T1_vals.getaxis('power'),T1_vals.data,1)
#T1_fit = lambda p: (m*p) + b
#print('T1(p) fit: m = %0.6f s/W b = %0.6f s'%(m,b))
#m1,b1 = np.polyfit(T1_vals.getaxis('power'),1/T1_vals.data,1)
#R1_fit = lambda p: (m1*p) + b1
#poweraxis = r_[0:powers.max():100j]
#print('R1(p) fit: m = %0.6f W/s b = %0.6f 1/s'%(m1,b1))
##{ Plot all T1's together
#fig1 = figure(figsize=(7,5))
#title('%s $T_{1}$(p)'%titl)
#plot(T1_vals,color='xkcd:dark magenta',ls='',marker='o',label='data')
#plot(poweraxis,T1_fit(poweraxis),'g--',label='fit')
#plt.ylim(0,T1_vals.data.max())
#plt.ylabel('$T_{1}$ (s)')
#plt.xlabel('power (W)')
##plt.text(powers[1]/4,T1_vals.data[0],'T1(s) = %0.4f(s/W) * p (W) + %0.4f(s)'%(m,b),fontsize='medium')
#plt.legend(**dict(bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.))
#plt.savefig('%s_T1s.png'%outname,overwrite=True,bbox_inches='tight')
#plt.show()
#plt.close()
#R1_vals = nddata(1/T1_vals.data,'power').labels('power',powers)
#fig2 = figure(figsize=(7,5))
#title('%s $R_{1}$(p)'%titl)
#plot(R1_vals,'xkcd:dark magenta',ls='',marker='o',label='data')
#plot(poweraxis,R1_fit(poweraxis),'g--',label='fit')
#plt.ylim(0,R1_vals.data.max())
#plt.ylabel('$R_{1}$ $s^{-1}$')
#plt.xlabel('power (W)')
##plt.text(powers[0],(1/T1_vals.data[0])*0.8,'R1(1/s) = %0.4f (1/Ws) * p (W) + %0.4f (1/s)'%(m1,b1),fontsize='medium')
#plt.legend(**dict(bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.))
#plt.savefig('%s_R1s.png'%outname,overwrite=True,bbox_inches='tight')
#plt.show()
#plt.close()
##}
##{ Save T1's as csv
#if export_csv:
#    import pandas as pd
#    powerdf = pd.DataFrame(data=T1_vals.getaxis('power'),columns=['power (W)'])
#    T1df = pd.DataFrame(data=T1_vals.data,columns=['$T_1$ (s)'])
#    df = powerdf.join(T1df)
#    df.to_csv('%s_T1s.csv'%outname,index=False)
#}
