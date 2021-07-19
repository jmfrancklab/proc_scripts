from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import postproc_dict
from pyspecProcScripts.third_level.process_IR import process_IR
import time
fl = fl_mod()
# {{{ input parameters
filename = '210707_Q183R1a_pR_DDM_ODNP'
powers = np.r_[0,2,2.5,3.16,3.98] #W
exp_type = 'odnp'
postproc = 'spincore_IR_v1'
export_csv = True

date = time.strftime('%y%m%d')
outname = '%s_%s'%(date,'_'.join(filename.split('_')[1:-1]))
titl = ' '.join(filename.split('_')[:-1])
T1_vals = nddata(np.zeros((len(powers))),'power').labels('power',powers)
T1w = 2.6 #s
coherence_pathway = {'ph1':0,'ph2':1}

for idx,filename,nodename,f_range,t_range,rep,clock_correction,IR,ILT in [
    (0,filename,'FIR_noPower_real_newvd2',(-200,25),(None,60e-3),3,
        False,False,False),
    (1,filename,'FIR_33dBm_real_newvd2',(-175,150),(None,83e-3),3,
        False,False,False),
    (2,filename,'FIR_34dBm_real_newvd2',(-200,150),(None,65e-3),3,
        False,False,False),
    (3,filename,'FIR_35dBm_real_newvd2',(-200,150),(None,65e-3),3,
        False,False,False),
    (4,filename,'FIR_36dBm_real_newvd2',(-200,175),(None,50e-3),3,
        False,False,False),
        ]:
#}}}
    fl.basename = filename
    s = find_file(filename,exp_type=exp_type,expno=nodename,
            postproc=postproc,lookup=postproc_dict,fl=fl)
    myslice = s['t2':f_range]
    mysgn = select_pathway(myslice,coherence_pathway).real.sum('t2').run(np.sign)
    fl.next('\nshowing coherence channel zoom')
    fl.image(s.C['ph1',0]['ph2',1]['t2':(-750,750)]*mysgn)
##    T1 = process_IR(s,W=rep,f_range=f_range,t_range=t_range,
##            clock_correction=clock_correction,IR=IR,flip=True,sgn=mysgn,fl=fl)
    T1 = process_IR(s,W=rep,f_range=f_range,t_range=t_range,clock_correction=clock_correction,IR=IR,
            flip=False,sgn=mysgn,fl=fl,hermitian_phasing=True,best_shift=0.004167) 
    T1_vals['power',idx] = T1
#    fl.show()
#{ Fit T1 curve
m,b = np.polyfit(T1_vals.getaxis('power'),T1_vals.data,1)
print('T1(p) fit: m = %0.6f s/W b = %0.6f s'%(m,b))
m1,b1 = np.polyfit(T1_vals.getaxis('power'),1/T1_vals.data,1)
print('R1(p) fit: m = %0.6f W/s b = %0.6f 1/s'%(m1,b1))
#{ Plot all T1's together
figure(figsize=(7,5))
title('%s $T_{1}$(p)'%titl)
plot(T1_vals,color='xkcd:dark magenta',ls='-',marker='o')
plt.ylabel('$T_{1}$ (s)')
plt.xlabel('power (W)')
plt.text(powers[1]/4,T1_vals.data[0],'T1(s) = %0.4f(s/W) * p (W) + %0.4f(s)'%(m,b),fontsize='medium')
plt.savefig('%s_T1s.png'%outname,overwrite=True,bbox_inches='tight')
plt.show()
plt.close()
figure(figsize=(7,5))
title('%s $R_{1}$(p)'%titl)
plot(powers,1/T1_vals.data,'xkcd:dark magenta',ls='-',marker='o')
plt.ylabel('$R_{1}$ $s^{-1}$')
plt.xlabel('power (W)')
plt.text(powers[0],(1/T1_vals.data[0])*0.8,'R1(1/s) = %0.4f (1/Ws) * p (W) + %0.4f (1/s)'%(m1,b1),fontsize='medium')
plt.savefig('%s_R1s.png'%outname,overwrite=True,bbox_inches='tight')
plt.show()
plt.close()
#}
#{ Save T1's as csv
if export_csv:
    import pandas as pd
    powerdf = pd.DataFrame(data=T1_vals.getaxis('power'),columns=['power (W)'])
    T1df = pd.DataFrame(data=T1_vals.data,columns=['$T_1$ (s)'])
    df = powerdf.join(T1df)
    df.to_csv('%s_T1s.csv'%outname,index=False)
#}
