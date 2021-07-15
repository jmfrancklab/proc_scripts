from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import postproc_dict
from pyspecProcScripts.third_level.process_IR import process_IR
fl = fl_mod()
# {{{ input parameters
save_npz = False
T1w = 2.6 #s
#}}}
coherence_pathway = {'ph1':0,'ph2':1}
exp_type = 'odnp'
postproc = 'spincore_IR_v1'
export_csv = True
outname = '210707_Q183R1a_pR_DDM'
powers = np.r_[0,2,2.5,3.16,3.98] #W
T1_vals = nddata(np.zeros((len(powers))),'power').labels('power',powers)
for idx,thisfile,nodename,f_range,t_range,rep,clock_correction,IR,ILT in [
        (0,'210707_Q183R1a_pR_DDM_ODNP','FIR_noPower_real_newvd2',
            (-200,25),(None,60e-3),3,False,False,False),
        (1,'210707_Q183R1a_pR_DDM_ODNP','FIR_33dBm_real_newvd2',
            (-175,150),(None,83e-3),3,False,False,False),
        (2,'210707_Q183R1a_pR_DDM_ODNP','FIR_34dBm_real_newvd2',
            (-200,150),(None,65e-3),3,False,False,False),
        (3,'210707_Q183R1a_pR_DDM_ODNP','FIR_35dBm_real_newvd2',
            (-200,150),(None,65e-3),3,False,False,False),
        (4,'210707_Q183R1a_pR_DDM_ODNP','FIR_36dBm_real_newvd2',
            (-200,175),(None,50e-3),3,False,False,False),
        ]:
    fl.basename = thisfile
    s = find_file(thisfile,exp_type=exp_type,expno=nodename,
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
#fitfxn = lambda p,x: (p[0]*x) + p[1]
#errfxn = lambda p_arg,x_arg,y_arg: abs(fitfxn(p_arg,x_arg) - y_arg)
m,b = np.polyfit(T1_vals.getaxis('power'),T1_vals.data,1)
print('T1(p) fit: $T_{1}=%0.6fp+%0.6f$'%(m,b))
print('R1(p): $R_{1} = %0.6fp+%0.6f$'%(1/m,1/b))
#{ Plot all T1's together
figure(figsize=(7,5))
title('$T_{1}$(p) for %s'%(' '.join(outname.split('_'))))
plot(T1_vals,color='xkcd:dark magenta',ls='-',marker='o')
plt.ylabel('$T_{1}$ (s)')
plt.xlabel('power (W)')
plt.text(1,powers[-2],'T1 = %0.4fp + %0.4f'%(m,b))
plt.savefig('%s_T1s.png'%outname,overwrite=True,bbox_inches='tight')
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
