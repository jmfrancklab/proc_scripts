from pyspecdata import *
import os
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
filename = '210707_Q183R1a_pR_DDM_ODNP' #'210707_Q183R1a_pR_DDM_ODNP' #'210714_150uM_TEMPOL_SMB_ODNP'
exp_type = 'odnp'
postproc = 'spincore_IR_v1'
date = time.strftime('%y%m%d')
coherence_pathway = {'ph1':0,'ph2':1}

for filename,nodename,f_range,t_range,rep,clock_correction,IR,ILT in [
    ('210811_S179R1a_pR_DDM_KH2PO4_ODNP','FIR_0dBm',(-200,175),(None,83e-3),4,True,False,False),
#    ('210811_S179R1a_pR_DDM_KH2PO4_ODNP','FIR_36dBm',(-250,250),(None,83e-3),4,True,False,False),
#    ('210811_S179R1a_pR_DDM_KI_ODNP','FIR_36dBm_2',(-250,250),(None,83e-3),4,True,False,False),
#    ('210811_S179R1a_pR_DDM_KI_ODNP','FIR_0dBm',(-250,150),(None,83e-3),4,True,False,False),
#    ('210811_S179R1a_pR_DDM_KI_ODNP','FIR_36dBm',(-250,250),(None,83e-3),4,True,False,False),
#    ('210810_S179R1a_pR_DHPC_ODNP','FIR_0dBm',(-300,100),(None,83e-3),4,True,False,False),
#    ('210810_S179R1a_pR_DHPC_ODNP','FIR_36dBm',(-250,250),(None,83e-3),4,True,False,False),
#    ('210810_S179R1a_pR_DDM_ODNP','FIR_0dBm',(-250,175),(None,83e-3),4,True,False,False),
#    ('210810_S179R1a_pR_DDM_ODNP','FIR_36dBm',(-250,175),(None,83e-3),4,False,False,False),
#    ('210806_129uM_TEMPOL_capProbe','FIR_noPower_1',(-250,250),(None,83e-3),6,False,False,False),
#    ('210806_100mM_TEMPOL_capProbe','FIR_noPower_3',(-200,200),(None,83e-3),3,False,True,False),
#    ('210714_150uM_TEMPOL_SMB_ODNP','FIR_0dBm',(-150,150),(None,83e-3),6.5,False,False,False),
#    ('210728_T177R1a_pR_DDM_ODNP','FIR_36dBm',(-200,75),(None,83e-3),4,False,False,False),
#    ('210729_T177R1a_pR_DHPC_ODNP','FIR_36dBm',(-200,125),(None,83e-3),4,True,False,False),
#    ('210714_A174R1a_pR_DDM_ODNP','FIR_0dBm',(-200,0),(None,83e-3),4,False,False,False),
#    ('210714_A174R1a_pR_DHPC_ODNP','FIR_0dBm',(-200,75),(None,83e-3),4,False,False,False),
#    ('210707_Q183R1a_pR_DDM_ODNP','FIR_noPower_real_newvd2',(-225,75),(None,60e-3),3,False,False,False),
#    ('210707_Q183R1a_pR_DHPC_ODNP','FIR_0dBm',(-200,50),(None,83e-3),3,False,False,False),
#    ('210616_F195R1a_pR_DDM_ODNP','FIR_0dBm',(-225,125),(None,83e-3),5,False,False,False),
#    ('210623_F195R1a_pR_DHPC_ODNP','FIR_0dBm',(-200,150),(None,83e-3),6,False,False,False),
        ]:
#}}}
    fl.basename = filename
    s = find_file(filename,exp_type=exp_type,expno=nodename,
            postproc=postproc,lookup=postproc_dict,fl=fl)
    myslice = s['t2':f_range]
    mysgn = [select_pathway(myslice,coherence_pathway).real['t2',np.argmax(i)].run(np.sign) for i in abs(select_pathway(myslice,coherence_pathway).real.data)][0]
#    mysgn = nddata(r_[-1.,-1.,-1.,-1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'vd').labels('vd',myslice.getaxis('vd')) #S179DDM 36 & DHPC 0 and 36 & KH2PO4 36
    mysgn = nddata(r_[-1.,-1.,-1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'vd').labels('vd',myslice.getaxis('vd')) #S179DDM 0 & KI 36 and 0 & KH2PO4 0
#    mysgn = nddata(r_[-1.,-1.,-1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'vd').labels('vd',myslice.getaxis('vd'))
#    mysgn = select_pathway(myslice,coherence_pathway).real.sum('t2').run(np.sign)
#    mysgn = nddata(np.ones((len(myslice.getaxis('vd')))),'vd').labels('vd',myslice.getaxis('vd'))
#    if len(s.getaxis('vd')) == 16:
#        mysgn = nddata(r_[-1.,-1.,-1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'vd').labels('vd',myslice.getaxis('vd'))
#    elif len(s.getaixs('vd')) == 12:
#        mysgn = nddata(r_[-1.,-1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'vd').labels('vd',myslice.getaxis('vd'))
#    else: 
#        mysgn = nddata(r_[-1.,-1.,ones((30))],'vd').labels('vd',myslice.getaxis('vd'))
 #   mysgn = nddata(r_[-1.,-1.,-1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'vd').labels('vd',myslice.getaxis('vd'))
    if filename == '210729_T177R1a_pR_DHPC_ODNP':
        mysgn = nddata(r_[-1.,-1.,-1.,-1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'vd').labels('vd',myslice.getaxis('vd'))
    fl.next('\nshowing coherence channel zoom')
    fl.image(s.C['ph1',0]['ph2',1]['t2':(-750,750)]*mysgn)
    T1 = process_IR(s,rd=rep,f_range=f_range,t_range=t_range,clock_correction=clock_correction,IR=IR,
            flip=False,sgn=mysgn,fl=fl,hermitian_phasing=True,best_shift=0.004167) 
    fl.show()
print(T1*1e6)
