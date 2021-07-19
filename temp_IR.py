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
for thisfile,exp_type,nodename,postproc,f_range,t_range,rep,clock_correction,IR,ILT in [
        ('210715_A174R1a_pR_KH2PO4_ODNP','odnp',
            'FIR_36dBm','spincore_IR_v1',
            (-300,300),(None,83e-3),4,False,False,False),
#        ('210715_A174R1a_pR_KI_ODNP','odnp',
#            'FIR_36dBm','spincore_IR_v1',
#            (-300,300),(None,83e-3),4,False,False,False),
#        ('210714_A174R1a_pR_DDM_ODNP','odnp',
#            'FIR_36dBm','spincore_IR_v1',
#            (-300,300),(None,83e-3),4,False,False,False),
#        ('210714_A174R1a_pR_DHPC_ODNP','odnp',
#            'FIR_36dBm','spincore_IR_v1',
#            (-200,150),(None,83e-3),4,False,False,False),
#        ('210714_150uM_TEMPOL_SMB_ODNP','odnp',
#            'FIR_0dBm','spincore_IR_v1',
#            (-300,300),(None,83e-3),6.5,False,False,False),
#        ('210708_Q183R1a_pR_KH2PO4_ODNP','odnp',
#            'FIR_0dBm','spincore_IR_v1',
#            (-300,300),(None,83e-3),2.5,False,False,False),
#        ('210707_Q183R1a_pR_DHPC_ODNP','odnp',
#            'FIR_36dBm','spincore_IR_v1',
#            (-300,300),(None,83e-3),3,False,False,False),
#        ('210707_Q183R1a_pR_DDM_ODNP','odnp',
#            'FIR_noPower_real_newvd2','spincore_IR_v1',
#            (-300,300),(None,83e-3),3,False,False,False),
#        ('210707_Q183R1a_pR_DDM_ODNP','odnp',
#            'FIR_36dBm_real','spincore_IR_v1',
#            (-300,300),(None,83e-3),6,
#            False,False,False)
#       ('210617_T177R1a_pR_DDM_ODNP','odnp',
#           'FIR_36dBm','spincore_IR_v1',
#           (-200,200),(None,83e-3),6,
#           False,False,False),
#        ('210507_TEMPOL_150uM_cap_probe_FIR_36dBm','odnp',
#            'signal','spincore_IR_v1',
#            (-200,250),(None,83e-3),4,
#            False,False,False),
        ]:
    fl.basename = thisfile
    s = find_file(thisfile,exp_type=exp_type,expno=nodename,
            postproc=postproc,lookup=postproc_dict,fl=fl)
    if thisfile == '210617_T177R1a_pR_DDM_ODNP' and nodename == 'FIR_36dBm':
        s = s['vd',:-1]
    myslice = s['t2':f_range]
    mysgn = select_pathway(myslice,coherence_pathway).real.sum('t2').run(np.sign)
    mysgn = nddata(r_[-1.,-1.,-1.,-1.,-1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'vd').labels('vd',s.getaxis('vd'))
    fl.next('\nshowing coherence channel zoom')
    fl.image(s.C['ph1',0]['ph2',1]['t2':(-750,750)]*mysgn)
##    T1 = process_IR(s,W=rep,f_range=f_range,t_range=t_range,
##            clock_correction=clock_correction,IR=IR,flip=True,sgn=mysgn,fl=fl)
    T1 = process_IR(s,W=rep,f_range=f_range,t_range=t_range,clock_correction=clock_correction,IR=IR,
            flip=True,sgn=mysgn,fl=fl,hermitian_phasing=True,best_shift=0.004167) 
    fl.show()
