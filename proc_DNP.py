from pyspecdata import *
from proc_scripts import *
from proc_scripts import postproc_dict
from proc_scripts.third_level.process_IR import process_IR
from proc_scripts.third_level.process_enhancement import process_enhancement
from sympy import symbols, Symbol, latex,limit,init_printing
fl = fl_mod()
# {{{ input parameters
save_npz = False
power_list = r_[0.001,0.5,1,1.5,2]
R1w = 1/2.172
C = 0.000007
signal_pathway = {'ph1':1,'ph2':-2}
excluded_pathways = [(0,3),(0,0)]
nPowers=25
#}}}
T1_list = []
for thisfile,exp_type,nodename,postproc,f_range,t_range,IR,ILT in [
       ('210525_TEMPOL7uM_cap_probe_DNP','ODNP_NMR_comp/test_equipment','0p001_Watts','spincore_IR_v1',
           (-0.5e3,0.8e3),(None,83e-3),True,False),
        ('210525_TEMPOL7uM_cap_probe_DNP','ODNP_NMR_comp/test_equipment','0p5_Watts','spincore_IR_v1',
           (-0.5e3,0.8e3),(None,83e-3),True,False),
        ('210525_TEMPOL7uM_cap_probe_DNP','ODNP_NMR_comp/test_equipment','1_Watts','spincore_IR_v1',
           (-0.5e3,0.8e3),(None,83e-3),True,False),
        ('210525_TEMPOL7uM_cap_probe_DNP','ODNP_NMR_comp/test_equipment','1p5_Watts','spincore_IR_v1',
           (-0.5e3,0.8e3),(None,83e-3),True,False),
        ('210525_TEMPOL7uM_cap_probe_DNP','ODNP_NMR_comp/test_equipment','2_Watts','spincore_IR_v1',
           (-0.5e3,0.8e3),(None,83e-3),True,False),
        ('210525_TEMPOL7uM_cap_probe_DNP','ODNP_NMR_comp/test_equipment','2p5_Watts','spincore_IR_v1',
           (-0.5e3,0.8e3),(None,83e-3),True,False),

        ]:
    s = find_file(thisfile,exp_type=exp_type,expno=nodename,
            postproc=postproc,lookup=postproc_dict,fl=fl)
   #fl.show();quit()
    T1 = process_IR(s,label=thisfile,W=6,f_range=f_range,IR=False,fl=fl) 
    T1_list.append(T1)
T1_list.pop(-1)
searchstr ='210525_TEMPOL7uM_cap_probe_DNP' 
d = find_file(searchstr,exp_type='ODNP_NMR_comp/test_equipment',
        expno='enhancement',postproc='spincore_ODNP_v1',lookup=postproc_dict,fl=fl)
enhancement,idx_maxpower = process_enhancement(d, searchstr = searchstr,freq_range=(-8e2,8e2),
        fl=fl)
T1p = nddata(T1_list,[-1],['power']).setaxis('power',power_list)
R1p = T1p**-1   
fl.next(r'$T_{1}$(p) vs power')
fl.plot(T1p,'o')
#{{{making Flinear and fitting
Flinear = ((R1p - R1p['power':0.001] + R1w)**-1)
polyorder = 1
coeff,_ = Flinear.polyfit('power',order=polyorder)
power = nddata(np.linspace(0,R1p.getaxis('power')[-1],nPowers),'power')
Flinear_fine = 0
for j in range(polyorder + 1):
    Flinear_fine += coeff[j] * power **j
fl.next('Flinear')
Flinear.set_units('power',enhancement.get_units('power'))
fl.plot(Flinear,'o',label='Flinear')
Flinear_fine.set_units('power',enhancement.get_units('power'))
fl.plot(Flinear_fine,label='Flinear_fine')
plt.title('polynomial fit of linear equation')
plt.ylabel("$F_{linear}$")
fl.next('R1p vs power')
R1p_fine = ((Flinear_fine)**-1) + R1p['power':0.001]-R1w
fl.plot(R1p,"x")
R1p_fine.set_units('power',R1p.get_units('power'))
fl.plot(R1p_fine)
plt.title("relaxation rates")
plt.ylabel("$R_1(p)$")
#{{{plotting with correcting for heating
ksigs_T=(0.0015167/C)*(1-enhancement['power',:idx_maxpower+1])*(R1p_fine)
fl.next('ksig_smax for %s'%searchstr)
ksigs_T.set_units('power','mW')
#}}}
#{{{plotting with correction for heating
x = enhancement['power',:idx_maxpower+1].fromaxis('power')
fitting_line = fitdata(ksigs_T)
k,p_half,power = symbols("k, p_half, power",real=True)
fitting_line.functional_form = (k*power)/(p_half+power)
fitting_line.fit()
fl.plot(ksigs_T,'o',label='with heating correction')
fl.plot(ksigs_T.imag,'o',label='imaginary')
fit = fitting_line.eval(25)
fit.set_units('power','mW')
#fl.plot(fit,label='fit')
#plt.text(0.75, 0.25, fitting_line.latex(), transform=plt.gca().transAxes,size='large',
#        horizontalalignment='center',color='k')
plt.title('ksigmas(p) vs Power')
plt.ylabel('ksigmas(p)')
    #}}}
fl.show();quit()

