from pyspecdata import *
date = '200303'
id_string = 'T1CPMG_AER'

np_dist = np.load('proc_'+date+'_'+id_string+'_1.npz')

nd_dist = nddata(np_dist['data'],['logT1','logT2'])

nd_dist.setaxis('logT1',np_dist['logT1'])
nd_dist.setaxis('logT2',np_dist['logT2'])

title('Dist')
image(nd_dist)
xlabel(r'$log(T_{1})$')
ylabel(r'$log(T_{2})$')
show()

