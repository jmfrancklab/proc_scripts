from pyspecdata import *
import tables as tb
#f = tb.open_file('C:/Users/saman/Research/Data/pR_ODNP/210611_S175R1a_pR_DDM_field_dep.h5')
f = tb.open_file('C:/Users/saman/Research/Data/pR_ODNP/210616_S175R1a_pR_DDM_ODNP.h5')
for node in f.walk_nodes():
    print(node)
