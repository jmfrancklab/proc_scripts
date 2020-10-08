from pyspecdata import *
with figlist_var() as fl:
    for exp, expno in [('w8.*200731',5),('w8.*201008',3)]:
        data = find_file(exp, exp_type='test_equip', expno=expno)
        print(f"for exp {exp} L20 is",data.get_prop('acq')['L'][20],"and td1 is",data.get_prop('acq')['TD1'],
                "td1/l20 is",data.get_prop('acq')['TD1']/data.get_prop('acq')['L'][20])
        fl.next(f'raw (t2,coh) data from {exp}')
        fl.image(data, interpolation='bilinear')
