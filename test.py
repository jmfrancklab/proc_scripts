from pyspecdata import *
with figlist_var() as fl:
    for exp, exp_type, nodename in [
            ('w8_200731','test_equip',5),
            ('w8_201008','test_equip',3)
            ]:
        data = find_file(exp, exp_type=exp_type, expno=nodename)
        print(f"for exp {exp} L20 is",data.get_prop('acq')['L'][20],"and td1 is",data.get_prop('acq')['TD1'],
                "td1/l20 is",data.get_prop('acq')['TD1']/data.get_prop('acq')['L'][20])
        data.chunk('indirect',['indirect','ph1','ph2'],[-1,2,4])
    #{{{removes CP aspect
    #s = s['ph2',[1,3]]
        data.setaxis('ph1',r_[0,2]/4).setaxis('ph2',r_[0:4]/4)
        data.ft(['ph1','ph2'])
        data.reorder(['ph1','ph2','indirect'])
        with figlist_var() as fl:
            fl.next('raw data(t2,coh) data from {exp}')
            fl.image(data,interpolation='bilinear')
fl.show();quit()
