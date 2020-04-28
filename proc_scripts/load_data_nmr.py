from pyspecdata import *
def load_data_nmr(searchstr,expno):
    files = search_filename(searchstr, 'test_equip')
    s = find_file(searchstr, exp_type='test_equip', dimname = 'indirect', expno=expno)
    return s
