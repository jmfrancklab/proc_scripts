from pyspecdata import *
#to use type s = load_data("nameoffile")
def load_data_nut(searchstr):
    files = search_filename(searchstr, 'test_equip')
    assert len(files)==1, "I found %d files matching the pattern %s"%(len(files),searchstr)
    dirname, filename = os.path.split(files[0])
    nodename = 'nutation'
    s = nddata_hdf5(filename+'/nutation',
            directory=dirname)
    return s
