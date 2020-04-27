from pyspecdata import *
#to use type s = load_data("nameoffile")
def load_data(searchstr):
    files = search_filename(searchstr, 'test_equip')
    assert len(files)==1, "I found %d files matching the pattern %s"%(len(files),searchstr)
    dirname, filename = os.path.split(files[0])
    nodename = 'signal'
    s = nddata_hdf5(filename+'/signal',
            directory=dirname)
    return s

