"""Validate test.h5 written by test_hdf5 with h5py.  Usage: python3 check_hdf5.py [dir]"""
import sys, h5py, numpy as np
d = sys.argv[1] if len(sys.argv) > 1 else '.'
f = h5py.File(d + '/test.h5', 'r')
sh = (12, 8, 16)
ok = True
def chk(name, cond):
    global ok
    print(('OK  ' if cond else 'FAIL') + ' ' + name)
    ok = ok and bool(cond)
raw = lambda n, dt, nc: np.fromfile(d + '/' + n, dtype=dt).reshape(sh + ((nc,) if nc > 1 else ()))
chk('s', np.array_equal(f['s'][...], raw('s.bin', '<f8', 1)))
chk('c', np.array_equal(f['c'][...], raw('c.bin', '<f8', 2)))
u = raw('u.bin', '<f8', 18)
chk('gauge/u', np.array_equal(f['gauge/u'][...], u))
chk('gauge/u_float', np.array_equal(f['gauge/u_float'][...], u.astype('<f4')))
z, y, x = np.meshgrid(np.arange(12), np.arange(8), np.arange(16), indexing='ij')
chk('ints/k', np.array_equal(f['ints/k'][...], x + 100 * y + 10000 * z))
chk('ints/cv', np.array_equal(f['ints/cv'][...], np.stack([x, y, z], axis=-1)))
chk('arrays/v', np.array_equal(f['arrays/v'][...], (np.arange(15) * 0.5).reshape(3, 5)))
chk('arrays/vc', np.array_equal(f['arrays/vc'][...], np.array([[i, -i] for i in range(4)], 'f4')))
chk('arrays/LLL', np.array_equal(f['arrays/' + 'L' * 300][...], (np.arange(15) * 0.5).astype('f4')))
a = f.attrs
chk('attrs', a['beta'] == 6.0 and a['ntraj'] == 100 and list(a['lattice size']) == [16, 8, 12]
    and a['name'] == b'hila test' and list(a['vec']) == [1.5, 2.5] and list(a['cplx']) == [1, 2]
    and a['i64'] == 1 << 40)
chk('dataset/group attrs', f['s'].attrs['what'] == b'gaussian' and f['gauge'].attrs['group attr'] == 3.25
    and f['empty/sub'].attrs['note'] == b'empty group')
print('ALL OK' if ok else 'SOME FAILED')
sys.exit(0 if ok else 1)
