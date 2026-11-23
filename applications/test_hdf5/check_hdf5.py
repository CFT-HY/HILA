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

# subvolumes, slices, projections; arrays are indexed [z, y, x, ...]
sa = raw('s.bin', '<f8', 1)
ca = raw('c.bin', '<f8', 2)
ka = x + 100 * y + 10000 * z
box = lambda a, lo, hi: a[lo[2]:hi[2] + 1, lo[1]:hi[1] + 1, lo[0]:hi[0] + 1]
chk('sub/s', np.array_equal(f['sub/s'][...], box(sa, (2, 1, 3), (13, 6, 9)))
    and list(f['sub/s'].attrs['cmin']) == [2, 1, 3] and list(f['sub/s'].attrs['cmax']) == [13, 6, 9])
chk('sub/k_small', np.array_equal(f['sub/k_small'][...], box(ka, (1, 1, 1), (5, 2, 3))))
chk('sub/u_float', f['sub/u_float'].dtype == np.float32
    and np.array_equal(f['sub/u_float'][...], u.astype('<f4')))
chk('slice/k_xz', np.array_equal(f['slice/k_xz'][...], ka[:, 3, :])
    and list(f['slice/k_xz'].attrs['slice']) == [-1, 3, -1])
chk('slice/c_y', np.array_equal(f['slice/c_y'][...], ca[7, :, 5]))
chk('slice/k_point', np.array_equal(f['slice/k_point'][...], [ka[6, 5, 4]]))
chk('proj/k_y', np.array_equal(f['proj/k_y'][...], ka.sum(axis=1))
    and list(f['proj/k_y'].attrs['projected directions']) == [1])
chk('proj/s_xz', np.allclose(f['proj/s_xz'][...], box(sa, (2, 1, 0), (9, 6, 10)).sum(axis=(0, 2)),
                             rtol=1e-12, atol=1e-12))
chk('proj/c_all', np.allclose(f['proj/c_all'][...], ca.sum(axis=(0, 1, 2)), rtol=1e-12, atol=1e-12)
    and f['proj/c_all'].shape == (2,))
chk('proj/u_z_float', f['proj/u_z_float'].dtype == np.float32
    and np.allclose(f['proj/u_z_float'][...], u.sum(axis=0), rtol=1e-5, atol=1e-5))

# VTKHDF files: datasets [nz, ny, nx(, ncomp)] in /VTKHDF/PointData, grid attributes
def vtk_grid(fn, extent, origin):
    g = h5py.File(d + '/' + fn, 'r')['VTKHDF']
    a = g.attrs
    return (list(a['Version']) == [1, 0] and a['Type'] == b'ImageData'
            and list(a['WholeExtent']) == extent and list(a['Origin']) == origin
            and list(a['Spacing']) == [1, 1, 1] and list(a['Direction']) == [1, 0, 0, 0, 1, 0, 0, 0, 1]), g['PointData']
okg, p = vtk_grid('vtk.hdf', [0, 15, 0, 7, 0, 11], [0, 0, 0])
chk('vtk.hdf', okg and sorted(p.keys()) == ['c', 'k', 's', 'u'] and np.array_equal(p['s'][...], sa)
    and np.array_equal(p['c'][...], ca) and np.array_equal(p['u'][...], u.astype('<f4'))
    and np.array_equal(p['k'][...], ka) and h5py.File(d + '/vtk.hdf', 'r').attrs['beta'] == 6.0)
okg, p = vtk_grid('vtk_box.vtkhdf', [0, 11, 0, 5, 0, 6], [2, 1, 3])
chk('vtk_box.vtkhdf', okg and np.array_equal(p['data'][...], box(sa, (2, 1, 3), (13, 6, 9))))
okg, p = vtk_grid('vtk_slice.hdf', [0, 15, 0, 11, 0, 0], [0, 0, 0])
chk('vtk_slice.hdf', okg and np.array_equal(p['data'][...], ca[:, 3, :].reshape(1, 12, 16, 2)))
okg, p = vtk_grid('vtk_proj.hdf', [0, 13, 0, 8, 0, 0], [1, 2, 0])
chk('vtk_proj.hdf', okg and np.array_equal(p['data'][...],
                                           box(ka, (1, 0, 2), (14, 7, 10)).sum(axis=1).reshape(1, 9, 14)))
okg, p = vtk_grid('vtk_k.vtkhdf', [0, 15, 0, 7, 0, 11], [0, 0, 0])
chk('vtk_k.vtkhdf', okg and np.array_equal(p['data'][...], ka))

print('ALL OK' if ok else 'SOME FAILED')
sys.exit(0 if ok else 1)
