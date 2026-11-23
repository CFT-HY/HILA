"""Read the VTKHDF files written by test_hdf5 with the VTK reader (vtkHDFReader, as ParaView).
Usage: python3 check_vtk.py [dir]
Only the 3-dimensional grids are read: the VTK 9.1 reader (and ParaView 5.11) crashes on
grids of size 1 in y or z (also for files written with h5py), so vtk_slice.hdf and
vtk_proj.hdf need a newer VTK/ParaView."""
import sys
try:
    import vtk
except ImportError:
    print('vtk python module not available, skipping')
    sys.exit(0)
d = sys.argv[1] if len(sys.argv) > 1 else '.'
ok = True
def chk(name, cond):
    global ok
    print(('OK  ' if cond else 'FAIL') + ' ' + name)
    ok = ok and bool(cond)
def read(fn):
    r = vtk.vtkHDFReader()
    r.SetFileName(d + '/' + fn)
    r.Update()
    return r.GetOutput()
kval = lambda x, y, z: x + 100 * y + 10000 * z

o = read('vtk.hdf')
pd = o.GetPointData()
names = sorted(pd.GetArrayName(i) for i in range(pd.GetNumberOfArrays()))
ncomp = {n: pd.GetArray(n).GetNumberOfComponents() for n in names}
pid = o.FindPoint(3, 5, 7)
chk('vtk.hdf', o.GetDimensions() == (16, 8, 12) and names == ['c', 'k', 's', 'u']
    and ncomp == {'c': 2, 'k': 1, 's': 1, 'u': 18} and pd.GetArray('u').GetDataTypeAsString() == 'float'
    and o.GetPoint(pid) == (3, 5, 7) and pd.GetArray('k').GetTuple1(pid) == kval(3, 5, 7))

o = read('vtk_box.vtkhdf')
ok_box = o.GetDimensions() == (12, 6, 7) and o.GetOrigin() == (2, 1, 3) and o.GetBounds() == (2, 13, 1, 6, 3, 9)
o = read('vtk_k.vtkhdf')
pid = o.FindPoint(11, 2, 9)
chk('vtk_box.vtkhdf', ok_box)
chk('vtk_k.vtkhdf', o.GetDimensions() == (16, 8, 12)
    and o.GetPointData().GetArray('data').GetTuple1(pid) == kval(11, 2, 9))

print('ALL OK' if ok else 'SOME FAILED')
sys.exit(0 if ok else 1)
