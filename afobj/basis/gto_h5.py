import numpy as np

def write_esh5_orbitals(cell, name, kpts):
  import h5py
  from pyscf.pbc.dft import numint
  from pyscf.pbc import tools
  def to_qmcpack_complex(array):
    shape = array.shape
    return array.view(np.float64).reshape(shape+(2,))
  nao = cell.nao_nr()

  fh5 = h5py.File(name, 'w')
  coords = cell.gen_uniform_grids(cell.mesh)
  
  kpts = np.asarray(kpts)
  nkpts = len(kpts)
  norbs = np.zeros((nkpts,),dtype=int)
  norbs[:] = nao

  grp = fh5.create_group("OrbsG")
  dset = grp.create_dataset("reciprocal_vectors", data=cell.reciprocal_vectors())
  dset = grp.create_dataset("number_of_kpoints", data=len(kpts))
  dset = grp.create_dataset("kpoints", data=kpts)
  dset = grp.create_dataset("number_of_orbitals", data=norbs)
  dset = grp.create_dataset("fft_grid", data=cell.mesh)
  dset = grp.create_dataset("grid_type", data=int(0))
  nnr = cell.mesh[0]*cell.mesh[1]*cell.mesh[2]
  # loop over kpoints later
  aol = []
  for (ik,k) in enumerate(kpts):
    ao = numint.KNumInt().eval_ao(cell, coords, k)[0]
    fac = np.exp(-1j * np.dot(coords, k))
    for i in range(norbs[ik]):
      aoi = fac * np.asarray(ao[:,i].T, order='C')
      aoi_G = tools.fft(aoi, cell.mesh)
      aoi_G = aoi_G.reshape(cell.mesh).transpose(2,1,0).reshape(nnr)
      aol.append(aoi_G)
      dset = grp.create_dataset('kp'+str(ik)+'_b'+str(i),
        data=to_qmcpack_complex(aoi_G)
      )
  fh5.close()
  return np.array(aol)

def gen_cell_no_basis(atoms):
  from pyscf.pbc.gto import Cell
  cell = Cell()
  cell.from_ase(atoms)
  return cell

def gen_cell(atoms, bset, x, mesh=None, prec=1e-12, verbose=0):
  from pyscf import gto as molgto
  from pyscf.pbc.gto import Cell
  # extract info from ase atoms
  axes = atoms.get_cell()
  elem = atoms.get_chemical_symbols()
  pos = atoms.get_positions()
  cell = Cell()
  cell.verbose = verbose  # shut the cell up
  cell.units = 'A'
  cell.precision = prec
  cell.a = axes
  atext = ''
  for name, r in zip(elem, pos):
    atext += '%s %.6f %.6f %.6f\n' % (name, *r)
  cell.atom = atext
  cell.pseudo = 'gthpbe'  # !!!! wut
  if mesh is not None:
    cell.mesh = mesh
  # build basis
  basis = {}
  for atm in elem:
    basis.update({atm: molgto.parse(
      bset.basis_str(atm, x)
    )})
  cell.basis = basis
  cell.build()
  return cell

def gen_qe_gto(atoms, bset, x, kpts, fname='pyscf.orbitals.h5',
  mesh=None, prec=1e-12, verbose=0):
  assert(len(x) == bset.number_of_params)
  cell = gen_cell(atoms, bset, x, mesh=mesh, prec=prec, verbose=verbose)
  nao = cell.nao_nr()
  aos = write_esh5_orbitals(cell, fname, kpts=kpts)
  #return aos
  return nao
