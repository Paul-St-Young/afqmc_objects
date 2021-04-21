#!/usr/bin/env python3
import os
import numpy as np

def parse_kline(line, ik=None):
  from qharv.reel import ascii_out
  assert 'k(' in line
  ikt, kvect, wkt = line.split('=')
  myik = int(ascii_out.lr_mark(ikt, '(', ')'))
  if ik is not None:  # check k index
    assert ik == myik-1  # fortran 1-based indexing
  wk = float(wkt)
  klist = ascii_out.lr_mark(kvect, '(', ')').split()
  kvec = np.array(klist, dtype=float)
  return kvec, wk

def read_kpoints(scf_out):
  from qharv.reel import ascii_out
  mm = ascii_out.read(scf_out)
  # get lattice units
  alat = ascii_out.name_sep_val(mm, 'lattice parameter (alat)')
  blat = 2*np.pi/alat
  # start parsing k points
  idx = mm.find(b'number of k points')
  mm.seek(idx)
  # read first line
  #  e.g. number of k points=    32  Fermi-Dirac smearing ...
  line = mm.readline().decode()
  nk = int(line.split('=')[1].split()[0])
  # confirm units in second line
  line = mm.readline().decode()
  assert '2pi/alat' in line
  # start parsing kvectors
  data = np.zeros([nk, 4])  # ik, kx, ky, kz, wk
  for ik in range(nk):
    line = mm.readline().decode()
    kvec, wk = parse_kline(line, ik=ik)
    data[ik, :3] = kvec*blat
    data[ik, 3] = wk
  mm.close()
  return data

def read_settings(qeinp, qeout):
  from ase import io
  from afobj.basis.opt import get_input, get_elem, get_params
  atoms = io.read(qeinp, format='espresso-in')
  outdir = get_input(qeinp, 'outdir')
  elem = get_elem(qeinp)
  params = get_params(qeout)
  data = read_kpoints(qeout)
  kpts = data[:, :3]
  params['kpts'] = kpts
  return atoms, outdir, elem, params

def write_orbs(fout, scf_inp, scf_out, lmax, x, kc=None):
  from afqmctools.utils.optimizable_basis_set import default_basis_set
  from afobj.basis.gto_h5 import gen_qe_gto
  # read parameters
  atoms, outdir, elem, params = read_settings(scf_inp, scf_out)
  _, basis_set = default_basis_set(lmax, elem)
  nao = gen_qe_gto(atoms, basis_set, x, params['kpts'],
    fname=fout, mesh=params['mesh'], kc=kc)
  return nao

if __name__ == '__main__':
  msg = 'example:\n'
  msg += 'python3 write_pyscf_orbitals.py orbitals.h5'
  msg += 'scf.inp scf.out 2 x0.dat\n'
  from argparse import ArgumentParser
  parser = ArgumentParser(description=msg)
  parser.add_argument('forb', type=str)
  parser.add_argument('scf_inp', type=str)
  parser.add_argument('scf_out', type=str)
  parser.add_argument('lmax', type=int)
  parser.add_argument('fx0', type=str)
  parser.add_argument('--kcut', '-kc', type=float, default=None)
  args = parser.parse_args()
  kc = args.kcut

  lmax = args.lmax
  try:
    x = np.loadtxt(args.fx0)
  except ValueError:
    x = args.fx0
  nao = write_orbs(args.forb, args.scf_inp, args.scf_out, lmax, x, kc=kc)
  print(nao)
# end __main__
