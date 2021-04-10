import numpy as np

def get_orbs(forb):
  import h5py
  fp = h5py.File(forb, 'r')
  alat = fp['OrbsG/alat'][()][0]
  kpts = fp['OrbsG/kpoints'][()]
  nbl = fp['OrbsG/number_of_orbitals'][()]
  cmatl = []
  kvl = []
  for ik, (kpt, nb) in enumerate(zip(kpts, nbl)):
    gvecs = fp['OrbsG/kp%d_g' % ik][()]
    npw = len(gvecs)
    kvecs = 2*np.pi/alat*(gvecs+kpt)
    kvl.append(kvecs)
    cmat = []
    for ib in range(nb):
      psigl = fp['OrbsG/kp%d_b%d' % (ik, ib)][()]
      #psig = psigl[:, 0] + 1j*psigl[:, 1]
      psig = psigl.view(complex).squeeze()
      assert np.allclose(psig[npw:], 0)
      cmat.append(psig[:npw])
    cmatl.append(np.array(cmat))
  fp.close()
  return kvl, cmatl

def get_pw_basis(forb):
  from scipy.linalg import block_diag
  from afobj.basis.pw_h5 import get_orbs
  kvl, cmatl = get_orbs(forb)
  kvecs = np.concatenate(kvl, axis=0)
  cmat = block_diag(*cmatl)  # (nbas, npw)
  return kvecs, cmat

def get_h1kin(kvecs, phik):
  npw, ndim = kvecs.shape
  nbas, npw1 = phik.shape
  assert npw1 == npw
  k2 = np.einsum('ki,ki->k', kvecs, kvecs)
  h1 = np.einsum('ik,k,kj->ij', phik.conj(), k2/2, phik.T)
  return h1
