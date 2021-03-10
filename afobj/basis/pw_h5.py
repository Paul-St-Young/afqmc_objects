import numpy as np

def get_orbs(forb):
  import h5py
  fp = h5py.File(forb, 'r')
  raxes = fp['OrbsG/reciprocal_vectors'][()]
  kpts = fp['OrbsG/kpoints'][()]
  nbl = fp['OrbsG/number_of_orbitals'][()]
  cmatl = []
  kvl = []
  for ik, (kpt, nb) in enumerate(zip(kpts, nbl)):
    gvecs = fp['OrbsG/kp%d_gvectors' % ik][()]
    npw = len(gvecs)
    kvecs = np.dot(gvecs+kpt, raxes)
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
