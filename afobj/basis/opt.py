import numpy as np

# ========================== level 0: setup =========================
def get_elem(qeinp):
  with open(qeinp, 'r') as f:
    lines = f.readlines()
  for i, line in enumerate(lines):
    if 'ATOMIC_SPECIES' in line:
      break
  elem = []
  for line in lines[i+1:]:
    tokens = line.split()
    if len(tokens) != 3:
      break
    name, mass, pseudo = tokens
    elem.append(name)
  assert len(elem) > 0
  return elem

def get_input(qeinp, name):
  with open(qeinp, 'r') as f:
    lines = f.readlines()
  for i, line in enumerate(lines):
    if name in line:
      out = line.split('=')[1]
      outdir = out.split("'")[1]
      return outdir

def get_params(qeout):
  with open(qeout, 'r') as f:
    lines = f.readlines()
  params = {}
  for i, line in enumerate(lines):
    if 'FFT dimensions' in line:
      ffts = line.split('(')[-1].replace(')', '')
      mesh = np.array(ffts.split(','), dtype=int)
      break
  params['mesh'] = mesh
  lines = lines[i+1:]
  for i, line in enumerate(lines):
    if 'occupation numbers' in line:
      break
  lines = lines[i+1:]
  occl = []
  for i, line in enumerate(lines):
    tokens = line.split()
    if len(tokens) < 1:
      break
    occl += list(map(float, tokens))
  nks = sum(occl)
  nks = int(round(nks))
  params['nks'] = nks
  return params

def calc_lmax(x, ml=10):
  nx = len(x)
  for lmax in range(2, ml):
    n2expect = (lmax+4)*(lmax-1)
    if n2expect == 2*nx:
      break
  if lmax >= ml-1:
    msg = 'increase ml=%d' % ml
    raise RuntimeError(msg)
  return lmax

def index_nl(n, l):
  return n*(n+5)//2+l

# ========================= level 0: gather =========================
def read_opt_log(flog):
  from qharv.reel import scalar_dat
  with open(flog, 'r') as f:
    lines = f.readlines()
  text = '\n'.join(lines[:-1])  # skip last line (cur. best)
  df = scalar_dat.parse(text)
  xcols = get_xcols(df.columns)
  sel = np.ones(len(df), dtype=bool)
  ## drop finite difference runs (from basinhopping)
  #nx = len(xcols)
  #sel = (df.icalc % (nx+1)) == 0
  return df.loc[sel].reset_index(drop=True)

def read_qe_out(qe_out):
  import pandas as pd
  from qharv.reel import ascii_out
  mm = ascii_out.read(qe_out)
  idx = ascii_out.all_lines_with_tag(mm, 'Exponents:')
  emp2s = []
  expos = []
  for i in idx:
    mm.seek(i)
    line = mm.readline().decode()
    if '\\n' in line: continue
    mm.seek(i)
    text = ascii_out.block_text(mm, 'Exponents:', 'Driver Energy')
    line = ascii_out.lr_mark(text, '[', ']')
    x = list(map(float, line.split()))
    j = mm.find(b'Driver Energy')
    mm.seek(j)
    line = mm.readline()
    emp2 = float(line.split()[2])
    emp2s.append(emp2)
    expos.append(x)
  ncalc = len(emp2s)
  assert len(expos) == ncalc
  entryl = []
  for icalc in range(ncalc):
    emp2 = emp2s[icalc]
    entry = {'icalc': icalc, 'emp2': emp2}
    x = expos[icalc]
    for ix, x1 in enumerate(x):
      entry['x%d' % ix] = x1
    entryl.append(entry)
  df = pd.DataFrame(entryl)
  return df

def get_xcols(cols):
  xcols = [col for col in cols if col.startswith('x')]
  return xcols

def get_descriptors(df, extra_cols=[]):
  xcols = get_xcols(df.columns)
  xarr = df[extra_cols+xcols].values
  return xarr

def define_pcs(xarr, ncomp=2):
  from sklearn.decomposition import PCA
  pca = PCA(n_components=ncomp)
  pca.fit(xarr)
  return pca

def get_min(df):
  xcols = get_xcols(df.columns)
  z = df['emp2'].values
  i = np.argmin(z)
  xi = df.iloc[i][xcols].values
  return z[i], xi

def default_contracted_basis(lmax, elems):
  from afqmctools.utils.optimizable_basis_set import OptimizableBasisSet
  from afqmctools.utils.optimizable_basis_set import EvenTemperedBasisBlock
  from afqmctools.utils.optimizable_basis_set import BasisBlock
  if (lmax < 2) or (lmax > 4):
    raise NotImplementedError('lmax=%d not implemented' % lmax)
  def contracted_shell(lmax, elem, expo_min=0.1, expo_max=10.0):
    bbl = []
    x0 = []  # initial guess
    xbl = []  # bounds
    for l, ang in zip(range(lmax+1), angs):
      ng = 5-lmax-l
      if ng > 1:
        bb = EvenTemperedBasisBlock(elem, ang, n_gauss=ng)
        x0 += [0.5+lmax] + [1.0]*ng
        xbl += [(expo_min, expo_max)] + [(-1, 1)]*ng
      else:
        bb = BasisBlock(elem, ang)
        x0 += [0.5+lmax]
        xbl += [(expo_min, expo_max)]
      bbl.append(bb)
    return bbl, x0, xbl
  x0 = []  # initial guess for parameters
  xbounds = []
  bset = OptimizableBasisSet(elems)
  angs = ['S', 'P', 'D', 'F', 'G', 'H']
  for elem in elems:
    for lm in range(2, lmax+1):
      bbl, x1, xb = contracted_shell(lm, elem)
      x0 += x1
      xbounds += xb
      for bb in bbl:
        bset.add(bb)
  return bset, x0, xbounds

# ========================= level 1: constraint =========================
class ExponentBounds(object):
  def __init__(self, lmax):
    self.lmax = lmax
    self.nl = 3  # first shell has 3 bases
  def __call__(self, **kwargs):
    x = kwargs['x_new']
    valid = True
    lmax = self.lmax
    for n in range(1, lmax-1):  # shell
      for l in range(self.nl):  # tigher than last shell
        i0 = index_nl(n-1, l)
        expo0 = x[i0]
        i1 = index_nl(n, l)
        expo1 = x[i1]
        if expo1 <= expo0:
          valid = False
          break
    return valid
  def make_valid(self, x):
    lmax = self.lmax
    for n in range(1, lmax-1):  # shell
      for l in range(self.nl):  # tigher than last shell
        i0 = index_nl(n-1, l)
        expo0 = x[i0]
        i1 = index_nl(n, l)
        expo1 = x[i1]
        if expo1 <= expo0:
          x[i1] = expo0
          x[i0] = expo1

# ========================= level 2: orbital =========================
def read_orbs(fh5):
  import h5py
  fp = h5py.File(fh5, 'r')
  norb = fp['OrbsG/number_of_orbitals'][()][0]
  psigl = []
  for iorb in range(norb):
    path = 'OrbsG/kp0_b%d' % iorb
    carr = fp[path][()]
    psig = carr[:, 0] + 1j*carr[:, 1]
    psigl.append(psig)
  fp.close()
  cmat = np.array(psigl)
  return cmat
