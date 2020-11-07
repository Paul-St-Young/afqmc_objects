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

# ========================= level 0: gather =========================
def read_opt_log(flog):
  from qharv.reel import scalar_dat
  with open(flog, 'r') as f:
    lines = f.readlines()
  text = '\n'.join(lines[:-1])  # skip last line (cur. best)
  df = scalar_dat.parse(text)
  xcols = get_xcols(df.columns)
  # drop finite difference runs (from basinhopping)
  nx = len(xcols)
  sel = (df.icalc % (nx+1)) == 0
  return df.loc[sel].reset_index(drop=True)

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

def et_to_x(ab):
  nab = len(ab)//2
  assert len(ab) == 2*nab
  xl = []
  for i in range(nab):
    a = ab[2*i]
    b = ab[2*i+1]
    nexpo = i+3
    xl += [a*np.exp(-b*i) for i in range(nexpo)]
  x = np.array(xl)
  return x

def x_to_et(x, ne, ml=10):
  nx = len(x)
  # find lmax
  for lm in range(2, ml+1):
    nx1 = (lm+4)*(lm-1)//2
    if nx1 == nx//ne:
      break
  if lm >= ml:
    msg = 'failed to find lmax'
    raise RuntimeError(msg)
  ab = []
  for ie in range(ne):
    for l in range(2, lm+1):
      i0 = ie*nx1+(l-2)*l
      x0 = x[i0]
      x1 = x[i0+1]
      a = x0
      b = np.log(a/x1)
      ab += [a, b]
  return np.array(ab)
