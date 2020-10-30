def read_opt_log(flog):
  from qharv.reel import scalar_dat
  df = scalar_dat.read(flog)
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
