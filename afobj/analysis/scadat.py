import os
import numpy as np

def calc_neq(finp, teq, min_block=8):
  from qharv.seed import xml
  doc = xml.read(finp)
  ex = doc.find('.//execute')
  s = ex.find('.//parameter[@name="steps"]')
  t = ex.find('.//parameter[@name="timestep"]')
  ts = float(t.text)
  nstep = int(s.text)
  # substep
  ss = ex.find('.//parameter[@name="substeps"]')
  nsub = 1
  if ss is not None:
    nsub = int(ss.text)
  # block time
  tblock = ts*nstep*nsub
  neq = int(np.ceil(teq/tblock))
  b = ex.find('.//parameter[@name="blocks"]')
  nblock = int(b.text)
  if neq >= nblock-min_block:
    msg = 'refusing to throw out %d/%d blocks' % (neq, nblock)
    raise RuntimeError(msg)
  return neq

def get_mdf(floc, nequil):
  from qharv.reel import mole, scalar_dat
  from qharv.sieve import mean_df
  df0 = scalar_dat.read(floc)
  mdf = mean_df.create(df0.iloc[nequil:])
  # add metadata
  path = os.path.dirname(floc)
  fdat = os.path.basename(floc)
  mdf['path'] = mole.clean_path(path)
  mdf['fdat'] = fdat
  return mdf
