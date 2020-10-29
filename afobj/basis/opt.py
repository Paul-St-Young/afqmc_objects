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
