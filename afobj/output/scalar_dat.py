import pandas as pd
class ScalarDat:
  def __init__(self):
    pass
  def parse(self,dat_fname):
    """ read the scalar.dat file, should be table format
    Args:
      dat_fname (str): name of input file
    Returns:
      df (pd.DataFrame): table of data, effect: self.df=df """
    df = pd.read_csv(dat_fname,sep='\s+')
    # remove first column name '#'
    columns = df.columns
    df.drop(columns[-1],axis=1,inplace=True)
    df.columns = columns[1:]
    # save df
    self.df = df
    return df
