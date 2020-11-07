import os
from ase.calculators.calculator import FileIOCalculator
from ase.calculators.calculator import CalculatorError, InputError

class QEMP2(FileIOCalculator):
  implemented_properties = ['energy']
  default_parameters = dict(
    eigcut = 1e-3,
    nextracut = 1e-6,
    verbose = '.true.',
  )  # automagically added to self.parameters

  def __init__(self, label='qemp2', **kwargs):
    super().__init__(label=label, **kwargs)

  def write_input(self, atoms, properties=None, system_changes=None):
    # extract needed parameters
    params = {}
    required_params = ['outdir', 'nks', 'gto_h5', 'ngto',
      'eigcut', 'nextracut', 'verbose']
    missing_params = []
    for key in required_params:
      if key in self.parameters:
        params[key] = self.parameters[key]
      else:
        missing_params.append(key)
    if len(missing_params) > 0:
      msg = '\nPlease provide missing parameters:\n%s' % missing_params
      msg += '\n e.g. using calc.set(key=val)'
      raise InputError(msg)
    # write input
    fin = '%s.pwi' % self.label
    super().write_input(atoms, properties, system_changes)
    text_fmt = '''&inputpp
  outdir = '{outdir:s}'
  run_type = 'mp2_driver'
  diag_type = 'keep_occ'
  number_of_orbitals = {nks:d}
  h5_add_orbs = '{gto_h5:s}'
  read_from_h5 = {ngto:d}
  eigcut = {eigcut:e}
  nextracut = {nextracut:e}
  verbose = {verbose:s}
/
'''
    text = text_fmt.format(**params)
    with open(fin, 'w') as f:
      f.write(text)

  def read_results(self):
    fout = '%s.pwo' % self.label
    with open(fout, 'r') as f:
      lines = f.readlines()
    for line in lines:
      if 'EMP2 (Ha)' in line:
        break
    et = line.split(':')[1]
    e1 = et.split(',')[0]
    emp2 = float(e1.replace('(', ''))
    self.results['energy'] = emp2
