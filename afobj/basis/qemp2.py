import os
from ase.calculators.calculator import FileIOCalculator
from ase.calculators.calculator import CalculatorError, InputError

class QEMP2(FileIOCalculator):
  implemented_properties = ['energy']
  default_parameters = dict(
    prefix = 'pwscf',
    out_prefix = 'pwscf',
    eigcut = 1e-3,
    nextracut = 1e-6,
    verbose = '.true.',
  )  # automagically added to self.parameters

  def __init__(self, label='qemp2', **kwargs):
    super().__init__(label=label, **kwargs)

  def write_input(self, atoms, properties=None, system_changes=None):
    # extract needed parameters
    params = {}
    required_params = ['prefix', 'out_prefix', 'outdir', 'gto_h5',
      'nks', 'ngto', 'eigcut', 'nextracut', 'verbose']
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
  prefix = '{prefix:s}'
  out_prefix = '{out_prefix:s}'
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

  def parse_results(self, lines):
    results = {}
    # read band
    evals = []
    iline = 0
    for line in lines:
      iline += 1
      if 'Eigenvalues for k-point' in line:
        break
    for line in lines[iline:]:
      iline += 1
      if 'Starting MP2' in line:
        break
      toks = line.split()
      if len(toks) != 2:
        break
      ib, ev = list(map(float, toks))
      evals.append(ev)
    if len(evals) > 0:
      results['evals'] = evals
    # read energy
    found_emp2 = False
    for line in lines[iline:]:
      if 'EMP2 (Ha)' in line:
        found_emp2 = True
        break
    if found_emp2:
      et = line.split(':')[1]
      e1 = et.split(',')[0]
      emp2 = float(e1.replace('(', ''))
      results['energy'] = emp2
    return results

  def read_results(self):
    fout = '%s.pwo' % self.label
    with open(fout, 'r') as f:
      lines = f.readlines()
    results = self.parse_results(lines)
    self.results.update(results)
