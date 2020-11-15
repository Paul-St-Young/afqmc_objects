import os
import h5py
import asyncio
from ase.calculators.calculator import Calculator, CalculatorSetupError
from ase.calculators.calculator import CalculationFailed
from afobj.basis.qemp2 import QEMP2
from afobj.basis.gto_h5 import gen_qe_gto

class AQEMP2(QEMP2):
  async def acalc(self, atoms=None, properties=['energy']):
    system_changes = ['positions', 'numbers', 'cell', 'pbc',
      'initial_charges', 'initial_magmoms']
    Calculator.calculate(self, atoms, properties, system_changes)
    self.write_input(self.atoms, properties, system_changes)
    if self.command is None:
      raise CalculatorSetupError(
        'Please set ${} environment variable '
        .format('ASE_' + self.name.upper() + '_COMMAND') +
        'or supply the command keyword')
    command = self.command
    if 'PREFIX' in command:
      command = command.replace('PREFIX', self.prefix)

    try:
      command = 'cd %s; %s' % (self.directory, command)
      proc = await asyncio.create_subprocess_shell(command)
    except OSError as err:
      # Actually this may never happen with shell=True, since
      # probably the shell launches successfully.  But we soon want
      # to allow calling the subprocess directly, and then this
      # distinction (failed to launch vs failed to run) is useful.
      msg = 'Failed to execute "{}"'.format(command)
      raise EnvironmentError(msg) from err

    await proc.communicate()
    errorcode = proc.returncode

    if errorcode:
      path = os.path.abspath(self.directory)
      msg = ('Calculator "{}" failed with command "{}" failed in '
             '{} with error code {}'.format(self.name, command,
                                            path, errorcode))
      raise CalculationFailed(msg)

    self.read_results()

class QEGTO:

  def __init__(self, atoms, basis_set, parameters, **kwargs):
    outdir = kwargs.pop('outdir', 'qeout/') # !!!! must terminante with '/'
    if not outdir.endswith('/'):
      raise RuntimeError('"%s" does not end with "/"' % outdir)
    prefix = kwargs.pop('prefix', 'pwscf')
    out_prefix = kwargs.pop('out_prefix', 'pwscf')
    self.verbose = kwargs.pop('verbose', True)
    self.tmpdir = kwargs.pop('tmpdir', './')  # may not end in '/'
    self.atoms = atoms
    self.basis_set = basis_set
    self.iteration = 0
    self.parameters = parameters
    self.default_parameters = dict(
      prefix = prefix,
      out_prefix = out_prefix,
      outdir = outdir,
      gto_h5 = 'orbitals.h5',
    )
    self.clean_gto_h5 = kwargs.pop('clean_gto_h5', False)
    self.fband_h5 = kwargs.pop('fband_h5', 'band.h5')
    self.calc = AQEMP2()

  def __del__(self):
    self.fbh5.close()

  async def get_mp2_energy(self, x, keep_qe_io=True):
    # get parameters
    params = self.default_parameters.copy()
    params.update(self.parameters)
    label = 'i%05d' % self.iteration
    # make run directory
    path = self.tmpdir + label
    if not os.path.isdir(path):
      os.mkdir(path)

    # step 1: generate GTO orbitals
    forb = os.path.join(path, params['gto_h5'])
    nao = gen_qe_gto(self.atoms, self.basis_set, x, params['kpts'],
      mesh=params['mesh'], fname=forb)
    params['ngto'] = nao


    # step 2: get MP2 energy
    outdir = os.path.relpath(
      os.path.abspath(params['outdir']),
      os.path.abspath(path)
    )
    params['outdir'] = outdir
    self.calc.set(**params)
    self.calc.directory = path
    await self.calc.acalc(self.atoms)
    emp2 = self.calc.get_potential_energy()
    if self.verbose:
      print(self.iteration, emp2, *x, flush=True)
      self.fbh5 = h5py.File(self.fband_h5, 'a')
      evals = self.calc.results['evals']
      dset = self.fbh5.create_dataset(label, data=evals)
      self.fbh5.close()
    if self.clean_gto_h5:
      import subprocess as sp
      sp.check_call(['rm', forb])
    if not keep_qe_io:
      import subprocess as sp
      sp.check_call(['rm', '-r', path])
    self.iteration += 1
    return emp2
