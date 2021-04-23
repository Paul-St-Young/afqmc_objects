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
    self.write_str = kwargs.pop('write_str', False)
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

  async def run(self, cmd):
    proc = await asyncio.create_subprocess_shell(
      cmd,
      stdout=asyncio.subprocess.PIPE,
      stderr=asyncio.subprocess.PIPE)
    stdout, stderr = await proc.communicate()
    return proc.returncode, stdout, stderr

  def wpo_command(self, path, fgto_h5):
    cmd = os.environ.get('WPO_COMMAND', None)
    if cmd is None:
      msg = 'please define WPO_COMMAND environment variable, e.g.\n'
      msg += 'cd PATH; write_pyscf_orbitals.py FORB'
      msg += '../scf.inp ../scf.out $lmax x0.dat\n'
      raise RuntimeError(cmd)
    for key in ['PATH', 'FORB']:
      if key not in cmd:
        msg = '%s must be in WPO_COMMAND' % key
        raise RuntimeError(msg)
    cmd = cmd.replace('PATH', path).replace('FORB', fgto_h5)
    return cmd

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

    ## begin 2021-02-05 write all orbitals.h5 simultaneously
    ##  call external program write_pyscf_orbitals.py in async subprocess
    #nao = gen_qe_gto(self.atoms, self.basis_set, x, params['kpts'],
    #  mesh=params['mesh'], fname=forb)
    #params['ngto'] = nao

    # execute write_pyscf_orbitals.py
    cmd = self.wpo_command(path, params['gto_h5'])
    # setup inputs for write_pyscf_orbitals.py
    import numpy as np
    fx0 = '%s/x0.dat' % path  # !!!! hard-coded filename
    if self.write_str:
      elem = np.unique(self.atoms.get_chemical_symbols())
      assert len(elem) == 1
      with open(fx0, 'w') as f:
        text = self.basis_set.basis_str(elem[0], x)
        f.write(text)
    else:
      np.savetxt(fx0, x)
    # go to calculator path
    ret = await self.run(cmd)
    iret, out, err = ret
    if iret != 0:
      msg = 'failed to execut:\n%s' % cmd
      raise RuntimeError(msg)
    nao = int(out.decode().strip('\n'))
    # end 2021-02-05 async write forb

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
      fbh5 = h5py.File(self.fband_h5, 'a')
      evals = self.calc.results['evals']
      dset = fbh5.create_dataset(label, data=evals)
      fbh5.close()
    if self.clean_gto_h5:
      import subprocess as sp
      sp.check_call(['rm', forb])
    if not keep_qe_io:
      import subprocess as sp
      sp.check_call(['rm', '-r', path])
    self.iteration += 1
    return emp2
