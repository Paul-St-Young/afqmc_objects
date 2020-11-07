import asyncio
from ase.calculators.calculator import Calculator, CalculatorSetupError
from afobj.basis.qemp2 import QEMP2

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
