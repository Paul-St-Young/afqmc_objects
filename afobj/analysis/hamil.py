def get_real(e1t):
  return float(e1t.split(',')[0].split('(')[1])

def read_ham(fham):
  from qharv.reel import ascii_out
  mm = ascii_out.read(fham)
  idx = mm.find(b'E0')
  mm.seek(idx)
  line = mm.readline().decode()
  toks = line.split(':')[1].split()
  e0t = toks[0]
  e1t = toks[1]
  e0 = float(e0t)
  e1 = get_real(e1t)
  # read Hartree energy
  idx = mm.find(b'EJ(1Det)')
  mm.seek(idx)
  line = mm.readline().decode()
  ejt = line.split(':')[1]
  ej = get_real(ejt)
  # read exchange energy
  idx = mm.find(b'EXX(1Det)')
  mm.seek(idx)
  line = mm.readline().decode()
  ext = line.split(':')[1]
  ex = get_real(ext)
  mm.close()
  entry = {'e0': e0, 'e1': e1, 'ej': ej, 'ex': ex}
  return entry
