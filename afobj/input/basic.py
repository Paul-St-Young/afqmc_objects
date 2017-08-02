from lxml import etree
class BasicInput(etree._ElementTree):
  def __init__(self):
    """ initialize BasicInput class, inherent all methods from lxml.etree.ElementTree """
    pass
 
  def read(self,fname):
    parser = etree.XMLParser(remove_blank_text=True)
    return self.parse(fname,parser)
  def write(self,fname):
    super(BasicInput,self).write(fname,pretty_print=True)

  def show(self,node):
    """ print text representation of an xml node
    Args:
      node (lxml.etree.Element): the xml node to print
    Returns:
      None (None): effct print xml node
    """
    raw_text = etree.tostring(node,pretty_print=True)
    print( raw_text.decode('ascii') )

  def get_xml_node(self,group_name,**kwargs):
    """ build <group_name> 
    Args:
      group_name (str): xml node tag, must be in "allowed_groups"
    Returns:
      node (etree.Element): <group_name> node
    allow inputs:
     AFQMCInfo
      nmo (int): number of molecular orbitals
      naea (int): number of valence up electrons 
      naeb (int): number of valence down electrons 
      nca (int): number of core up electrons 
      ncb (int): number of core down electrons 
    """
    allowed_groups = set(['project','AFQMCInfo','Hamiltonian','Wavefunction',
      'WalkerSet','Propagator','execute'])
    allowed_inputs_map = {
      'project':set(['id','series']),
      'AFQMCInfo':set(['nmo','naea','naeb','nca','ncb']),
      'Hamiltonian':set(['filetype','filename','cutoff_1bar','cutoff_2bar',
         'cutoff_decomposition','hdf_write_file']),
      'WalkerSet':set(['min_weight','max_weight','reset_weight','extra_spaces']),
      'Propagator':set(['cutoff_propg','hdf_write_file','parallel_factorization']),
      'execute':set(['timestep','blocks','steps','substeps','nWalkers'])
    }

    # initialize xml node
    if group_name not in allowed_groups:
      raise RuntimeError('unknwon group name "%s" to AFQMC'%group_name)
    node = etree.Element(group_name)
    attrib = kwargs.pop('attrib',None)
    if attrib is not None:
      node.attrib.update(attrib)

    # for a simple node, fill parameters
    allowed_inputs = allowed_inputs_map[group_name]
    for name,val in kwargs.items():
      if name not in allowed_inputs:
        raise RuntimeError('unknown input "%s" to %s'%(name,group_name))
      entry = etree.Element('parameter',{'name':name})
      entry.text = str(val)
      node.append(entry)

    return node

  def example_n2_vdz(self):
    root = etree.Element('simulation',{'method':'afqmc'})

    proj_node = self.get_xml_node('project',attrib={'id':'N2','series':'0'})
    info_node = self.get_xml_node('AFQMCInfo',nmo=28,naea=14,naeb=14,nca=0,ncb=0,attrib={'name':'info0'})
    ham_node  = self.get_xml_node('Hamiltonian',filetype='fcidump',filename='FCIDUMP',cutoff_1bar=1e-6,cutoff_2bar=1e-6,cutoff_decomposition=1e-5,hdf_write_file='ham.h5',attrib={'name':'ham0','type':'SparseGeneral','info':'info0'})

    # build <Wavefunction>
    wf_node  = etree.Element('WaveFunction',{'name':'wfn0','info':'info0'})
    det_node = etree.Element('ImpSamp',{'name':'impsamp0','type':'MultiPureSD'})
    name_val_map = {
      'filetype':'ascii',
      'filename':'ifort.100',
      'diagHam':'no',
      'cutoff':1e-6,
      'hdf_write_file':'wfn.h5'
    }
    for name,val in name_val_map.items():
      param_node = etree.Element('parameter',{'name':name})
      param_node.text = str(val)
      det_node.append(param_node)
    wf_node.append(det_node)

    wset_node = self.get_xml_node('WalkerSet',min_weight=0.05,max_weight=4,reset_weight=1,extra_spaces=10,attrib={'name':'wset0','type':'distributed'})
    prop_node = self.get_xml_node('Propagator',cutoff_propg=1e-6,hdf_write_file='prop.h5',parallel_factorization='yes',attrib={'name':'prop0','phaseless':'yes','localenergy':'yes','drift':'yes','info':'info0'})
    exe_node  = self.get_xml_node('execute',timestep=0.005,blocks=100,steps=4,substeps=4,nWalkers=100,attrib={'wset':'wset0','ham':'ham0','wfn':'wfn0','prop':'prop0','info':'info0'})

    for child in [proj_node,info_node,ham_node,wf_node,wset_node,prop_node,exe_node]:
      root.append(child)
    return root

  def hdf5_input(self,proj_id,nmo,naea,naeb,fcidump_fname):
    template = self.example_n2_vdz()

    proj_node = template.find('.//project')
    proj_node.set('id',proj_id)

    info_node = template.find('.//AFQMCInfo')
    for name,val in zip(['nmo','naea','naeb'],[nmo,naea,naeb]):
      param_node = info_node.find('.//parameter[@name="%s"]'%name)
      param_node.text = str(val)

    ham_node = template.find('.//Hamiltonian')
    for name,val in zip(['filetype','filename'],['hdf5',fcidump_fname]):
      param_node = ham_node.find('.//parameter[@name="%s"]'%name)
      param_node.text = str(val)

    return template
