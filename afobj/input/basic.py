from lxml import etree
class BasicInput(etree._ElementTree):
  def __init__(self):
    """ initialize BasicInput class """
    pass

  def show(self,node):
    """ print text representation of an xml node
    Args:
      node (lxml.etree.Element): the xml node to print
    Returns:
      None (None): effct print xml node
    """
    raw_text = etree.tostring(node,pretty_print=True)
    print( raw_text.decode('ascii') )
