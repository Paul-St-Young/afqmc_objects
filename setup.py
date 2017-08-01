from distutils.core import setup

setup(
  name           = 'afqmc_objects',
  version        = '0.0',
  description    = 'Modular Python classes representing AFQMC objects.',
  author         = 'Yubo "Paul" Yang',
  author_email   = 'yyang173@illinois.edu',
  url            = 'http://publish.illinois.edu/yubo-paul-yang/',
  packages       = ['afobj','afobj.input','afobj.output'],
  install_requires = ['numpy','lxml']
)
