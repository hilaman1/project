from setuptools import setup, Extension

setup(name='mykmeanssp',
      version='0.1.0',
      description="kmeans++ C-API",
      ext_modules=[Extension('mykmeanssp', sources=['spkmeansmodule.c'])])
