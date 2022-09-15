from setuptools import setup, Extension

setup(name="myspkmeans",
      version="1",
      description="Implementation of spkmeans",
      ext_modules=[Extension("myspkmeans", sources=["spkmeansmodule.c", "spkmeans.c"])])
