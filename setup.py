from setuptools import setup , find_namespace_packages

setup(name='PeriPyDIC',
      version='0.3',
      description='Peridynamics (PD) computations for state-based PD in 1D, 2D for elastic and viscoelastic materials. Also possible to import Digital Image Correlation results and compute PD forces for each pixel as a node.',
      author='Patrick Diehl, Rolland Delorme, Ilyass Tabiai',
      author_email='patrickdiehl@lsu.edu, rolland.delorme@polymtl.ca, ilyass.tabiai@polymtl.ca',
      url='https://github.com/lm2-poly/PeriPyDIC',
      keywords='material science, peridynamics, digital image correlation',
      license='GPL-3.0',
      packages=find_namespace_packages(where="./"),
      zip_safe=False)
