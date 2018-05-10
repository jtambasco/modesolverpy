from setuptools import setup

setup(name='modesolverpy',
      version='0.2.2',
      description='Photonic mode solver.',
      url='https://github.com/jtambasco/modesolverpy',
      author='Jean-Luc Tambasco',
      author_email='an.obscurity@gmail.com',
      license='MIT',
      install_requires=[
          'tqdm',
          'scipy',
          'numpy',
          'gnuplotpy',
          'opticalmaterialspy'
      ],
      packages=['modesolverpy'],
      include_package_data=True,
      zip_safe=False)
