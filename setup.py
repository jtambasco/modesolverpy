from setuptools import setup

setup(name='modesolverpy',
      version='0.4.4',
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
          'opticalmaterialspy',
          'matplotlib',
          'future'
      ],
      packages=['modesolverpy'],
      include_package_data=True,
      zip_safe=False)
