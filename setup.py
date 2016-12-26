from setuptools import setup

setup(name='modesolverpy',
      version='0.1',
      description='Photonic mode solver with a simple interface.',
      url='https://github.com/jtambasco/modesolverpy',
      author='Jean-Luc Tambasco',
      author_email='an.obscurity@gmail.com',
      license='MIT',
      install_requires=['tqdm', 'scipy', 'numpy'],
      packages=['modesolverpy'],
      zip_safe=False)
