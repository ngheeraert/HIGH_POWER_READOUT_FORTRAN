from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='dynamics_analysis',
      version='0.1',
      description='This program loads simulation data and allows further analysis of the results.',
      long_description=readme(),
      author='Nicolas Gheeraert',
      author_email='n.gheeraert@physics.iitm.ac.in',
      license='',
      packages=['dynamics_analysis'],
      install_requires=[],
      include_package_data=True,
      zip_safe=False,
      test_suite='nose.collector',
      tests_require=['nose'],)
