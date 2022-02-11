from setuptools import setup

###############################################################################
# TODO
# 04/10/2018 : detrend_median_seasonal: periods option to be implemented 
#        think how seasonal terms are then removed for the part not used in parametric seasonal signal 
###############################################################################

#from sphinx.setup_command import BuildDoc
#cmdclass = {'build_sphinx': BuildDoc}

name = 'pyacs'
version = '0.65'
release = '0.65.66'

setup(name=name,
      version=release,
      description='PYACS: Geodetic analysis and modeling tools for Tectonics',
      long_description='Geodetic analysis and modeling tools for Tectonics.',
#      cmdclass=cmdclass,
#      include_package_data=True,
      classifiers=[
          'Development Status :: 3 - Alpha',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3.6',
          'Topic :: Geodesy :: Geophysics',
      ],
      keywords='geodesy GPS earthquake elastic dislocation tectonics time series',
      url='',
      author='Jean-Mathieu Nocquet (Geoazur, IRD, CNRS, OCA, Cote d Azur University, France)',
      author_email='nocquet@geoazur.unice.fr',
      license='NA',
      packages=['pyacs',
                'pyacs.message',
                'pyacs.lib',
                'pyacs.lib.faultslip',
                # gts
                'pyacs.gts',
                'pyacs.gts.lib',
#                'pyacs.gts.lib.non_linear_model',
                'pyacs.gts.lib.filters',
                'pyacs.gts.lib.format',
                'pyacs.gts.lib.primitive',
                'pyacs.gts.lib.plot',
                'pyacs.gts.lib.outliers',
                'pyacs.gts.lib.model',
                'pyacs.gts.lib.tensor_ts',
                # sol & sinex
                'pyacs.sol',
                'pyacs.sinex',
                # Sgts
                'pyacs.gts.Sgts_methods',
                # refactoring 0.63.8
                'pyacs.glinalg',
                ],


      scripts=[
          # IPYACS
          'pyacs/scripts/ipyacs.py',
          # PYACS MAKE TIME SERIES
          'pyacs/scripts/pyacs_make_time_series.py',
          #                # QGIS TRANSLATERS
          'pyacs/scripts/pyacs_qgis_psvelo_2_shapefile.py',
          #                # PYACS GAMIT
          #                # PYACS GVEL POLE & STRAIN
          'pyacs/scripts/pyacs_gvel_pole.py',
          'pyacs/scripts/pyacs_gvel_estimate_pole.py',
          'pyacs/scripts/pyacs_gvel_prep_comb.py',
          'pyacs/scripts/pyacs_gvel_comb.py',
          'pyacs/scripts/pyacs_gvel_strain.py',
          #
      ],
#       install_requires=['ipython',
#                         'numpy',
#                         'scipy',
#                         'matplotlib',
#                         'argparse',
#                         'pyshp',
#                         'pyaml',
#                         'ansicolors',
#                         'pyshp>=2.0.1',
#                         'contextily',
#                         'geopandas',
#                         #'cartopy',
#                         'descartes',
#                         #'cvxpy',
#                         'sphinx_jekyll_builder',
#                         'sphinx_rtd_theme'
# #                        'l1tf @ git+https://github.com/ivannz/l1_tf.git'
#                         #                        'pwlf' , \
#                         #                        'prox_tv'
#                         ],
#       #                        'hdf5',  \
#       #                        'netCDF4'],

# strict import version
#       install_requires=['ipython==7.23.1',
#                         'numpy==1.20.3',
#                         'scipy==1.6.3',
#                         'matplotlib==3.3.4',
#                         'argparse==1.4.0',
#                         'pyshp==2.1.3',
#                         'pyaml==20.4.0',
#                         'ansicolors==1.1.8',
#                         'contextily==1.1.0',
#                         'geopandas==0.9.0',
#                         'descartes==1.1.0',
#                         'netcdf4==1.5.7'],

# relax import version
      install_requires=['ipython',
                        'numpy',
                        'scipy',
                        'matplotlib',
                        'argparse',
                        'pyshp',
                        'pyaml',
                        'ansicolors',
                        'contextily',
                        'geopandas',
                        'descartes',
                        'netcdf4'],

# old stuff
#                        'Polygon3',
#                        'cartopy==0.18.0b1',
# 'cvxpy',
#                        'sphinx_jekyll_builder==0.3.0',
#                        'sphinx_rtd_theme==0.5.2'
#                        'l1tf @ git+https://github.com/ivannz/l1_tf.git'
#                        'pwlf' , \
#                        'prox_tv'
#                        ],
      #                        'hdf5',  \
      #                        'netCDF4'],



      zip_safe=False,

      test_suite='nose.collector',
      tests_require=['nose'],

      # for building documentation using sphinx
      # these are optional and override conf.py settings
      command_options={
          'build_sphinx': {
              'project': ('setup.py', name),
              'version': ('setup.py', release),
              'release': ('setup.py', release),
              'source_dir': ('setup.py', 'documentation/source'),
              'build_dir': ('setup.py', 'documentation/build')}}
      )
