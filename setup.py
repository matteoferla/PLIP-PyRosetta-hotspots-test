from setuptools import setup, find_packages
import os

this_directory = os.path.abspath(os.path.dirname(__file__))

if os.path.exists(os.path.join(this_directory, 'README.md')):
    with open(os.path.join(this_directory, 'README.md'), 'r') as f:
        long_description = f.read()
else:
    long_description = '''
    Scaffold hopping between bound compounds by stitching them together like a reanimated corpse.
    <img src="https://github.com/matteoferla/Fragmenstein/blob/master/images/fragmenstein.jpg?raw=true" width="300px">
    Documentation in [GitHub](https://github.com/matteoferla/Fragmenstein).
    [![colab demo](https://img.shields.io/badge/Run--demo--in--colab-colab_fragmenstein.ipynb-f9ab00?logo=googlecolab)](https://colab.research.google.com/github/matteoferla/Fragmenstein/blob/master/colab_fragmenstein.ipynb)
    ![Ox](https://upload.wikimedia.org/wikipedia/en/thumb/2/2f/University_of_Oxford.svg/132px-University_of_Oxford.svg.png)
    '''

if os.path.exists(os.path.join(this_directory, 'requirements.txt')):
    with open(os.path.join(this_directory, 'requirements.txt'), 'r') as f:
        requirements = [line.split('#')[0].strip() for line in f.readlines()]
        requirements = [line for line in requirements if line]
else:
    requirements = []

setup(
    name='plipspots_docking',
    version='0.1',
    packages=find_packages(),
    url='https://github.com/matteoferla/PLIP-PyRosetta-hotspots-test',
    license='MIT',
    author='Matteo Ferla',
    python_requires='>=3.7',  # tested on 3.9 only
    include_package_data=True,
    install_requires=requirements,
    author_email='matteo.ferla@gmail.com',
    description='PLIP-PyRosetta-hotspots-test aka Plipspots',
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[  # https://pypi.org/classifiers/
        'Development Status :: 2 - Pre-Alpha',  # Development Status :: 5 - Production/Stable
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)
