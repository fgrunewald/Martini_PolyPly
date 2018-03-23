from setuptools import find_packages, setup

setup(
    name='PolyPly',
    version='0.1.0',
    author='Fabian Grünewald',
    author_email='f.grunewald@student.rug.nl',
    packages=find_packages(),
    include_package_data=True,
    scripts=['bin/polyply', ],
    url='https://github.com/fgrunewald/Martini_PolyPly',
    license='GPLv3',
    description='Tool for generating MARTINI Polymer itps and structures',
    long_description=open('README.md').read(),
    install_requires=['numpy', 'scipy','networkx','tqdm'],
)
