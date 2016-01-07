from distutils.core import setup

version = '0.1'
setup(
    name='aakbar',
    version=version,
    url='http://github.com/GenerisBio/aakbar',
    download_url='http://github.com/GenerisBio/aakbar/tarball/'+version,
    keywords=['science', 'biology', 'bioinformatics', 'phylogenomics', 'peptide', 'signatures'],
    platforms=['Linux', 'Mac OSX', 'Windows', 'Unix'],
    license='GenerisBio',
    description='Amino-Acid k-mer phylogenetic signature tools',
    long_description='Creating, searching, and analyzing k-mer signatures in space',
    author='Joel Berendzen',
    author_email='joel@generisbio.com',
    packages=['aakbar'],
    zip_safe=True,
    install_requires=['click'],
    entry_points={
                 'console_scripts': [
                                     'aakbar = aakbar:aakbar'
                                    ]
                 },
    classifiers=[
                        'Development Status :: 3 - Alpha',
                        'Environment :: Console',
                        'Environment :: MacOS X',
                        'Environment :: Win32 (MS Windows)',
                        'Intended Audience :: Science/Research',
                        'License :: Other/Proprietary License ',
                        'Operating System :: OS Independent',
                        'Programming Language :: Python',
                        'Topic :: Scientific/Engineering :: Bio-Informatics',
                        ]
)