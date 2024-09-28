from setuptools import setup

setup(
    name='redspec',
    version='0.1',
    author='Sebastian Gomez',
    author_email='sgomez@cfa.harvard.edu',
    description='Functions to reduce single slit spectroscopic data.',
    url='https://github.com/gmzsebastian/redspec',
    license='MIT License',
    python_requires='>=3.6',
    packages=['redspec'],
    include_package_data=True,
    package_data={'redspec': ['ref_data/*']},
    install_requires=[
        'numpy',
        'matplotlib',
        'astropy',
        'scipy',
        'emcee'
    ]
)
