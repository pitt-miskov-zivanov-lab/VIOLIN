from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(
    name='violin',
    version='1.1',

    author='Casey Hansen and Haomiao Luo',
    description='VIOLIN',
    long_description='Verifying Interactions Of Likely Importance to the Network, built by the Mechanisms and Logic of Dynamics Lab at the University of Pittsburgh',
    keywords='dynamic system boolean logical qualitative modeling simulation',
    package_dir={'': 'src'},
    packages=['violin'],
    include_package_data=True,

    install_requires=[
        'openpyxl',
        'networkx',
        'numpy',
        'pandas',
        'httplib2',
        "matplotlib",
        "requests",
        "tqdm",
        "nltk",
    ],
    zip_safe=False # install as directory
    )
