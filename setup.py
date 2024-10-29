from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(
    name='violin',
    version='1.0',

    author='Casey Hansen',
    #author_email='',
    description='VIOLIN',
    long_description='Verifying Interactions Of Likely Importance to the Network, built by the Mechanisms and Logic of Dynamics Lab at the University of Pittsburgh',
    #license='',
    keywords='dynamic system boolean logical qualitative modeling simulation',
    package_dir={'': 'src'},
    packages=['violin'],
    include_package_data=True,

    install_requires=[
        'networkx',
        'numpy',
        'pandas>=1.5.3',
        'tornado', # to not interfere with jupyter
        'httplib2',
    ],
    zip_safe=False # install as directory
    )
