from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(
    name='VIOLIN',
    version='1.0',
    
    author='Casey Hansen',
    #author_email='',
    description='VIOLIN',
    long_description='Verifying Interactions Of Likely Importance to the Network, built by the Mechanisms and Logic of Dynamics Lab at the University of Pittsburgh',
    #license='',
    keywords='dynamic system boolean logical qualitative modeling simulation',

    packages=['VIOLIN','violin_tutorial'],
    include_package_data=True,

    install_requires=[
        'networkx',
        'numpy',
        'pandas',
        'tornado==4.5.3' # to not interfere with jupyter
    ],
    zip_safe=False # install as directory
    )