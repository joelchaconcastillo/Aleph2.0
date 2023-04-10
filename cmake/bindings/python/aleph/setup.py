from setuptools import setup

setup(
    name='Aleph',
    version='',  # TODO: might want to use commit ID here
    packages=['aleph'],
    package_dir={
        '': '/home/jcc/Aleph2.0/cmake/bindings/python'
    },
    package_data={
        '': ['aleph.so']
    },
    zip_safe=False
)
