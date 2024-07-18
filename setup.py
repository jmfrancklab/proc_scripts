from setuptools import setup

setup(
    name="pyspecProcScripts",
    version="0.1",
    author="A Guinness, A A Beaton, J M Franck",
    packages=[
        "pyspecProcScripts",
    ],
    license=open("LICENSE.md").read(),
    long_description=open("README.rst").read(),
    install_requires=[
        "pyspecdata",
        "pint",
    ],
)
