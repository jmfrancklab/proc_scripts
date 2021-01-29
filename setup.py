from setuptools import setup

setup(
    name='proc_scripts',
    version='0.1',
    author="A Guinness, A A Beaton, J M Franck",
    packages=['proc_scripts',],
    license=open('LICENSE.md').read(),
    long_description=open('README.rst').read(),
    # I leave these commented out b/c they're good examples
    #package_data={'pyDiffTools':['diff-doc.js','xml2xlsx.vbs']},
    #entry_points=dict(
    #    console_scripts=["pydifft = pydifftools.command_line:main",])
)
