import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Lindel",
    version="1.0",
    author="Wei Chen",
    author_email="wchen108@uw.edu",
    description="A package for Lindel",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/shendurelab/CRISPR_NHEJ_prediction/tree/master/scripts",
    packages=['Lindel'],
    package_dir={'Lidel': 'Lindel'},
    package_data={'Lindel': ['data/*.pkl']},
    install_requires=['biopython','numpy','scipy'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,zip_safe=False)