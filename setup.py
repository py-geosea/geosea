import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="geosea", # Replace with your own username
    version="1.21",
    author="Florian Petersen and Katrin Hannemann",
    author_email="flpetersen@geomar.de",
    description="A processing package for seafloor geodesy",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/flp-geo/geosea",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
