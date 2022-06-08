"""A setuptools based setup module.
See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/ygidtu/scATS
Modified by Madoshakalaka@Github (dependency links added)
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name="scats",  # Required
    version="0.0.0",  # Required
    description="A sample Python project",  # Optional
    long_description=long_description,  # Optional
    long_description_content_type="text/markdown",  # Optional (see note above)
    url="https://github.com/ygidtu/scATS",  # Optional
    author="Zhang Yiming; chengl7",  # Optional
    author_email="ygidtu@gmail.com",  # Optional
    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 3 - Alpha",
        # Indicate who your project is intended for
        "Intended Audience :: Developers",
        "Topic :: Software Development :: Build Tools",
        # Pick your license as you wish
        "License :: OSI Approved :: MIT License",

        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    keywords="scRNA-seq; ATS",  # Optional
    packages=find_packages(),  # Required
    include_package_data=True,
    zip_safe=False,
    python_requires=">=3.6, <4",
    entry_points={
        'console_scripts':
            [
                'scats = cli.cli:cli'
            ]
    },
    install_requires=[
        "click==8.0.1",
        "colorama==0.4.4; python_version >= '2.7' and python_version not in '3.0, 3.1, 3.2, 3.3, 3.4'",
        "commonmark==0.9.1",
        "importlib-metadata==4.4.0; python_version < '3.8'",
        "numpy>=1.21",
        "pygments==2.9.0; python_version >= '3.5'",
        "pysam==0.16.0.1",
        "rich==10.2.2",
        "scipy==1.6.3",
        "typing-extensions==3.10.0.0; python_version < '3.8'",
        "zipp==3.4.1; python_version >= '3.6'",
    ],  # Optional
    project_urls={  # Optional
        "Bug Reports": "https://github.com/ygidtu/scATS/issues",
        "Source": "https://github.com/ygidtu/scATS/",
    },
)
