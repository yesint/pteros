from setuptools import setup

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="pteros",
    version="0.0.1",
    author="Semen Yesylevskyy",
    author_email="yesint4@gmail.com",
    description="Library for molecular modeling",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yesint/pteros",
    scripts=['install_pteros.sh'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: POSIX :: Linux",
        'License :: OSI Approved :: MIT License'
    ],
    python_requires='>=3.6'
)