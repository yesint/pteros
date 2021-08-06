import os
import pathlib

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as build_ext_orig
from shutil import copytree, copy as copyfile


class CMakeExtension(Extension):

    def __init__(self, name):
        # don't invoke the original build_ext for this special extension
        super().__init__(name, sources=[])


class build_ext(build_ext_orig):

    def run(self):
        for ext in self.extensions:
            self.build_cmake(ext)
        super().run()

    def build_cmake(self, ext):
        cwd = pathlib.Path().absolute()

        # these dirs will be created in build_py, so if you don't have
        # any python sources to bundle, the dirs will be missing        
        build_temp = pathlib.Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)

        extdir = pathlib.Path(self.get_ext_fullpath(ext.name))
        extdir.mkdir(parents=True, exist_ok=True)

        print(f'build_temp: {build_temp}')
        print(f'extdir: {extdir}')

        # example of cmake args
        config = 'Debug' if self.debug else 'Release'
        cmake_args = [
            '-DDOWNLOAD_DEPENDENCIES=ON',
            '-DTRY_SYSTEM_DEPENDENCIES=OFF',
            #'-DWITH_GROMACS=OFF',
            #'-DWITH_OPENBABEL=OFF',
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + str(extdir.parent.absolute()),
            '-DCMAKE_BUILD_TYPE=' + config,
            #'-DCMAKE_INSTALL_DIRECTORY=' + str(cwd / 'pteros')
        ]

        # example of build args
        build_args = [
            '--config', config,
            '--', '-j8'
        ]

        os.chdir(str(build_temp))
        self.spawn(['cmake', str(cwd)] + cmake_args)
        if not self.dry_run:
            #self.spawn(['cmake', '--build', '.', '--target', 'install'] + build_args)
            self.spawn(['cmake', '--build', '.'] + build_args)
        # Troubleshooting: if fail on line above then delete all possible
        # temporary CMake files including "CMakeCache.txt" in top level dir.
        os.chdir(str(cwd))


# Before running setup we need to form a correct file structure for setup.py
# We create a pteros directory and copy necessary *py files there
cwd = pathlib.Path().absolute()
package_root = cwd / 'pteros'
copytree(cwd / 'src' / 'python' / 'pteros', cwd / 'pteros', dirs_exist_ok=True)
copyfile(cwd / 'src' / 'python' / 'scripts' / 'pteros_analysis.py', cwd / 'pteros')

setup(
    name='pteros',
    version='3.0',
    author='Semen Yesylevskyy',
    author_email='yesint4@yahoo.com',
    description='Pteros molecular modeling library',
    long_description=open("./README.md", 'r').read(),
    long_description_content_type="text/markdown",
    keywords="molecular modeling, extension, molecular dyanamics",
    classifiers=["Intended Audience :: Science/Research",
                 "Intended Audience :: Information Technology",
                 "License :: OSI Approved :: Artistic License",
                 "Natural Language :: English",
                 "Programming Language :: C",
                 "Programming Language :: C++",
                 "Programming Language :: Python",
                ],
    license='Artistic License',
    packages=['pteros','pteros.extras'],
    ext_modules=[CMakeExtension('pteros')],
    cmdclass={
        'build_ext': build_ext,
    }
)
