#!/bin/bash

## Script to compile a nightly build of OpenStructure from the git repo
##
## This script requires that you can clone the OST repo and it uses the sciCORE 
## software stack for dependencies so it will only work in the sciCORE login nodes
##
## Author: Pablo Escobar <pablo.escobarlopez@unibas.ch>


# print help if no arguments supplied or -h or --help
if [[ $# -eq 0 ]] || [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
	echo "Usage: $0 branch_to_clone"
	exit 1
fi	

function main {

CURRENT_DATE=$(date +"%Y%m%dT%H%M%SZ")
GIT_BRANCH=$1

# create a temporary directory to download the sources and do the build
TEMP_DIR="${TMPDIR}/${USER}/ost-builds/$CURRENT_DATE"
echo -e "\ncreating temporary dir to clone and build in $TEMP_DIR\n"
mkdir -p $TEMP_DIR

# clone the development branch of OST to $TEMP_DIR
git_cmd="git clone -b $GIT_BRANCH ssh://git@git.scicore.unibas.ch:2222/schwede/openstructure.git $TEMP_DIR/openstructure"
echo "cloning the OST code: $GIT_BRANCH branch"
echo -e "running $git_cmd"
$git_cmd

# sanity check to verify that the cloned branch exists
cd $TEMP_DIR/openstructure
git checkout $GIT_BRANCH
if [[ $? != 0 ]]; then
	echo -e "\ngit branch $GIT_BRANCH not found in cloned repository. Are you sure this is an existing branch?" 
	echo -e "Aborting.... \n" 
	exit 1
fi

# fetch the latest git commit ID. We will use this in the module name together with the branch name
cd $TEMP_DIR/openstructure
LATEST_COMMIT=`git rev-parse HEAD | cut -c1-8`

# create a tarball with the sources. We will pass this to easybuild
TARBALL_NAME="OpenStructure-nightly-${CURRENT_DATE}_${GIT_BRANCH}_commit${LATEST_COMMIT}.tar.bz2"
echo -e "\nCreating a tarball with the source code to feed EasyBuild"
echo -e "running tar -cjf $TARBALL_NAME openstructure"
cd $TEMP_DIR
tar -cjf $TARBALL_NAME openstructure

# create the installation directory if it doesn't exists
# by default $INSTALL_DIR is in the cluster shared scratch so the module
# can be used from any compute node
INSTALL_DIR="/scicore/scratch/$USER/ost-nightly-builds"
echo -e "\nCreating install dir $INSTALL_DIR if it doesn't exists"
if [ ! -d "$INSTALL_DIR" ]; then
	mkdir -p $INSTALL_DIR
fi

# create a easyconfig to build OST
echo -e "\nCreating an easyconfig in $TEMP_DIR/ost-easyconfig.eb based on the provided template"
generate_easyconfig_template  # this function is defined at the bottom of the script to improve readability

echo -e "\nSetting up the modules environment by executing 'source /etc/profile.d/lmod.sh'"
source /etc/profile.d/lmod.sh

echo -e "\nLoading the Easybuild module executing 'module load EasyBuild'"
module load EasyBuild

echo -e "\nRunning EasyBuild"
eb_command="eb --optarch=GENERIC --sourcepath=$TEMP_DIR --buildpath=$TEMP_DIR/easybuild_build_dir --installpath=$INSTALL_DIR $TEMP_DIR/ost-easyconfig.eb"
echo $eb_command
$eb_command

if [[ $? != 0 ]]; then
	echo -e "\nBuild failed. Not deleting temp dir $TEMP_DIR"
	exit 1
else
        echo -e "\nBuild completed succesfully. Deleting temp dir $TEMP_DIR"
	rm -fr $TEMP_DIR
	if [[ $? != 0 ]]; then
		echo -e "\nError deleting temp dir in $TEMP_DIR"
		exit 1
	else
		echo -e "\nrun 'module use $INSTALL_DIR/modules/all/' to add the new modules to your \$MODULEPATH"
		echo -e "If you don't see the new module try deleting your Lmod's cache doing 'rm ~/.lmod.d/.cache/*'\n"
	fi

fi

} # end of main function

function generate_easyconfig_template {

cat > $TEMP_DIR/ost-easyconfig.eb << EOF
# This file is an EasyBuild reciPY as per https://github.com/hpcugent/easybuild
# Author: Pablo Escobar Lopez
# sciCORE - University of Basel
# SIB Swiss Institute of Bioinformatics

easyblock = "CMakeMake"

name = 'OpenStructure-nightly'
version = '${CURRENT_DATE}'
versionsuffix = '_${GIT_BRANCH}_commit${LATEST_COMMIT}'

homepage = 'http://www.openstructure.org'
description = "Open-Source Computational Structural Biology Framework"

toolchain = {'name': 'goolf', 'version': '1.4.10'}

sources = ['$TARBALL_NAME']

builddependencies = [('CMake', '3.10.2')]

pythonversion = "3.6.6"
pythonshortversion = ".".join(pythonversion.split(".")[:-1])

dependencies = [
    ('Python', pythonversion),
    ('Boost', '1.53.0', '-Python-%s' % (pythonversion)),
    ('numpy', '1.9.1', '-Python-%s' % (pythonversion)),
    ('Qt', '4.8.4'),
    ('Eigen', '3.2.1'),
#   ('FFTW', '3.3.3'),  # this dependency is already included with toolchain goolf-1.4.10
    ('zlib', '1.2.8'),
    ('bzip2', '1.0.6'),
    ('libpng', '1.6.17'), 
    ('LibTIFF', '4.0.3'), 
    ('libjpeg-turbo', '1.4.0'), 
    ('OpenMM', '6.1-Linux64', '', True), 
    ('FFTW', '3.3.3', '-dynamic', ('gompi', '1.4.10')), 
]

configopts = " -DOPTIMIZE=1 -DCOMPILE_TMTOOLS=1 -DENABLE_GFX=OFF -DENABLE_GUI=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS_RELEASE='-O3 -DNDEBUG -Wno-unused-local-typedefs' -DUSE_NUMPY=1 -DUSE_RPATH=1"

configopts += " -DCOMPOUND_LIB=/scicore/home/schwede/GROUP/OpenStructure/ChemLib/1.6/compounds.chemlib"

configopts += " -DPython_ROOT_DIR=\$EBROOTPYTHON -DEIGEN3_INCLUDE_DIR=\$EBROOTEIGEN -DFFTW_LIBRARY=\$EBROOTFFTW/lib/libfftw3f.a -DFFTW_INCLUDE_DIR=\$EBROOTFFTW/include -DBOOST_ROOT=\$EBROOTBOOST -DQT_QMAKE_EXECUTABLE=\$EBROOTQT/bin/qmake -DBoost_USE_MULTITHREADED=OFF"

configopts += " -DENABLE_MM=1 -DOPEN_MM_LIBRARY=\$EBROOTOPENMM/lib/libOpenMM.so -DOPEN_MM_PLUGIN_DIR=\$EBROOTOPENMM/lib/plugins -DOPEN_MM_INCLUDE_DIR=\$EBROOTOPENMM/include"

runtest = 'check'

sanity_check_paths = { 
    'files': ["bin/ost", "bin/chemdict_tool"],
    'dirs': ["include/ost"],
}

# add \$INSTALLDIR/lib64/python3.6/site-packages to PYTHONPATH
modextrapaths = {
    'PYTHONPATH': 'lib64/python%s/site-packages' % (pythonshortversion)
}

# OST_ROOT env var points to install directory
modextravars = {
    'OST_ROOT': '%(installdir)s'
}

postinstallcmds = ["rm -fr %(installdir)s/share/openstructure/compounds.chemlib",
                   "cd %(installdir)s/share/openstructure/ && ln -s /scicore/home/schwede/GROUP/OpenStructure/ChemLib/1.6/compounds.chemlib compounds.chemlib",]

moduleclass = 'bio'
EOF

}

main "$@"
