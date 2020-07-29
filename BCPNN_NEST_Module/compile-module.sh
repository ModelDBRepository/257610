#!/bin/sh
#Input: Take module-DATE as input


# on the cluster:
# compile with the GNU compiler
#module swap PrgEnv-cray PrgEnv-gnu

#For the correct directory
CURR=$(pwd)
echo "current dir: $CURR"

#Start time watch
START=$(date +%s)

#Get number of processors on the system
#noProcs=1
noProcs=$(grep -c 'model name' /proc/cpuinfo) 

#Source directory
srcDir="$CURR/$1"

#Bootstrap directory, used temoprally. Removed at end of script.
bootstrapDir="$CURR/bootstrap-$1"

#Build directory
buildDir="$CURR/build-$1"

echo "Source dir: $srcDir"
echo $bootstrapDir
echo "Build dir: $buildDir"

#Copy source to bootstrap directory
cp -r $srcDir $bootstrapDir

echo $(pwd)
#Go into bootstrap dir and run bootstrap
cd $bootstrapDir
"$bootstrapDir/bootstrap.sh"

#Move out
cd ..

#Remove old build
echo "Removing previous build directory" 
rm -r "$CURR/build*"

#Make new build directory, configure and run make, make install and make installcheck
mkdir $buildDir
echo "Entering $buildDir"
cd $buildDir

export NEST_INSTALL_DIR=/home/fiebig/Programs/nest-2.4.2-build
#export NEST_INSTALL_DIR=/home/bernhard/Downloads/nest-2.4.2-build/
#export NEST_INSTALL_DIR=/pdc/vol/nest/2.2.2/
#$bootstrapDir"configure" --with-nest=${NEST_INSTALL_DIR}/bin/nest-config --prefix=${NEST_INSTALL_DIR}/ 2>&1 | tee "$CURR/mymodule-configure.log"
$bootstrapDir"configure" --with-nest=${NEST_INSTALL_DIR}/bin/nest-config --prefix=$CURR/ 2>&1 | tee "$CURR/mymodule-configure.log"
make -j $noProcs 2>&1 | tee "$CURR/mymodule-make.log"
make -j $noProcs install  

#Stop time watch
END=$(date +%s)
DIFF=$(( $END - $START ))

# Move out
cd ..

#ln -s /home/bernhard/workspace/BCPNN-Module/share/nest/sli/pt_module.sli share/nest/sli/pt_module.sli
ln -s $CURR/share/nest/sli/pt_module.sli $NEST_INSTALL_DIR/share/nest/sli/pt_module.sli

#sudo ln -s $buildDir/pt_module /usr/bin/pt_module


#Display script execution time
echo "It took $DIFF seconds"
