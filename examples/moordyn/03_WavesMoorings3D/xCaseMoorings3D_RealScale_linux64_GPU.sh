#!/bin/bash 

fail () { 
 echo Execution aborted. 
 read -n1 -r -p "Press any key to continue..." key 
 exit 1 
}

# "name" and "dirout" are named according to the testcase

export name=CaseMoorings3D_RealScale
export dirout=${name}_out
export diroutdata=${dirout}/data

# "executables" are renamed and called from their directory

export dirbin=../../../bin/linux
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${dirbin}
export gencase="${dirbin}/GenCase_linux64"
export dualsphysicscpu="${dirbin}/DualSPHysics5.0CPU_linux64"
export dualsphysicsgpu="${dirbin}/DualSPHysics5.0_linux64"
export boundaryvtk="${dirbin}/BoundaryVTK_linux64"
export partvtk="${dirbin}/PartVTK_linux64"
export partvtkout="${dirbin}/PartVTKOut_linux64"
export measuretool="${dirbin}/MeasureTool_linux64"
export computeforces="${dirbin}/ComputeForces_linux64"
export isosurface="${dirbin}/IsoSurface_linux64"
export flowtool="${dirbin}/FlowTool_linux64"
export floatinginfo="${dirbin}/FloatingInfo_linux64"

option=-1
 if [ -e $dirout ]; then
 while [ "$option" != 1 -a "$option" != 2 -a "$option" != 3 ] 
 do 

	echo -e "The folder "${dirout}" already exists. Choose an option.
  [1]- Delete it and continue.
  [2]- Execute post-processing.
  [3]- Abort and exit.
"
 read -n 1 option 
 done 
  else 
   option=1 
fi 

if [ $option -eq 1 ]; then
# "dirout" to store results is removed if it already exists
if [ -e ${dirout} ]; then rm -r ${dirout}; fi

# CODES are executed according the selected parameters of execution in this testcase
${gencase} ${name}_Def ${dirout}/${name} -save:all 
if [ $? -ne 0 ] ; then fail; fi

${dualsphysicsgpu} -gpu ${dirout}/${name} ${dirout} -dirdataout data -svres -stable
if [ $? -ne 0 ] ; then fail; fi

fi

if [ $option -eq 2 -o $option -eq 1 ]; then
export dirout2=${dirout}/floatinginfo
${floatinginfo} -dirin ${diroutdata} -onlymk:60 -savemotion -savedata ${dirout2}/FloatingMotion 
if [ $? -ne 0 ] ; then fail; fi

export dirout2=${dirout}/particles
${partvtk} -dirin ${diroutdata} -savevtk ${dirout2}/PartFluid -onlytype:-all,+fluid 
if [ $? -ne 0 ] ; then fail; fi
${partvtk} -dirin ${diroutdata} -savevtk ${dirout2}/PartMoving -onlytype:-all,+moving 
if [ $? -ne 0 ] ; then fail; fi
${partvtk} -dirin ${diroutdata} -savevtk ${dirout2}/PartFloating -onlytype:-all,+floating 
if [ $? -ne 0 ] ; then fail; fi

export dirout2=${dirout}/boundary
${boundaryvtk} -loadvtk ${dirout}/${name}_Tank_Dp.vtk -motiondata ${diroutdata} -savevtkdata ${dirout2}/Tank.vtk -onlymk:10
if [ $? -ne 0 ] ; then fail; fi
${boundaryvtk} -loadvtk ${dirout}/${name}_Floater_Dp.vtk -motiondata ${diroutdata} -savevtkdata ${dirout2}/MotionFloating -onlymk:60
if [ $? -ne 0 ] ; then fail; fi
${boundaryvtk} -loadvtk ${dirout}/${name}_Piston_Dp.vtk -motiondata ${diroutdata} -savevtkdata ${dirout2}/MotionPiston -onlymk:20
if [ $? -ne 0 ] ; then fail; fi

export dirout2=${dirout}/surface
${isosurface} -dirin ${diroutdata} -saveiso ${dirout2}/Surface 
if [ $? -ne 0 ] ; then fail; fi

fi
if [ $option != 3 ];then
 echo All done
 else
 echo Execution aborted
fi

read -n1 -r -p "Press any key to continue..." key
