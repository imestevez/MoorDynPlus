//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2016, Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/).

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics.

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License, along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>.
*/

/// \file DSphMoorDyn.cpp \brief Implements basic interface functions between MoorDyn+ and DualSPHysics.

#include "DSphMoorDyn.h"
#include "MoorDyn.h"

MoorDyn moordyn;

//==============================================================================
/// Initializes MoorDyn (returns true in case of error).
//==============================================================================
bool MoorDyn_LinesInit(const std::string filexml,const std::string nodexml
                      ,const std::string dirout,const unsigned numFt
                      ,const unsigned ftmkbound[],const tdouble3 vellin[]
                      ,const tdouble3 velang[],const tfloat3 gravity
                      ,const double tmax,const double dtout) {
  // double pos[6]={center.x,center.y,center.z,angles.x,angles.y,angles.z};
  // double vel[6]={vellin.x,vellin.y,vellin.z,velang.x,velang.y,velang.z};
  double pos[6]={ 0,0,0,0,0,0 };
  double vel[6]={ 0,0,0,0,0,0 };
  return (moordyn.LinesInit(pos,vel,filexml,nodexml,dirout.c_str(),gravity,tmax,dtout,ftmkbound,numFt) != 0);
}

//==============================================================================
/// Deallocates the variables used by MoorDyn (returns true in case of error).
//==============================================================================
bool MoorDyn_LinesClose() {
  return (moordyn.LinesClose() != 0);
}

//==============================================================================
/// Force calculation of moorings by MoorDyn (returns true in case of error).
//==============================================================================
bool MoorDyn_LinesCalc(const unsigned ftid,tdouble3 center,tdouble3 angles,tdouble3 vellin,tdouble3 velang,double t,double dt,tdouble3 &forcelin,tdouble3 &forceang) {
  angles=TDouble3(0);
  //printf("LinesCalc_%u> time:%f  dt:%f \n",id,t,dt);
  //printf("LinesCalc_%u> center:(%f,%f,%f)  angles:(%f,%f,%f) \n",id,center.x,center.y,center.z,angles.x,angles.y,angles.z);
  //printf("LinesCalc_%u> fvel..:(%f,%f,%f)  fomega:(%f,%f,%f) \n",id,vellin.x,vellin.y,vellin.z,velang.x,velang.y,velang.z);
  double pos[6]={ center.x,center.y,center.z,angles.x,angles.y,angles.z };
  double vel[6]={ vellin.x,vellin.y,vellin.z,velang.x,velang.y,velang.z };
  double flines[6]={ 0,0,0,0,0,0 };
  bool error=false;//(LinesCalc(pos,vel,flines,&t,&dt) != 0);
  forcelin=TDouble3(flines[0],flines[1],flines[2]);
  forceang=TDouble3(flines[3],flines[4],flines[5]);
  return(error);
}

//==============================================================================
/// Force calculation of moorings by MoorDyn (returns true in case of error).
//==============================================================================
bool MoorDyn_FairleadsCalc(const unsigned numFts,double*** fairpos,double*** fairvel,double*** fairforce,double t,double dt) {
  bool error=(moordyn.FairleadsCalc(numFts,fairpos,fairvel,fairforce,&t,&dt) != 0);
  return(error);
}

//==============================================================================
/// Returns the tension at the fairlead of a given line. (line=0...)
//==============================================================================
double MoorDyn_GetFairTen(unsigned line) {
  return (moordyn.GetFairTen(line));
}

//==============================================================================
/// Returns number of fairlead connections.
//==============================================================================
unsigned MoorDyn_FairsCount(const unsigned ftid) {
  return (unsigned(moordyn.GetNFairs(ftid)));
}

//==============================================================================
/// Returns number of lines.
//==============================================================================
unsigned MoorDyn_LinesCount() {
  return (unsigned(moordyn.GetNLines()));
}

//==============================================================================
/// Returns number of segments in indicated line. (line=0...)
//==============================================================================
unsigned MoorDyn_SegsCount(const unsigned line) {
  return (unsigned(moordyn.GetSegsCount(line)));
}

//==============================================================================
/// Returns number of segments in indicated line. (line=0...)
//==============================================================================
unsigned MoorDyn_SegsCount(const unsigned ftid,const unsigned line) {
  return (unsigned(moordyn.GetSegsCount(ftid,line)));
}

//==============================================================================
/// Returns position of node. (line=0... node=0...)
//==============================================================================
tdouble3 MoorDyn_GetNodePos(const unsigned line,const unsigned node) {
  double pos[3];
  if(moordyn.GetNodePos(line,node,pos) == 0)return(TDouble3(pos[0],pos[1],pos[2]));
  return(TDouble3(0));
}

//==============================================================================
/// Returns position of node. (line=0... node=0...)
//==============================================================================
tdouble3 MoorDyn_GetNodePosLink(const unsigned ftid,const unsigned line) {
  double pos[3];
  if(moordyn.GetNodePosLink(ftid,line,pos) == 0)return(TDouble3(pos[0],pos[1],pos[2]));
  return(TDouble3(0));
}

//==============================================================================
/// Returns number of Moorings created
//==============================================================================
unsigned MoorDyn_MooringsCount() {
  return(moordyn.GetNBodies());
}

//==============================================================================
/// Returns the mkbound of the Mooring
//==============================================================================
unsigned MoorDyn_GetMooringReference(const unsigned ftid) {
  return(moordyn.GetBodyReference(ftid));
}

//==============================================================================
/// Sends to MoorDyn the DualSPHysics log to store and print messages
//==============================================================================
void MoorDyn_LogInit(JLog2 * log) {
  return(moordyn.LogInit(log));
}

