/*===================================================================================
<MOORDYN+> Copyright (c) 2021
Ivan Martinez Estevez and Jose M. Dominguez (Universidade de Vigo, Spain)
Matt Hall (github.com/mattEhall)

This file is part of MoorDyn+.  MoorDyn+ is free software: you can redistribute
it and/or modify it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.

Linking the MoorDyn+ library statically or dynamically with other modules is
making a combined work based on this library. Thus, the terms and conditions
of the GNU General Public License cover the whole combination. As a special
exception, the copyright holders of MoorDyn+ give you permission to dynamically
link this library with the program DualSPHysics to produce a combined model
featuring the capabilities of both DualSPHysics and MoorDyn+. This exception
is strictly limited to linking between the compiled MoorDyn+ library and
DualSPHysics. It does not extend to other programs or the use of the MoorDyn+
source code beyond the stipulations of the GPL. When the exception is used,
this paragraph must be included in the copyright notice.

MoorDyn+ is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for details.

You should have received a copy of the GNU General Public License along with
MoorDyn+. If not, see <http://www.gnu.org/licenses/>.
===================================================================================*/

/// \file FunMoorDyn.cpp \brief Implements the common funtions for MoorDyn+.

#include "FunMoorDyn.h"
#include "Functions.h"
#include <cstdlib>
#include <cstring>

namespace funmd{
  //=====================================================================
  /// 2D double array destruction functions.
  //=====================================================================
  void FreeArray2D(double** v,unsigned sizex){
    for(unsigned c=0;c<sizex;c++) delete[] v[c];
    delete[] v;
  }

  //==============================================================================
  /// 3D double array destruction functions.
  //==============================================================================
  void FreeArray3D(double*** v,unsigned sizex) {
    for(unsigned c=0;c<sizex;c++) {
      unsigned sizey=sizeof(v[c])/sizeof(double);
      for(unsigned f=0;f<sizey;f++)delete[] v[c][f];
      delete[] v[c];
    }
    delete[] v;
  }

  //==============================================================================
  /// Return a Pointer to Pointer od doubles
  //==============================================================================
  double **GetArray2D(const unsigned sizex,const unsigned sizey) {
    double **array2d;
    if(sizex) {
      try {
        array2d=new double*[sizex];
      }
      catch (const std::bad_alloc){ fun::Run_ExceptioonFun(fun::PrintStr("Could not allocate the requested memory.")); }
      for(unsigned c=0;c<sizex;c++) {
        array2d[c]=new double[sizey];
        for(unsigned e=0;e<sizey;e++) {
          array2d[c][e]=0;
        }
      }
    }
    return array2d;
  }

  //==============================================================================
  /// Return a Pointer to Pointer od doubles
  //==============================================================================
  double ***GetArray3D(const unsigned sizex,const unsigned sizey, const unsigned sizez) {
    double*** array3d;
    if(sizex) {
      try {
        array3d=new double**[sizex];
      }
      catch (const std::bad_alloc){ fun::Run_ExceptioonFun(fun::PrintStr("Could not allocate the requested memory.")); }
      for(unsigned c=0;c<sizex;c++) {  
        array3d[c]=new double*[sizey];
        for(unsigned e=0;e<sizey;e++) {
          array3d[c][e]=new double[sizez];
          for(unsigned f=0; f<sizey; f++) { array3d[c][e][f]=0; }
        }
      }
    }
    return array3d;
  }

  //==============================================================================
  /// Return a Pointer to Pointer of int
  //==============================================================================
  int *GetArrayInt2D(const unsigned size) {
    int *pToInt;
    if(size) {
      try {
        pToInt=new int[size];
      }
      catch (const std::bad_alloc){ fun::Run_ExceptioonFun(fun::PrintStr("Could not allocate the requested memory.")); }
      memset(pToInt,0,sizeof(int)*size);
    }
    return pToInt;
  }

  //==============================================================================
  /// Return a Pointer to Pointer of doubles
  //==============================================================================
  double *GetArray(const unsigned size) {
    double *array;
    if(size) {
      try {
        array=new double[size];
      }
      catch (const std::bad_alloc){ fun::Run_ExceptioonFun(fun::PrintStr("Could not allocate the requested memory.")); }
      memset(array,0,sizeof(double)*size);
    }
    return array;
  }

};