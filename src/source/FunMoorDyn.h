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

/// \file FunMoorDyn.cpp \brief Defines the common funtions for MoorDyn+.

#ifndef _FunMoorDyn_
#define _FunMoorDyn_

#include <string>
#include <vector>
#include <sstream>

namespace funmd{
  //==============================================================================
  /// Converts numeric values to string and returns it
  //==============================================================================
  template <class T> std::string ToString(T value){
    std::string ret="";
    std::stringstream ss;
    ss<<value;
    ret=(ss.str());
    return ret;
  };  
 
  //==============================================================================
  /// Frees memory vector
  //==============================================================================
  template <class T> void FreeVector(std::vector<T *> v){
    for(unsigned i=0;i<v.size();i++){delete v[i];}
    v.clear();
  }

  //=====================================================================
  /// 2D double array destruction functions.
  //=====================================================================
  void FreeArray2D(double** v,unsigned sizex);

  //==============================================================================
  /// 3D double array destruction functions.
  //==============================================================================
  void FreeArray3D(double*** v,unsigned sizex);

  //==============================================================================
  /// Return a Pointer to Pointer od doubles
  //==============================================================================
  double **GetArray2D(const unsigned sizex,const unsigned sizey);

  //==============================================================================
  /// Return a Pointer to Pointer od doubles
  //==============================================================================
  double ***GetArray3D(const unsigned sizex,const unsigned sizey,const unsigned sizez);

  //==============================================================================
  /// Return a Pointer to Pointer of int
  //==============================================================================
  int *GetArrayInt2D(const unsigned size);

  //==============================================================================
  /// Return a Pointer to Pointer of doubles
  //==============================================================================
  double *GetArray(const unsigned size);


};
#endif // !_FunMoorDyn_
