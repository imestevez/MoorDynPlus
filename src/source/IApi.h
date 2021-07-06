/*===================================================================================
<MOORDYN+> Copyright (c) 2020
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


/// \file IApi.h \brief Defines the information of the application.

#include <string>
#include <cstdlib>
#include <iostream>
#include <iomanip>      // std::setfill,std::setw -> For refill
#include <sstream>

#ifndef _IApi_
#define _IApi_
class IApi{

#define MOORDYN_VERSION "v2.023 (05-07-2021)"

public:
static std::string GetLicense() {
  std::string tx="";
  tx=tx+"\n\n<MOORDYN+>  Copyright (c) 2020";
  tx=tx+"\nIvan Martinez Estevez and Jose M. Dominguez (Universidade de Vigo, Spain)";
  tx=tx+"\nMatt Hall (github.com/mattEhall)\n";
  tx=tx+"\nThis file is part of MoorDyn+.  MoorDyn+ is free software: you can redistribute";
  tx=tx+"\nit and/or modify it under the terms of the GNU General Public License as";
  tx=tx+"\npublished by the Free Software Foundation, either version 3 of the License,";
  tx=tx+"\nor (at your option) any later version.\n";
  tx=tx+"\nLinking the MoorDyn+ library statically or dynamically with other modules is";
  tx=tx+"\nmaking a combined work based on this library. Thus, the terms and conditions";
  tx=tx+"\nof the GNU General Public License cover the whole combination. As a special";
  tx=tx+"\nexception, the copyright holders of MoorDyn+ give you permission to dynamically";
  tx=tx+"\nlink this library with the program DualSPHysics to produce a combined model";
  tx=tx+"\nfeaturing the capabilities of both DualSPHysics and MoorDyn+ .This exception";
  tx=tx+"\nis strictly limited to linking between the compiled MoorDyn+ library and";
  tx=tx+"\nDualSPHysics.It does not extend to other programs or the use of the MoorDyn+";
  tx=tx+"\nsource code beyond the stipulations of the GPL.When the exception is used,";
  tx=tx+"\nthis paragraph must be included in the copyright notice.\n";
  tx=tx+"\nMoorDyn+ is distributed in the hope that it will be useful, but WITHOUT ANY";
  tx=tx+"\nWARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS";
  tx=tx+"\nFOR A PARTICULAR PURPOSE.See the GNU General Public License for details.\n";
  tx=tx+"\nYou should have received a copy of the GNU General Public License along with ";
  tx=tx+"\nMoorDyn+.  If not, see <http://www.gnu.org/licenses/>.\n";
  tx=tx+"\n\nRunning MoorDyn+ "+MOORDYN_VERSION+" GPL v3";
  tx=tx+"\n===================================================";
  return tx;
}


//==============================================================================
/// Creates the progress bar during the dynamic relaxation calculation
//==============================================================================
static std::string BuildProgressdBar(const unsigned p=100) {
  std::stringstream ss;
  /*ss << "Progress: "  << std::setw(3) << percent << "% [";
  for(unsigned i=0; i<100; i+=2) {
    if(i <= percent)ss << '#';
    else ss << '-';
  }
  ss << ']';*/
  ss << "Computing IC generation [ Progress: "  << std::setw(3) << p << " % ]";
  return ss.str();
}

//==============================================================================
/// Prints in console the progress bar
//==============================================================================
static void ShowProgressBar(const unsigned p=100) {
  std::string bar=BuildProgressdBar(p);
  std::string endLine="\r";
  const std::string backspace(bar.size(),'\b');
  if(p >= 100)endLine="\n";
  std::cout << backspace << bar << endLine;
}

};
#endif