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

#include <vector>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include "JObject.h"
#include "JLog2.h"

using namespace std;

class QSlines: protected JObject
{
public:
	QSlines(JLog2 * log):Log(log){}; /// Constructor
	int Catenary(double XF, double ZF, double L, double EA, double W, double CB, double Tol,
		double* HFout, double* VFout, double* HAout, double* VAout, int Nnodes, vector<double> & s, 
		vector<double> & X, vector<double> & Z, vector<double> & Te, const unsigned numLine);/// This function return the positions and tensions of a single mooring line
private:
	static const int longwinded=0;	///< Switch to turn on excessive output for locating crashes
	JLog2 * Log;					///< Shows and stores messages
};
