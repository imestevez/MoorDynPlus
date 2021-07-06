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

/// \file Environment.h \brief Defines the class \ref Environment.

#ifndef _Environment_
#define _Environment_

#include "JXml.h"
#include "TypesDef.h"
#include "JObject.h"
#include "JLog2.h"
#include <algorithm>
#include <cmath>
#include <iostream>

class Environment : protected JObject
{
private:
	JLog2 * Log;
	double TimeMax; ///< Time of simulation
	double G; ///< Gravity (m/s^2)
	double WtrDpth; ///< Depth of the water (m)
	double Rho_w; ///< Density
	double Kb; ///< Bottom stiffness (Pa/m)
	double Cb; ///< Bottom damping   (Pa/m/s)
	int WaveKin; ///< Wave kinematics flag (0=off, gt than 0=on)
	int WriteUnits; ///< Global switch for whether to show the units line in the output files (1, default), or skip it (0)
	double FrictionCoefficient; ///< General bottom friction coefficient, as a start
	double FricDamp; ///< Damping coefficient used to model the friction at speeds near zero
	double StatDynFricScale; ///< Ratio of static to dynamic friction ( = mu_static/mu_dynamic)
	double ICDfac; ///< Factor by which to boost drag coefficients during dynamic relaxation IC generation
	double ICdt; ///< Convergence analysis time step for IC generation
	double ICTmax; ///< Max time for IC generation
	double ICthresh; ///< Threshold for relative change in tensions to call it converged
	double DtM0; ///< Default value for desired mooring model time step
	double FreeSurface; ///< Z position of water free surface	
	void ReadXml(JXml *sxml, TiXmlElement* lis); /// Reads the Xml file

public:
	Environment(JLog2 * log); /// Constructor
	~Environment(); /// Destructor
	void Reset(); /// Restore attributes
	void LoadXml(JXml *sxml, const std::string &place);	/// Loads Xml file
	double GetTimeMax() { return TimeMax; }; /// Time of simulation
	double GetG() { return G; }; /// Returns the Gravity		
	double * GetWtrDpth() { return &WtrDpth; }; /// Returns the depth of the water
	double GetRho_w() { return Rho_w; }; /// Returns the density	
	double GetKb() { return Kb; }; /// Returns the bottom stiffnes  		
	double GetCb(){ return Cb;}; /// Returns the botom damping     			
	int GetWaveKin(){ return WaveKin;}; /// Returns the waves kinematics flag	
	int GetWriteUnits(){ return WriteUnits;}; /// Returns the Global switch for whether to show the units		
	double GetFrictionCoefficient() { return FrictionCoefficient; }; /// Returns the General bottom friction coefficient
	double GetFricDamp(){ return FricDamp;}; /// Returns the Damping coefficient 
	double GetStatDynFricScale(){ return StatDynFricScale;}; /// Returns the Ratio of static to dynamic friction 
	double GetICDfac(){ return ICDfac;}; /// Returns the factor by which to boost drag coefficients during dynamic relaxation IC generation
	double GetICdt(){ return ICdt;}; /// Returns the Convergence analysis time step for IC generation
	double GetICTmax(){ return ICTmax;}; /// Returns the Max time for IC generation
	double GetICthresh(){ return ICthresh;}; /// Returns the Threshold for relative change in tensions to call it converged
	double GetDtM0(){ return DtM0;}; /// Returns the Default value for desired mooring model time step	
	void SetDtM0(const double dtm) { DtM0=dtm; }; /// Sets a new value of Dtm0
  void SetG(const tfloat3 g) {G=std::abs(g.z);}; /// Sets the Gravity		
  void SetTimeMax(const double t) {TimeMax=t;}; /// Sets Time of simulation	
	double GetFreeSurface(){ return FreeSurface;}; /// Returns the Z position of water free surface	
	void SetFreeSurface(const double fs) { FreeSurface=fs; }; /// Sets the Z position of water free surface
  void VisuConfig()const; /// Shows object configuration using Log.
};


#endif