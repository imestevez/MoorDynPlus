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

#ifndef _Connection_
#define _Connection_

#include "Misc.h"
#include "JXml.h"
#include "JObject.h"
#include "Types_moordyn.h"

using namespace std;

class Environment;
class Line;

class Connection : protected JObject
{
private:
	bool Disabled;
	unsigned Node;					///< Indicates which is its position in the line
	bool HasMass;					///< Indicates if the Mass was introduced
	bool HasVol;					///< Indicates if the Volumn was introduced
	unsigned NumExec;				///< Indicates the number of executions of CalculateMassVol()
	unsigned Ref;					///-TODO: Use string or double 10.1 , 11.2 ....
	std::vector<Line*> Attached; 	///< Pointers to lines attached to this Connection Node
	std::vector<unsigned> Top; 		///< Which end of line are we attached to? 1=top/Fairlead, 0=bottom/anchor
	unsigned Number;				///< Identifier
	int Type;  						///< Defining whether fixed 0, vessel 1, or Connect 2 	
	unsigned nAttached; 			///< Mumber of attached lines
	double ConX;					///< Position at X
	double ConY;					///< Position at Y
	double ConZ;					///< Position at Z
	double ConM;					///< Initial Mass
	double ConV;					///< Initial Velocity
	double ConFX;					///< Force at X
	double ConFY;					///< Force at X
	double ConFZ;					///< Force at X
	double ConCdA; 					///< Product of drag coefficient and frontal area
	double ConCa; 					///< Added mass coefficient
	Environment * Env; 				///< Ponter to environmental settings
	//-Common properties with line internal Nodes 	
	vector< double > R; 			///< Node position [x/y/z] // double* to continuous state or input state or parameter/other-state
	vector< double > Rd;			///< Node velocity[x/y/z]  // double* to continuous state or input state or parameter/other-state
	vector< double > Q;      		///< Unit tangent vector for end Node (unused)
	vector< double > R_ves; 		///< Fairlead position for vessel Node types [x/y/z]
	vector< double > Rd_ves;		///< Fairlead velocity for vessel Node  types [x/y/z]
	vector< double > Fnet;			///< Total force on Node
	vector< double > Fnet_i;		///< 
	vector< double > RHS;			///< RHS of state-space equation (Forces divided by mass matrix)
	vector< vector< double > > S;  	///< Inverse mass matrices (3x3) for each Node
	vector< vector< double > > M; 	///< Node mass + added mass matrices
	vector< vector< double > > M_i;	///< New Mass calculated 
	double t; 						///< Simulation time
	double t0; 						///< Simulation time current integration was started at (used for BC function)
	double tlast;					///< 
	
	void ReadXml(JXml *sxml, TiXmlElement* lis, const int id_connection); /// Reads input file
	void AllocateMemoryVectors();	/// Reserves memory for vectors
	void Reset();					/// Restores attributes
	void CalculateMassVol();		/// This function calculates the Mass and Volumn of the connect Connections if the user didn't introduce value for them


public:
	Connection();												/// Constructor
	~Connection();												/// Destructor
	void Setup();												/// Makes the setup of object
	unsigned GetNumber() { return Number; };					/// Returns the number of Connection
	unsigned GetNode() { return Node; };						/// Returns the node position
	void SetNumber(const unsigned number) { Number=number; };	/// Stores the number id
	void SetNode(const unsigned node) { Node=node; };			/// Stores the node position
	unsigned GetRef() { return Ref; };							/// Returns the reference of connect
	void SetRef(unsigned ref) { Ref=ref; };						/// Stores the reference of connect
	//unsigned GetNode() { return Node; };						/// Returns the node of connect
	//void SetNode(unsigned node) { Node=node; };				/// Stores the node of connect
	int GetType() { return Type; };								/// Returns the connection Type
	void GetConnectState(vector<double> &r_out, vector<double> &rd_out);/// Writes the postion and velocities of Connection
	void GetFnet(double Fnet_out[]);							/// Writes on Fnet_out the tension
	void SetEnv(Environment * env_in) { Env=env_in; };			/// Stores a pointer to Environment settings
	void InitializeFairlead(double pX[], double TransMat[]);	/// Initialize the positions for Connection (Type Vessel)
	void InitializeConnect(double* X);							/// Initialize the positions for the Connection
	void GetNetForceAndMass();									/// Stores the Tension and Mass in their variables (Fnet & M)
	void DoRHS(const double* X, double* Xd, const double time);	/// Function for updates the States. Also includes hydrodynamic forces
	void InitiateStep(double FairIn[3], double rdFairIn[3], double time);/// Updates the boundary conditions for each coupling step at te begining
	void UpdateFairlead(const double time);						/// Updates the positions and velocities for the tine passed
	void LoadXml(JXml *sxml, const std::string &place, TiXmlElement * elec, const int id_connection);/// Load the Xml file
	void AddLineToConnect(Line * line);							/// This function handles assigning a line to a connection node
	double GetX() { return ConX; };								/// Returns the initial poistion in X
	double GetY() { return ConY; };								/// Returns the initial poistion in Y
	double GetZ() { return ConZ; };								/// Returns the initial poistion in Z
	tdouble3 GetPositions();									/// Returns a tdouble3 with the current positions [X Y Z]
	void SetPositions(tdouble3 pos);							/// Stores the new positions [X Y Z]
	bool operator==( Connection* conn ) const { return ((conn->GetNumber()== Number)&&(conn->GetType()==Type)); }
	bool GetDisabled() {return Disabled; };						/// Returns true if the line is broken, else returns false
	void SetDisabled(const bool b) { Disabled = b; };			/// Stores the new value for Disabled variable
	void ResetTension();										/// Restores the tension for the Connection when it is disabled

#ifdef USEGL
	void drawGL(void);
#endif
};

#endif
