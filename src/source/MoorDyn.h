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

/* 
* ==================================================================================================
* ================================= MoorDyn+ versions ==============================================
* ==================================================================================================
*
* ==== v1.001 ====
* - C++ reimplementation.
* ==== v1.002 ====
* - New Xml.
* ==== v1.003 ====
* - <connect> connections.
* ==== v1.004 ====
* - Calculation Mass and Volumn of the connect connections.
* - Updated Moorings.cpp and JMooredFloatings to allow more than one moored floating.
* - Two Fairleads for each line.
* ==== v1.005 ====
* - Output files will not be created when user doesn't choose save data.
* ==== v1.006 ====
* - Moored floatings to different water depth.
* ==== v1.007 ====
* - Log system to print messages by console and store it in Run.out.
* - Warning System to show it by console.
* ==== v1.008 ====
* - Class QSlines. It allows Log and warning system.
* ==== v1.009 ====
* - Added Progress bar
* ==== v1.010 ====
* - Allowed different set different water free surface
*
* ==================================================================================================
* ==================================================================================================
* ==================================================================================================
*/
#ifndef _Moordyn_
#define _Moordyn_

#include "TypesDef.h"
#include "Functions.h"
#include "JXml.h"
#include "JObject.h"
#include "Misc.h"
#include "Environment.h"
#include "Body.h"
#include "Output.h"
#include "JLog2.h"

#include <stdio.h>
#include <string>


 class MoorDyn : protected JObject
{
private:
	std::string FileIN; 		///< Name of input file.
	std::string Dirout;			///< Default Output directory (Is replaced when an external software calls to LinesInit() )

	bool Coupling;							///< Used to know if MoorDyn is coupling with external software
	std::string RootNodeXml;				///< Root element of Xml for read MoorDyn configutation
	Environment * Env;						///< Environment Object
	LineDefault * LineDef;					///< Pointer to default parameters of the lines
	std::vector<Body *> Bodies;				///< Array of pointers to Bodies
	std::vector<Line *> Lines;				///< Array of pointers to Lines
	std::vector<Line *> OutLines;			///< Array of pointers to Lines which will be stored into the csv file
	std::vector<Connection *> Connections;	///< Array of pointers to Connections
	std::vector<unsigned> References;		///< Stores all references of Bodies and connects
	
	unsigned NBodies; 						///< Number of Bodies
	unsigned NLines; 						///< Number of Lines
	unsigned NConnections; 					///< Number of Connection objects
	unsigned NFairs; 						///< Number of fairlead connections
	unsigned NAnchs; 						///< Number of anchor connections
	unsigned NConns; 						///< Number of "connect" connections
	bool XmlReaded;							///< Used for know if Xml already was readed
	bool HeaderProg;						///< Used to know if header progres was print
	bool FairArrays;						///< Indicate if the fairlead arrays was initializated.
	bool ItArrays;							///< Indicate if the integration arrays to store the stats were initialised.
	bool DSPH_Log;							///< True if DualSPHysics sends a JLog2 object, else false
	Output * OutP; 							///< Pointer to Properties of the outpus
	JLog2 *Log;								///< Shows and stores messages

	double*** FairleadPos;					///<Stores link positions.  
	double*** FairleadVel;					///<Stores link velocities. 
	double*** FairleadForce;				///<Stores link forces.     

	//-Static vectors for fairleads
	vector<double> FlinesS;					///< Net line force vector (6-DOF) - retains last solution for use when inputted dt is zero (such as during FAST predictor steps) when the model does not time step
	double** rFairtS;						///< Fairlead locations in turbine/platform coordinates
	vector< vector< double > > rFairRel;	///< Fairlead locations relative to platform ref point but in inertial orientation
	double** rFairi;						///< Fairlead locations in inertial reference frame
	double** rdFairi;						///< Fairlead velocities in inertial reference frame
	vector< int > FairIs;  					///< Vector of fairlead connection indices in ConnectList vector
	vector< int > ConnIs;  					///< Vector of connect connection indices in ConnectList vector
	int NStateVariables; 					///< Used for add six state variables for each "connect"-type Connection 
	double** Ffair;							///< Pointer to 2-d array holding fairlead forces
	unsigned NX; 							///< Size of state vector array
	double* States;							///< Pointer to array comprising global state vector
	//-State vectors things for Rk2/rk4 integration 
	double* Xt;								///< Temporal states of positions for RHSmaster
	double* F0;								///< Temporal state of velocities for RHSmaster
	double* F1;								///< Temporal state of velocities for RHSmaster
	//double* F2;
	//double* F3;	
	double TransMat[9];						///< Calculate TransMat 

	void Reset();													/// Restores attributes
	template <class T>	void FreeVector(std::vector<T *> v);		/// Frees memory vector	
	void LoadXml(const std::string FileXml); 						/// Loads the Xml file
	void ReadXml(JXml *xml); 										/// Reads the Xml file
	void CheckDepth(Environment * env, Body * body);				/// Checks if depth was introduced bay user in body section. If not, uses waterDepth SolverOptions
	void FreeArray3D(double*** v, unsigned sizex);
	void AllocateMemoryIntegration(const unsigned size); 			/// Reserves memory for vectors used for Integration
	double ** GetPtrToPtr(const unsigned  sizex, const unsigned sizey); /// Allocates memory for an Pointer to Pointer for doubles an Returns it [used for temporaly variables]
	double *** GetPtrToPtrToPtr(const unsigned  sizex, const unsigned sizey, const unsigned sizez); /// Allocates memory for an Pointer to Pointer to Pointer for doubles an Returns it [used for temporaly variables]
	int * GetPtrToInt(const unsigned size); 						/// Allocates memory for an Pointer to int and returns it [used for temporaly variables]
	double * GetPtrToDouble(const unsigned size); 					/// Allocates memory for an Pointer to int and returns it [used for temporaly variables]
	void CreateOutDirectory(); 										/// Creates the ouptuf directory when there isn't a coupling with external software
	void PrepareState(); 											/// Inits the state for calculate
	void StoreConnPositions();										/// Stores the connetions in a array with their positions
	void CheckReferences(const unsigned ref); 						/// Checks if exists references equals
	void CheckConnReferences(const unsigned ref, std::vector<Connection *> conns); /// Checks if exists references equals
	void Setup(); 													/// Initialize objects
	void InitializeSystem(double X[], double XD[]);					/// Initialize system, including trying catenary IC gen of Lines
	void DoDynamicRelaxIC();										/// Finalizing ICs using dynamic relaxation
	void AllOutput(double t, double dtC);							/// Calls to writer function of the Body passed
	void UpdatePos(const double divisor);							/// Funtion for update the positions of simulated floating
	void PrintProgress(double time);								/// Prints the execution progress in the console
	void PrintHeader();												/// Prints the header for the execution progress
	double PrintPercent(double time);								/// Prints the completed perecent of the simulation
	void ConfigBodies(); 											/// Searches Nulls and reorder the Bodies by Ref
	void ConfigConnections(std::vector<Connection *> connects); 	/// Configures the Connections
	void CheckFtIds(const unsigned ids[], const unsigned numFt); 	/// Checks if floating exists in MoorDyn the floating ids passed by parameters
	bool CheckBody(const unsigned ftid);							/// Checks if floating id and mooring exists. Returns true in error case
	bool CheckLine(const unsigned numLine);							/// Checks if line id exists. Returns true in error case
	bool CheckLine(const unsigned ftid, const unsigned numLine);	/// Checks if line id exists. Returns true in error case
	std::vector<Line *> GetLines(){return Lines;}; 					/// Returns a pointer to all lines
	std::vector<Line *> GetLines(Body * body); 						/// Returns a pointer to lines of a body
	std::vector<Line *> GetLines(const unsigned ref); 				/// Returns a pointer to lines of a body
	std::vector<Line *> GetLinesNoBody(); 					    	/// Return a pointer to array of lines without body
	Body * GetBody(const unsigned ref); 							/// Returns a pointer to Body with reference passed
	Body * GetBody(Line * line); 									/// Returns a pointer to Body of the line passed
	std::vector<Connection *> GetConnections(){return Connections;};/// Returns a pointer to pointer of all connections
	std::vector<Connection *> GetConnections(const unsigned type);	/// Returns a pointer to pointer of Connections of a type [0=anchor, 1=vessel, 2=connects]
	void AddAttachedLines();										/// Searches the connected lines with a connect Connection
	void RHSmaster(const double X[], double Xd[], const double t, const double dt);	/// Master function to handle time stepping (updated in v1.0.1 to follow MoorDyn F)
	void Rk2(double *t0, const double dt);							/// Runge-Kutta 2 integration routine  (integrates States and time)
	void CheckBodyRefs();											/// Checks if there are Fair conection with Body refence which doesn't exist
	void RemoveLine(Line* line);									/// Removes the broken line from the system
	void PrintLicense() {Log->Print(GetLicense());}					/// Shows by console the License and MoorDyn+ version

public:
	MoorDyn();																/// Constructor
	~MoorDyn();																/// Destructor
	void Run(); 															/// Executes the funtions for to run the program. Used when there isn't coupling with external sowtware
	int LinesInit(double X[], double XD[], const std::string filexml, const std::string nodexml, const char *dir,const tfloat3 gravity, const unsigned ftmkbound[]=NULL, const unsigned numFts=0);///Funtion to Init the system when there is a coupling with externan software
	int FairleadsCalc(const unsigned numFts, double*** rFairIn, double*** rdFairIn, double*** fFairIn, double* t_in, double *dt_in);
	unsigned GetNFairs(const unsigned ftid);  								///Returns the number of fairleads (Vessel connections) of a Mooring created
	unsigned GetNLines(){return NLines;}; 									///Returns the number of lines created of a Mooring
	unsigned GetSegsCount(const unsigned numLine);							///Returns the number of segmentos of one Line and Mooring selected
	unsigned GetSegsCount(const unsigned ftid, const unsigned numLine);		///Returns the number of segmentos of one Line and Mooring selected
	double GetFairTen(const unsigned numLine);								///Returns the Tension of one Line and Mooring selecteds
	int GetNodePos(const unsigned numLine, const int NodeNum, double pos[3]);///Returns the position of one node, Line and Mooring selected
	int GetNodePosLink(const unsigned ftid, const unsigned numLine, double pos[3]);///Returns the position of one node, Line and Mooring selected
	unsigned GetNBodies() { return NBodies; };								///Returns the number of Moorings
	unsigned GetBodyReference(const unsigned ftid=0);						///Returns the mooring reference of the ftid sent	
	int LinesClose();														///Ends the simulation. Used for free memory and close files output
	void LogInit(JLog2 * log);												///Receives the log to store (in run.out) and print messages
	void CreateLog();														/// Creates a JLog2 system if DualSPHysics doesn't send its JLog2 object

};
#endif