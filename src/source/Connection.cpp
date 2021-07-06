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

/// \file Connection.cpp \brief Implements the class \ref Connection.

#include "Connection.h"
#include "Environment.h"
#include "Line.h"

 //==============================================================================
 /// Constructor.
 //==============================================================================
Connection::Connection()  {
	Reset();
	ClassName="Connection";
}

//==============================================================================
/// Destructor.
//==============================================================================
Connection::~Connection(){
	Reset();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void Connection::Reset() {
	Ref=0;    nAttached=0;
	ConX=0; 	ConY=0;		ConZ=0;
	ConFX=0;	ConFY=0;	ConFZ=0;
	ConM=0;		ConV=0;	  NumExec=0;
	ConCdA=0;	ConCa=0;
	R.clear();		R_ves.clear();	
  Rd.clear();   Rd_ves.clear();
	Fnet.clear();	Fnet_i.clear();
	RHS.clear();	S.clear();
	M.clear();		M_i.clear();
  Attached.clear();	Top.clear();
	Q.clear();
	HasMass=false;	HasVol=false;
	Disabled=false;
}
//==============================================================================
/// Restores the tension for the Connection when it is disabled
//==============================================================================
void Connection::ResetTension(){
	Fnet[0]=0.0;
	Fnet[1]=0.0;
	Fnet[2]=0.0;
}
//==============================================================================
/// Allocates memory of vectors
//==============================================================================
void Connection::AllocateMemoryVectors() {
	// size vectors (could change a lot of these to arrays)
	R.resize(3,0.0); // node positions [i][x/y/z]
	Rd.resize(3,0.0); // node velocities [i][x/y/z]
	Q.resize(3,0.0); // unit tangent vectors for each node

	R_ves.resize(3,0.0);
	Rd_ves.resize(3,0.0);

	Fnet.resize(3,0.0); // total force on node
	Fnet_i.resize(3,0.0);
	RHS.resize(3,0.0);	// RHS of state-space equation (Forces divided by mass matrix)

	S.resize(3,vector< double >(3,0.0));  // inverse mass matrices (3x3) for each node
	M.resize(3,vector< double >(3,0.0)); // node mass + added mass matrix
	M_i.resize(3,vector< double >(3,0.0));
}
//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void Connection::LoadXml(JXml *sxml,const std::string &place,TiXmlElement* elec,const int id_connection) {
	std::string function="LoadXml";
	//TiXmlNode* node=sxml->GetNode(place,false);
	//if(!node)Run_Exceptioon("Cannot find the element \'" + place + "\'.");
	if(!elec)Run_Exceptioon("Cannot find the element \'" + place + "\'.");
	ReadXml(sxml,elec,id_connection);
}


//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void Connection::ReadXml(JXml *sxml,TiXmlElement* eleConn,const int id_connection) {
	std::string function="ReadXml";
	//-Loads connections zones.
	std::string value=eleConn->Value(); // Text of connection tag
	//-Connect connections
	if(value=="connects"){
		eleConn=eleConn->FirstChildElement("connect");
		value=eleConn->Value(); 
	}

	//-If the current connection tag exists
	if(eleConn) {
		//-Stores the Connection Type [Anch=0; Fair=1 ; Con=2]
		if(value=="fixconnection") {		 	Type=0;}
		else if(value=="vesselconnection") { 	Type=1; Ref=sxml->GetAttributeUnsigned(eleConn,"bodyref",false);}
		else if(value=="connect") { 			Type=2; Ref=sxml->GetAttributeUnsigned(eleConn,"conref",false);}
		else {
			string tex="Cannot find the element \'"; tex=tex + sxml->ErrGetFileRow(eleConn) + "\'.";
			Run_Exceptioon(tex);
		}

		Number=id_connection; //Stores a secuencial Number of connections in the program for each connection
		//PRINT: Look number of connection stored
		if(Debug) {printf("\t\tnumber of connection: %d\n",Number);}
		
		//-Gets positions
		tdouble3 position=sxml->GetAttributeDouble3(eleConn); 

		//-Stores positions
		ConX=position.x; 
		ConY=position.y;
		ConZ=position.z;
	 	//-Gets Mass and volum attributes
		ConM=sxml->GetAttributeDouble(eleConn,"M",true,0);
		ConV=sxml->GetAttributeDouble(eleConn,"V",true,0);
		if(ConM!=0)HasMass=true;
		if(ConV!=0)HasVol=true;
		//-Stores Forces and Coefficients (Only type=2 [Connects])
		if(Type==2){
			//-Gets forces
			if(sxml->ExistsElement(eleConn,"initialForce")){
				TiXmlElement* eleConnF=eleConn->FirstChildElement("initialForce");
				tdouble3 forces=sxml->GetAttributeDouble3(eleConnF);
				ConFX=forces.x;
				ConFY=forces.y;
				ConFZ=forces.z;
			}
			//-Gets Coefficients
			if(sxml->ExistsElement(eleConn,"initialCoeffs")){
				TiXmlElement* eleCoeff=eleConn->FirstChildElement("initialCoeffs");
				ConCdA=sxml->GetAttributeDouble(eleCoeff,"cda",false);
				ConCa=sxml->GetAttributeDouble(eleCoeff,"ca",false);
			}
		}
	}
	return;
}
//==============================================================================
///  Makes the setup for this Connection. Initializes the position of this
//==============================================================================
void Connection::Setup() {
	//PRINT: look when Connection::Setup() is executed
	if(Debug) { printf("\t\tConnection::setup() - Number: %d\n",Number); }
	t=0.;
	//beta=0.0;
	nAttached=0;  // start off with zero connections
	AllocateMemoryVectors();
	R[0]=ConX;  // start off position at that specified in input file 
	R[1]=ConY;  //  (will be starting point for connect connections
	R[2]=ConZ;  //   and the permanent location of anchor connections.)		
	return;
}

//==============================================================================
// This function handles assigning a line to a connection node
//==============================================================================
void Connection::AddLineToConnect(Line * line) {
	if(misc::wordy>0) { printf("L %u ->N %u \n",line->GetNumber(),Number); }
	if(std::find(Attached.begin(),Attached.end(),line)==Attached.end()) {
		//-Attached doesn't contain line
		unsigned top=line->GetN(); // Gets the number of segments
		Connection * anch=line->GetAnchConnect();
		if(anch==this){top=0;} //If the connect is an "anchor", stores the segment=0
		if(nAttached<10) { // this is currently just a maximum imposed by a fixed array size.  could be improved.
			Attached.push_back(line);
			Top.push_back(top);
			nAttached += 1;
		}
	}
	return;
}
//==============================================================================
/// Function to return connection position and velocity to Line object
//==============================================================================
void Connection::GetConnectState(vector<double> &r_out,vector<double> &rd_out) {
	for(int J=0; J<3; J++) {
		r_out[J]=R[J];
		rd_out[J]=Rd[J];
	}
	return;
}
//==============================================================================
/// Return a TDouble3 with th positions of connection
//==============================================================================
tdouble3 Connection::GetPositions() {
	return TDouble3(R[0],R[1],R[2]);
}
//==============================================================================
/// Stores the new positions [X Y Z]
//==============================================================================
void Connection::SetPositions(tdouble3 pos) {
	R[0]=pos.x;
	R[1]=pos.y;
	R[2]=pos.z;
}
//==============================================================================
/// Function to return net force on fairlead (just to allow public reading of Fnet)
//==============================================================================
void Connection::GetFnet(double Fnet_out[]) {
	for(int I=0; I<3; I++) { Fnet_out[I]=Fnet[I];}
	//Fnet_out[2] += M[0][0]*(-env.g); // add weight  NO this is alread in Fnet !!! (removed Oct 20)
	// should it include inertial "force"?  i.e.	for(int J=0; J<3; J++)  Fnet_out += M[I][J]*acceleration[J] 	
	return;
}


//==============================================================================
/// Initialize Positions and velocities of fairlead. 
/// pX -> Positions sent by LinesInit()
//==============================================================================
void Connection::InitializeFairlead(double pX[],double TransMat[]) {
	std::string function="InitializeFairlead";

	if(Type==1) { //If is vessel Type
		R[0]=TransMat[0] * ConX + TransMat[1] * ConY + TransMat[2] * ConZ + pX[0];	// x
		R[1]=TransMat[3] * ConX + TransMat[4] * ConY + TransMat[5] * ConZ + pX[1];	// y
		R[2]=TransMat[6] * ConX + TransMat[7] * ConY + TransMat[8] * ConZ + pX[2];	// z

		for(int I=0; I<3; I++) { Rd[I]=0.0; }
		// also specify the above in the "vessel" arrays that prescribe the kinematics over the following time steps, for IC gen
		for(int J=0; J<3; J++) {
			R_ves[J]=R[J];
			Rd_ves[J]=0.0;
		}
	}
	else { 
		string tex="Error: wrong connection Type given. Something's not right.";
		Run_Exceptioon(tex); 
	}
	// TODO: should handle.
	return;
};

//==============================================================================
/// Initialize the positions and velocities of this Connection
//==============================================================================
void Connection::InitializeConnect(double* X) {
	std::string function="InitializeConnect";

	if(Type==2) {  // If is connect Type
		// assign initial node kinematics to state vector
		for(int I=0; I<3; I++) {
			X[3 + I]=R[I];
			X[I]=Rd[I];
		}
	}else { 
		string tex="Error: wrong connection Type given. Something's not right.";
		Run_Exceptioon(tex); 
	}
	// TODO: should handle.
	return;
};

//==============================================================================
/// Helper function to sum forces and mass from attached lines 
/// Used for connect dynamics and fair/anch tensions
//==============================================================================
void Connection::GetNetForceAndMass() {
	// loop through each connected line,summing to get the final result
	unsigned Nl=nAttached;				// Number of attached line segments
	// clear before re-summing	
	for(int I=0; I<3; I++) {
		Fnet[I]=0;
		for(int J=0; J<3; J++) { M[I][J]=0.0; }
	}
	// loop through attached lines
	for(unsigned l=0; l<Nl; l++) {	// get quantities
		if(Top[l]==0) { (Attached[l])->GetAnchStuff(Fnet_i,M_i);}		// if attached to bottom/anchor of a line...
		else {	(Attached[l])->GetFairStuff(Fnet_i,M_i);}// attached to top/fairlead
		// sum quantitites
		for(int I=0; I<3; I++) {
			Fnet[I] += Fnet_i[I];
			for(int J=0; J<3; J++) { M[I][J] += M_i[I][J]; }
		}
	}
	//-Checks if the Mass and the Volum were initialized and calculates it (only connects)
	if(Type==2){CalculateMassVol(); } 

	// add Constant quantities for Connection if applicable from input file
	Fnet[0] += ConFX;
	Fnet[1] += ConFY;
	Fnet[2] += ConFZ + ConV * (Env->GetRho_w()) * (Env->GetG()) - ConM * (Env->GetG());
	for(int I=0; I<3; I++) { M[I][I] += ConM; }

	return;
}
//==============================================================================
/// This function calculates the Mass and Volumn of the Connect Connections if
/// the user didn't introduce value for them
//==============================================================================
void Connection::CalculateMassVol(){
	//-Calculates Mass
	if(!HasMass){
		ConM=0;
		unsigned numF=0; //Number of connects that are "Fairs"
		for(unsigned lm=0; lm<nAttached; lm++){
			unsigned posm=Top[lm]; //Gets the node of this connection in the current line
			if(posm>0){	ConM+=Attached[lm]->GetM(posm); numF++;} //If the connect is an "fairlead"
		}
		ConM /= numF; //Divides by number of lines that have Connect<------>Anchor Connections
	}
	//-Calculates Volumn
	if(!HasVol){
		ConV=0;
		for(unsigned lv=0; lv<nAttached; lv++){
			unsigned posv=Top[lv];//Gets the node of this connection in the current line
			if(posv>0){ ConV+=Attached[lv]->GetV(0)*posv*Attached[lv]->GetLength();}
		}
	}
	if(NumExec>2){HasMass=true;	HasVol=true; }
	else{NumExec++; }
}
//==============================================================================
/// This is the function that updates the States. Also includes hydrodynamic forces
//==============================================================================
void Connection::DoRHS(const double* X,double* Xd,const double time) {
	std::string function="DoRHS";

	t=time;

	// below now called by master RHS function, for all connection types
	//GetNetForceAndMass(); // in all cases, call to get net force and mass of connection (especially for returning fairlead forces to FAST)

	// ------ behavior dependant on connect Type -------

 //	if(Type==0){ // fixed Type
 //		R[0]=ConX;
 //		R[1]=ConY;
 //		R[2]=ConZ;
 //		for(int I=0; I<3; I++){ Rd[I]=0; }
 //	}else if(Type==1) {// vessel (moves with platform)					
 //		// set fairlead position and velocity based on BCs (linear model for now)
 //		for(int J=0; J<3; J++)  {
 //			R[J]=R_ves[J] + Rd_ves[J]*(t-t0);
 //			Rd[J]=Rd_ves[J];
 //		}		
 //		
 //		// assign States
 //		for(int I=0; I<3; I++)  {
 //			Xd[3+I]=Rd[I];  // velocities - these are unused in integration
 //			Xd[I]=0.;     // accelerations - these are unused in integration
 //		}
 //	}			
	if(Type==2) { // "connect" Type
		if(t==0) {   // with current IC gen approach, we skip the first call to the line objects, because they're set AFTER the call to the connects
			// above is no longer true!!! <<<
			for(int I=0; I<3; I++) {
				Xd[3 + I]=X[I];  	// velocities - these are unused in integration
				Xd[I]=0.;     		// accelerations - these are unused in integration
			}
		}else {
			// from state values, get R and rdot values 
			for(int J=0; J<3; J++) {
				R[J]=X[3 + J]; // get positions
				Rd[J]=X[J]; // get velocities
			}
			//cout << "ConRHS: m: " << M[0][0] << ", f: " << Fnet[0] << " " << Fnet[1] << " " << Fnet[2] << endl;
			// add dynamic quantities for connection as specified in input file (feature added 2015/01/15)
			Fnet[0] -= 0.5* (Env->GetRho_w()) *Rd[0] * abs(Rd[0])*ConCdA;
			Fnet[1] -= 0.5* (Env->GetRho_w()) *Rd[1] * abs(Rd[1])*ConCdA;
			Fnet[2] -= 0.5* (Env->GetRho_w()) *Rd[2] * abs(Rd[2])*ConCdA;
			for(int I=0; I<3; I++) 	M[I][I] += ConV*(Env->GetRho_w())*ConCa;
			// invert node mass matrix
			misc::inverse3by3(S,M);
			// RHS constant - (premultiplying force vector by inverse of mass matrix  ... i.e. rhs=S*Forces
			for(int I=0; I<3; I++) {
				double RHSI=0.0; // temporary accumulator 
				for(int J=0; J<3; J++) { RHSI += S[I][J] * Fnet[J]; } //  matrix multiplication [S i]{Forces i}
				// update States
				Xd[3 + I]=X[I];    // dxdt=V    (velocities)
				Xd[I]=RHSI;      // dVdt=RHS * A  (accelerations)
			}
		}
	}else { 
		string tex="Error: wrong connection Type given. Something's not right.";
		Run_Exceptioon(tex); 
	}
	return;
}

//==============================================================================
/// Called at the beginning of each coupling step to update the boundary conditions 
/// (fairlead kinematics) for the proceeding line time steps
//==============================================================================
void Connection::InitiateStep(double rFairIn[3],double rdFairIn[3],double time) {
	t0=time; // set start time for BC functions
	if(Type==1) {  // if is vessel Type 
		// update values to fairlead position and velocity functions (fn of time)
		for(int J=0; J<3; J++) {
			R_ves[J]=rFairIn[J];
			Rd_ves[J]=rdFairIn[J];
		}
	}
	// do I want to get precalculated values here at each FAST time step or at each line time step?
	return;
}

//==============================================================================
/// Updates the positions and velocities for the tine passed
//==============================================================================
void Connection::UpdateFairlead(const double time) {
	std::string function="UpdateFairlead";
	t=time;
	if(Type==1) { // vessel (moves with platform)					
		// set fairlead position and velocity based on BCs (linear model for now)
		for(int J=0; J<3; J++) {
			R[J]=R_ves[J] + Rd_ves[J] * (time - t0);
			Rd[J]=Rd_ves[J];
		}
	}else { 
		string tex="Error: wrong connection Type given. Something's not right.";
		Run_Exceptioon(tex); 
	}
	return;
}


// new function to draw instantaneous line positions in openGL context
#ifdef USEGL
void Connection::drawGL(void)
{
	double radius=pow(ConV / (4 / 3 * PI),0.33333);  //ConV
	Sphere(R[0],R[1],R[2],radius);
};
#endif