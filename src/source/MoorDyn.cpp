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

#include "MoorDyn.h"
#include "JLog2.h"

#include <string>
//#include <omp.h>
#include <time.h>
//==============================================================================
/// Constructor.
//==============================================================================
MoorDyn::MoorDyn() {
	FairArrays=false;
	ItArrays=false;
	NStateVariables=6;
	FileIN="moordyn.xml";
	RootNodeXml="moordyn";
	Env=NULL;
	OutP=NULL;
	LineDef=NULL;
	Log=NULL;
	Coupling=true;
	DSPH_Log=false;
#ifdef WIN32
	Dirout="\\MoorDyn_out";			///< Default Output directory (Is replaced when an external software calls to LinesInit() )
#else
	Dirout="/MoorDyn_out";			///< Default Output directory (Is replaced when an external software calls to LinesInit() )
#endif // WIN32
	Reset();
	ClassName="MoorDyn";
}

//==============================================================================
/// Destructor.
//==============================================================================
MoorDyn::~MoorDyn() {
	FreeVector(Lines);
	FreeVector(Connections); 
	FreeVector(Bodies); 
	free2Darray(Ffair, NFairs);
	free2Darray(rFairi, NFairs);
	free2Darray(rdFairi, NFairs);
	Reset();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void MoorDyn::Reset() {
	delete Env; 	Env=NULL;
	delete OutP; 	OutP=NULL;
	delete LineDef; LineDef=NULL;
	if(!DSPH_Log){delete Log;	Log=NULL;};
	if (FairArrays) {
		if (FairleadPos)   FreeArray3D(FairleadPos, NBodies);  	FairleadPos=NULL;
		if (FairleadVel)   FreeArray3D(FairleadVel, NBodies);  	FairleadVel=NULL;
		if (FairleadForce) FreeArray3D(FairleadForce, NBodies); FairleadForce=NULL;
	}

	FairArrays=false;

	NBodies=0;	NLines=0;	NConnections=0;
	NFairs=0;	NAnchs=0;	NConns=0;
	XmlReaded=false;		HeaderProg=false;

	if(ItArrays){
		delete States;			States=NULL;
		delete F0;				F0=NULL;
		delete F1;  			F1=NULL;
		//delete F2;  			F2=NULL;
		//delete F3;  			F3=NULL;
		delete Xt;  			Xt=NULL;
	}
	ItArrays=false;
	Debug=false;
	NX=0;
}

//==============================================================================
/// Frees memory vector
//==============================================================================
template <class T>
void MoorDyn::FreeVector(std::vector<T *> v){
	std::string function="FreeVector";
	for(unsigned i=0; i<v.size(); i++){delete v[i];}
	v.clear();
}

//==============================================================================
/// Creates the output directory for store LineX_BodyX.txt files.
//==============================================================================
void MoorDyn::CreateOutDirectory() {
	int check; //Used for know if output directory was created correctly
	Dirout=fun::GetCurrentDir()+Dirout;
	check=fun::Mkdir(Dirout);

#ifdef WIN32
	Dirout=Dirout + "\\";
#else
	Dirout=Dirout + "/";
#endif
	//PRINT: look the output directory
	if (Debug) { printf("Output directory: %s\n", Dirout.c_str()); }
		
	// check if directory wasn't created
	//TODO: if the folder is empty throws the exception
	// if (!check) {Run_Exceptioon("CreateOutDirectory", std::string("Cannot create the output directory \'") + Dirout + "\'.");} 
}
//==============================================================================
/// Creates a JLog2 system if DualSPHysics doesn't send its JLog2 object
//==============================================================================
void MoorDyn::CreateLog(){
	Log=new JLog2;
	Log->Init(Dirout+"Run_MoorDyn.out");
	Log->AddFileInfo(Dirout+"Run_MoorDyn.out","Log file of the MoorDyn simulation.");
}
//==============================================================================
/// Allocates memory for Bodies.
//==============================================================================
/* void MoorDyn::AllocateMemory() {
std::string function="AllocateMemory"; 
delete[] OutP;			OutP=NULL;
try {
OutP=new Output;
}
catch (const std::bad_alloc) {
Run_Exceptioon("Could not allocate the requested memory.");
}

}*/
//==============================================================================
/// 3D double array destruction functions.
//==============================================================================
void MoorDyn::FreeArray3D(double*** v, unsigned sizex) {
	for (unsigned c=0;c<sizex;c++) {
		unsigned sizey=sizeof(v[c]) / sizeof(double);
		for (unsigned f=0;f<sizey;f++)delete[] v[c][f];
		delete[] v[c];
	}
	delete[] v;
}
//==============================================================================
/// Return a Pointer to Pointer od doubles
//==============================================================================
double ** MoorDyn::GetPtrToPtr(const unsigned sizex, const unsigned sizey) {
	std::string function="GetPtrToPtr"; 
	double **pToPtr;
	if (sizex) {
		try {
			pToPtr=new double*[sizex];
		}
		catch (const std::bad_alloc) {
			Run_Exceptioon("Could not allocate the requested memory.");
		}
		for (unsigned c=0;c<sizex;c++) {
			pToPtr[c]=new double[sizey];
			for (unsigned e=0;e<sizey;e++) {
				pToPtr[c][e]=0;
			}
		}
	}
	return pToPtr;
}
//==============================================================================
/// Return a Pointer to Pointer od doubles
//==============================================================================
double*** MoorDyn::GetPtrToPtrToPtr(const unsigned sizex, const unsigned sizey, const unsigned sizez) {
	std::string function="GetPtrToPtr";
	double*** pToPtr;
	if (sizex) {
		try {
			pToPtr=new double**[sizex];
		}
		catch (const std::bad_alloc) {
			Run_Exceptioon("Could not allocate the requested memory.");
		}
		for (unsigned c=0;c<sizex;c++) {	
			pToPtr[c]=new double*[sizey];
			for (unsigned e=0;e<sizey;e++) {
				pToPtr[c][e]=new double[sizez];
				for (unsigned f=0; f<sizey; f++) { pToPtr[c][e][f]=0; }
			}
		}
	}
	return pToPtr;
}
//==============================================================================
/// Return a Pointer to Pointer of int
//==============================================================================
int * MoorDyn::GetPtrToInt(const unsigned size) {
	std::string function="GetPtrToInt"; 
	int *pToInt;
	if (size) {
		try {
			pToInt=new int[size];
		}
		catch (const std::bad_alloc) {
			Run_Exceptioon("Could not allocate the requested memory.");
		}
		std::memset(pToInt, 0, sizeof(int)*size);
	}
	return pToInt;
}
//==============================================================================
/// Return a Pointer to Pointer of doubles
//==============================================================================
double * MoorDyn::GetPtrToDouble(const unsigned size) {
	std::string function="GetPtrToDouble"; 
	double *pToDbl;
	if (size) {
		try {
			pToDbl=new double[size];
		}
		catch (const std::bad_alloc) {
			Run_Exceptioon("Could not allocate the requested memory.");
		}
		std::memset(pToDbl, 0, sizeof(double)*size);
	}
	return pToDbl;
}
//==============================================================================
/// Allocates memory for arrays of integration
//==============================================================================
void MoorDyn::AllocateMemoryIntegration(const unsigned size) {
	std::string function="AllocateMemoryIntegration"; 

	if (size) {
		try {
			States=new double[size];
			F0=new double[size];
			F1=new double[size];
			//F2=new double[size];	
			//F3=new double[size];																				
			Xt=new double[size];
		}
		catch (const std::bad_alloc) {
			Run_Exceptioon("Could not allocate the requested memory.");
		}
		ItArrays=true;
		std::memset(States, 0, sizeof(double)*size);
		std::memset(F0, 0, sizeof(double)*size);
		std::memset(F1, 0, sizeof(double)*size);
		//std::memset(F2,0,sizeof(double)*size);
		//std::memset(F3,0,sizeof(double)*size);
		std::memset(Xt, 0, sizeof(double)*size);
	}
}
//==============================================================================
/// Allocates memory for persistent vectors
//==============================================================================
/*
void MoorDyn::AllocateMemoryPersistent(){
FlinesS.resize(nStateVariables);
rFairtS=GetPtrToPtr(NFairs,3);
rFairRel.GetPtrToPtr(NFairs,3);


//rFairi.resize  (NFairs);  // after applying platform DOF ICs, should eventually pass this rather than rFairt to Line.setup()
//	rdFairi.resize (NFairs);
}*/
//==============================================================================
///Receives the log to store (in run.out) and print messages
//==============================================================================
void MoorDyn::LogInit(JLog2 * log) {
	DSPH_Log=true;
	if(log)Log=log;
	else Run_Exceptioon("ERROR: JLog2 Class wasn't initialized.");
}
//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void MoorDyn::LoadXml(const std::string FileXml) {
	std::string function="LoadXml"; 
	XmlReaded=true;
	FileIN=FileXml;
	if (!fun::FileExists(FileXml))Run_ExceptioonFile("MoorDyn configuration was not found.", FileXml);
	JXml xml; xml.LoadFile(FileXml);
	TiXmlNode* node=xml.GetNode(RootNodeXml, false);
	if (!node)Run_Exceptioon("Cannot find the element \'" + RootNodeXml + "\'.");
	ReadXml(&xml);
}

//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void MoorDyn::ReadXml(JXml *sxml) {
	std::string function="ReadXml";
	//AllocateMemory();
	//-Gets execution mode
	std::string place = RootNodeXml+".mode";
	if (sxml->GetNode(place, false)) {
		TiXmlElement* eleMode = (sxml->GetNode(RootNodeXml, false))->ToElement();
		std::string debug_str=sxml->ReadElementStr(eleMode, "mode", "value", true, "release");
		if (debug_str == "debug") { Debug = true; }
	}
	//-Configuration of Environment
	place=RootNodeXml + ".solverOptions";
	if (sxml->GetNode(place, false)) {
		Env=new Environment(Log);
		Env->LoadXml(sxml, RootNodeXml);
		//PRINT: look when reads solverOptions tag
		if (Debug) { printf("%s() - Environment\n", function.c_str()); }
	}

	//-Configuration of Bodies
	place=RootNodeXml + ".bodies";
	if (sxml->GetNode(place, false)) {
		unsigned nBodiesTag=sxml->CountElements(sxml->GetNode(RootNodeXml, false), "bodies");
		if(nBodiesTag>1){Run_Exceptioon("Error: There are more \"bodies\" section than one");}
		TiXmlElement* eleBodies=(sxml->GetNode(place, false))->ToElement();
		//PRINT: look when reads MoorDyn tag
		if (Debug) { printf("%s() - Body - Num Bodies: %d\n", function.c_str(), NBodies); }
		NBodies=sxml->CountElements(sxml->GetNode(place, false), "body");	

		TiXmlElement* eleBody=eleBodies->FirstChildElement("body");
		for(unsigned b=0; b<NBodies; b++){
			Body * body=new Body();
			body->LoadXml(sxml, place, eleBody, (unsigned)Bodies.size());
			References.push_back(body->GetNumber());
			CheckReferences(body->GetRef()); //check if exist someone reference equals than this
			CheckDepth(Env, body); //Checks if the user inserts a depth for a body (if not, stores a Environment detph)
			Bodies.push_back(body); //stores the Body
			if(b+1<NBodies)eleBody=sxml->GetNextElement(eleBody, "body");
		}
		ConfigBodies();
	}

	//-Configuration of connects
	place=RootNodeXml + ".connects";
	std::vector<Connection *> connects;
	if (sxml->GetNode(place, false)) {	
		//-Count "connects" blocks	
		unsigned nConnectsTag=sxml->CountElements(sxml->GetNode(RootNodeXml, false), "connects");
		if(nConnectsTag>1){Run_Exceptioon("Error: There are more \"connects\" section than one");}
		//-Count connect tags	
		TiXmlElement* eleConnects=(sxml->GetNode(place, false))->ToElement();
		unsigned NConnections=sxml->CountElements(sxml->GetNode(place, false), "connect");
		//PRINT: look when reads MoorDyn tag
		if (Debug) { printf("%s() - Num Connects: %d\n", function.c_str() , NConnections); }
		TiXmlElement* eleConn=eleConnects->FirstChildElement("connect");
		for(unsigned c=0; c<NConnections; c++){
			Connection * conn=new Connection();
			conn->LoadXml(sxml, place, eleConn, (unsigned)connects.size());
			CheckConnReferences(conn->GetRef(), connects); //check if exist someone reference equals than this
			connects.push_back(conn); //stores the connect
			if(c+1<NConnections)eleConn=sxml->GetNextElement(eleConn, "connect");
		}
	}

	//-Configuration of lines
	place=RootNodeXml + ".lines";
	if (!sxml->GetNode(place, false)) {Run_Exceptioon("Cannot find the element \'" + place + "\'.");	}
	else{
		unsigned nLinesTag=sxml->CountElements(sxml->GetNode(RootNodeXml, false), "lines");
		if(nLinesTag>1){Run_Exceptioon("Error: There are more \"lines\" section than one");}
		TiXmlElement* eleLines=(sxml->GetNode(place, false))->ToElement();
		//-Configuration of LineDefault
		std::string placeLinDef=place + ".linedefault";
		LineDef=NULL;
		if (sxml->GetNode(placeLinDef, false)) {
			unsigned nLinesDefTag=sxml->CountElements(sxml->GetNode(place, false), "linedefault");
			if(nLinesDefTag>1){Run_Exceptioon("Error: There are more \"linedefault\" section than one");}
			LineDef=new LineDefault();
			LineDef->LoadXml(sxml, place);
		}

		//-Configuration of Line
		std::string placeLine=place + ".line";
		unsigned NLines=sxml->CountElements(sxml->GetNode(place, false), "line");
		//PRINT: look when reads MoorDyn tag
		if (Debug) { printf("%s() - Num Lines: %d\n", function.c_str(), NLines); }
		TiXmlElement* eleLine=eleLines->FirstChildElement("line");
		for(unsigned l=0; l<NLines; l++){
			Line * line=new Line(Log);
			line->LoadXml(sxml, placeLine, eleLine, (unsigned)Lines.size(), connects, LineDef);
			Lines.push_back(line); 
			if(line->HasOutput()){OutLines.push_back(line);} //Stores the lines which will store the data in csv file
			Connections.push_back(line->GetAnchConnect());
			Connections.push_back(line->GetFairConnect());
			if(l+1<NLines)eleLine=sxml->GetNextElement(eleLine, "line");
		}
	}
	ConfigConnections(connects);
	CheckBodyRefs();
	//-Configuration of output
	place=RootNodeXml + ".output";
	if (sxml->GetNode(place, false)) {	
		OutP=new Output();
		OutP->LoadXml(sxml, RootNodeXml, OutLines);
	}else{OutP=NULL;}
	return;
}
//==================================================================================
// Checks if there are Fair conections with Body refence which doesn't exist
//==================================================================================
void MoorDyn::CheckBodyRefs(){
	std::string function="CheckBodyRefs";	
	std::vector<Connection *> fairs = GetConnections(1);
	for(unsigned c=0; c<fairs.size(); c++){
		unsigned bRef=fairs[c]->GetRef(); //Gets the body reference for each fairlead 
		if(!GetBody(bRef)){ Run_Exceptioon("Error: Cannot find someone Body with id " + ToString<unsigned>(bRef) + "."); }
	}
	return;
}
//==================================================================================
/// Checks if all ids of floatings are matching with ids of Bodies
//==================================================================================
void MoorDyn::CheckFtIds(const unsigned ids[], const unsigned numFt) {
	std::string function="CheckFtIds";
	for (unsigned i=0; i<numFt; i++) {
		Body * body=GetBody(ids[i]);
		if (!body) { Run_Exceptioon("Error: Cannot find someone Body with id " + ToString<unsigned>(ids[i]) + "."); }
	}
	return;
}
//==============================================================================
/// Searches if Bodies arrays contais some NULL 
/// ( When user doesn't write lines for any Body )
/// If founds some NULL deletes it and updates number of Bodies
//==============================================================================
void MoorDyn::ConfigBodies() {
	//-Searches Null and deletes it
	unsigned numBodies=(unsigned)Bodies.size();
	for (unsigned m=0; m<numBodies; m++) {
		if (Bodies[m]==NULL) {
			for (unsigned mn=m; mn<numBodies-1; mn++) {Bodies[mn]=Bodies[mn + 1];}
			NBodies--;
		}
	}
	//-Order Bodies by id
	numBodies=(unsigned)Bodies.size(); //Updates if was modify
	for (unsigned c=0;c<numBodies - 1;c++) {
		for (unsigned c2=c + 1;c2<numBodies;c2++) {
			if (Bodies[c]->GetRef()>Bodies[c2]->GetRef()) { //Sort by ftid and changes Number
				Body* mo=Bodies[c];
				Bodies[c]=Bodies[c2];
				Bodies[c]->SetNumber(c2);
				Bodies[c2]=mo;
				Bodies[c2]->SetNumber(c);
			}
		}
	}
	NBodies=(unsigned)Bodies.size();
}
//==============================================================================
/// Configures the Connections
//==============================================================================
void MoorDyn::ConfigConnections(std::vector<Connection *> connects){
	std::string function="ConfigConnections";
	NLines=(unsigned)Lines.size();
	//-Configures connect Connections
	std::vector<Connection*> conns_temp;
	for(unsigned c=0; c<connects.size(); c++){
		Connection * conn=connects[c]; //Gets the curren connect
		unsigned ref=conn->GetRef();
		for(unsigned l=0; l<NLines; l++){
			Line * line=Lines[l]; //gets the current line
			if(line->HasConnect()){ 
				std::vector<Connection*> connsLine=line->GetConnections(2); //Gets connects
				for(unsigned c=0; c<line->GetNConns(); c++){
					if(connsLine[c]->GetRef()==ref){ 
						if(std::find(conns_temp.begin(), conns_temp.end(), connsLine[c])==conns_temp.end()) {
							//-Attached doesn't contain line
							conns_temp.push_back(conn);
							break; // Goes out of the loop
						}
					}
				}
			}
		}
	}

	//-Checks if all connects are used
	if(connects.size()!=connects.size()){Run_Exceptioon("Error: There are \"connects\" not used by the lines.");}
	//-Stores the number of connects used
	NConns=(unsigned)conns_temp.size();
	//-Configures anch and fair Connections
	std::vector<Connection*> connections_temp; 
	for(unsigned l=0; l<NLines; l++){
		Line * line=Lines[l]; //gets the current line
		std::vector<Connection*> connections=line->GetConnections(); //Gets all connections
		for(unsigned c=0; c<line->GetNConnections(); c++){
			Connection * connection=connections[c]; //Gets the curren connect
			if(connection->GetType()!=2){ //If the connection isn't a connect
				connections_temp.push_back(connection);
			}
		}
	}

	//-Stores the number of all connections
	NConnections=(unsigned)connections_temp.size() + NConns;
	Connections.clear();
	//-Orders the connections
	//Stores connect Connections
	for(unsigned c=0; c<NConns; c++){ 
		Connections.push_back(conns_temp[c]);
		Connections[c]->SetNumber(c);
	}
	//Stores anch and fair Connections
	for(unsigned c=NConns; c<NConnections; c++){
		Connections.push_back(connections_temp[c-NConns]);
		Connections[c]->SetNumber(c);
		if(Connections[c]->GetType()==0){NAnchs++;}
		else if(Connections[c]->GetType()==1){NFairs++;}
	}
	NConnections=(unsigned)Connections.size();
}
//==============================================================================
/// Checks if exists equals references
//==============================================================================
void MoorDyn::CheckReferences(const unsigned ref) {
	std::string function="CheckReferences";
	for (unsigned m=0; m<Bodies.size(); m++) {
		if (Bodies[m]) {	//Check if isnt't NULL
			if (Bodies[m]->GetRef()==ref) {
				string tex="Error: The ref="+ToString(ref) +" of the bodies is duplicated.";
				Run_Exceptioon(tex);
			}
		}
	}
}
//==============================================================================
/// Checks if exists equals references
//==============================================================================
void MoorDyn::CheckConnReferences(const unsigned ref, std::vector<Connection *> conns) {
	std::string function="CheckConnReferences";
	for (unsigned m=0; m<conns.size(); m++) {
		if (conns[m]->GetRef()==ref) {
			string tex="Error: The conref=" + ToString(ref) + " of the connects is duplicated.";
			Run_Exceptioon(tex);
		}
	}
}
//==============================================================================
/// Checks depth value of body and changes it if doesn't initialized
//==============================================================================
void MoorDyn::CheckDepth(Environment * env, Body * body) {
	if (*body->GetDepth()==DBL_MAX) {
		body->SetDepth(*env->GetWtrDpth());
		//PRINT: look if depth is changed because this tag doesn't exists
		if (Debug) { printf("\n**** The wather depth of Body %d is %.2f ( WaterDepth choosen in solverOptions ) \n", body->GetRef(), *body->GetDepth()); }
	}
}
//==============================================================================
/// Return a pointer to array of Body
//==============================================================================
Body * MoorDyn::GetBody(const unsigned ref) {
	for (unsigned m=0; m<NBodies; m++) {
		Body * body=Bodies[m];
		if (body->GetRef()==ref) { return body; }
	}
	return NULL;
}
//==============================================================================
/// Returns a pointer to Body with reference passed
//==============================================================================
Body * MoorDyn::GetBody(Line * line) {
	std::vector<Connection *> fairs=line->GetConnections(1); //Gets Fairs
	for(unsigned c=0; c<fairs.size(); c++){ //For each fair
		unsigned ref=fairs[c]->GetRef(); //Gets ref of the current fairlead
		for (unsigned m=0; m<NBodies; m++) { 
			Body * body=Bodies[m];//Searches the body with the same ref
			if (body->GetRef()==ref) { return body; }
		}
	}
	return NULL;
}

//==============================================================================
/// Return a pointer to array of lines of Body 
//==============================================================================
std::vector<Line *> MoorDyn::GetLines(Body * body_in) {
	std::vector<Line *> lines;
	unsigned nL=body_in->GetNLines(); // num lines

	Body * body=body_in; //Current body
	for (unsigned l=0; l<nL; l++) {
		lines.push_back(body->GetLines()[l]);
	}
	return lines;
}

//==============================================================================
/// Return a pointer to array of lines of Body by id
//==============================================================================
std::vector<Line *> MoorDyn::GetLines(const unsigned ref) {
	std::vector<Line *> lines;
	for (unsigned l=0; l<Lines.size(); l++) {
		Line * line=Lines[l];
		if(line->HasFair()){
			std::vector<Connection *> connsF=line->GetConnections(1);
			for(unsigned c=0; c<connsF.size(); c++){
				if(connsF[c]->GetRef()==ref){lines.push_back(line);break;}
			}
		}
	}
	return lines;
}
//==============================================================================
/// Return a pointer to array of lines without body
//==============================================================================
std::vector<Line *> MoorDyn::GetLinesNoBody() {
	std::vector<Line *> lines;
	for (unsigned l=0; l<Lines.size(); l++) {
		Line * line=Lines[l];
		if(!line->HasFair()){lines.push_back(line);	}
	}
	return lines;
}

//==============================================================================
/// Return a pointer to array of connections of a Type
//==============================================================================
std::vector<Connection *> MoorDyn::GetConnections(const unsigned type) {
	std::vector<Connection *> connections;//Stores the connections by type
	for (unsigned c=0; c<NConnections; c++) {
		Connection * conn=Connections[c];
		if(conn->GetType()==type){connections.push_back(conn);}
	}
	return connections;
}
//==============================================================================
/// Removes the broken line from the system
//==============================================================================
void MoorDyn::RemoveLine(Line *line){
	if(NLines){
		Lines.erase(std::remove(Lines.begin(), Lines.end(), line), Lines.end());
		NLines--;
	}
}
//==============================================================================
/// Searches the connected lines with a connect Connection
//==============================================================================
void MoorDyn::AddAttachedLines(){
	std::vector<Connection *> conns=GetConnections(2); //-Catches the connects connections
	for(unsigned cc=0; cc<NConns; cc++){
		Connection * conn=conns[cc]; //For each connect connection
		unsigned ref=conn->GetRef();
		for(unsigned l=0; l<NLines; l++){ // Searches the connected lines to this connection
			Line * line=Lines[l];
			//-If there are connects connections with the same reference, add the line to this connect
			if(line->GetNConns(ref)>0) {conn->AddLineToConnect(line);	}
		}
	}
	return;
}

//==============================================================================
/// Calls to setups of all objects.
//==============================================================================
void MoorDyn::Setup() {
	std::string function="Setup";
	//PRINT:
	if (Debug) { printf("\n    ---------- START Setup --------\n\n"); }
	double * Ucurrent=GetPtrToDouble(3);//TODO: checks the meanig of this variable

	//-Setup of Connections
	for (unsigned c=0; c<NConnections; c++) {Connections[c]->Setup();}

	//-Setup of lines
	for (unsigned l=0; l<NLines; l++) {
		Line * line=Lines[l]; //current line		
		if (Debug) { printf("\tLine: %d\n", l); }
		Body * body=NULL;
		if(line->HasFair()){body=GetBody(line);} //-Gets its body
		line->Setup(&Dirout, Env, body);
		line->SetupWaves(Ucurrent, 0.0);
		for (unsigned c=0; c<line->GetNConnections(); c++) {
			Connection * connection=line->GetConnections()[c]; //current Connection
			connection->SetEnv(Env);
			if((connection->GetType()<0)||(connection->GetType()>2)){Run_Exceptioon("Error: Wrong Connection type found: "+ ToString(connection->GetType()) + " in Line: " + ToString(line->GetNumber())+ ".");}
			connection->AddLineToConnect(line);
		}
	}

	//-Setup of Bodies	
	for (unsigned m=0; m<NBodies; m++) {
		Body * body=Bodies[m]; //current body
		std::vector<Line *> lines=GetLines(body->GetRef());
		body->Setup(&Dirout, Env, lines);
		//PRINT: look num body which is computed
		if (Debug) { printf("\nBody: %d\n", m); }
	}

	//-Setup of Outputs
	if(OutP){OutP->Setup(&Dirout);}

	//-Attaches the lines with connects
	AddAttachedLines();
	return;
}

//==================================================================================
/// Prepares the States
//==================================================================================
void MoorDyn::PrepareState() {
	unsigned n=NConns*NStateVariables;
	for (unsigned l=0; l<NLines; l++) {
		Line * line=Lines[l];
		line->SetIndex(n);  		// assign start index of each line
		n +=(line->GetN()-1)*6;	// add 6 state variables for each internal node
	}
	// make state vector	
	NX=n;  // size of state vector array
	if (wordy>1) { printf("\n   Creating state vectors of size %d \n", NX);};

	AllocateMemoryIntegration(NX);  //reserves memory for integration arrays arrays with size NX

									// make array used for passing fairlead kinematics and forces between fairlead- and platform-centric interface functions
	Ffair=make2Darray(NFairs, 3);
	rFairi=make2Darray(NFairs, 3);
	rdFairi=make2Darray(NFairs, 3);
}

//==================================================================================
/// Inserts in a vector the connections positions
//==================================================================================
void MoorDyn::StoreConnPositions() {
	unsigned type=1; // Fair connections
	unsigned nConn_in=0;// Number of connections to store in array
						//Select the Number of connections of the type introduced
	if (type==0) { nConn_in=NAnchs; } //Fixed
	else if (type==1) { nConn_in=NFairs; } //Vessel
	else { nConn_in=NConns; } // Connects

	rFairtS=GetPtrToPtr(nConn_in, 3); //returns a double **
	std::vector<Connection *> connections=GetConnections(type); // returns a pointer to connections of type selected

	for (unsigned c=0; c<nConn_in; c++) {
		Connection * connection=connections[c];
		rFairtS[c][0]=connection->GetX();
		rFairtS[c][1]=connection->GetY();
		rFairtS[c][2]=connection->GetZ();
	}
}
//==================================================================================
/// Initialize system, including trying catenary IC gen of Lines
//==================================================================================
void MoorDyn::InitializeSystem(double X[], double XD[]) {
	//PRINT:
	if (Debug) { printf("\n    ---------- START InitializeSystem --------\n\n"); }
	Log->Printf("\nCreating mooring system:\n Bodies: %d\n Lines: %d\n Connections: %d [ Fairleads: %d, Anchors: %d, Connects: %d ]", NBodies, NLines, NConnections, NFairs, NAnchs, NConns);
	double * x=X;
	double * xd=XD;

	RotMat(x[3], x[4], x[5], TransMat);
	std::vector<Connection *> connAnch=GetConnections(0);
	std::vector<Connection *> connFair=GetConnections(1);
	std::vector<Connection *> connConns=GetConnections(2);

	// set positions of fairleads based on inputted platform position
	for (unsigned l=0; l<NFairs; l++) {
		connFair[l]->InitializeFairlead(x, TransMat); //FairLeads
	}

	// for connect types, write the coordinates to the state vector
	for (unsigned l=0; l<NConns; l++) { 
		connConns[l]->InitializeConnect( States + 6*l ); //Connects // 6 for each connect
	}
	// go through lines and Initialize internal node positions using quasi-static model
	for (unsigned l=0; l<NLines; l++) {
		Line * line=Lines[l];
		line->Initialize(States + line->GetIndex());   // lines
	}

	// write t=-1 output line for troubleshooting preliminary ICs
	//AllOutput(-1.0);
	//cout << "outputting ICs for troubleshooting" << endl;

	//PRINT:
	if (Debug) { printf("\n    ---------- END InitializeSystem --------\n\n"); }
	return;
}

//==================================================================================
/// Finalizing ICs using dynamic relaxation
//==================================================================================
void MoorDyn::DoDynamicRelaxIC() {
	//PRINT:
	if (Debug) { printf("\n    ---------- START DoDynamicRelaxIC --------\n\n"); }
	Log->Printf("\nFinalizing ICs using dynamic relaxation ( %.2f X normal drag ).", Env->GetICDfac());

	for (unsigned l=0; l<NLines; l++) { Lines[l]->scaleDrag(Env->GetICDfac()); }; // boost drag coefficient

	int niic=(int)round(Env->GetICTmax()/Env->GetICdt());	// max Number of IC gen time steps
	double Ffair[3];										// array to temporarily store fairlead force components
	double * Tensions=GetPtrToDouble(NFairs * 3 * niic);	// pointer to double for store tensions for analyzing convergence
	double * FairTens=GetPtrToDouble(NFairs); 				// pointer to double for store tensions for analyzing convergence
	double * FairTensLast=GetPtrToDouble(NFairs); 			// pointer to double for store tensions for analyzing convergence
	double * FairTensLast2=GetPtrToDouble(NFairs); 			// pointer to double for store tensions for analyzing convergence

	unsigned lf=0;		// Fairlead index
	float percent = 0.0;// Percent of Fairlead tensions convergence
	unsigned progress;	// Used to store the current progress of the simulation
	double t = 0.0;		// Convergence current time
	unsigned exp = 0;	// Exponent to round the Percent
	// round to get appropriate mooring model time step
	int NdtM=(int)(ceil(Env->GetICdt()/Env->GetDtM0()));   // Number of mooring model time steps per outer time step (10%) to acelerate
	unsigned step = (unsigned)(NdtM >= 1e6 ? NdtM*1e-1: NdtM*25e-2); //Adjust the step to print progress during the loop for IC generation
	double dtM=Env->GetICdt()/NdtM;	// mooring model time step size (s)
	int iic;						//Main loop iterator
	std::vector< std::vector< double > > fairTension; //Store [time, FairTension] for each step
	// loop through IC generation time analysis time steps
	for (iic=0; iic<niic; iic++) {
		progress = iic * 100 / niic;	// Get current generation progress
		ShowProgressBar(progress);		// Show progres by console
		t=iic*Env->GetICdt();// IC gen time (s).  << is this a robust way to handle time progression?
		//#ifdef _OPENMP
		//#pragma omp parallel for schedule(static)
		//#endif
		// loop through line integration time steps
		for (int its = 0; its < NdtM; its++) {
			if (its % step == 0) { ShowProgressBar((unsigned) (progress + ((its*(100.f/niic)/(float)NdtM))));}
			Rk2(&t, dtM);// call RK2 time integrator (which calls the model)	
		} 	

		// store previous fairlead tensions for comparison
		for (lf=0; lf<NFairs; lf++) {
			FairTensLast2[lf]=FairTensLast[lf];
			FairTensLast[lf]=FairTens[lf];
		}

		std::vector<Connection *> conFair=GetConnections(1); //Returns Fairlead Connections
		// go through connections to get fairlead forces
		for (lf=0; lf<NFairs; lf++) {
			conFair[lf]->GetFnet(Ffair);
			FairTens[lf]=0.0;
			for (int j=0; j<3; j++) { FairTens[lf] +=Ffair[j] * Ffair[j]; }
			FairTens[lf]=sqrt(FairTens[lf]);
		}
		//Store fairlead tensions for print by console
		int size = (int)fairTension.size();
		fairTension.resize(size+1, std::vector<double>(2, 0.0));
		fairTension[size][0]=t;
		fairTension[size][1]=FairTens[0];

		 // check for convergence (compare current tension at each fairlead with previous two values)
		if (iic>2) {
			for (lf=0; lf<NFairs; lf++) {
				if ((abs(FairTens[lf]/FairTensLast[lf] - 1.0)>Env->GetICthresh()) || (abs(FairTens[lf]/FairTensLast2[lf] - 1.0)>Env->GetICthresh()))break;
			}
			if (lf==NFairs) {  // if we made it with all cases satisfying the threshold
				percent=(float)(100.0f*Env->GetICthresh()); 
				exp=(unsigned)(round(log10(percent)));
				//tConvergence = t;
				break; // break out of the time stepping loop
			}
		}
	}
	ShowProgressBar(100);
	if (lf = NFairs) {	Log->Printf("\nFairlead tensions converged to %.*f %% after %.f seconds.", (exp*(-1)), percent, t);	}
	for (unsigned l=0; l<NLines; l++) {
		Lines[l]->scaleDrag(1.0/Env->GetICDfac()); // restore drag coefficients
		Lines[l]->setTime(0.0);		// reset time to zero so first output line doesn't start at>0 s
	}
	for (int i = 0; i < fairTension.size(); i++) {
		Log->Printf("  t=%.f s, tension at first fairlead is %.4f N.", fairTension[i][0], fairTension[i][1]);   // write status update and send cursor back to start of line
	}
	fairTension.clear();

	// @mth: new approach to be implemented
	//	// ------------------------- calculate wave time series if needed -------------------
	//	if (env.WaveKin==2)
	//	{
	//		for (unsigned l=0; l<NLines; l++) 
	//			LineList[l].makeWaveKinematics( 0.0 );
	//	}
	//PRINT:
	if (Debug) { printf("\n    ---------- END DoDynamicRelaxIC --------\n\n"); }
	AllOutput(0, 1.0);   // writes outputs for t=0
	return;
}

//==============================================================================
/// Runge-Kutta 2 integration routine  (integrates States and time)
//==============================================================================
void MoorDyn::Rk2(double *t0, const double dt){
	std::string function="Rk2";
	//CHANGED: x0 by States
	RHSmaster(States, F0, *t0, dt);    // gets derivatives at t0.      F0=f ( t0, x0 );

	for (unsigned i=0; i<NX; i++) {
		Xt[i]=States[i] + 0.5*dt*F0[i];  // integrates to t0  + dt/2.        x1=x0 + dt*F0/2.0;
	}
	for (unsigned i=0; i<NX; i++) {
		if (std::isnan(States[i])) { 
			string tex="Error: 1 NaN value detected in MoorDyn state at dynamic relaxation time " 
				+ ToString(*t0) + " s.";
			Run_Exceptioon(tex); 
		}
	}
	RHSmaster(Xt, F1, *t0 + 0.5*dt, dt);  // gets derivatives at t0  + dt/2.	F1=f ( t1, x1 );

	for (unsigned i=0; i<NX; i++) {
		States[i]=States[i] + dt*F1[i]; // integrates States to t0 + dt
	}
	*t0=*t0 + dt;						// update time			
										// check for NaNs
	for (unsigned i=0; i<NX; i++) {
		if (std::isnan(States[i])) { 
			string tex="Error: 2 NaN value detected in MoorDyn state at dynamic relaxation time " 
				+ ToString(*t0) + " s.";
			Run_Exceptioon(tex); 
		}
	}
	return;
}
//==============================================================================
/// Master function to handle time stepping (updated in v1.0.1 to follow MoorDyn F)
//==============================================================================
void MoorDyn::RHSmaster(const double X[], double Xd[], const double t, const double dt){
	//std::vector<Connection *> connAnch=GetConnections(0); ///< Anchors connections //This kind of connection doesn't used here
	std::vector<Connection *> connFairs=GetConnections(1); 	///< Fairs connections
	std::vector<Connection *> connConns=GetConnections(2); 	///< Connects connections
	std::vector<Connection *> connections=GetConnections(); ///< All connections
															// extrapolate instaneous fairlead positions
	for (unsigned l=0; l<NFairs; l++) { connFairs[l]->UpdateFairlead(t); }
	// calculate forces on all connection objects	
	for (unsigned l=0; l<NConnections; l++) { connections[l]->GetNetForceAndMass(); }
	// calculate connect dynamics (including contributions from latest line dynamics, above, as well as hydrodynamic forces)
	for (unsigned l=0; l<NConns; l++) { connConns[l]->DoRHS((X + 6 * l), (Xd + 6 * l), t); }
	// calculate line dynamics
	for (unsigned l=0; l<NLines; l++) {
		Line * line=Lines[l];
		if(!line->GetDisabled())line->DoRHS((X + line->GetIndex()), (Xd + line->GetIndex()), t, dt);
	}
	return;
}
//==================================================================================
///Funtion to Configure the initial boundary conditions.
/// Used by external software coupling
//==================================================================================
int MoorDyn::LinesInit(double X[], double XD[], const std::string filexml, const std::string nodexml, const char *dir,const tfloat3 gravity,const unsigned ftmkbound[],const unsigned numFts) {
	std::string function="LinesInit";

	string dir_in=dir;
	if ((numFts==0) || (!ftmkbound)) { Run_Exceptioon("Error: The number of driven-object floatings doesn't be 0."); }
	if (dir_in!="") { Dirout=dir; }
	else {CreateOutDirectory();}
	if(!Log){CreateLog();}
	PrintLicense();
	FileIN=filexml;
	RootNodeXml=nodexml;
	if (!XmlReaded) { LoadXml(FileIN); } //Check if Xml was readed before
	Env->SetG(gravity);
	Setup();
	CheckFtIds(ftmkbound, numFts);
	PrepareState();
	StoreConnPositions();
	InitializeSystem(X, XD);
	Env->VisuConfig();
	DoDynamicRelaxIC();
	return 0;
}
//==================================================================================
/// This function now handles the assignment of fairlead boundary conditions, time stepping, and collection of resulting forces at fairleads
/// It can also be called externally for fairlead-centric coupling.
//==================================================================================
int MoorDyn::FairleadsCalc(const unsigned numFts, double*** rFairIn, double*** rdFairIn, double*** fFairIn, double* t_in, double *dt_in) {
	std::string function="FairleadsCalc";

	double t=*t_in;		// this is the current time
	double dtC=*dt_in;	// this is the coupling time step

	if (dtC>0) { // if DT>0, do simulation, otherwise leave passed fFairs unadjusted.
				 // send latest fairlead kinematics to fairlead objects

		for(unsigned b=0; b<numFts; b++){
			if(CheckBody(b)){
				string tex="There isn't a mooring with id " + ToString(b) + ".";
				Run_Exceptioon(tex);
			}
			Body * body=Bodies[b];
			std::vector<Connection *> connFairs=body->GetFairs(); //Get Fair connections of Body passed
			const unsigned nFairs_m=body->GetNFairs(); //Get num of Fairs connections of this body
			for (unsigned l=0; l<nFairs_m; l++) { if(!connFairs[l]->GetDisabled())connFairs[l]->InitiateStep(rFairIn[b][l], rdFairIn[b][l], t); 	}
		}
		// round to get appropriate mooring model time step
		int NdtM=(int)ceil(dtC/Env->GetDtM0());   // Number of mooring model time steps per outer time step
		if (NdtM<1) {
			string tex="Error: dtC is less than dtM. ("
				+ ToString(dtC) + "<" + ToString(Env->GetDtM0()) + ").";
			Run_Exceptioon(tex);
		}

		double dtM=dtC/NdtM;// mooring model time step size (s)

		// loop through line integration time steps (integrate solution forward by dtC)
		for (int its=0; its<NdtM; its++) { Rk2(&t, dtM); } // call RK2 time integrator (which calls the model)
		//for	(unsigned nl=0;nl<NLines;nl++){ if(!Lines[nl]->GetDisabled()&&Lines[nl]->GetBreakTension())Lines[nl]->CheckTension();}											   
		for(unsigned b=0; b<numFts; b++){
			// go through connections to get fairlead forces	
			Body * body=Bodies[b];
			std::vector<Connection *> connFairs=body->GetFairs(); //Get Fair connections of Body passed
			const unsigned nFairs_m=body->GetNFairs(); //Get num of Fairs connections of this body	
			for (unsigned l=0; l<nFairs_m; l++) { if(!connFairs[l]->GetDisabled())connFairs[l]->GetFnet(fFairIn[b][l]);	 }
		}
		AllOutput(t,dtC);   // write outputs

	}
	return 0;
}
//==================================================================================
/// Checks if floating id and body exists. Returns true in error case
//==================================================================================
bool MoorDyn::CheckBody(const unsigned ftid) {
	std::string function="CheckBody";

	bool error=false;
	if ((ftid<0) || (ftid>=NBodies)) {error=true;}
	//-If the body doesn't exist
	if (!Bodies[ftid]) {error=true;}
	return error;
}
//==================================================================================
/// Checks if floating id and body exists. Returns true in error case
//==================================================================================
bool MoorDyn::CheckLine(const unsigned ftid, const unsigned numLine) {
	std::string function="CheckLine";
	bool error=false;

	if(CheckBody(ftid)){
		string tex="There isn't a body with id " + ToString(ftid) + ".";
		Run_Exceptioon(tex);
	}
	Body * body=Bodies[ftid]; 	
	unsigned nLines_b=body->GetNLines();
	if ((numLine<0) || (numLine>=nLines_b)) {error=true;}
	Line * line=Lines[numLine]; 
	//-If the line doesn't exist
	if (!line){ error=true;}
	return error;
}

//==================================================================================
/// Checks if floating id and body exists. Returns true in error case
//==================================================================================
bool MoorDyn::CheckLine(const unsigned numLine) {
	std::string function="CheckLine";
	bool error=false;
	if ((numLine<0) || (numLine>=NLines)) {
		string tex="Error: The range of Lines is ["
			+ ToString(0) + "," + ToString(NLines - 1)
			+ "]. Sent: " + ToString(numLine);
		Run_Exceptioon(tex);
	}
	Line * line=Lines[numLine]; 
	//-If the line doesn't exist
	if (!line){ error=true;}
	return error;
}
//==================================================================================
///Returns the number of fairleads (Vessel connections) of a Mooring created
//==================================================================================
unsigned MoorDyn::GetNFairs(const unsigned ftid) {
	std::string function="GetNFairs";

	if(CheckBody(ftid)){
		string tex="There isn't a body with id " + ToString(ftid) + ".";
		Run_Exceptioon(tex);
	}
	unsigned numF=0;
	Body * body=Bodies[ftid];
	numF=body->GetNFairs();

	return numF;
}

//==================================================================================
/// Returns the Number of segments of the line of floating selected
//==================================================================================
unsigned MoorDyn::GetSegsCount(const unsigned numLine) {
	std::string function="GetSegsCount";
	if(CheckLine(numLine)){
		string tex="There isn't a line with number " + ToString(numLine) + ".";
		Run_Exceptioon(tex);
	}
	unsigned numSeg=0;
	Line * line=Lines[numLine]; 
	numSeg=line->GetN();

	return numSeg;
}
//==================================================================================
/// Returns the Number of segments of the line of floating selected
//==================================================================================
unsigned MoorDyn::GetSegsCount(const unsigned ftid, const unsigned numLine) {
	std::string function="GetSegsCount";

	if(CheckLine(ftid,numLine)){
		string tex="There isn't a line with number " + ToString(numLine) + " for the Body "+ToString(ftid)+".";
		Run_Exceptioon(tex);
	}
	unsigned numSeg=0;
	Body * body=Bodies[ftid]; 	
	Line * line=body->GetLines()[numLine]; 
	numSeg=line->GetN();

	return numSeg;
}
//==================================================================================
/// Returns the Tension of the Fair Selected
//==================================================================================
double MoorDyn::GetFairTen(const unsigned numLine) {
	std::string function="GetFairTen";
	if(CheckLine(numLine)){
		string tex="There isn't a line with number " + ToString(numLine) + ".";
		Run_Exceptioon(tex);
	}
	Line * line=Lines[numLine]; 
	int nNodes=line->GetN();
	double fairTen=line->GetNodeTen(nNodes);

	return fairTen;
}
//==================================================================================
/// Returns the position of Node selected
//==================================================================================
int MoorDyn::GetNodePos(const unsigned numLine, const int NodeNum, double pos[3]) {
	std::string function="GetNodePos";
	if(CheckLine(numLine)){
		string tex="There isn't a line with number " + ToString(numLine) + ".";
		Run_Exceptioon(tex);
	}
	Line * line=Lines[numLine]; 
	int nodePos=line->GetNodePos(NodeNum, pos);	// call line member function to fill in coordinates of the node of interest 
	if (nodePos<0) { //if is -1, an error was produced
		string tex="Error: The range of the Nodes is ["
			+ ToString(0) + "," + ToString(line->GetN())
			+ "]. Sent: " + ToString(NodeNum)
			+ " of Line " + ToString(numLine) + ".";
		Run_Exceptioon(tex);
	}
	return 0;
}

//==================================================================================
/// Returns the position of Node selected
//==================================================================================
int MoorDyn::GetNodePosLink(const unsigned ftid, const unsigned numLine, double pos[3]) {
	std::string function="GetNodePosLink";

	if(CheckLine(ftid,numLine)){
		string tex="There isn't a line with number " + ToString(numLine) + " for the Body "+ToString(ftid)+".";
		Run_Exceptioon(tex);
	}
	Body * body=Bodies[ftid]; 	
	Line * line=body->GetLines()[numLine]; 
	std::vector<Connection *> fairs=line->GetFairs(body->GetRef()); //Gets the fairleads connected to this body
	unsigned nodeNum=fairs[0]->GetNode();//Gets the node number of the fairlead in the line
	int nodePos=line->GetNodePos(nodeNum, pos);	// call line member function to fill in coordinates of the node of interest 
	if (nodePos<0) { //if is -1, an error was produced
		string tex="Error: The range of the Nodes is ["
			+ ToString(0) + "," + ToString(line->GetN())
			+ "]. Sent: " + ToString(nodeNum)
			+ " of Line " + ToString(numLine) + ".";
		Run_Exceptioon(tex);
	}
	return 0;
}
//=================================================================================
///Returns the body reference of the ftid passed	
//=================================================================================
unsigned MoorDyn::GetBodyReference(const unsigned ftid) {
	std::string function="GetBodyReference";

	if(CheckBody(ftid)){
		string tex="There isn't a body with id " + ToString(ftid) + ".";
		Run_Exceptioon(tex);
	}
	Body * body=Bodies[ftid]; 	
	unsigned ref=body->GetRef();

	return ref;
}
//==================================================================================
/// Calls to writer function of the Body passed
//==================================================================================
void MoorDyn::AllOutput(double t, double dtC) {
	//-Writes Connection output of lines
	if(OutP){OutP->WriteOutput(t, dtC);};
	//-Writes individual line output files (of Nodes)
	for (unsigned l=0; l<NLines; l++) { Lines[l]->WriteOutput(t,dtC); }
	return;
}
//==================================================================================
/// Function for finish the program and free the memory
//==================================================================================
int MoorDyn::LinesClose() {
	//PRINT
	Log->Print("\n Lines Close.");
	return 0;
}
//==================================================================================
/// Main function of MoorDyn. Order the necessaries executions. 
/// This function is executed when this library isn't coupling to external software.
//==================================================================================
void MoorDyn::Run() {
	//Only calls CreateOutDirectory() funtion if executes this library how module. 
	//When is coupling with DSPH, use the DSPH output directory
	//PRINT:
	if (Debug) { printf("\n ---------- START RUN --------\n"); }

	Coupling=false;
	CreateOutDirectory();
	if(!Log){CreateLog();}
	PrintLicense();
	if (!XmlReaded) { LoadXml(FileIN); }
	Setup();
	PrepareState();
	StoreConnPositions();
	//****** For init the lines ******
	double * X=new double[6]; //Positions
	double * XD=new double[6]; //Velocities
	std::memset(X, 0, sizeof(double) * 6);
	std::memset(XD, 0, sizeof(double) * 6);
	//********************************
	InitializeSystem(X, XD);
	Env->VisuConfig();
	DoDynamicRelaxIC();
	//*********************************** FAIRLEADS CALC ***********************************
	double t=0; 		///< Initial time
	double * dtm; 		///< Time step
	double dtOut=0.5; 	///< Out time step
	dtm=new double; 	///< Reserves memory
	std::memset(dtm, 0, sizeof(double));
	*dtm=Env->GetDtM0();
	const double divisor=0.00001; ///< Used for increment positions and to make the movement

	while (t<Env->GetTimeMax() + Env->GetDtM0()) {//loop time
		if (!FairArrays) {
			FairleadPos=new double**[NBodies];
			FairleadVel=new double**[NBodies];
			FairleadForce=new double**[NBodies];
		}
		for (unsigned c=0;c<NBodies;c++) {
			Body* b=Bodies[c];
			unsigned nFairs=b->GetNFairs();
			if (!FairArrays) {
				FairleadPos[c]=new double*[nFairs];
				FairleadVel[c]=new double*[nFairs];
				FairleadForce[c]=new double*[nFairs];
			}
			std::vector<Connection *> fairs=b->GetFairs();
			for (unsigned cp=0;cp<nFairs;cp++) {
				const word ptid=fairs[cp]->GetNumber();
				const tdouble3 pos=fairs[cp]->GetPositions();
				const tfloat3  vel=TFloat3(0);
				if (!FairArrays) {
					FairleadPos[c][cp]=new double[3];
					FairleadVel[c][cp]=new double[3];
					FairleadForce[c][cp]=new double[3];
				}
				FairleadPos[c][cp][0]=pos.x;  FairleadPos[c][cp][1]=pos.y;  FairleadPos[c][cp][2]=pos.z;
				FairleadVel[c][cp][0]=vel.x;  FairleadVel[c][cp][1]=vel.y;  FairleadVel[c][cp][2]=vel.z;
			}
		}
		unsigned out=FairleadsCalc(NBodies, FairleadPos, FairleadVel, FairleadForce, &t, dtm);
		UpdatePos(divisor); //Update positions

		if( (t>(floor((t - Env->GetDtM0())/dtOut) + 1.0)* dtOut) || (t==Env->GetTimeMax()) ){ PrintProgress(t); } // if output should print to console
		t +=Env->GetDtM0();
		FairArrays=true;
	}
	if (Debug) { printf("\n ----------- END RUN ---------"); }
	LinesClose();
	Log->Print("\n Execution completed successfully. Press any key to continue...");
	getchar();
	return;
}

//==================================================================================
/// Funtion for simulated the driven-object floating
//==================================================================================
void MoorDyn::UpdatePos(const double divisor) {
	for (unsigned b=0; b<NBodies; b++) {
		Body * body=Bodies[b];
		for (unsigned l=0; l<body->GetNFairs(); l++) {
			Connection * fairlead=body->GetFairs()[l];
			tdouble3 currentPos=fairlead->GetPositions();
			fairlead->SetPositions(TDouble3(currentPos.x, currentPos.y, (currentPos.z+divisor)));
		}
	}
	return;
}
//==============================================================================
/// Prints the execution progress in the console
//==============================================================================
void MoorDyn::PrintProgress(double time) {
	if (!HeaderProg) {
		PrintHeader();
		HeaderProg=true;
	}
	Log->Printf("%.2f\t\t%.2f\t\t%.2f %%", time, Env->GetTimeMax(), PrintPercent(time));
}
//==============================================================================
/// Prints the header for the execution progress
//==============================================================================
void MoorDyn::PrintHeader() {
	Log->Print("\n==================================================");
	Log->Print("Time [s]\tTimeMax [s]\tProgress [%]");
	Log->Print("==================================================");

}
//==============================================================================
/// Prints the completed perecent of the simulation
//==============================================================================
double MoorDyn::PrintPercent(double time) {
	return (100 * time)/Env->GetTimeMax();
}
