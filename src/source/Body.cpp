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

/// \file Body.cpp \brief Implements the class \ref Body.

#include "Body.h"

//==============================================================================
/// Constructor.
//==============================================================================

Body::Body() {
	Reset();
	ClassName="Body";
}

//==============================================================================
/// Destructor.
//==============================================================================
Body::~Body() {
	Reset();
	Lines.clear();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void Body::Reset() {
	NLines=0;
	NFairs=0;
	Env=NULL;
}

//==============================================================================
/// Allocates memory.
//==============================================================================
void Body::AllocateMemory() {
	std::string function="AllocateMemory";

	/*try {
		Env=new Environment;
	}
	catch (const std::bad_alloc) {
		Run_Exceptioon("Could not allocate the requested memory.");
	}
	std::memset(Env, 0, sizeof(Environment));*/
}


//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void Body::LoadXml(JXml *sxml, const std::string &place, TiXmlElement* eleb, const unsigned numBody) {
	std::string function="LoadXml";
	TiXmlNode* node=sxml->GetNode(place, false);
	if(!node)Run_Exceptioon("Cannot find the element \'" + place + "\'.");
	if(!eleb)Run_Exceptioon("Cannot find the element \'" + place + "."+eleb->Value()+"\'.");	
	ReadXml(sxml, place, eleb,  numBody);
}

//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void Body::ReadXml(JXml *sxml, const std::string &place_in, TiXmlElement* eleBody, const unsigned numBody) {
	//-Loads mooring zones.
	//PRINT: look the Number of Mooting tag
	sxml->CheckElementNames(eleBody,false,"depth");

	if(Debug) { printf("numBody - Body: %d\n", numBody); }

	if(eleBody) { 
		Number=numBody; //ID of mooring
		Ref=sxml->GetAttributeInt(eleBody, "ref", true, Number); //If "ref" doesn't exist, stores Number value of creation
		Depth=sxml->ReadElementDouble(eleBody, "depth", "value", true, DBL_MAX); //If depth doesn't exist, stores more later the WaterDepth of SolverOptions
	}
	return;
}

//==============================================================================
/// Returns a pointer of connections of this Body for type
//==============================================================================
Line * Body::GetLine(unsigned number) {
	for(unsigned l=0; l<NLines; l++) {
		Line * line=Lines[l];
		if(line->GetNumber() == number) { return line; }
	}
	return NULL;
}

//==========================================================================================
/// Makes the setup for this Body. Initializes the Outptfile if exists this option in XML
//==========================================================================================
void Body::Setup(std::string &dir, Environment * env_in, std::vector<Line *> lines_in) {
	std::string function="Setup";
	AllocateMemory();
	Env=env_in;
	Lines=lines_in;
	NLines=(unsigned)Lines.size();
	//-Stores the fairs connected to this body
	for(unsigned l=0; l<NLines; l++) {
		Line * line=Lines[l];
		if(line->HasFair()){
			std::vector<Connection *> fairs=line->GetFairs(Ref);
			for(unsigned c=0; c<fairs.size(); c++){Fairs.push_back(fairs[c]);}
		}
	}
	//-Stores the number of fairs
	NFairs=(unsigned)Fairs.size();

	return;
}




