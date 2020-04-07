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

#ifndef _Body_
#define _Body_

#include "JXml.h"
#include "Line.h"
#include "Environment.h"
#include "Connection.h"
#include "JObject.h"
#include "Types_moordyn.h"

class Body : protected JObject
{
private:	
	unsigned Ref;					///< Identifier of body introduced in xml file
	unsigned Number; 				///< Identifier of body. Secuential by creation [0, 1, 2...]
	unsigned NFairs;				///< Number of fairleads
	unsigned NLines; 				///< Number of lines of this Body
	std::vector<Line *> Lines; 		///< Vector of pointers to lines of this Body
	std::vector<Connection *> Fairs;///< Vector of pointers to Fairs of its lines
	double Depth; 					///< Stores the Depth of the water (m)
	Environment * Env;				///< Pointer to environmental settings

	void Reset();					/// Restores the attributes	
	void AllocateMemory(); 			/// Reserves memory for arrays
	void ReadXml(JXml *sxml, const std::string &place, TiXmlElement* eleb, const unsigned num_tag);/// Reads Xml file


public:
	Body();														/// Constructor
	~Body();													/// Destructor
	void LoadXml(JXml *sxml, const std::string &place, TiXmlElement* eleb, const unsigned num_tag);	/// Loads Xml file
	void Setup(string * dir, Environment * env_in, std::vector<Line *> lines); /// Makes a setup of this Body
	std::vector<Line *> GetLines() { return Lines; }			/// Returns a vector of Lines
	unsigned GetNFairs() { return NFairs; };					/// Returns the number of Fairleads
	unsigned GetNumber() { return Number; };					/// Returns the number of Body
	unsigned GetNLines() { return NLines; };					/// Returns the number of Body
	void SetNumber(const unsigned num) { Number=num; };			/// Stores the new number
	unsigned GetRef() { return Ref; };							/// Returns the floating id of Body
	unsigned * GetPtrRef() { return &Ref; };					/// Returns a pointer to floating id of Body
	Line * GetLine(unsigned number);							/// Returns a pointer of a line selected
	double * GetDepth() { return &Depth; }						/// Returns a pointer of the depth
	void SetDepth(double value) { Depth=value; }				/// Sets depth value
	unsigned GetFtMk() { return Ref; }							/// Returns the mk value
	std::vector<Connection *> GetFairs(){return Fairs;}; 		/// Returns a pointer of connections of this body
	
};

#endif