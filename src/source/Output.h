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

/// \file Output.h \brief Defines the class \ref Output.

#ifndef _Output_
#define _Output_

#include "JObject.h"
#include "JXml.h"
#include "Line.h"
#include "Connection.h"
#include "TypesMoorDyn.h"


class OutProperties : protected JObject
{
private:
  Magnitude Mag; ///< Output magnitude [tension, force, velocity, position]
  TypeConnection Type_con; ///< Type of output connection[fixed, vessel, all]
  TypeObject Type_obj; ///< Type object output [line, connection]
  std::vector<Line*> Lines; ///< Array of pointers of lines selected
  unsigned NLines; ///< Num of lines selected
  std::string Units; ///< Pointer to units charcters
  unsigned NodeID; ///< Node identifier [0, 1, 2,..., N]
  void Reset(); /// Restores attributes
public:
  OutProperties(); /// Constructor
  ~OutProperties(); /// Destructor
  unsigned GetNodeID() { return NodeID; }; /// Returns the Id of Node
  Magnitude GetMagnitude() { return Mag; }; /// Returns the Magnitude 
  TypeConnection GetTypeCon() { return Type_con;}; /// Returns the type of Connection
  TypeObject GetTypeObj() { return Type_obj; }; /// Returns the type of Object
  std::vector<Line*> GetLines() { return Lines; }; /// Returns an array of Lines
  std::string GetUnits() { return Units; }; /// Returns a pointer to Units
  unsigned GetNLines() { return NLines; }; /// Returns the number of lines
  void SetNodeID(const unsigned id) { NodeID=id; }; /// Sets the Id of Node
  void SetMagnitude(const Magnitude mag) { Mag=mag;}; /// Sets the Magnitude  
  void SetTypeObj(const TypeObject type) { Type_obj=type;}; /// Sets the type of object
  void SetTypeCon(const TypeConnection type) { Type_con=type;}; /// Sets the type of connection
  void SetLines(std::vector<Line*> lines) { Lines=lines;  }; /// Sets the Lines
  void SetNLines(const unsigned num) { NLines=num; }; /// Sets the number of lines
  void SetUnits(const std::string units) { Units=std::string(units); }  /// Sets the units
};
class Output : protected JObject
{
private:
  unsigned NLines; ///< Number of lines
  unsigned NTypes; ///< Number of types for write
  double DtOut; ///< Time step for write
  double TimeStart; ///< Start time for begin to write
  double TimeMax; ///< Time of simulation
  std::string FileName;

  std::vector<OutProperties *> OutProps; ///< Pointer to output properties
  std::string DataDir; ///< string with path 
  void Reset(); /// Restores attributes
  void ReadXml(JXml *sxml, TiXmlElement *ele, std::vector<Line*> lines); /// Reads Xml file
  void SetLinesSelected(std::string str, OutProperties * oProps, std::vector<Line*> lines); /// Checks the lines selected by user and assinngs it to Output class

public:
  Output(double timemax,double dtout); /// Constructor
  ~Output(); /// Destructor
  void LoadXml(JXml *sxml, const std::string &place, std::vector<Line*> lines); /// Loads Xml file
  unsigned GetNLines() { return NLines; }; /// Returns the number of lines                    
  unsigned GetNTypes() { return NTypes; }; /// Returns the nuber of types    
  double GetDtOut() { return DtOut; }; /// Returns the DtOut
  std::vector<OutProperties *> GetOutProps() { return OutProps; }; /// Returns an array of OutProperties  
  double GetStartTime() { return TimeStart; }; /// Returns the StartTime
  double GetEndTime() { return TimeMax; }; /// Returns the EndTime
  void SetTypeOut(std::string str, OutProperties *oProps); /// Sets the connection type selected
  void SetNLines(const unsigned numLines) { NLines=numLines; }; /// Sets the number of Lines    
  void Setup(std::string dir); /// Initialize objects
  void SaveCsv(double timestep, double dt);
};
#endif //!Output
