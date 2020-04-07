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

#include "Output.h"

//==============================================================================
/// Constructor.
//==============================================================================
OutProperties::OutProperties() {
	Reset();
	ClassName="OutProperties";
}

//==============================================================================
/// Destructor.
//==============================================================================
OutProperties::~OutProperties() {
	Reset();
	Lines.clear();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void OutProperties::Reset() {
	NLines=0;
	Units=NULL;
}

//==============================================================================
/// Constructor.
//==============================================================================
Output::Output() {
	Reset();
	ClassName="Output";
}

//==============================================================================
/// Destructor.
//==============================================================================
Output::~Output() {
	FreeVector(OutProps);
	Reset();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void Output::Reset() {
	NLines=0;
	NTypes=0;
	StartTime=0;
	EndTime=0;
	Outfile=NULL;
	FileName=NULL;
}
//==============================================================================
/// Frees memory vector
//==============================================================================
void Output::FreeVector(std::vector<OutProperties *> v){
	std::string function="FreeVector";
	for(unsigned i=0; i<v.size(); i++){delete v[i];}
	v.clear();
}
//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void Output::LoadXml(JXml *sxml, const std::string &place, std::vector<Line *> lines) {
	std::string function="LoadXml";

	TiXmlNode* node=sxml->GetNode(place, false);
	if (!node)Run_Exceptioon("Cannot find the element \'" + place + "\'.");
	ReadXml(sxml,node->ToElement(), lines);
}

//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void Output::ReadXml(JXml *sxml, TiXmlElement * ele, std::vector<Line *> lines) {
	NLines=(unsigned)lines.size();
	TiXmlElement* eleOut=ele->FirstChildElement("output");
	
	if (eleOut) {
		TiXmlElement* eleTime=eleOut->FirstChildElement("time");
		if (eleTime) {
			StartTime=sxml->GetAttributeDouble(eleTime, "starTime", true, 0);
			EndTime=sxml->GetAttributeDouble(eleTime, "endTime", true, 0);
			DtOut=sxml->GetAttributeDouble(eleTime, "dtOut", true, 0.01);
		}
		std::string type;
		std::string lines_str;
		NTypes=0;
		//-Tensions
		if (sxml->ExistsElement(eleOut, "tension")) {
			TiXmlElement* elet=eleOut->FirstChildElement("tension");
			OutProps.push_back(new OutProperties());
			OutProps[NTypes]->SetMagnitude(Magnitude::tension);
			type=sxml->GetAttributeStr(elet, "type", true, "all");
			SetTypeOut(type, OutProps[NTypes]);
			OutProps[NTypes]->SetTypeObj(TypeObject::line);
			//TODO: check the number of lines selected
			lines_str="all";
			SetLinesSelected(lines_str, OutProps[NTypes], lines);
			OutProps[NTypes]->SetUnits("[N]");
			NTypes++;
		}
		//-Forces
		if (sxml->ExistsElement(eleOut, "force")) {
			TiXmlElement* elef=eleOut->FirstChildElement("force");
			OutProps.push_back(new OutProperties());
			OutProps[NTypes]->SetMagnitude(Magnitude::force);
			type=sxml->GetAttributeStr(elef, "type", true, "all");
			SetTypeOut(type, OutProps[NTypes]);
			OutProps[NTypes]->SetTypeObj(TypeObject::line);
			//TODO: check the number of lines selected
			lines_str="all";
			SetLinesSelected(lines_str, OutProps[NTypes], lines);
			OutProps[NTypes]->SetUnits("[N]");
			NTypes++;
		}
		//-Velocities
		if (sxml->ExistsElement(eleOut, "velocity")) {
			TiXmlElement* elev=eleOut->FirstChildElement("velocity");
			OutProps.push_back(new OutProperties());
			OutProps[NTypes]->SetMagnitude(Magnitude::velocity);
			type=sxml->GetAttributeStr(elev, "type", true, "all");
			SetTypeOut(type, OutProps[NTypes]);
			OutProps[NTypes]->SetTypeObj(TypeObject::line);
			//TODO: check the number of lines selected
			lines_str="all";
			SetLinesSelected(lines_str, OutProps[NTypes], lines);
			OutProps[NTypes]->SetUnits("[m/s]");
			NTypes++;
		}
		//-Positions
		if (sxml->ExistsElement(eleOut, "position")) {
			TiXmlElement* elep=eleOut->FirstChildElement("position");
			OutProps.push_back(new OutProperties());
			OutProps[NTypes]->SetMagnitude(Magnitude::position);
			type=sxml->GetAttributeStr(elep, "type", true, "all");
			SetTypeOut(type, OutProps[NTypes]);
			//TODO:Check is is line or connection
			OutProps[NTypes]->SetTypeObj(TypeObject::line);
			//TODO: check the number of lines selected
			lines_str="all";
			SetLinesSelected(lines_str, OutProps[NTypes], lines);
			OutProps[NTypes]->SetUnits("[m]");
			NTypes++;
		}
		NTypes=(unsigned)OutProps.size();
		eleOut=eleOut->NextSiblingElement();
	}
}
//==============================================================================
/// Initializes output file
//==============================================================================
void Output::Setup(std::string * dir){
	std::string function = "Setup";
	stringstream  oname;
	//oname << *dir << "MoorDyn" << std::setw(2) << std::setfill('0') << Ref << ".csv"; //Format BodyXX.csv
	oname << *dir << "MoorDyn.csv"; //Format MoorDyn.csv
	FileName=new string(oname.str());
	Outfile=new ofstream(*FileName, std::ofstream::out);
	Outfile->clear();

	if (Outfile) { // check it's not null.  Null signals no individual line output files
		if (Outfile->is_open()) {
			// --- channel titles ---
			*Outfile << "Time [s]" << ";";
			// output all LINE fairlead (top end) tensions
			for (unsigned lf=0; lf<GetNTypes(); lf++) { InitializeHeader(GetOutProps()[lf]); }
			*Outfile << "\n";
			Outfile->close();
		}
		else { Run_Exceptioon("Cannot open file ouput \'" + *FileName + "\'\n"); }
	}
	else { Run_Exceptioon("Cannot create file ouput \'" + *FileName + "\'\n"); }
	
}
//==============================================================================
/// Checks the string passed and store in OutProperties the lines selected
//==============================================================================
void Output::SetLinesSelected(string str, OutProperties * oProps, std::vector<Line *> lines) {
	//TODO: Complete for when the selection isn't "all" -> F.Ex: 1,2,3......
	for(unsigned nl=0;nl<NLines;nl++){lines[nl]->SetDtOut(DtOut);}
	if (str=="all") {
		oProps->SetNLines(NLines);
		oProps->SetLines(lines);
	}
}

//==============================================================================
/// Sets the connection type selected
//==============================================================================
void Output::SetTypeOut(std::string str, OutProperties * oProps) {
	if (str=="all") { oProps->SetTypeCon(TypeConnection::all); }
	else if (str=="fixed") { oProps->SetTypeCon(TypeConnection::fixed); }
	else if (str=="vessel") { oProps->SetTypeCon(TypeConnection::vessel); }
}
//==============================================================================
/// Creates the Header of Output file with the options selected
//==============================================================================
void Output::InitializeHeader(OutProperties * out_p) {
	string type=""; ///< Title for each column
	string units;	///< Units fo each magnitude

	switch (out_p->GetMagnitude()) { //Checks the magnitude 
	case Magnitude::tension: //If wants tension
		for (unsigned l=0; l<out_p->GetNLines() ; l++) {
			if ((out_p->GetTypeCon()==TypeConnection::all) || (out_p->GetTypeCon()==TypeConnection::fixed)) {
				type="AnchTen_L";
				units=*out_p->GetUnits();
				*Outfile << type << out_p->GetLines()[l]->GetNumber() << " " << units << ";";
			}
			if ((out_p->GetTypeCon()==TypeConnection::all) || (out_p->GetTypeCon()==TypeConnection::vessel)) {
				type="FairTen_L";
				units=*out_p->GetUnits();
				*Outfile << type << out_p->GetLines()[l]->GetNumber() << " " << units << ";";
			}
		}
		break;
	case Magnitude::force: //If wants forces
		for (unsigned l=0; l<out_p->GetNLines(); l++) {
			if ((out_p->GetTypeCon()==TypeConnection::all) || (out_p->GetTypeCon()==TypeConnection::fixed)) {
				type="ForceAnch";
				units=*out_p->GetUnits();
					*Outfile << type << "X_L" << out_p->GetLines()[l]->GetNumber() << " " << units << ";"
						<< type << "Y_L" << out_p->GetLines()[l]->GetNumber() << " " << units << ";"
						<< type << "Z_L" << out_p->GetLines()[l]->GetNumber() << " " << units << ";";
			}
			if ((out_p->GetTypeCon()==TypeConnection::all) || (out_p->GetTypeCon()==TypeConnection::vessel)) {
				type="ForceFair";
				units=*out_p->GetUnits();
					*Outfile << type << "X_L" << out_p->GetLines()[l]->GetNumber() << " " << units << ";"
						<< type << "Y_L" << out_p->GetLines()[l]->GetNumber() << " " << units << ";"
						<< type << "Z_L" << out_p->GetLines()[l]->GetNumber() << " " << units << ";";
			}	
		}
		break;
	case Magnitude::velocity: //If wants velocities
		for (unsigned l=0; l<out_p->GetNLines(); l++) {
			if ((out_p->GetTypeCon()==TypeConnection::all) || (out_p->GetTypeCon()==TypeConnection::fixed)) {
				type="VelAnch";
				units=*out_p->GetUnits();
				*Outfile << type << "X_L" << out_p->GetLines()[l]->GetNumber()  << " " << units << ";"
					<< type << "Y_L" << out_p->GetLines()[l]->GetNumber() << " " << units << ";"
					<< type << "Z_L" << out_p->GetLines()[l]->GetNumber() << " " << units << ";";
			}
			if ((out_p->GetTypeCon()==TypeConnection::all) || (out_p->GetTypeCon()==TypeConnection::vessel)) {
				type="VelFair";
				units=*out_p->GetUnits();
				*Outfile << type << "X_L" << out_p->GetLines()[l]->GetNumber()  << " " << units << ";"
					<< type << "Y_L" << out_p->GetLines()[l]->GetNumber() << " " << units << ";"
					<< type << "Z_L" << out_p->GetLines()[l]->GetNumber() << " " << units << ";";
			}
		}
		break;
	case Magnitude::position://If wants positions
		for (unsigned l=0; l<out_p->GetNLines(); l++) {
			if ((out_p->GetTypeCon()==TypeConnection::all) || (out_p->GetTypeCon()==TypeConnection::fixed)) {
				type="PosAnch";
				units=*out_p->GetUnits();
					*Outfile << type << "X_L" << out_p->GetLines()[l]->GetNumber()  << " " << units << ";"
						<< type << "Y_L" << out_p->GetLines()[l]->GetNumber() << " " << units << ";"
						<< type << "Z_L" << out_p->GetLines()[l]->GetNumber() << " " << units << ";";
			}
			if ((out_p->GetTypeCon()==TypeConnection::all) || (out_p->GetTypeCon()==TypeConnection::vessel)) {
				type="PosFair";
				units=*out_p->GetUnits();
				*Outfile << type << "X_L" << out_p->GetLines()[l]->GetNumber() << " " << units << ";"
					<< type << "Y_L" << out_p->GetLines()[l]->GetNumber() << " " << units << ";"
					<< type << "Z_L" << out_p->GetLines()[l]->GetNumber() << " " << units << ";";
			}
		}
		break;
	default:
		break;
	}
	return;
}
//==============================================================================
/// Writes the outputs of all lines
//==============================================================================
void Output::WriteOutput(double t, double dtC) {
	std::string function="WriteOutput";

	if (DtOut>0) { //If out exists and dtOut is greater than 0
		if (t<(floor((t - dtC) / DtOut) + 1.0)*DtOut) {return; } // if output should occur over the course of this time step, then do it!
		if ((StartTime>t) || (EndTime<t)) {return; } //if t isn't into range of start and end times
	}
	else {return; }// If not, ends the write function

	 //PRINT:
	if (Debug) { printf("\n       ---------- START AllOutput t= %f--------\n",t); }
	if (Outfile) {
		Outfile=new ofstream(*FileName, std::ofstream::app); //Concat the new content		
		if (Outfile->is_open()) {
			*Outfile << t << ";"; // output time
			// output all LINE fairlead (top end) tensions
			//for (unsigned l=0; l<NLines; l++) outfileMain << 0.001*(LineList[l].GetNodeTen(LineList[l].GetN())) << "; ";
			for (unsigned lf=0; lf<GetNTypes(); lf++) {
				OutProperties * oP=GetOutProps()[lf];
				for (unsigned l=0; l<oP->GetNLines(); l++) {
					Line * line=oP->GetLines()[l];
					WriteOutputLine(oP, line);	// output each channel's value
				}
			}
			*Outfile << "\n";
		}
		else { Run_Exceptioon("Cannot open file ouput \'" + *FileName + "\'\n"); }
	}
	else { Run_Exceptioon("Cannot create file ouput \'" + *FileName + "\'\n"); }

	Outfile->close();

	//PRINT:
	if (Debug) { printf("\n       ---------- END AllOutput --------\n"); }
	return;
}
//==============================================================================
/// Write the outputs for each line
//==============================================================================
void Output::WriteOutputLine(OutProperties * out_p, Line * line) {
	switch (out_p->GetMagnitude()) {
	case Magnitude::tension:
		if ((out_p->GetTypeCon()==TypeConnection::all) || (out_p->GetTypeCon()==TypeConnection::fixed)) {
			out_p->SetNodeID(0) ;
			*Outfile << line->GetTensionOutput(out_p) << ";";
		}
		if ((out_p->GetTypeCon()==TypeConnection::all) || (out_p->GetTypeCon()==TypeConnection::vessel)) {
			out_p->SetNodeID(line->GetN()) ;
			*Outfile << line->GetTensionOutput(out_p) << ";";
		}
		break;
	case Magnitude::force:
		if ((out_p->GetTypeCon()==TypeConnection::all) || (out_p->GetTypeCon()==TypeConnection::fixed)) {
			out_p->SetNodeID(0);
			tdouble3 out_val=line->GetForceOutput(out_p);
			*Outfile << out_val.x << ";" << out_val.y << ";" << out_val.z << ";";
		}
		if ((out_p->GetTypeCon()==TypeConnection::all) || (out_p->GetTypeCon()==TypeConnection::vessel)) {
			out_p->SetNodeID(line->GetN()) ;
			tdouble3 out_val=line->GetForceOutput(out_p);
			*Outfile << out_val.x << ";" << out_val.y << ";" << out_val.z << ";";
		}
		break;
	case Magnitude::velocity:
		if ((out_p->GetTypeCon()==TypeConnection::all) || (out_p->GetTypeCon()==TypeConnection::fixed)) {
			out_p->SetNodeID(0) ;
			tdouble3 out_val=line->GetVelocityOutput(out_p);
			*Outfile << out_val.x << ";" << out_val.y << ";" << out_val.z << ";";
		}
		if ((out_p->GetTypeCon()==TypeConnection::all) || (out_p->GetTypeCon()==TypeConnection::vessel)) {
			out_p->SetNodeID(line->GetN()) ;
			tdouble3 out_val=line->GetVelocityOutput(out_p);
			*Outfile << out_val.x << ";" << out_val.y << ";" << out_val.z << ";";
		}
		break;
	case Magnitude::position:
		if ((out_p->GetTypeCon()==TypeConnection::all) || (out_p->GetTypeCon()==TypeConnection::fixed)) {
			out_p->SetNodeID(0);
			Connection * conn=line->GetAnchConnect();
			tdouble3 out_val=conn->GetPositions();
			*Outfile << out_val.x << ";" << out_val.y << ";" << out_val.z << ";";
		}
		if ((out_p->GetTypeCon()==TypeConnection::all) || (out_p->GetTypeCon()==TypeConnection::vessel)) {
			out_p->SetNodeID(line->GetN());
			Connection * conn=line->GetFairConnect();
			tdouble3 out_val=conn->GetPositions();
			*Outfile << out_val.x << ";" << out_val.y << ";" << out_val.z << ";";
		}
		break;
	default:
		break;
	}
	return;
}
