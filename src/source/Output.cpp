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

/// \file Output.cpp \brief Implements the class \ref Output.

#include "Output.h"
#include "FunMoorDyn.h"
#include "Functions.h"
#include "JSaveCsv2.h"

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
  Units="";
  Type_con=TypeConnection::all;
}

//==============================================================================
/// Constructor.
//==============================================================================
Output::Output(double timemax,double dtout) {
  Reset();
  ClassName="Output";
  TimeMax=timemax;
  DtOut=dtout;
}

//==============================================================================
/// Destructor.
//==============================================================================
Output::~Output() {
  funmd::FreeVector(OutProps);
  Reset();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void Output::Reset() {
  NLines=0;
  NTypes=0;
  TimeStart=0;
  TimeMax=0;
  FileName="";
}

//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void Output::LoadXml(JXml *sxml,const std::string &place,std::vector<Line*> lines) {
  std::string function="LoadXml";

  TiXmlNode* node=sxml->GetNode(place,false);
  if(!node)Run_Exceptioon("Cannot find the element \'" + place + "\'.");
  ReadXml(sxml,node->ToElement(),lines);
}

//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void Output::ReadXml(JXml *sxml,TiXmlElement *ele,std::vector<Line*> lines) {
  NLines=(unsigned)lines.size();

  TiXmlElement* eleOut;
  if     (sxml->ExistsElement(ele,"savedata"))eleOut=ele->FirstChildElement("savedata");
  else if(sxml->ExistsElement(ele,"output")  )eleOut=ele->FirstChildElement("output");//-For compatibility with old versions

  if(eleOut) {
    TiXmlElement* eleTime=eleOut->FirstChildElement("time");
    if(eleTime) {
      TimeStart=sxml->GetAttributeDouble(eleTime,"starTime",true,0);
      TimeMax=sxml->GetAttributeDouble(eleTime,"endTime",true,TimeMax);
      DtOut=sxml->GetAttributeDouble(eleTime,"dtOut",true,DtOut);
    }
    std::string type;
    std::string lines_str;
    NTypes=0;
    //-Tensions
    if(sxml->ExistsElement(eleOut,"tension")) {
      TiXmlElement* elet=eleOut->FirstChildElement("tension");
      if(sxml->GetAttributeBool(elet,"value",true,true)){
        OutProps.push_back(new OutProperties());
        OutProps[NTypes]->SetMagnitude(Magnitude::tension);
        type=sxml->GetAttributeStr(elet,"type",true,"all");
        SetTypeOut(type,OutProps[NTypes]);
        OutProps[NTypes]->SetTypeObj(TypeObject::line);
        //TODO: check the number of lines selected
        lines_str="all";
        SetLinesSelected(lines_str,OutProps[NTypes],lines);
        OutProps[NTypes]->SetUnits("[N]");
        NTypes++;
      }
    }
    //-Forces
    if(sxml->ExistsElement(eleOut,"force")) {
      TiXmlElement* elef=eleOut->FirstChildElement("force");
      if(sxml->GetAttributeBool(elef,"value",true,true)){
        OutProps.push_back(new OutProperties());
        OutProps[NTypes]->SetMagnitude(Magnitude::force);
        type=sxml->GetAttributeStr(elef,"type",true,"all");
        SetTypeOut(type,OutProps[NTypes]);
        OutProps[NTypes]->SetTypeObj(TypeObject::line);
        //TODO: check the number of lines selected
        lines_str="all";
        SetLinesSelected(lines_str,OutProps[NTypes],lines);
        OutProps[NTypes]->SetUnits("[N]");
        NTypes++;
      }
    }
    //-Velocities
    if(sxml->ExistsElement(eleOut,"velocity")) {
      TiXmlElement* elev=eleOut->FirstChildElement("velocity");
      if(sxml->GetAttributeBool(elev,"value",true,true)){
        OutProps.push_back(new OutProperties());
        OutProps[NTypes]->SetMagnitude(Magnitude::velocity);
        type=sxml->GetAttributeStr(elev,"type",true,"all");
        SetTypeOut(type,OutProps[NTypes]);
        OutProps[NTypes]->SetTypeObj(TypeObject::line);
        //TODO: check the number of lines selected
        lines_str="all";
        SetLinesSelected(lines_str,OutProps[NTypes],lines);
        OutProps[NTypes]->SetUnits("[m/s]");
        NTypes++;
      }
    }
    //-Positions
    if(sxml->ExistsElement(eleOut,"position")) {
      TiXmlElement* elep=eleOut->FirstChildElement("position");
      if(sxml->GetAttributeBool(elep,"value",true,true)){
        OutProps.push_back(new OutProperties());
        OutProps[NTypes]->SetMagnitude(Magnitude::position);
        type=sxml->GetAttributeStr(elep,"type",true,"all");
        SetTypeOut(type,OutProps[NTypes]);
        //TODO:Check is is line or connection
        OutProps[NTypes]->SetTypeObj(TypeObject::line);
        //TODO: check the number of lines selected
        lines_str="all";
        SetLinesSelected(lines_str,OutProps[NTypes],lines);
        OutProps[NTypes]->SetUnits("[m]");
        NTypes++;
      }
    }
    NTypes=(unsigned)OutProps.size();
    eleOut=eleOut->NextSiblingElement();
  }
}
//==============================================================================
/// Initializes output file
//==============================================================================
void Output::Setup(std::string dir){
  std::string function="Setup";
  DataDir=dir;
}

//==============================================================================
/// Checks the string passed and store in OutProperties the lines selected
//==============================================================================
void Output::SetLinesSelected(string str,OutProperties *oProps,std::vector<Line*> lines) {
  //TODO: Complete for when the selection isn't "all" -> F.Ex: 1,2,3......
  for(unsigned nl=0;nl<NLines;nl++){ lines[nl]->SetDtOut(DtOut); }
  if(str=="all") {
    oProps->SetNLines(NLines);
    oProps->SetLines(lines);
  }
}

//==============================================================================
/// Sets the connection type selected
//==============================================================================
void Output::SetTypeOut(std::string str,OutProperties *oProps) {
  if(str=="all") { oProps->SetTypeCon(TypeConnection::all); }
  else if(str=="fixed") { oProps->SetTypeCon(TypeConnection::fixed); }
  else if(str=="vessel") { oProps->SetTypeCon(TypeConnection::vessel); }
}

//==============================================================================
/// Write the outputs for each line
//==============================================================================
void Output::SaveCsv(double timestep,double dt) {
  const double start=TimeStart;  //-Start time. | Tiempo de inicio.
  const double finish=TimeMax;   //-End time. | Tiempo de finalizacion.
  const double dtout=DtOut;      //-Save frequency. | Frecuencia de guardado.
  bool savedata=false;
  if(start<=timestep && timestep<=finish) savedata=!(timestep<(floor((timestep-dt)/dtout)+1.0)*dtout);

  if(savedata){
    std::string units="";
    const std::string path=DataDir+"/";
    if(!fun::DirExists(path))fun::MkdirPath(path);
    for(unsigned lf=0; lf<GetNTypes(); lf++){
      OutProperties *out_p=GetOutProps()[lf];
      //-Tensions
      if(out_p->GetMagnitude()==Magnitude::tension){
        const std::string file=path+fun::PrintStr("MoorDyn_tension.csv");
        jcsv::JSaveCsv2 scsv(file,true,0);
        if(!scsv.GetAppendMode()){
          //-Saves head.
          //-Guarda cabecera.
          scsv.SetHead();
          scsv << "time [s];dt [s]";
          for(unsigned l=0; l<out_p->GetNLines(); l++) {
            units=out_p->GetUnits();
            scsv << fun::PrintStr("Anchor_L%u %s",out_p->GetLines()[l]->GetNumber(),units.c_str());
            scsv << fun::PrintStr("Fairlead_L%u %s",out_p->GetLines()[l]->GetNumber(),units.c_str());
          }
          scsv << jcsv::Endl();
        }
        //-Saves data.
        //-Guarda datos.
        scsv.SetData();
        scsv << timestep << dt;
        for(unsigned l=0; l<out_p->GetNLines(); l++) {
          scsv << out_p->GetLines()[l]->GetTensionOutput(0);
          scsv << out_p->GetLines()[l]->GetTensionOutput(out_p->GetLines()[l]->GetN());
        }
        scsv << jcsv::Endl();
        scsv.SaveData();
      }

      //-Forces
      if(out_p->GetMagnitude()==Magnitude::force){
        const std::string file=path+fun::PrintStr("MoorDyn_force.csv");
        jcsv::JSaveCsv2 scsv(file,true,0);
        if(!scsv.GetAppendMode()){
          //-Saves head.
          //-Guarda cabecera.
          scsv.SetHead();
          scsv << "time [s];dt [s]";
          for(unsigned l=0; l<out_p->GetNLines(); l++) {
            units=out_p->GetUnits();
            scsv << fun::PrintStr("Anchor.x_L%u %s",out_p->GetLines()[l]->GetNumber(),units.c_str())
              << fun::PrintStr("Anchor.y_L%u %s",out_p->GetLines()[l]->GetNumber(),units.c_str())
              << fun::PrintStr("Anchor.z_L%u %s",out_p->GetLines()[l]->GetNumber(),units.c_str());
            scsv << fun::PrintStr("Fairlead.x_L%u %s",out_p->GetLines()[l]->GetNumber(),units.c_str())
              << fun::PrintStr("Fairlead.y_L%u %s",out_p->GetLines()[l]->GetNumber(),units.c_str())
              << fun::PrintStr("Fairlead.z_L%u %s",out_p->GetLines()[l]->GetNumber(),units.c_str());
          }
          scsv << jcsv::Endl();
        }
        //-Saves data.
        //-Guarda datos.
        scsv.SetData();
        scsv << timestep << dt;
        for(unsigned l=0; l<out_p->GetNLines(); l++) {
          tdouble3 out_val=out_p->GetLines()[l]->GetForceOutput(0);
          scsv << out_val.x << out_val.y << out_val.z;
          out_val=out_p->GetLines()[l]->GetForceOutput(out_p->GetLines()[l]->GetN());
          scsv << out_val.x << out_val.y << out_val.z;
        }
        scsv << jcsv::Endl();
        scsv.SaveData();
      }

      //-Velocities
      if(out_p->GetMagnitude()==Magnitude::velocity){
        const std::string file=path+fun::PrintStr("MoorDyn_velocity.csv");
        jcsv::JSaveCsv2 scsv(file,true,0);
        if(!scsv.GetAppendMode()){
          //-Saves head.
          //-Guarda cabecera.
          scsv.SetHead();
          scsv << "time [s];dt [s]";
          for(unsigned l=0; l<out_p->GetNLines(); l++) {
            units=out_p->GetUnits();
            scsv << fun::PrintStr("Anchor.x_L%u %s",out_p->GetLines()[l]->GetNumber(),units.c_str())
              << fun::PrintStr("Anchor.y_L%u %s",out_p->GetLines()[l]->GetNumber(),units.c_str())
              << fun::PrintStr("Anchor.z_L%u %s",out_p->GetLines()[l]->GetNumber(),units.c_str());
            scsv << fun::PrintStr("Fairlead.x_L%u %s",out_p->GetLines()[l]->GetNumber(),units.c_str())
              << fun::PrintStr("Fairlead.y_L%u %s",out_p->GetLines()[l]->GetNumber(),units.c_str())
              << fun::PrintStr("Fairlead.z_L%u %s",out_p->GetLines()[l]->GetNumber(),units.c_str());
          }
          scsv << jcsv::Endl();
        }
        //-Saves data.
        //-Guarda datos.
        scsv.SetData();
        scsv << timestep << dt;
        for(unsigned l=0; l<out_p->GetNLines(); l++) {
          tdouble3 out_val=out_p->GetLines()[l]->GetVelocityOutput(0);
          scsv << out_val.x  << out_val.y << out_val.z;
          out_val=out_p->GetLines()[l]->GetVelocityOutput(out_p->GetLines()[l]->GetN());
          scsv << out_val.x << out_val.y << out_val.z;
        }
        scsv << jcsv::Endl();
        scsv.SaveData();
      }

      //-Position
      if(out_p->GetMagnitude()==Magnitude::position){
        const std::string file=path+fun::PrintStr("MoorDyn_position.csv");
        jcsv::JSaveCsv2 scsv(file,true,0);
        if(!scsv.GetAppendMode()){
          //-Saves head.
          //-Guarda cabecera.
          scsv.SetHead();
          scsv << "time [s];dt [s]";
          for(unsigned l=0; l<out_p->GetNLines(); l++) {
            units=out_p->GetUnits();
            scsv << fun::PrintStr("Anchor.x_L%u %s",out_p->GetLines()[l]->GetNumber(),units.c_str())
              << fun::PrintStr("Anchor.y_L%u %s",out_p->GetLines()[l]->GetNumber(),units.c_str())
              << fun::PrintStr("Anchor.z_L%u %s",out_p->GetLines()[l]->GetNumber(),units.c_str());
            scsv << fun::PrintStr("Fairlead.x_L%u %s",out_p->GetLines()[l]->GetNumber(),units.c_str())
              << fun::PrintStr("Fairlead.y_L%u %s",out_p->GetLines()[l]->GetNumber(),units.c_str())
              << fun::PrintStr("Fairlead.z_L%u %s",out_p->GetLines()[l]->GetNumber(),units.c_str());
          }
          scsv << jcsv::Endl();
        }
        //-Saves data.
        //-Guarda datos.
        scsv.SetData();
        scsv << timestep << dt;
        for(unsigned l=0; l<out_p->GetNLines(); l++) {
          tdouble3 out_val=out_p->GetLines()[l]->GetPositionOutput(0);
          scsv << out_val.x << out_val.y  << out_val.z;
          out_val=out_p->GetLines()[l]->GetPositionOutput(out_p->GetLines()[l]->GetN());
          scsv << out_val.x  << out_val.y  << out_val.z;
        }
        scsv << jcsv::Endl();
        scsv.SaveData();
      }
    }
  }
}