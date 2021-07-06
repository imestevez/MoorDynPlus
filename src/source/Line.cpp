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

/// \file Line.cpp \brief Implements the class \ref Line.

#include "Line.h"
#include "Environment.h"
#include "Connection.h"
#include "Body.h"
#include "QSlines.h" // the c++ version of quasi-static model Catenary
#include "FunMoorDyn.h"

#include <string>


//*******************************************************************************
//*********************** LineDefault Section ***********************************
//*******************************************************************************

//==============================================================================
/// Constructor
//==============================================================================
LineDefault::LineDefault() {
  Reset();
  ClassName="LineDefault";
}
//==============================================================================
/// Destructor
//==============================================================================
LineDefault::~LineDefault() {
  //Reset();
  Ea=0;
  Diameter=0;
  Can=0;
  Cat=0;
  Cdn=0;
  Cdt=0;
  W=0;
}
//==============================================================================
/// Initialization of variables.
//==============================================================================
void LineDefault::Reset() {
  Ea=DBL_MAX;
  Diameter=DBL_MAX;
  Can=DBL_MAX;
  Cat=DBL_MAX;
  Cdn=DBL_MAX;
  Cdt=DBL_MAX;
  W=DBL_MAX;
  C=DBL_MAX;
  Channels="";
  BreakTension=0.0;
}

//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void LineDefault::LoadXml(JXml *sxml,const std::string &place) {
  std::string function="LoadXml";
  TiXmlNode* node=sxml->GetNode(place,false);
  if(!node )Run_Exceptioon("Cannot find the element \'"+place+"\'.");
  ReadXml(sxml,node->ToElement(),place);
}
//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void LineDefault::ReadXml(JXml *sxml,TiXmlElement* lis,const std::string &place) {
  //-Loads linedefault zone.
  TiXmlElement* elelDef=lis->FirstChildElement("linedefault");
    
  sxml->CheckElementNames(elelDef,false,"ea diameter massDenInAir ba can cat cdn cdt breaktension outputFlags");

  // If wants overwrite default parameters
  Ea=sxml->ReadElementDouble(elelDef,"ea","value",false);
  Diameter=sxml->ReadElementDouble(elelDef,"diameter","value",false);
  W=sxml->ReadElementDouble(elelDef,"massDenInAir","value",false);
  C=sxml->ReadElementDouble(elelDef,"ba","value",true,-0.8);
  Can=sxml->ReadElementDouble(elelDef,"can","value",true,1.0);
  Cat=sxml->ReadElementDouble(elelDef,"cat","value",true,0.0);
  Cdn=sxml->ReadElementDouble(elelDef,"cdn","value",true,1.6);
  Cdt=sxml->ReadElementDouble(elelDef,"cdt","value",true,0.05);
  BreakTension=sxml->ReadElementDouble(elelDef,"breaktension","value",true,0.0); 
  Channels=sxml->ReadElementStr(elelDef,"outputFlags","value",true,"-");
  CheckChannels(sxml,elelDef);
  
}
//==============================================================================
/// Check if the output flag is in the list of posibilities.
//==============================================================================
void LineDefault::CheckChannels(JXml *sxml,TiXmlElement *eleLine) {
  std::string function="CheckChannels";

  bool equals=false;
  const char *value=Channels.c_str(); //stores the Channels read into xml file
  const unsigned nChannels_in=(unsigned)strlen(value)/sizeof(char);
  const unsigned nChannels=(unsigned)sizeof(outputChannels)/sizeof(char);
  for(unsigned i=0; i<nChannels_in; i++) { //for each input channes insert by user
    for(unsigned j=0; j<nChannels; j++) { //for each char of nodes
      if(outputChannels[j]==value[i]) { //If both of them are equals
        equals=true;
        break;
      }
    }
    if(!equals) {
      string tex="Value of ouputsFlags tag isn't correctly. Checks the options. \'"; 
      tex=tex+sxml->ErrGetFileRow(eleLine)+"\'.";
      Run_Exceptioon(tex);
    }
    equals=false;
  }
  return;
}
//*******************************************************************************
//****************************** Line Section ***********************************
//*******************************************************************************

// here is the new numbering scheme (N segments per line)
//   [connect (node 0)]  --- segment 0 --- [ node 1 ] --- seg 1 --- [node2] --- ... --- seg n-2 --- [node n-1] --- seg n-1 ---  [connect (node N)]

//==============================================================================
/// Constructor
//==============================================================================
Line::Line(JLog2 *log):Log(log)  {
  Reset();
  Qslines=new QSlines(log);
  ClassName="Line";
}
//==============================================================================
/// Destructor.
//==============================================================================
Line::~Line() {
  Reset();  
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void Line::Reset() {
  Connections.clear();
  NConnections=0;   // total Number of Connection objects
  NFairs=0;      // Number of fairlead connections
  NAnchs=0;      // Number of anchor connections
  NConns=0;     // Number of "connect" connections
  WaveKin=0;    // start off with wave kinematics disabled.  Can be enabled after initial conditions are found and wave kinematics are calculated  
  StartIndex=0;
  ConnectFlag=false;
  FairFlag=false;
  OutputLine=false;
  OutputNode=false;
  BodyId=NULL;
  FileName=NULL;
  Outfile=NULL;
  Depth=NULL;
  Env=NULL;
  AnchConnect=NULL;
  FairConnect=NULL;
  Qslines=NULL;
  BreakTension=0.0;
  Disabled=false;
  DtOut=0.0;
}
//==============================================================================
/// Allocates memory.
//==============================================================================
void Line::AllocateMemory() {
  std::string function="AllocateMemory";
  try {
    Outfile=new ofstream;
  //  Env=new Environment;
    AnchConnect=new Connection;
    FairConnect=new Connection;
    Depth=new double;
    BodyId=new unsigned;
  }
  catch (const std::bad_alloc) {
    Run_Exceptioon("Could not allocate the requested memory.");
  }
  std::memset(Outfile,0,sizeof(ofstream));
  //std::memset(Env,0,sizeof(Environment));
  std::memset(AnchConnect,0,sizeof(Connection));
  std::memset(FairConnect,0,sizeof(Connection));
  std::memset(Depth,0,sizeof(double));
  std::memset(BodyId,0,sizeof(unsigned));
  
}

//==============================================================================
/// Allocates memory of vectors
//==============================================================================
void Line::AllocateMemoryVectors() {
  // =============== size vectors =========================
  R.resize(N+1,vector<double>(3,0.0));    // node positions [i][x/y/z]
  Rd.resize(N+1,vector<double>(3,0.0));    // node velocities [i][x/y/z]
  Q.resize(N+1,vector<double>(3,0.0));   //< unit tangent vectors for each node

  // forces 
  T.resize(N,vector<double>(3,0.0));  // line tensions
  Td.resize(N,vector<double>(3,0.0));   // line damping forces
    //  Tmag.resize(N,0.0);        // segment tension magnitudes << hardly used
  W.resize(N+1,vector<double>(3,0.0));  // node weights

  Dp.resize(N+1,vector<double>(3,0.0));    // node drag (transverse)
  Dq.resize(N+1,vector<double>(3,0.0));    // node drag (axial)
  Ap.resize(N+1,vector<double>(3,0.0));    // node added mass forcing (transverse)
  Aq.resize(N+1,vector<double>(3,0.0));    // node added mass forcing (axial)
  B.resize(N+1,vector<double>(3,0.0));    // node bottom contact force
  Fnet.resize(N+1,vector<double>(3,0.0));  // total force on node

  S.resize(N+1,vector<vector< double>>(3,vector<double>(3,0.0)));  // inverse mass matrices (3x3) for each node
  M.resize(N+1,vector<vector< double>>(3,vector<double>(3,0.0)));  // mass matrices (3x3) for each node

  L.resize(N,0.0);       // line unstretched segment lengths
  Lstr.resize(N,0.0);     // stretched lengths
  Ldstr.resize(N,0.0);     // rate of stretch
  V.resize(N,0.0);      // volume?

  Zeta.resize(N+1,0.0);          // wave elevation above each node
  F.resize(N+1,0.0);   // fixed 2014-12-07  // VOF scalar for each NODE (mean of two half adjacent segments) (1=fully submerged,0=out of water)
  U.resize(N+1,vector<double>(3,0.));     // wave velocities
  Ud.resize(N+1,vector<double>(3,0.));;   // wave accelerations
}
//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void Line::LoadXml(JXml *sxml,const std::string &place,TiXmlElement* eleL,const unsigned numLine,std::vector<Connection*> connects,LineDefault *lineDef) {
  std::string function="LoadXml";  
  //TiXmlNode* node=sxml->GetNode(place,false);
  //if(!node)Run_Exceptioon("Cannot find the element \'"+place+"\'.");
  if(!eleL)Run_Exceptioon("Cannot find the element \'"+place+"\'.");
  ReadXml(sxml,place,eleL,numLine,connects,lineDef);
}
//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void Line::ReadXml(JXml *sxml,const std::string &place_in,TiXmlElement* eleLine,const unsigned numLine,std::vector<Connection*> connects,LineDefault *lineDef) {
  std::string function="ReadXml";
  std::string place=place_in;
  sxml->CheckElementNames(eleLine,false,"ea diameter massDenInAir ba can cat cdn cdt breaktension outputFlags length fixconnection vesselconnection segments connect");

  //-If there is line tag
  if(eleLine) { 
    //-If there are default parameters
    if(lineDef) {
      D=lineDef->Diameter;
      E=lineDef->Ea;
      C=lineDef->C;
      Can=lineDef->Can;
      Cat=lineDef->Cat;
      Cdn=lineDef->Cdn;
      Cdt=lineDef->Cdt;
      Channels=lineDef->Channels;
      MassDenInAir=lineDef->W;
      BreakTension=lineDef->BreakTension;
    }
    //-If exists,then overwrites the default parameters
    if(sxml->ExistsElement(eleLine,"ea"))           { E=sxml->ReadElementDouble(eleLine,"ea","value",false); }
    if(sxml->ExistsElement(eleLine,"massDenInAir")) { MassDenInAir=sxml->ReadElementDouble(eleLine,"massDenInAir","value",false); }
    if(sxml->ExistsElement(eleLine,"diameter"))     { D=sxml->ReadElementDouble(eleLine,"diameter","value",false); }
    if(sxml->ExistsElement(eleLine,"ba"))           { C=sxml->ReadElementDouble(eleLine,"ba","value",false); }
    if(sxml->ExistsElement(eleLine,"can"))          { Can=sxml->ReadElementDouble(eleLine,"can","value",false); }
    if(sxml->ExistsElement(eleLine,"cat"))          { Cat=sxml->ReadElementDouble(eleLine,"cat","value",false); }
    if(sxml->ExistsElement(eleLine,"cdn"))          { Cdn=sxml->ReadElementDouble(eleLine,"cdn","value",false); }
    if(sxml->ExistsElement(eleLine,"cdt"))          { Cdt=sxml->ReadElementDouble(eleLine,"cdt","value",false); }
    if(sxml->ExistsElement(eleLine,"outputFlags"))  { Channels=sxml->ReadElementStr(eleLine,"outputFlags","value",false); }
    if(sxml->ExistsElement(eleLine,"breaktension")) { BreakTension=sxml->ReadElementDouble(eleLine,"breaktension","value",false); }

    //-Checks if the output Channels are rigth
    CheckChannels(sxml,eleLine);

    //-Stores the elements of the line
    OutputLine=sxml->GetAttributeBool(eleLine,"savedata",true,true);
    Number=numLine; // stores a secuencial Number of lines in the program for each body
    UnstrLen=sxml->ReadElementDouble(eleLine,"length","value",false);
    N=sxml->ReadElementInt(eleLine,"segments","value",false);

    //PRINT: look num of line executed
    if(Debug) { printf("\tReadLine()-num: %d\n",Number); }
    AllocateMemory();//Calls to memory reserve  
    //-Reads the connections of the Line
    vector<Connection*> connections_temp;
    //-Fixed connection
    place=place_in+".fixconnection";
    TiXmlElement* elecF=NULL;
    if(sxml->ExistsElement(eleLine,"fixconnection")) {
      elecF=eleLine->FirstChildElement("fixconnection");
      if(elecF) {
        connections_temp.push_back(new Connection()); //Pointer to Fix connection
        connections_temp[NConnections]->LoadXml(sxml,place,elecF,(unsigned)connections_temp.size());//Build Fix Connection
        NAnchs++;
        NConnections++;
      }
    }

    //-Vessel connection
    place=place_in+".vesselconnection";
    TiXmlElement* elecV=NULL;
    if(sxml->ExistsElement(eleLine,"vesselconnection")) {  
      elecV=eleLine->FirstChildElement("vesselconnection");
      if(elecV) {
        connections_temp.push_back(new Connection()); //Pointer to Vessel connection
        connections_temp[NConnections]->LoadXml(sxml,place,elecV,(unsigned)connections_temp.size()); //Build VesselConnection
        FairFlag=true;
        NFairs++;
        NConnections++;
      }
    }
    //-Searches if there is some connect
    place=place_in+".connect";
    if(sxml->ExistsElement(eleLine,"connect")) {  
      //-Connect connection
      TiXmlElement* elecC=eleLine->FirstChildElement("connect");
      if(elecC) {
        if(connects.size()>0) {
          unsigned ref=sxml->GetAttributeInt(elecC,"conref",false);
          Connection *conn_temp=SearchConnect(ref,connects);
          connections_temp.push_back(conn_temp);
          ConnectFlag=true;
          NConns++;
          NConnections++;
        }else{Run_Exceptioon("Error: There aren't any connect Connections declarated into the \"connnects\" block.");}
      }
    }
    //-Searches if there is the 2nd vessel
    unsigned nVessel=sxml->CountElements(eleLine,"vesselconnection");
    if(nVessel>1){
      place=place_in+".vesselconnection";
      if(sxml->ExistsElement(eleLine,"vesselconnection")) {  
        TiXmlElement *elecV2=sxml->GetNextElement(elecV,"vesselconnection"); //- CHECK -TODO: Two Vessel
        //-2nd Vessel connection
        if(elecV2) {
          connections_temp.push_back(new Connection()); //Pointer to Vessel connection
          connections_temp[NConnections]->LoadXml(sxml,place,elecV2,(unsigned)connections_temp.size()); //Build VesselConnection
          FairFlag=true;
          NFairs++;
          NConnections++;
        }
      }
    }
    //-Stores the number of connections readed
    NConnections=(unsigned)connections_temp.size();
    if(NConnections!= 2) { Run_Exceptioon("Error: Just 2 connections are allowed. Found: "+funmd::ToString(NConnections)+"."); }
    ConfigLine(connections_temp);
    connections_temp.clear();

    //PRINT: look Number and types of Connections
    if(Debug) {
      for(unsigned i=0; i<NConnections; i++) {
        printf("\t\tconnection: %d Type: %d\n",connections_temp[i]->GetNumber(),connections_temp[i]->GetType());
        printf("\t\t\tx= %.4f ; y= %.4f ; z= %.4f\n",connections_temp[i]->GetX(),connections_temp[i]->GetY(),connections_temp[i]->GetZ());
      }
    }
    CheckInputParameters();
  }
}
//==============================================================================
/// Configures the start and end Connection
//==============================================================================
void Line::ConfigLine(vector<Connection*> connections){
  std::string function="ConfigLine";
  if((connections[0]->GetType()==0) 
  || (connections[0]->GetType()==1 && connections[1]->GetType()==1)  ){ //If is an Anchor or Two fairs
    AnchConnect=connections[0];
    FairConnect=connections[1];
  }else if(connections[0]->GetType()==1){//If is a Fairlead and Connect
    AnchConnect=connections[1];
    FairConnect=connections[0];
  }
  else{Run_Exceptioon("Error: The First Connection of the Line "+funmd::ToString(Number)+" must be fixed of vessel.");}
  
  //-Sets the id number of connection
  //AnchConnect->SetNumber(0);  
  //FairConnect->SetNumber(1);  

  //-Sets the node position
  AnchConnect->SetNode(0);
  FairConnect->SetNode(N);
  //-Stores the Connections
  Connections.clear();
  Connections.push_back(AnchConnect); //assign Start Connection to the array
  Connections.push_back(FairConnect); //assign End Connection to the array

}
//====================================================================================
/// Checks if the connect Connection selected in Lines was configurated previously
//====================================================================================
Connection *Line::SearchConnect(const unsigned ref,std::vector<Connection*> connects){
  std::string function="SearchConnect";
  
  for(unsigned c=0; c<connects.size(); c++){
    if(connects[c]->GetRef()==ref){
      return connects[c];
    }
  }
  Run_Exceptioon("Error: There aren't connconnection with ref: "+funmd::ToString(ref)+" in linedefault block." );
  return NULL  ;
}
//==============================================================================
/// Check if the output flag is in the list of posibilities.
//==============================================================================
void Line::CheckChannels(JXml *sxml,TiXmlElement *eleLine) {
  std::string function="CheckChannels";
  OutputNode=true;
  bool equals=false;
  const char *value=Channels.c_str(); //stores the Channels read into xml file
  const unsigned nChannels_in=(unsigned)strlen(value)/sizeof(char);
  const unsigned nChannels=(unsigned)sizeof(outputChannels)/sizeof(char);
  for(unsigned i=0; i<nChannels_in; i++) { //for each input channes insert by user
    for(unsigned j=0; j<nChannels; j++) { //for each char of nodes
      if(outputChannels[j]==value[i]) { //If both of them are equals
        equals=true;
        break;
      }
    }
    if(!equals) {
      string tex="Value of ouputsFlags tag isn't correctly. Checks the options. \'"; 
      tex=tex+sxml->ErrGetFileRow(eleLine)+"\'.";
      Run_Exceptioon(tex);
    }
    equals=false;
  }

  if(nChannels_in==1&&value[0]=='-'){OutputNode=false; }

  return;
}
//==============================================================================
/// Checks if the input parameters are into the range
//==============================================================================
void Line::CheckInputParameters() {
  std::string function="CheckInputParameters";
  
  //-Checks distance between 2 points
  double distance=sqrt((pow((abs(AnchConnect->GetX())-abs(FairConnect->GetX())),2))
   +(pow((abs(AnchConnect->GetY())-abs(FairConnect->GetY())),2))
   +(pow((abs(AnchConnect->GetZ())-abs(FairConnect->GetZ())),2))
  ); //distance between nodeAnch and nodeFair
  
  if(UnstrLen<distance) { 
    string tex="The unstretched length of the line "+funmd::ToString(Number) 
     +" should be at least of "+funmd::ToString(distance)+"m.";
    //Run_Exceptioon(tex); 
    Log->PrintfWarning("%s",tex.c_str());
  }
  return;
}
//==============================================================================
/// Returns a pointer to Anch connect
//==============================================================================
Connection* Line::GetAnchConnect() {
  return AnchConnect;
}
//==============================================================================
/// Returns a pointer to Fair connect
//==============================================================================
Connection* Line::GetFairConnect() {
  return FairConnect;
}
//==============================================================================
/// Returns a pointer to Connections by type
//==============================================================================
std::vector<Connection*> Line::GetConnections(const unsigned type) {
  std::vector<Connection*> conns;
  for(unsigned c=0; c<NConnections; c++) {
    if(Connections[c]->GetType()==type) {
      conns.push_back(Connections[c]);
    }
  }
  return conns;
}
//==============================================================================
/// Returns the array of Connections by type
//==============================================================================
std::vector<Connection*> Line::GetFairs(const unsigned ref) {
  const unsigned type=1;//Fairs
  std::vector<Connection*> conns;
  for(unsigned c=0; c<NConnections; c++) {
    if( (Connections[c]->GetType()==type)
      &&(Connections[c]->GetRef()==ref) ) {
      conns.push_back(Connections[c]);
    }
  }
  return conns;
}
//==============================================================================
/// Returns a pointer to Connecs connect
//==============================================================================
std::vector<Connection*> Line::GetConnConnects() {
  unsigned type=2;
  std::vector<Connection*> conns;
  for(unsigned c=0; c<NConnections; c++) {
    if(Connections[c]->GetType()==type) {
      conns.push_back(Connections[c]);
    }
  }
  return conns;
}

//==============================================================================
/// Returns the number of Connects by reference
//=============================================================================
unsigned Line::GetNConns(const unsigned ref){
  std::vector<Connection*> conns=GetConnConnects();
  int nc=0;
  for(unsigned c=0; c<conns.size(); c++) {
    if(conns[c]->GetRef()==ref) {nc++;}
  }
  return nc;
}
//==============================================================================
/// Returns the "Connects" by reference
//==============================================================================
std::vector<Connection*> Line::GetConnConnects(const unsigned ref) {
  std::vector<Connection*> conns;
  for(unsigned c=0; c<NConnections; c++) {
    if(Connections[c]->GetRef()==ref) {conns.push_back(Connections[c]);}
  }
  return conns;
}
//==============================================================================
/// Check if the relationship between dtM, N and EA is rigth
//==============================================================================
void Line::CheckDtMRelationship() {
  //-Gets Parameters
  double NdtM=(int)ceil(Env->GetICdt()/Env->GetDtM0());   // Number of body model time steps per outer time step
  double dtM=Env->GetICdt()/NdtM;    // body model time step size (s)
  double weight=((Rho-Env->GetRho_w())*(PI/4.*D*D))*9.81; //Weigth per unit length
  double frecuency=(2*N/UnstrLen*sqrt(E/weight)*0.159155)/10; //Converts rad/s to Hz and divided by 10 for that the relationship was rigth
  double stepLimit=(1/dtM); //Used for knows the limit of the frecuency.
  //-Checks if the DtMRelationship is rigth
  if(frecuency>stepLimit) {
    dtM=(UnstrLen/(2*N *sqrt(E/weight)*0.159155))*10; // Multiply by 10 for the rigth relationship
    NdtM=(int)ceil(Env->GetICdt()/dtM);
    double dtm0_temp=Env->GetICdt()/NdtM;  // Number of body model time steps per outer time step
    //-Round the dtm to 0.0{0,n}1 values
    int exp=(int)(round(log10(dtm0_temp))); //Used for know the value of exponent 10^exp
    dtm0_temp=pow(10,exp)*0.1; //Multiply to make more small the dtm Value. 
    if(dtm0_temp<Env->GetDtM0()) { //If the new dtM calculated is more small than the initial
      Env->SetDtM0(dtm0_temp); //Store the new dtm value
      if(Debug) { Log->PrintfWarning("Adjusted dtM to %g for fairlead calcs.",Env->GetDtM0()); }
    }
  }
}
//==============================================================================
///  Makes the setup for this Line. Initializes this and creates output file
//==============================================================================
void Line::Setup(std::string &dir,Environment *env_in,Body *body_in) {
  Env=env_in;
  Depth=env_in->GetWtrDpth();
  BodyId=NULL;
  //-If this line is connected to a body
  if(body_in){BodyId=body_in->GetPtrRef(); Depth=body_in->GetDepth();}

  //PRINT: look when Line::Setup() is executed
  if(Debug) {printf("\tLine::setup()-Number: %d\n",Number);}

  //-Configures Output file
  if(OutputNode){
    stringstream oname;
    oname << dir;
    //-Checks if this line is conneected to a body
    if(BodyId){  oname << "Body"  << std::setw(2) << std::setfill('0') << *BodyId ;}
    else{ oname<<"NoBody"; }
    oname << "_Line" << std::setw(2) << std::setfill('0') << Number << ".csv"; //Stores in Dir directory. format BodyXX_LineXX.csv
    FileName=new string(oname.str());
  }else Outfile=NULL;
  // ================== set up properties ===========  

  WaveKin=0;  // start off with wave kinematics disabled.  Can be enabled after initial conditions are found and wave kinematics are calculated
  
  //-Configures Line parameters
  Rho=MassDenInAir/(PI/4.*D*D);
  E=E/(PI/4.*D*D);

  double temp_c=C; //Stores BA/-Zeta – the line internal damping (Ns)
  C=C/(PI/4.*D*D);

  CheckDtMRelationship();
  AllocateMemoryVectors();

  //- Automatic internal damping option (if negative BA provided,as damping ratio)
  if(temp_c<0) {
    double Zeta=-temp_c; // desired damping ratio
    C=Zeta*UnstrLen/N*sqrt(E*Rho);   // Rho=w/A
    if(misc::wordy>1) { printf( "\n   Line %u damping set to %.f Ns.\n",Number,C);} //wordy is in Misc.h
  }

  for(unsigned i=0; i<N; i++) {
    L[i]=UnstrLen/double(N);  // distribute line length evenly over segments
    V[i]=L[i]*0.25*PI*D;
  }

  return;
};

//==============================================================================
/// Initialize wave parameters for no waves situation
//==============================================================================
void Line::SetupWaves(double *Ucurrent_in,float dt_in) {
  Ucurrent=Ucurrent_in;
  WaveDT=dt_in; // new variable for wave time step (should be same as WaveDT I think...)  

  if(Env->GetWaveKin()>0) {  Log->PrintWarning("Dummy function called when waves are supposed to be enabled!");  }

  if(misc::wordy>1) { printf("\n    Setting up wave variables for Line %u !  ---------------------\n",Number);  }
  if(misc::wordy>1) { printf("\n   Nt= %d, and WaveDT= %.f, Depth= %.f\n",Nt,WaveDT,*Depth );}

  WGNC_Fact=1.0;
  S2Sd_Fact=1.0;

  Nw=0;    // use no components since not worrying about wave kinematics
  Nt=2; // this is a new variable containing the Number of wave time steps to be calculated

  WGNC_Fact=1.0;
  S2Sd_Fact=1.0;
  // resize the new time series vectors
  zetaTS.resize(N+1,vector<double>(Nt,0.));
  FTS.resize(N+1,vector<double>(Nt,0.));
  UTS.resize(N+1,vector<vector< double>>(Nt,vector<double>(3,0.)));
  UdTS.resize(N+1,vector<vector< double>>(Nt,vector<double>(3,0.)));
  tTS.resize(Nt,0.);
  return;
}
//==============================================================================
/// Get ICs for line using quasi-static approach
//==============================================================================
void Line::Initialize(double* X) {
  std::string function="Initialize";
  //-Creates new file if the outputFlags was selected into xml
  if(OutputNode){
    Outfile=new ofstream(*FileName,std::ofstream::out);
    Outfile->clear();
  }
  //-Checks it's not null.  Null signals no individual line output files
  if(Outfile) { 
    if(Outfile->is_open()) {
      // ------------- write channel names line --------------------
      // output time
      unsigned nChannels=(unsigned)Channels.length();
      if((nChannels==1) && Channels[0]=='-') Outfile->close(); //If there aren't output Channels selected
      else{
        *Outfile << "Time;";
        for(unsigned i=0; i<nChannels; i++) {
          char ch=Channels[i]; //takes a character
          if(ch!='-') { InitializeOutputs(ch); } //calls to store funtion
        }
        *Outfile << "\n";
        Outfile->close();
      }
    }
    else { Run_Exceptioon("Cannot open file ouput \'"+*FileName+"\'\n"); }
  }


  // set end node positions and velocities from connect objects
  //TODO: check these functions
  AnchConnect->GetConnectState(R[0],Rd[0]);
  FairConnect->GetConnectState(R[N],Rd[N]);

  //  cout << "IT: " << zmin << ";  X= "<< R[0][0]<< " ; Y= "<< R[0][1] << " ; Z= "<< R[0][2] << "-I: " <<0 <<endl ;

  double depth=*Depth;

  if(-(depth)>R[0][2]) { 
    string tex="Error: water depth is shallower than the anchor heigth in Line "+funmd::ToString(Number)+".";
    Run_Exceptioon(tex); 
  }

  // try to calculate initial line profile using catenary routine (from FAST v.7)
  // note: much of this function is adapted from the FAST source code
  //*
  // input variables for the Catenary function
  double XF=sqrt(pow((R[N][0]-R[0][0]),2.0)+pow((R[N][1]-R[0][1]),2.0)); // quasi-static body line coordinate system (vertical plane with corners at anchor and fairlead) 
  double ZF=R[N][2]-R[0][2];
  double W=((Rho-Env->GetRho_w())*(PI/4.*D*D))*9.81;
  double CB=0.;
  double Tol=0.00001;

  vector<double> snodes(N+1,0.0);             // locations of line nodes along line length-evenly distributed here 
  for(unsigned i=1; i<=N; i++) { snodes[i]=snodes[i-1]+L[i-1]; }
  snodes[N]=UnstrLen;                 // double check to ensure the last node does not surpass the line length

  // output variables
  double HF,VF,HA,VA,COSPhi,SINPhi;
  vector<double> Xl(N+1,0.0); // x location of line nodes
  vector<double> Zl(N+1,0.0);
  vector<double> Te(N+1,0.0);

  if(XF==0.0) {  // if the current body line is exactly vertical; thus, the solution below is ill-conditioned because the orientation is undefined; so set it such that the tensions and nodal positions are only vertical
    COSPhi=0.0;   SINPhi=0.0;
  }
  else {   // The current body line must not be vertical; use simple trigonometry
    COSPhi=(R[N][0]-R[0][0])/XF;
    SINPhi=(R[N][1]-R[0][1])/XF;
  }

  int success=Qslines->Catenary(XF,ZF,UnstrLen,E*PI/4.*D*D,W,CB,Tol,&HF,&VF,&HA,&VA,N,snodes,Xl,Zl,Te,Number);

  if(success >= 0) {  // assign the resulting line positions to the model
    for(unsigned i=1; i<N; i++) {
      R[i][0]=R[0][0]+Xl[i]*COSPhi;
      R[i][1]=R[0][1]+Xl[i]*SINPhi;
      R[i][2]=R[0][2]+Zl[i];
    }
  }
  else {  // otherwise just stretch the nodes between the endpoints linearly and hope for the best
    if(misc::wordy>0){  printf("\n   Catenary IC gen failed for Line %u, so using linear node spacing.",Number);}
    for(unsigned i=1; i<N; i++) {
      R[i][0]=R[0][0]+(R[N][0]-R[0][0])*(float(i)/float(N));
      R[i][1]=R[0][1]+(R[N][1]-R[0][1])*(float(i)/float(N));
      R[i][2]=R[0][2]+(R[N][2]-R[0][2])*(float(i)/float(N));
    }
  }
  // also assign the resulting internal node positions to the integrator initial state vector! (velocities leave at 0)
  for(unsigned i=1; i<N; i++) {
    for(int J=0; J<3; J++) {
      X[3*N-3+3*i-3+J]=R[i][J];  // positions
      X[3*i-3+J]=0.0;       // velocities=0
    }
  }
  // now we need to return to the integrator for the dynamic relaxation stuff
  return;
};
//==============================================================================
/// Smart (selective) function to get tension at any node including fairlead or anchor 
/// (accounting for weight in these latter cases) 
//==============================================================================
double Line::GetNodeTen(int i) {
  double NodeTen=0.0;

  if       (i==0) {  NodeTen=sqrt(Fnet[i][0]*Fnet[i][0]+Fnet[i][1]*Fnet[i][1]+(Fnet[i][2]+M[i][0][0]*(-Env->GetG()))*(Fnet[i][2]+M[i][0][0]*(-Env->GetG())));}
  else if(i==N) {  NodeTen=sqrt(Fnet[i][0]*Fnet[i][0]+Fnet[i][1]*Fnet[i][1]+(Fnet[i][2]+M[i][0][0]*(-Env->GetG()))*(Fnet[i][2]+M[i][0][0]*(-Env->GetG())));}
  else {
    double Tmag_squared=0.;
    for(int J=0; J<3; J++)  Tmag_squared+=0.25*(T[i][J]+T[i-1][J])*(T[i][J]+T[i-1][J]);  // take average of tension in adjacent segments 
    NodeTen=sqrt(Tmag_squared);    //     previously used: NodeTen=0.5*(Tmag[i-1]+Tmag[i]); // should add damping in here too <<<<<<<<<<<<<
  }
  return NodeTen;
}
//==============================================================================
/// Function to get position of any node along the line
//==============================================================================
unsigned Line::GetNodePos(unsigned NodeNum,double pos[3]) {
  if((NodeNum >= 0) && (NodeNum<=N)) {
    for(int i=0; i<3; i++) { pos[i]=R[NodeNum][i]; }
    return 0;
  }
  else { return -1; } // indicate an error
}
//==============================================================================
/// FASTv7 style line tension outputs
//==============================================================================
void Line::GetFASTtens(float* FairHTen,float* FairVTen,float* AnchHTen,float* AnchVTen) {
  *FairHTen=(float)sqrt(Fnet[N][0]*Fnet[N][0]+Fnet[N][1]*Fnet[N][1]);
  *FairVTen=(float)(Fnet[N][2]+M[N][0][0]*(-Env->GetG()));
  *AnchHTen=(float)sqrt(Fnet[0][0]*Fnet[0][0]+Fnet[0][1]*Fnet[0][1]);
  *AnchVTen=(float)(Fnet[0][2]+M[0][0][0]*(-Env->GetG()));
  return;
}
//==============================================================================
/// Writes Forces and Mass into the vectors, which are pased how parameters
//==============================================================================
void Line::GetAnchStuff(vector<double> &Fnet_out,vector<vector<double>>&M_out) {
  for(int I=0; I<3; I++) {
    Fnet_out[I]=Fnet[0][I];
    for(int J=0; J<3; J++)   M_out[I][J]=M[0][I][J];
  }
  return;
}
//==============================================================================
/// Writes Forces and Mass into the vectors, which are pased how parameters
//==============================================================================
void Line::GetFairStuff(vector<double> &Fnet_out,vector<vector<double>>&M_out) {
  for(int I=0; I<3; I++) {
    Fnet_out[I]=Fnet[N][I];
    for(int J=0; J<3; J++)   M_out[I][J]=M[N][I][J];
  }
  return;
}
//==============================================================================
/// Returns the position (X,Y,Z) of one Connection of this Line
//==============================================================================
tdouble3 Line::GetPositionOutput(const unsigned n) {
  return TDouble3(R[n][0],R[n][1],R[n][2]);
}
//==============================================================================
/// Returns the velocity (X,Y,Z)
//==============================================================================
tdouble3 Line::GetVelocityOutput(const unsigned n) {
  return TDouble3(Rd[n][0],Rd[n][1],Rd[n][2]);
}
//==============================================================================
/// Returns the force (X,Y,Z)
//==============================================================================
tdouble3 Line::GetForceOutput(const unsigned n) {
  return TDouble3(Fnet[n][0],Fnet[n][1],Fnet[n][2]);
}
//==============================================================================
/// Returns the Tension
//==============================================================================
double Line::GetTensionOutput(const unsigned n) {
  return GetNodeTen(n);
}
//==============================================================================
/// Function for boosting drag coefficients during IC generation  
//==============================================================================
void Line::scaleDrag(double scaler) {
  Cdn=Cdn*scaler;
  Cdt=Cdt*scaler;
  return;
}
//==============================================================================
/// Function to reset time after IC generation
//==============================================================================
void Line::setTime(double time) {
  t=time;
  return;
}
//==============================================================================
// This is the big function that updates the States
//==============================================================================
void Line::DoRHS(const double* X,double* Xd,const double time,const double dt) {
  std::string function="DoRHS";

  t=time;

  AnchConnect->GetConnectState(R[0],Rd[0]);
  FairConnect->GetConnectState(R[N],Rd[N]);

  //-Set interior node positions and velocities
  for(unsigned i=1; i<N; i++) {
    for(int J=0; J<3; J++) {
      R[i][J]=X[3*N-3+3*i-3+J]; // get positions
      Rd[i][J]=X[3*i-3+J]; // get velocities
    }
  }

  //-Calculate current (Stretched) segment lengths
  for(unsigned i=0; i<N; i++) {
    double lstr_squared=0.0;
    for(int J=0; J<3; J++) {lstr_squared+=(R[i+1][J]-R[i][J])*(R[i+1][J]-R[i][J]);  }
    Lstr[i]=sqrt(lstr_squared);   // stretched segment length

    double ldstr_top=0.0;

    for(int J=0; J<3; J++) {
      ldstr_top+=(R[i+1][J]-R[i][J])*(Rd[i+1][J]-Rd[i][J]);
      if(std::isnan(ldstr_top)) { 
        string tex="Error:  NaN value detected in MoorDyn state at dynamic relaxation time "
         +funmd::ToString(time)+" s. Try to reduce the EA or massDenInAir or increase the segments in line "+funmd::ToString(Number)+".";
        Run_Exceptioon(tex); 
      }

    }
    if(Lstr[i]==double(0)) { 
      string tex="Error: Lstr["+funmd::ToString(i)+"]=="         
       +funmd::ToString(Lstr[i])+". lstr_squared= "+funmd::ToString(lstr_squared) 
       +" sqrt= "+funmd::ToString(sqrt(lstr_squared))+" Cannot divide by 0. Check line Number " 
       +funmd::ToString(Number)+" of body "+funmd::ToString(*BodyId)+".";
      Run_Exceptioon(tex); 
    }
    Ldstr[i]=ldstr_top/Lstr[i];   // strain rate of segment
    V[i]=PI/4. *(D*D*L[i]);    // volume attributed to segment
  }
  //-Calculate unit tangent vectors (Q) for each node (including ends)  note: I think these are pointing toward 0 rather than N!
  for(unsigned i=0; i<=N; i++) {
    if      (i==0) { misc::unitvector(Q[i],R[i+1],R[i]);   } // compute unit vector Q
    else if(i==N) { misc::unitvector(Q[i],R[i],R[i-1]);   } // compute unit vector Q
    else           { misc::unitvector(Q[i],R[i+1],R[i-1]); } // compute unit vector Q ... using adjacent two nodes!
  }
  //============================================================================================
  // --------------------------------- apply wave kinematics ------------------------------------
  if(WaveKin==0) {   // if Wave Kinematics haven't been calculated.   ...this is a local Line switch (so wave kinematics can be enabled/disabled for individual lines)
    for(unsigned i=0; i<=N; i++) {
      Zeta[i]=0.0;  F[i]=1.0;
      for(int J=0; J<3; J++) { U[i][J]=0.0;  Ud[i][J]=0.0;}
    }
  }else {  // if wave kinematics time series have been precalculated.
     // =========== obtain (precalculated) wave kinematics at current time instant ============
     // get precalculated wave kinematics at previously-defined node positions for time instant t
     // get interpolation constant and wave time step index
    double frac;
    for(int ts=0; ts<Nt-1; ts++) {      // loop through precalculated wave time series time steps  (start ts at Ts0 to save time)
      if(tTS[ts+1]>t) {         // moving precalculated time bracked "up".  Stop once upper end of bracket is greater than current time t.
        Ts0=ts;
        frac=(t-tTS[ts])/(tTS[ts+1]-tTS[ts]);
        break;
      }
    }

    //-Loop through nodes 
    for(unsigned i=0; i<=N; i++) {
      Zeta[i]=zetaTS[i][Ts0]+frac*(zetaTS[i][Ts0+1]-zetaTS[i][Ts0]);
      F[i]=1.0;   // FTS[i][Ts0]+frac*(FTS[i][Ts0+1]-FTS[i][Ts0]);
      for(int J=0; J<3; J++) {
        U[i][J]=UTS[i][Ts0][J]+frac*(UTS[i][Ts0+1][J]-UTS[i][Ts0][J]);
        Ud[i][J]=UdTS[i][Ts0][J]+frac*(UdTS[i][Ts0+1][J]-UdTS[i][Ts0][J]);
      }
    }
    if(misc::wordy>2) {if(Number==1) { printf("\n t=%.f, U[4][0]= %.f\n",t,U[4][0]); }  }
  }

  CalculateMassMatrix(); //Calls calculate Mass matrix

  // ============  CALCULATE FORCES ON EACH NODE ===============================
  // loop through the segments
  for(unsigned i=0; i<N; i++) {
  // line tension
    if(Lstr[i]/L[i]>1.0) {  for(int J=0; J<3; J++) { T[i][J]=E*PI/4.*D*D* (1./L[i]-1./Lstr[i])*(R[i+1][J]-R[i][J]); }}
    else                     {  for(int J=0; J<3; J++) { T[i][J]=0.; } }// cable can't "push"
    // line internal damping force
    for(int J=0; J<3; J++) {   Td[i][J]=C*PI/4.*D*D* (Ldstr[i]/L[i])*(R[i+1][J]-R[i][J])/Lstr[i]; }
  }

  // loop through the nodes
  for(unsigned i=0; i<=N; i++) {
  // submerged weight (including buoyancy)
    if      (i==0) {  W[i][2]=PI/8.*(D*D*L[i]*(Rho-F[i]*Env->GetRho_w()))*(-Env->GetG()); }
    else if(i==N) {   W[i][2]=PI/8.*(D*D*L[i-1]*(Rho-F[i-1]*Env->GetRho_w()))*(-Env->GetG()); } // missing the "W[i][2] =" previously!
    else           {  W[i][2]=PI/8.*(D*D*L[i]*(Rho-F[i]*Env->GetRho_w())+D*D*L[i-1]*(Rho-F[i-1]*Env->GetRho_w()))*(-Env->GetG()); }
    // flow velocity calculations       
    double vq_squared=0.;
    double vp_squared=0.;
    for(int J=0; J<3; J++) { vi[J]=U[i][J]-Rd[i][J]; }// relative flow velocity over node
    for(int J=0; J<3; J++) {
      vq[J]=misc::dotprod(vi,Q[i])*Q[i][J];   // tangential relative flow component
      vp[J]=vi[J]-vq[J];          // transverse relative flow component
      vq_squared+=vq[J]*vq[J];
      vp_squared+=vp[J]*vp[J];
    }
    double vp_mag=sqrt(vp_squared);
    double vq_mag=sqrt(vp_squared);

    // transverse drag    
    if      (i==0)  for(int J=0; J<3; J++) { Dp[i][J]=1./2.*Env->GetRho_w()*Cdn* (F[i]*D*L[i])/2.*vp_mag*vp[J]; }
    else if(i==N)   for(int J=0; J<3; J++) { Dp[i][J]=1./2.*Env->GetRho_w()*Cdn* (F[i-1]*D*L[i-1])/2.*vp_mag*vp[J]; }
    else            for(int J=0; J<3; J++) { Dp[i][J]=1./2.*Env->GetRho_w()*Cdn* (F[i]*D*L[i]+F[i-1]*D*L[i-1])/2.*vp_mag*vp[J]; }

    // tangential drag    
    if      (i==0)  for(int J=0; J<3; J++) { Dq[i][J]=1./2.*Env->GetRho_w()*Cdt* PI*(F[i]*D*L[i])/2.*vq_mag*vq[J]; }
    else if(i==N)   for(int J=0; J<3; J++) { Dq[i][J]=1./2.*Env->GetRho_w()*Cdt* PI*(F[i-1]*D*L[i-1])/2.*vq_mag*vq[J]; }
    else            for(int J=0; J<3; J++) { Dq[i][J]=1./2.*Env->GetRho_w()*Cdt* PI*(F[i]*D*L[i]+F[i-1]*D*L[i-1])/2.*vq_mag*vq[J]; }
    // acceleration calculations          
    for(int J=0; J<3; J++) {
      aq[J]=misc::dotprod(Ud[i],Q[i])*Q[i][J]; // tangential component of fluid acceleration
      ap[J]=Ud[i][J]-aq[J];       // normal component of fluid acceleration
    }

    // transverse Froude-Krylov force
    if      (i==0)  for(int J=0; J<3; J++) { Ap[i][J]=Env->GetRho_w()*(1.+Can)*0.5*(V[i])*ap[J]; }
    else if(i==N)   for(int J=0; J<3; J++) { Ap[i][J]=Env->GetRho_w()*(1.+Can)*0.5*(V[i-1])*ap[J]; }
    else             for(int J=0; J<3; J++) { Ap[i][J]=Env->GetRho_w()*(1.+Can)*0.5*(V[i]+V[i-1])*ap[J]; }

    // tangential Froude-Krylov force          
    if      (i==0)  for(int J=0; J<3; J++) { Aq[i][J]=Env->GetRho_w()*(1.+Cat)*0.5*(V[i])*aq[J]; }
    else if(i==N)  for(int J=0; J<3; J++) { Aq[i][J]=Env->GetRho_w()*(1.+Cat)*0.5*(V[i-1])*aq[J]; }
    else             for(int J=0; J<3; J++) { Aq[i][J]=Env->GetRho_w()*(1.+Cat)*0.5*(V[i]+V[i-1])*aq[J]; }

    const double zmin=(Env->GetFreeSurface()-(*Depth));  //Gets the minimun position at z for contacts with ground
    // bottom contact (stiffness and damping, vertical-only for now)-updated for general case of potentially anchor or fairlead end in contact
    if(R[i][2]<zmin) {
      if      (i==0) {  B[i][2]=((zmin-R[i][2])*Env->GetKb()-Rd[i][2]*Env->GetCb())*0.5*(L[i]); }//B[i][2]=((zmin-R[i][2])*Env->GetKb()-Rd[i][2]*Env->GetCb())*0.5*(D*L[i-1]); //-ERROR if i==0 => L[i-1] -> segmentation fault
      else if(i==N) {   B[i][2]=((zmin-R[i][2])*Env->GetKb()-Rd[i][2]*Env->GetCb())*0.5*(D*L[i]); }
      else           {  B[i][2]=((zmin-R[i][2])*Env->GetKb()-Rd[i][2]*Env->GetCb())*0.5*(D*L[i]+D*L[i-1]); }

      // new rough-draft addition of seabed friction
      //double FrictionCoefficient=0.5;                  // just using one coefficient to start with          
      double FrictionMax=abs(B[i][2])*(Env->GetFrictionCoefficient()); // dynamic friction force saturation level based on bottom contact force
      // saturated damping approach to applying friction, for now
      double BottomVel=sqrt(Rd[i][0]*Rd[i][0]+Rd[i][1]*Rd[i][1]); // velocity of node along sea bed

      double FrictionForce=BottomVel*(Env->GetFrictionCoefficient())*(Env->GetFricDamp()); // some arbitrary damping scaling thing at end
      if(FrictionForce>Env->GetStatDynFricScale()*FrictionMax) { FrictionForce=FrictionMax; }    // saturate (quickly) to static/dynamic friction force level 
      // apply force in correct directions -- opposing direction of motion
      // could add ifs in here to handle end nodes
      if(BottomVel==double(0)) {BottomVel=DBL_MIN;}  
      if(std::isinf(BottomVel)) { 
        string tex="Error: BottomVel==+/- inf. Cannot divide by infinite. Check line Number " 
         +funmd::ToString(Number)+".";
        Run_Exceptioon(tex); 
      }
      B[i][0]=-FrictionForce*Rd[i][0]/BottomVel;
      B[i][1]=-FrictionForce*Rd[i][1]/BottomVel;
    }
    else {
      B[i][0]=0.;
      B[i][1]=0.;
      B[i][2]=0.;
    }
    // total forces
    if(i==0)       for(int J=0; J<3; J++) { Fnet[i][J]=T[i][J]+Td[i][J]+W[i][J]+(Dp[i][J]+Dq[i][J]+Ap[i][J]+Aq[i][J])+B[i][J];  }
    else if(i==N)     for(int J=0; J<3; J++) { Fnet[i][J]=-T[i-1][J]-Td[i-1][J]+W[i][J]+(Dp[i][J]+Dq[i][J]+Ap[i][J]+Aq[i][J])+B[i][J]; }
    else         for(int J=0; J<3; J++) { Fnet[i][J]=T[i][J]-T[i-1][J]+Td[i][J]-Td[i-1][J]+W[i][J]+(Dp[i][J]+Dq[i][J]+Ap[i][J]+Aq[i][J])+B[i][J]; }
    for(int J=0; J<3; J+=2) {if(Fnet[i][J]==double(0)) {  Fnet[i][J]=-DBL_MIN;} }

  }

  // loop through internal nodes and update their States
  for(unsigned i=1; i<N; i++) {
    // calculate RHS constant (premultiplying force vector by inverse of mass matrix  ... i.e. rhs=S*Forces)  
    for(int I=0; I<3; I++) {
      double RHSiI=0.0; // temporary accumulator 
      for(int J=0; J<3; J++) {
        RHSiI+=S[i][I][J]*Fnet[i][J];   //  matrix multiplication [S i]{Forces i}
        if(std::isnan(Fnet[i][J])) {
          string tex="Error: Fnet["+funmd::ToString(i)+"]["+funmd::ToString(J) 
           +"] NaN value detected in MoorDyn state at dynamic relaxation time "
           +funmd::ToString(time)+" s.";
          Run_Exceptioon(tex); 
        }
      }
      // update States
      Xd[3*N-3+3*i-3+I]=X[3*i-3+I];      // dxdt=V  (velocities)
      Xd[3*i-3+I]=RHSiI;          // dVdt=RHS*A  (accelerations)
    }
  }
  //-Checks the tension in the connections
  if(BreakTension&&!Disabled)CheckTension(); 

  return;
}
//============================================================================================
/// Checks if some Connection reached the maximun allowed value of tension
//============================================================================================
void Line::CheckTension(){
  double tenF=GetNodeTen(N);
  double tenA=GetNodeTen(0);
  if(std::abs(tenF)>=BreakTension||std::abs(tenA)>=BreakTension) BreakLine();
}
//============================================================================================
/// Disables the Line and its Connections when the maximun value of tension is reached. 
/// Throws a warning.
//============================================================================================
void Line::BreakLine(){
  FairConnect->SetDisabled(true);  FairConnect->ResetTension();
  AnchConnect->SetDisabled(true);  AnchConnect->ResetTension();
  for(unsigned i=0; i<N; i++) {Fnet[i][0]=0.0;Fnet[i][1]=0.0;Fnet[i][2]=0.0;}
  if(!Disabled)Log->PrintfWarning("The line %d reached the maximun value of tension (%fN). It will be disabled.",Number,BreakTension);
  Disabled=true;
}
//============================================================================================
/// Store in M the mass of nodes
//============================================================================================
void Line::CalculateMassMatrix() {
  // calculate mass matrix 
  for(unsigned i=0; i<=N; i++) {
    double m_i; // node mass
    double v_i; // node submerged volume 

    if(i==0) {      m_i=PI/8.*D*D*L[0]*Rho;        v_i=1./2. *F[i]*V[i];  }
    else if(i==N) {  m_i=PI/8.*D*D*L[N-2]*Rho;      v_i=1./2. *F[i-1]*V[i-1];  }
    else {        m_i=PI/8.*(D*D*Rho*(L[i]+L[i-1]));  v_i=1./2. *(F[i-1]*V[i-1]+F[i]*V[i]);  }

    // make node mass matrix
    for(int I=0; I<3; I++) {
      for(int J=0; J<3; J++) {M[i][I][J]=m_i*misc::eye(I,J)+Env->GetRho_w()*v_i *(Can*(misc::eye(I,J)-Q[i][I]*Q[i][J])+Cat*Q[i][I]*Q[i][J]);    }
    }
    misc::inverse3by3(S[i],M[i]);  // invert node mass matrix (written to S[i][:][:])  
  }
  return;
}
//============================================================================================
/// Write output file for line  (accepts time parameter since retained time value (t) 
/// will be behind by one line time step
//============================================================================================
void Line::WriteOutput(double t,double dtC) {
  // run through output flags
  // if channel is flagged for output, write to file.
  // Flags changed to just be one character (case sensitive) per output flag.  To match FASTv8 version.
  // create output file for writing output (and write channel header and units lines) if applicable
  std::string function="WriteOutput";
  bool savedata=!DtOut?true:!(t<(floor((t-dtC)/DtOut)+1.0)*DtOut);

  if(savedata && Outfile) { // check it's not null.  Null signals no individual line output files
    Outfile=new ofstream(*FileName,std::ofstream::app); //Concat the new content
    if(Outfile->is_open()) {
      // ------------- write channel names line --------------------
      unsigned nChannels=(unsigned)Channels.length();
      if((nChannels==1) && Channels[0]=='-') Outfile->close(); //If there aren't output Channels selected
      else {
        *Outfile << t << ";";   // output time
        for(unsigned i=0; i<nChannels; i++) {
          char ch=Channels[i]; //takes a character
          if(ch!='-') { WriteOutputValues(ch); } //calls to store funtion
        }
        *Outfile << "\n";
        Outfile->close();
      }
    }
    else { Run_Exceptioon("Cannot open file ouput \'"+*FileName+"\'\n"); }
  }

  return;
}

//==================================================================================
/// Write the ouptus selected in Channels tag
//==================================================================================
void Line::WriteOutputValues(const char ch) {
  bool writeUnits=(Env->GetWriteUnits()>0);
  switch (ch) {
  case 'p':// output positions
    for(unsigned i=0; i<=N; i++){  for(int J=0; J<3; J++) { *Outfile << R[i][J] << ";"; }  }
    break;
  case 'v': // output velocities
    for(unsigned i=0; i<=N; i++){  for(int J=0; J<3; J++) { *Outfile << Rd[i][J] << ";"; }}
    break;
  case 'U':// output wave velocities
    for(unsigned i=0; i<=N; i++) {  for(int J=0; J<3; J++) { *Outfile << U[i][J] << ";"; }  }
    break;
  case 'D':// output hydro force
    for(unsigned i=0; i<=N; i++) {  for(int J=0; J<3; J++) { *Outfile << Dp[i][J]+Dq[i][J]+Ap[i][J]+Aq[i][J] << ";"; }  }
    break;
  case 'c':// output internal damping force
    for(unsigned i=0; i<N; i++) {  for(int J=0; J<3; J++) { *Outfile << Td[i][J]+Td[i][J]+Td[i][J] << ";"; }  }
    break;
  case 't':  // output segment tensions
    for(unsigned i=0; i<N; i++) {
      double Tmag_squared=0.;
      for(int J=0; J<3; J++) { Tmag_squared+=T[i][J]*T[i][J]; } // doing this calculation here only,for the sake of speed
      *Outfile << sqrt(Tmag_squared) << ";";
    }
    break;
  case 's': //output segment strains
    for(unsigned i=0; i<N; i++) {  *Outfile << Lstr[i]/L[i]-1.0 << ";";  }
    break;
  case 'd': // output segment strain rates
    for(unsigned i=0; i<N; i++) {  *Outfile << Ldstr[i]/L[i] << ";";  }
    break;
  case 'b':
    for(unsigned i=0; i<=N; i++) {  for(int J=0; J<3; J++) { *Outfile << B[i][J] << ";"; }  }
    break;
  default:
    break;
  }
  return;
}

//==================================================================================
/// Initialize Node outputs file writing the header
//==================================================================================
void Line::InitializeOutputs(const char ch) {
  bool writeUnits=(Env->GetWriteUnits()>0);
  switch (ch) {
  case 'p':// output positions
    for(unsigned i=0; i<=N; i++) {  //loop through nodes
      if(!writeUnits) { *Outfile << "Node" << i << "px;Node" << i << "py;Node" << i << "pz;"; }
      else { *Outfile << "Node" << i << "px [m];Node" << i << "py [m];Node" << i << "pz [m];"; }
    }
    break;
  case 'v': // output velocities
    for(unsigned i=0; i<=N; i++) {
      if(!writeUnits) { *Outfile << "Node" << i << "vx;Node" << i << "vy;Node" << i << "vz;"; }
      else { *Outfile << "Node" << i << "vx [m/s];Node" << i << "vy [m/s];Node" << i << "vz [m/s];";; }
    }
    break;
  case 'U':// output wave velocities
    for(unsigned i=0; i<=N; i++) {
      if(!writeUnits) { *Outfile << "Node" << i << "Ux;Node" << i << "Uy;Node" << i << "Uz;"; }
      else { *Outfile << "Node" << i << "Ux [m/s];Node" << i << "Uy [m/s];Node" << i << "Uz [m/s];"; }
    }
    break;
  case 'D':// output hydro force
    for(unsigned i=0; i<=N; i++) {
      if(!writeUnits) { *Outfile << "Node" << i << "Dx;Node" << i << "Dy;Node" << i << "Dz;"; }
      else { *Outfile << "Node" << i << "Dx [N];Node" << i << "Dy [N];Node" << i << "Dz [N];"; }
    }
    break;
  case 'c':// output internal damping force
    for(unsigned i=0; i<N; i++) {
      if(!writeUnits) { *Outfile << "Seg" << i << "cx;Node" << i << "cy;Node" << i << "cz;"; }
      else { *Outfile << "Seg" << i << "cx [N];Node" << i << "cy [N];Node" << i << "cz [N];"; }
    }
    break;
  case 't':  // output segment tensions
    for(unsigned i=0; i<N; i++) {
      if(!writeUnits) { *Outfile << "Seg" << i << "Te;"; }
      else { *Outfile << "Seg" << i << "Te [N];"; }
    }
    break;
  case 's': //output segment strains
    for(unsigned i=0; i<N; i++) {
      if(!writeUnits) { *Outfile << "Seg" << i << "St;"; }
      else { *Outfile << "Seg" << i << "St [-];"; ; }
    }
    break;
  case 'd': // output segment strain rates
    for(unsigned i=0; i<N; i++) {
      if(!writeUnits) { *Outfile << "Seg" << i << "dSt;"; }
      else { *Outfile << "Seg" << i << "dSt [-/s];"; }
    }
    break;
  case 'b':
    for(unsigned i=0; i<=N; i++) {
      if(!writeUnits) { *Outfile << "Node" << i << "bx;Node" << i << "by;Node" << i << "bz;"; }
      else { *Outfile << "Node" << i << "bx [N];Node" << i << "by [N];Node" << i << "bz [N];"; }
    }
    break;
  default:
    break;
  }
  return;
}

// new function to draw instantaneous line positions in openGL context
#ifdef USEGL
void Line::drawGL(void) {
  double maxTen=0.0;
  double normTen;
  double rgb[3];
  for(unsigned i=0; i<=N; i++) {
    double newTen=GetNodeTen(i);
    if(newTen>maxTen) { maxTen=newTen; }
  }

  glColor3f(0.5,0.5,1.0);
  glBegin(GL_LINE_STRIP);
  for(unsigned i=0; i<=N; i++) {
    glVertex3d(R[i][0],R[i][1],R[i][2]);
    if(i<N) {
      normTen=GetNodeTen(i)/maxTen;
      ColorMap(normTen,rgb);
      glColor3d(rgb[0],rgb[1],rgb[2]);
    }
  }
  glEnd();
  return;
}



void Line::drawGL2(void) {
  double maxTen=0.0;
  double normTen;
  double rgb[3];
  for(unsigned i=0; i<=N; i++) {
    double newTen=GetNodeTen(i);
    if(newTen>maxTen) { maxTen=newTen; }
  }
  // line
  for(unsigned i=0; i<N; i++) {
    normTen=0.2+0.8*pow(GetNodeTen(i)/maxTen,4.0);
    ColorMap(normTen,rgb);
    glColor3d(rgb[0],rgb[1],rgb[2]);
    Cylinder(R[i][0],R[i][1],R[i][2],R[i+1][0],R[i+1][1],R[i+1][2],27,0.5);
  }
  // velocity vectors
  for(unsigned i=0; i<=N; i++) {
    glColor3d(0.0,0.2,0.8);
    double vscal=5.0;
    Arrow(R[i][0],R[i][1],R[i][2],vscal*Rd[i][0],vscal*Rd[i][1],vscal*Rd[i][2],0.1,0.7);
  }
  return;
}
#endif