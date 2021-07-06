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

/// \file Line.h \brief Defines the class \ref Line.

#ifndef _Line_
#define _Line_

#include "Misc.h"
#include "JXml.h"
#include "JObject.h"
#include "JLog2.h"
#include "TypesMoorDyn.h"

#include <limits.h>
#include <float.h>

class Body;
class Connection;
class Environment;
class QSlines;


//********************************************************************************
/// Class LineDefault
//********************************************************************************
class LineDefault : protected JObject
{

private:
  void Reset(); /// Restores the attributes of the class
  void ReadXml(JXml *sxml, TiXmlElement* lis, const std::string &place); ///Reads the input file
public:
  double Ea; ///< Stiffness of Line
  double Diameter; ///< Diameter of Line
  double Can; ///< Default=1.0
  double Cat; ///< Default=0.0
  double Cdn; ///< Default=1.6
  double Cdt; ///< Default=0.05
  double W; ///< Linear weight in air
  double C; ///< The line internal damping (Ns)
  double BreakTension; ///< Maximun value of tension allowed
  std::string Channels; ///< Node output properties
  LineDefault(); /// Constructor
  ~LineDefault(); /// Destructor
  void LoadXml(JXml *sxml, const std::string &place); /// Loads Xml file
  void CheckChannels(JXml *sxml, TiXmlElement * elel); /// Checks if the output flag is in the list of posibilities
};

//********************************************************************************
/// Class Line
//********************************************************************************
class Line : protected JObject
{
  // here is the numbering scheme (N segments per line):
  // [connect (node 0)]  --- segment 0 --- [ node 1 ] --- seg 1 --- [node2] --- ... --- seg n-2 --- [node n-1] --- seg n-1 ---  [connect (node N)]
private:
  double DtOut; ///< Time step for write
  double BreakTension; ///< Maximun value of tension allowed
  bool Disabled; ///< Indicates if tje line was broken
  bool ConnectFlag; ///< Indicates if the line has a connect connection
  bool FairFlag; ///< Indicates if the line has a fair connection
  bool OutputLine; ///< Indticates if this line will be stored in the csv file. default=true
  bool OutputNode; ///< Indticates if the nodes of this line will be stored in the csv files. default=true
  JLog2 *Log;
  QSlines *Qslines;
  double *Depth; ///< Pointer to Body Depth
  std::string *FileName; ///< Pointer to file output name
  Environment *Env; ///< Pointer to environmental settings
  unsigned N; ///< Number of line segments 
  double UnstrLen; ///< Length of line
  unsigned *BodyId; ///< Id of body and Floating object (used when there is a coupling with external software)
  std::vector<std::vector<double>> R; ///< Node positions [i][x/y/z]
  std::vector<std::vector<double>> Rd; ///< Node velocities [i][x/y/z]
  std::vector<std::vector<double>> Q; ///< Unit tangent vectors for each segment

  std::vector<Connection*> Connections; ///< Array of Connections
  Connection *AnchConnect; ///< Pointer to anchor connection
  Connection *FairConnect; ///< pointer to fairlead connection

  unsigned Number; ///< Line "number" id
  unsigned NConnections; ///< Total number of Connection objects
  unsigned NFairs; ///< Number of fairlead connections
  unsigned NAnchs; ///< Number of anchor connections
  unsigned NConns; ///< Number of "connect" connections
  unsigned StartIndex; ///< Start index for each node for States array 
  int WaveKin; ///< Flag indicating whether wave kinematics will be considered for this line

  double t; ///< Simulation time
  double t0; ///< Simulation time current integration was started at (used for BC function)
  double tlast;

  // declaring these here as they caused problems if declared in the loop
  double vi[3]; ///< relative velocity
  double vp[3]; ///< transverse component of relative velocity
  double vq[3]; ///< axial component of relative velocity
  double ap[3]; ///< transverse component of absolute acceleration-HAD TO MOVE THIS UP FROM LOWER DOWN (USED TO CRASH AFTER I=3)
  double aq[3]; ///< axial component of absolute acceleration

  // forces 
  std::vector<std::vector<double>> T; ///<
  std::vector<std::vector<double>> Td; ///<
  std::vector<double> Tmag; ///< Segment tension magnitude 
  std::vector<std::vector<double>> W; ///< Node weight 

  std::vector<std::vector<double>> Dp; ///< Node drag (transverse)
  std::vector<std::vector<double>> Dq; ///< Node drag (axial)
  std::vector<std::vector<double>> Ap; ///< Node added mass forcing (transverse)
  std::vector<std::vector<double>> Aq; ///< Node added mass forcing (axial)
  std::vector<std::vector<double>> B; ///< Node bottom contact force

  std::vector<std::vector<double>> Fnet; ///< total force on node

  std::vector<std::vector<std::vector<double>> > S; ///< Inverse mass matrices (3x3) for each node
  std::vector<std::vector<std::vector<double>> > M; ///< Node mass+added mass matrix

  std::vector<double> F; ///< VOF scalar for each segment (1=fully submerged, 0=out of water)
  std::vector<double> L; ///< Line unstretched segment lengths
  std::vector<double> Lstr; ///< Stretched lengths (m)
  std::vector<double> Ldstr; ///< Rate of stretch (m/s)

  double D; ///< line diameter
  double Rho; ///< line density
  double E; ///< line elasticity modulus [N]
  double C; ///< line axial internal damping coefficient [Ns]
  double Can; ///<
  double Cat; ///<
  double Cdn; ///<
  double Cdt; ///<
  double ReFac; ///< 
  double MassDenInAir; ///< Linear weight in air [ w in LineDefault ]

  std::vector<double> V; ///< Line segment volume

  //-set up output arrays, at each node i:
  std::vector<std::vector<double>>  U; ///< Wave velocities  
  std::vector<std::vector<double>>  Ud; ///< Wave accelerations
  std::vector<double> Zeta; ///< Free surface elevation

  //-file stuff
  std::ofstream *Outfile; ///< Pointer to a output file
  std::string Channels; ///< Channels selected by user
  //-new additions for handling waves in-object and precalculating them  (not necessarily used right now)
  int WaveMod; ///<   
  int WaveStMod; ///< 
  double Hs; ///<   
  double Tp; ///< 
  double Gamma; ///< 
  float Beta; ///< Wave heading

  int Nw; ///< Number of wave frequency components considered    //OK AS INT???
  std::vector<float> w; ///< 
  std::vector<float> k; ///< 
  float Dw; ///< FAST's dw (really small typically)

  std::vector<misc::floatC> ZetaC0; ///< Fourier transform of wave elevation at origin
  std::vector<misc::floatC> ZetaC; ///< Fourier transform of wave elevation
  std::vector<std::vector<misc::floatC>> UC; ///< Fourier transform of wave velocities
  std::vector<std::vector<misc::floatC>> UdC; ///< Fourier transform of wave accelerations

  std::vector<misc::doubleC> WGNC; ///< 
  std::vector<double> WaveS2Sdd; ///< 
  misc::doubleC WGNC_Fact; ///< This factor is needed by the discrete time inverse Fourier transform to ensure that the time series WGN process has unit variance // sqrt( PI/(dw*WaveDT) );   
  double S2Sd_Fact; ///< This factor is also needed by the discrete time inverse Fourier transform // 1.0/WaveDT;
  double * Ucurrent; ///< Constant uniform current to add (three components)
  // new additions for precalculating wave quantities
  std::vector<std::vector<double>> zetaTS; ///< Time series of wave elevations above each node
  std::vector<std::vector<double>> FTS; ///< 
  std::vector<std::vector< std::vector<double>> > UTS; ///< 
  std::vector<std::vector< std::vector<double>> > UdTS; ///< 
  int Nt; ///< Number of wave time steps
  double WaveDT; ///< Wave time step size (s)
  std::vector<double> tTS; ///< Time step vector
  int Ts0; ///< Time step index used for interpolating wave kinematics time series data (put here so it's persistent)

  void Reset(); /// Restore the variables
  void ReadXml(JXml *sxml, const std::string &place,TiXmlElement *eleL,const unsigned numLine,std::vector<Connection*> connects,LineDefault *lineDef=NULL); ///Reads the input file
  void AllocateMemory(); /// Reserves memory for the arrays
  void AllocateMemoryVectors(); /// Reserves memory for the vectors
  void InitializeOutputs(const char ch); /// Writes the header for output files
  void WriteOutputValues(const char ch); /// Writes the output values
  void CalculateMassMatrix(); /// Stores in M the mass of nodes
  void CheckInputParameters(); /// Checks if the output flag is in the array of posibilities (outputChannels)
  void CheckDtMRelationship(); /// Check if the relationship between dtM, N and EA is rigth
  void CheckChannels(JXml *sxml, TiXmlElement * elel); /// Checks if the output flag is in the list of posibilities
  void ConfigLine(std::vector<Connection*> connections); ///Configures the start and end Connection
  Connection * SearchConnect(const unsigned ref, std::vector<Connection*> connects); /// Checks if the connect Connection selected in Lines was configurated in linedefault


public:
  Line(JLog2 * log); /// Constructor
  ~Line(); /// Destructor
  void LoadXml(JXml *sxml,const std::string &place,TiXmlElement *eleL,const unsigned num_tag, std::vector<Connection*> connects, LineDefault *lineDef=NULL); /// Loads Xml file
  void Setup(std::string &dir,Environment *env_in,Body *body=NULL); /// Stores initial properties and creates output directory
  void SetupWaves(double *Ucurrent_in,float dt_in); /// Initializes  wave parameters
  void Initialize(double *X); /// Gets ICs for line using quasi-static approach
  double GetRho() { return Rho; }; /// Returns the density
  unsigned GetNConnections() { return NConnections; }; /// Returns the number of connections
  unsigned GetNFairs() { return NFairs; }; /// Returns the number of Fairleads
  unsigned GetNAnchs() { return NAnchs; }; /// Returns the number of Anchors
  unsigned GetNConns() { return NConns; }; /// Returns the number of Connects
  unsigned GetNConns(const unsigned ref); /// Returns the number of Connects by reference
  unsigned GetNumber() { return Number; }; /// Returns the number of Line
  unsigned GetN() { return N; }; /// Returns N (number of segments)  
  unsigned GetIndex() { return StartIndex; }; /// Returns the Start Index
  double GetLength() { return UnstrLen; }; /// Returns the length
  bool GetDisabled() {return Disabled; }; /// Returns the break tension value
  double GetBreakTension() {return BreakTension; }; /// Returns true if the line is broken, else returns false
  double GetM(const unsigned seg) { return M[seg][0][0]; }; /// Returns the Mass of one segment
  double GetV(const unsigned seg) { return V[seg]; }; /// Returns the Volumn of one segment
  void SetIndex(const unsigned index) { StartIndex=index; }; /// Stores the new Index
  void SetDisabled(const bool b) { Disabled=b; }; /// Stores the new value for Disabled variable
  void SetDtOut(const double d) { DtOut=d; }; /// Stores the new value for time step to write node properties
  bool HasConnect(){return ConnectFlag;}; /// Returns true if this line has a connect connection
  bool HasFair(){return FairFlag;}; /// Returns true if this line has a fair connection
  bool HasOutput(){return OutputLine;}; /// Returns true if will store data of this line
  double GetNodeTen(int i); /// Returns the tension of a node
  unsigned GetNodePos(unsigned i, double pos[3]); /// Returns the position of node number i
  double GetTensionOutput(const unsigned n); /// Returns the tension to write in output file
  tdouble3 GetPositionOutput(const unsigned n); /// Returns the force to write it in output file
  tdouble3 GetVelocityOutput(const unsigned n); /// Returns the force for write it in output file
  tdouble3 GetForceOutput(const unsigned n); /// Returns the force for write it in output file
  Connection * GetAnchConnect(); /// Returns the "Anchor" connections
  Connection * GetFairConnect(); /// Returns the "Fairlead" connections
  std::vector<Connection*> GetConnections(){return Connections;}; /// Returns the array of Connections
  std::vector<Connection*> GetConnections(const unsigned type); /// Returns the array of Connections by type
  std::vector<Connection*> GetFairs(const unsigned ref); /// Returns a vector of fairs ref
  std::vector<Connection*> GetConnConnects(); /// Returns the "Connects" connections
  std::vector<Connection*> GetConnConnects(const unsigned ref); /// Returns the "Connects" by reference
  void scaleDrag(double scaler); /// Function for boosting drag coefficients during IC generation  
  void setTime(double time); /// Function to reset time after IC generation
  void DoRHS(const double* X, double* Xd, const double time, const double dt); /// This is the big function that updates the States
  void WriteOutput(double t, double dtC); /// Write output file for line
  void GetFASTtens(float* FairHTen, float* FairVTen, float* AnchHTen, float* AnchVTen); /// FASTv7 style line tension outputs
  void GetAnchStuff(std::vector<double> &Fnet_out, std::vector<std::vector<double>> &M_out); /// Writes Forces and Mass into the vectors, which are pased how parameters. For Anchors
  void GetFairStuff(std::vector<double> &Fnet_out, std::vector<std::vector<double>> &M_out); /// Writes Forces and Mass into the vectors, which are pased how parameters. For Fairleads
  void BreakLine();
  void CheckTension();
#ifdef USEGL  
  void drawGL(void);
  void drawGL2(void);
#endif
};

#endif


