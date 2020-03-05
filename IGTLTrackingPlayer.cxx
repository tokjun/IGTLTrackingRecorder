/*=========================================================================

  Program:   
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

=========================================================================*/

#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include <cstring>

#include "igtlOSUtil.h"
#include "igtlMath.h"
#include "igtlTransformMessage.h"
#include "igtlServerSocket.h"

typedef struct {
  int index;
  std::string name;
  double ts;
  std::string pname0;
  double x0;
  double y0;
  double z0;
  std::string pname1;
  double x1;
  double y1;
  double z1;
  std::string pname2;
  double x2;
  double y2;
  double z2;
  std::string pname3;
  double x3;
  double y3;
  double z3;
} TrackingData;

typedef struct {
  double ts;
  igtl::Matrix4x4 matrix;
} MatrixData;

typedef std::vector<TrackingData> TrackingDataList;
typedef std::vector<MatrixData> MatrixList;


void ReadFile(std::string filename, TrackingDataList& coordinates);
void ConvertTrackingData(TrackingDataList& coordinates, MatrixList& matrices);
void GetRandomTestMatrix(igtl::Matrix4x4& matrix);

int main(int argc, char* argv[])
{
  //------------------------------------------------------------
  // Parse Arguments

  if (argc != 4) // check number of arguments
    {
    // If not correct, print usage
    std::cerr << "Usage: " << argv[0] << " <port> <fps> <file>"    << std::endl;
    std::cerr << "    <port>     : Port # (18944 in Slicer default)"   << std::endl;
    std::cerr << "    <fps>      : Frequency (fps) to send coordinate" << std::endl;
    std::cerr << "    <file>     : Tracking file" << std::endl;
    exit(0);
    }

  int    port     = atoi(argv[1]);
  double fps      = atof(argv[2]);
  std::string filename = argv[3];
  
  int    interval = (int) (1000.0 / fps);

  TrackingDataList coordinates;
  MatrixList matrices;
  ReadFile(filename, coordinates);
  ConvertTrackingData(coordinates, matrices);
  
  igtl::TransformMessage::Pointer transMsg;
  transMsg = igtl::TransformMessage::New();
  transMsg->SetDeviceName("Tracker");

  igtl::ServerSocket::Pointer serverSocket;
  serverSocket = igtl::ServerSocket::New();
  int r = serverSocket->CreateServer(port);

  if (r < 0)
    {
    std::cerr << "Cannot create a server socket." << std::endl;
    exit(0);
    }

  igtl::Socket::Pointer socket;
  
  while (1)
    {
    //------------------------------------------------------------
    // Waiting for Connection
    socket = serverSocket->WaitForConnection(1000);
    
    if (socket.IsNotNull()) // if client connected
      {
      //------------------------------------------------------------
      // loop
      MatrixList::iterator iter;
      for (iter = matrices.begin(); iter != matrices.end(); iter ++)
        {
        transMsg->InitPack();
        transMsg->SetMatrix(iter->matrix);
        transMsg->Pack();
        socket->Send(transMsg->GetPackPointer(), transMsg->GetPackSize());
        igtl::Sleep(interval); // wait
        }
      }
    }
    
  //------------------------------------------------------------
  // Close connection (The example code never reachs to this section ...)
  
  socket->CloseSocket();

}


void ReadFile(std::string filename, TrackingDataList& coordinates)
{
  std::ifstream ifs;

  ifs.open(filename);
  if (!ifs.is_open())
    {
    std::cerr << "Error opening file: " << filename << std::endl;
    exit(0);
    }
  
  while (!ifs.eof())
    {
    double ots; // time stamp at origin
    std::string type;
    std::string group0, group1, group2, group3;
    std::string owner0, owner1, owner2, owner3;
    int r0, g0, b0, a0, r1, g1, b1, a1, r2, g2, b2, a2, r3, g3, b3, a3;
    double rad0, rad1, rad2, rad3;
    TrackingData pt;
    //if (ifs >> pt.ts >> pt.x >> pt.y >> pt.z >> pt.vx >> pt.vy >> pt.vz)
    if (ifs >> pt.index >> type >> pt.name >> ots >> pt.ts
        >> pt.pname0 >> group0 >> owner0 >> r0 >> g0 >> b0 >> a0 >> pt.x0 >> pt.y0 >> pt.z0 >> rad0
        >> pt.pname1 >> group1 >> owner1 >> r1 >> g1 >> b1 >> a1 >> pt.x1 >> pt.y1 >> pt.z1 >> rad1
        >> pt.pname2 >> group2 >> owner2 >> r2 >> g2 >> b2 >> a2 >> pt.x2 >> pt.y2 >> pt.z2 >> rad2
        >> pt.pname3 >> group3 >> owner3 >> r3 >> g3 >> b3 >> a3 >> pt.x3 >> pt.y3 >> pt.z3 >> rad3)
      {
      // TODO: Need LPS to RAS conversion?
      coordinates.push_back(pt);
      }
    else
      {
      std::cerr << "Read: " << coordinates.size() << " coordinates from file." << std::endl;
      }
    }

  ifs.close();
  
}


void ConvertTrackingData(TrackingDataList& coordinates, MatrixList& matrices)
{
  
  TrackingDataList::iterator iter;

  for (iter = coordinates.begin(); iter != coordinates.end(); iter ++)
    {
    MatrixData matrixData;
    igtl::Matrix4x4& matrix = matrixData.matrix;
    float t[3];
    float s[3];
    float n[3];
    float nlen;

    matrixData.ts = iter->ts;
    
    n[0] = iter->vx;
    n[1] = iter->vy;
    n[2] = iter->vz;

    nlen = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    if (nlen > 0)
      {
      n[0] /= nlen;
      n[1] /= nlen;
      n[2] /= nlen;
      }
    else
      {
      n[0] = 0.0;
      n[1] = 0.0;
      n[2] = 1.0;
      }
    
    // Check if <n> is not parallel to <s>=(0.0, 1.0, 0.0)
    if (n[1] < 1.0)
      {
      s[0] = 0.0;
      s[1] = 1.0;
      s[2] = 0.0;
      igtl::Cross(t, s, n);
      igtl::Cross(s, n, t);
      }
    else
      {
      t[0] = 1.0;
      t[1] = 0.0;
      t[2] = 0.0;
      igtl::Cross(s, n, t);
      igtl::Cross(t, s, n);
      }
    
    matrix[0][0] = t[0];
    matrix[1][0] = t[1];
    matrix[2][0] = t[2];
    matrix[0][1] = s[0];
    matrix[1][1] = s[1];
    matrix[2][1] = s[2];
    matrix[0][2] = n[0];
    matrix[1][2] = n[1];
    matrix[2][2] = n[2];
    matrix[0][3] = iter->x;
    matrix[1][3] = iter->y;
    matrix[2][3] = iter->z;

    igtl::PrintMatrix(matrix);
    matrices.push_back(matrixData);
    }
}






