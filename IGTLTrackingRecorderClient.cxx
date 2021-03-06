/*=========================================================================

  Program:   OpenIGTLink -- Example for Data Receiving Server Program
  Language:  C++

  Copyright (c) Insight Software Consortium. All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <cstdlib>
#include <cstring>

#include "igtlOSUtil.h"
#include "igtlMessageHeader.h"
#include "igtlTransformMessage.h"
#include "igtlImageMessage.h"
#include "igtlClientSocket.h"
#include "igtlStatusMessage.h"
#include "igtlPositionMessage.h"

#if OpenIGTLink_PROTOCOL_VERSION >= 2
#include "igtlPointMessage.h"
#include "igtlTrajectoryMessage.h"
#include "igtlStringMessage.h"
#include "igtlBindMessage.h"
#include "igtlCapabilityMessage.h"
#include "igtlTrackingDataMessage.h"
#endif //OpenIGTLink_PROTOCOL_VERSION >= 2


int ReceiveTransform(igtl::Socket * socket, igtl::MessageHeader * header, std::ofstream& ofs);
int ReceivePosition(igtl::Socket * socket, igtl::MessageHeader * header);
int ReceiveStatus(igtl::Socket * socket, igtl::MessageHeader * header);

#if OpenIGTLink_PROTOCOL_VERSION >= 2
int ReceivePoint(igtl::Socket * socket, igtl::MessageHeader * header, std::ofstream& ofs);
int ReceiveTrajectory(igtl::Socket * socket, igtl::MessageHeader::Pointer& header);
int ReceiveString(igtl::Socket * socket, igtl::MessageHeader * header);
int ReceiveBind(igtl::Socket * socket, igtl::MessageHeader * header);
int ReceiveCapability(igtl::Socket * socket, igtl::MessageHeader * header);
int ReceiveTrackingData(igtl::Socket * socket, igtl::MessageHeader * header, std::ofstream& ofs);
#endif //OpenIGTLink_PROTOCOL_VERSION >= 2

int main(int argc, char* argv[])
{
  //------------------------------------------------------------
  // Parse Arguments

  if (argc != 4) // check number of arguments
    {
    // If not correct, print usage
    std::cerr << "    <hostname> : IP or host name"                    << std::endl;
    std::cerr << "    <port>     : Port # (18944 in Slicer default)"   << std::endl;
    std::cerr << "    <file>     : Output file name"   << std::endl;
    exit(0);
    }

  char*  hostname = argv[1];
  int    port     = atoi(argv[2]);
  std::string fname = argv[3];

  std::ofstream ofs (fname, std::ofstream::out);

  //------------------------------------------------------------
  // Establish Connection

  igtl::ClientSocket::Pointer socket;
  socket = igtl::ClientSocket::New();
  int r = socket->ConnectToServer(hostname, port);

  int msgCounter = 0;
  
  if (r != 0)
    {
    std::cerr << "Cannot connect to the server." << std::endl;
    exit(0);
    }

  //------------------------------------------------------------
  // Create a message buffer to receive header
  igtl::MessageHeader::Pointer headerMsg;
  headerMsg = igtl::MessageHeader::New();
  
  //------------------------------------------------------------
  // Allocate a time stamp 
  // Message Time Stamp
  igtl::TimeStamp::Pointer mts;
  mts = igtl::TimeStamp::New();
  
  // Arrival Time Stamp
  igtl::TimeStamp::Pointer ats;
  ats = igtl::TimeStamp::New();
  

  while (1)
    {
    //------------------------------------------------------------
    // loop
    for (int i = 0; i < 100; i ++)
      {
      
      // Initialize receive buffer
      headerMsg->InitPack();
      
      // Receive generic header from the socket
      int r = socket->Receive(headerMsg->GetPackPointer(), headerMsg->GetPackSize());
      if (r == 0)
        {
        socket->CloseSocket();
        exit(0);
        }
      if (r != headerMsg->GetPackSize())
        {
        continue;
        }
      
      // Record time of arraival
      ats->GetTime();
      
      // Deserialize the header
      headerMsg->Unpack();

      // Message type,name, and channel
      ofs << msgCounter << "," << headerMsg->GetDeviceType() << "," << headerMsg->GetDeviceName() << ",";
      
      // Get time stamp
      igtlUint32 sec;
      igtlUint32 nanosec;

      headerMsg->GetTimeStamp(mts);
      mts->GetTimeStamp(&sec, &nanosec);

      std::cerr << "Time stamp: "
                << sec << "." << std::setw(9) << std::setfill('0') 
                << nanosec << std::endl;
      
      ofs << sec << "." << std::setw(9) << std::setfill('0') << nanosec << ",";
        
      // Arrival time stamp
      ats->GetTimeStamp(&sec, &nanosec);
      ofs << sec << "." << std::setw(9) << std::setfill('0') << nanosec << ",";
      
      // Check data type and receive data body
      if (strcmp(headerMsg->GetDeviceType(), "TRANSFORM") == 0)
        {
	ReceiveTransform(socket, headerMsg, ofs);
        }
      else if (strcmp(headerMsg->GetDeviceType(), "POSITION") == 0)
        {
          ReceivePosition(socket, headerMsg);
        }
      else if (strcmp(headerMsg->GetDeviceType(), "STATUS") == 0)
        {
        ReceiveStatus(socket, headerMsg);
        }
#if OpenIGTLink_PROTOCOL_VERSION >= 2
      else if (strcmp(headerMsg->GetDeviceType(), "POINT") == 0)
        {
        ReceivePoint(socket, headerMsg, ofs);
        }
      else if (strcmp(headerMsg->GetDeviceType(), "TRAJ") == 0)
        {
        ReceiveTrajectory(socket, headerMsg);
        }
      else if (strcmp(headerMsg->GetDeviceType(), "STRING") == 0)
        {
        ReceiveString(socket, headerMsg);
        }
      else if (strcmp(headerMsg->GetDeviceType(), "CAPABILITY") == 0)
        {
        ReceiveCapability(socket, headerMsg);;
        }
      else if (strcmp(headerMsg->GetDeviceType(), "TDATA") == 0)
        {
          ReceiveTrackingData(socket, headerMsg, ofs);
        }
#endif //OpenIGTLink_PROTOCOL_VERSION >= 2
      else
        {
        std::cerr << "Receiving : " << headerMsg->GetDeviceType() << std::endl;
        socket->Skip(headerMsg->GetBodySizeToRead(), 0);
        }
      ofs << std::endl;
      msgCounter ++;
      }
    }

  //------------------------------------------------------------
  // Close connection (The example code never reaches this section ...)
  
  socket->CloseSocket();

}

int ReceiveTransform(igtl::Socket * socket, igtl::MessageHeader * header, std::ofstream& ofs)
{
  std::cerr << "Receiving TRANSFORM data type." << std::endl;
  
  // Create a message buffer to receive transform data
  igtl::TransformMessage::Pointer transMsg;
  transMsg = igtl::TransformMessage::New();
  transMsg->SetMessageHeader(header);
  transMsg->AllocatePack();
  
  // Receive transform data from the socket
  socket->Receive(transMsg->GetPackBodyPointer(), transMsg->GetPackBodySize());
  
  // Deserialize the transform data
  // If you want to skip CRC check, call Unpack() without argument.
  int c = transMsg->Unpack(1);
  
  if (c & igtl::MessageHeader::UNPACK_BODY) // if CRC check is OK
    {
    // Retrive the transform data
    igtl::Matrix4x4 matrix;
    transMsg->GetMatrix(matrix);
    //igtl::PrintMatrix(matrix);
    ofs << matrix[0][0] << "," << matrix[0][1] << "," << matrix[0][2] << "," << matrix[0][3] << ",";
    ofs << matrix[1][0] << "," << matrix[1][1] << "," << matrix[1][2] << "," << matrix[1][3] << ",";
    ofs << matrix[2][0] << "," << matrix[2][1] << "," << matrix[2][2] << "," << matrix[2][3] << ",";
    ofs << matrix[3][0] << "," << matrix[3][1] << "," << matrix[3][2] << "," << matrix[3][3] << ",";
    return 1;
    }

  return 0;

}


int ReceivePosition(igtl::Socket * socket, igtl::MessageHeader * header)
{
  std::cerr << "Receiving POSITION data type." << std::endl;
  
  // Create a message buffer to receive transform data
  igtl::PositionMessage::Pointer positionMsg;
  positionMsg = igtl::PositionMessage::New();
  positionMsg->SetMessageHeader(header);
  positionMsg->AllocatePack();
  
  // Receive position position data from the socket
  socket->Receive(positionMsg->GetPackBodyPointer(), positionMsg->GetPackBodySize());
  
  // Deserialize the transform data
  // If you want to skip CRC check, call Unpack() without argument.
  int c = positionMsg->Unpack(1);
  
  if (c & igtl::MessageHeader::UNPACK_BODY) // if CRC check is OK
    {
    // Retrive the transform data
    float position[3];
    float quaternion[4];

    positionMsg->GetPosition(position);
    positionMsg->GetQuaternion(quaternion);

    std::cerr << "position   = (" << position[0] << ", " << position[1] << ", " << position[2] << ")" << std::endl;
    std::cerr << "quaternion = (" << quaternion[0] << ", " << quaternion[1] << ", "
              << quaternion[2] << ", " << quaternion[3] << ")" << std::endl << std::endl;

    return 1;
    }

  return 0;

}

int ReceiveStatus(igtl::Socket * socket, igtl::MessageHeader * header)
{

  std::cerr << "Receiving STATUS data type." << std::endl;

  // Create a message buffer to receive transform data
  igtl::StatusMessage::Pointer statusMsg;
  statusMsg = igtl::StatusMessage::New();
  statusMsg->SetMessageHeader(header);
  statusMsg->AllocatePack();
  
  // Receive transform data from the socket
  socket->Receive(statusMsg->GetPackBodyPointer(), statusMsg->GetPackBodySize());
  
  // Deserialize the transform data
  // If you want to skip CRC check, call Unpack() without argument.
  int c = statusMsg->Unpack(1);
  
  if (c & igtl::MessageHeader::UNPACK_BODY) // if CRC check is OK
    {
    std::cerr << "========== STATUS ==========" << std::endl;
    std::cerr << " Code      : " << statusMsg->GetCode() << std::endl;
    std::cerr << " SubCode   : " << statusMsg->GetSubCode() << std::endl;
    std::cerr << " Error Name: " << statusMsg->GetErrorName() << std::endl;
    std::cerr << " Status    : " << statusMsg->GetStatusString() << std::endl;
    std::cerr << "============================" << std::endl;
    }

  return 0;

}


#if OpenIGTLink_PROTOCOL_VERSION >= 2
int ReceivePoint(igtl::Socket * socket, igtl::MessageHeader * header, std::ofstream& ofs)
{

  std::cerr << "Receiving POINT data type." << std::endl;

  // Create a message buffer to receive transform data
  igtl::PointMessage::Pointer pointMsg;
  pointMsg = igtl::PointMessage::New();
  pointMsg->SetMessageHeader(header);
  pointMsg->AllocatePack();

  // Receive transform data from the socket
  socket->Receive(pointMsg->GetPackBodyPointer(), pointMsg->GetPackBodySize());

  // Deserialize the transform data
  // If you want to skip CRC check, call Unpack() without argument.
  int c = pointMsg->Unpack(1);

  if (c & igtl::MessageHeader::UNPACK_BODY) // if CRC check is OK
    {
    int nElements = pointMsg->GetNumberOfPointElement();
    for (int i = 0; i < nElements; i ++)
      {
      igtl::PointElement::Pointer pointElement;
      pointMsg->GetPointElement(i, pointElement);

      igtlUint8 rgba[4];
      pointElement->GetRGBA(rgba);

      igtlFloat32 pos[3];
      pointElement->GetPosition(pos);

      // Format: element#, group name, owner, R, G, B, A, x, y, z, rad, element#, group name, owner, R, G, B, A, x, y, z, rad, ...
      ofs << pointElement->GetName() << ","
          << pointElement->GetGroupName() << ","
          << pointElement->GetOwner() << ","
          << (int)rgba[0] << "," << (int)rgba[1] << "," << (int)rgba[2] << "," << (int)rgba[3] << ","
          << std::fixed << pos[0] << "," << pos[1] << "," << pos[2] << ","
          << std::fixed << pointElement->GetRadius() << ",";
      
      //std::cerr << "========== Element #" << i << " ==========" << std::endl;
      //std::cerr << " Name      : " << pointElement->GetName() << std::endl;
      //std::cerr << " GroupName : " << pointElement->GetGroupName() << std::endl;
      //std::cerr << " RGBA      : ( " << (int)rgba[0] << ", " << (int)rgba[1] << ", " << (int)rgba[2] << ", " << (int)rgba[3] << " )" << std::endl;
      //std::cerr << " Position  : ( " << std::fixed << pos[0] << ", " << pos[1] << ", " << pos[2] << " )" << std::endl;
      //std::cerr << " Radius    : " << std::fixed << pointElement->GetRadius() << std::endl;
      //std::cerr << " Owner     : " << pointElement->GetOwner() << std::endl;
      //std::cerr << "================================" << std::endl;
      }
    }

  return 1;
}

int ReceiveTrajectory(igtl::Socket * socket, igtl::MessageHeader::Pointer& header)
{

  std::cerr << "Receiving TRAJECTORY data type." << std::endl;

  // Create a message buffer to receive transform data
  igtl::TrajectoryMessage::Pointer trajectoryMsg;
  trajectoryMsg = igtl::TrajectoryMessage::New();
  trajectoryMsg->SetMessageHeader(header);
  trajectoryMsg->AllocatePack();

  // Receive transform data from the socket
  socket->Receive(trajectoryMsg->GetPackBodyPointer(), trajectoryMsg->GetPackBodySize());

  // Deserialize the transform data
  // If you want to skip CRC check, call Unpack() without argument.
  int c = trajectoryMsg->Unpack(1);

  if (c & igtl::MessageHeader::UNPACK_BODY) // if CRC check is OK
    {
    int nElements = trajectoryMsg->GetNumberOfTrajectoryElement();
    for (int i = 0; i < nElements; i ++)
      {
      igtl::TrajectoryElement::Pointer trajectoryElement;
      trajectoryMsg->GetTrajectoryElement(i, trajectoryElement);

      igtlUint8 rgba[4];
      trajectoryElement->GetRGBA(rgba);

      igtlFloat32 entry[3];
      igtlFloat32 target[3];
      trajectoryElement->GetEntryPosition(entry);
      trajectoryElement->GetTargetPosition(target);

      std::cerr << "========== Element #" << i << " ==========" << std::endl;
      std::cerr << " Name      : " << trajectoryElement->GetName() << std::endl;
      std::cerr << " GroupName : " << trajectoryElement->GetGroupName() << std::endl;
      std::cerr << " RGBA      : ( " << (int)rgba[0] << ", " << (int)rgba[1] << ", " << (int)rgba[2] << ", " << (int)rgba[3] << " )" << std::endl;
      std::cerr << " Entry Pt  : ( " << std::fixed << entry[0] << ", " << entry[1] << ", " << entry[2] << " )" << std::endl;
      std::cerr << " Target Pt : ( " << std::fixed << target[0] << ", " << target[1] << ", " << target[2] << " )" << std::endl;
      std::cerr << " Radius    : " << std::fixed << trajectoryElement->GetRadius() << std::endl;
      std::cerr << " Owner     : " << trajectoryElement->GetOwner() << std::endl;
      std::cerr << "================================" << std::endl << std::endl;
      }
    }

  return 1;
}


int ReceiveString(igtl::Socket * socket, igtl::MessageHeader * header)
{

  std::cerr << "Receiving STRING data type." << std::endl;

  // Create a message buffer to receive transform data
  igtl::StringMessage::Pointer stringMsg;
  stringMsg = igtl::StringMessage::New();
  stringMsg->SetMessageHeader(header);
  stringMsg->AllocatePack();

  // Receive transform data from the socket
  socket->Receive(stringMsg->GetPackBodyPointer(), stringMsg->GetPackBodySize());

  // Deserialize the transform data
  // If you want to skip CRC check, call Unpack() without argument.
  int c = stringMsg->Unpack(1);

  if (c & igtl::MessageHeader::UNPACK_BODY) // if CRC check is OK
    {
    std::cerr << "Encoding: " << stringMsg->GetEncoding() << "; "
              << "String: " << stringMsg->GetString() << std::endl;
    }

  return 1;
}


int ReceiveBind(igtl::Socket * socket, igtl::MessageHeader * header)
{

  std::cerr << "Receiving BIND data type." << std::endl;

  // Create a message buffer to receive transform data
  igtl::BindMessage::Pointer bindMsg;
  bindMsg = igtl::BindMessage::New();
  bindMsg->SetMessageHeader(header);
  bindMsg->AllocatePack();

  // Receive transform data from the socket
  socket->Receive(bindMsg->GetPackBodyPointer(), bindMsg->GetPackBodySize());

  // Deserialize the transform data
  // If you want to skip CRC check, call Unpack() without argument.
  int c = bindMsg->Unpack(1);

  if (c & igtl::MessageHeader::UNPACK_BODY) // if CRC check is OK
    {
    int n = bindMsg->GetNumberOfChildMessages();

    for (int i = 0; i < n; i ++)
      {
      if (strcmp(bindMsg->GetChildMessageType(i), "STRING") == 0)
        {
        igtl::StringMessage::Pointer stringMsg;
        stringMsg = igtl::StringMessage::New();
        bindMsg->GetChildMessage(i, stringMsg);
        stringMsg->Unpack(0);
        std::cerr << "Message type: STRING" << std::endl;
        std::cerr << "Message name: " << stringMsg->GetDeviceName() << std::endl;
        std::cerr << "Encoding: " << stringMsg->GetEncoding() << "; "
                  << "String: " << stringMsg->GetString() << std::endl;
        }
      else if (strcmp(bindMsg->GetChildMessageType(i), "TRANSFORM") == 0)
        {
        igtl::TransformMessage::Pointer transMsg;
        transMsg = igtl::TransformMessage::New();
        bindMsg->GetChildMessage(i, transMsg);
        transMsg->Unpack(0);
        std::cerr << "Message type: TRANSFORM" << std::endl;
        std::cerr << "Message name: " << transMsg->GetDeviceName() << std::endl;
        igtl::Matrix4x4 matrix;
        transMsg->GetMatrix(matrix);
        igtl::PrintMatrix(matrix);
        }
      }
    }

  return 1;
}


int ReceiveCapability(igtl::Socket * socket, igtl::MessageHeader * header)
{
  
  std::cerr << "Receiving CAPABILITY data type." << std::endl;

  // Create a message buffer to receive transform data
  igtl::CapabilityMessage::Pointer capabilMsg;
  capabilMsg = igtl::CapabilityMessage::New();
  capabilMsg->SetMessageHeader(header);
  capabilMsg->AllocatePack();

  // Receive transform data from the socket
  socket->Receive(capabilMsg->GetPackBodyPointer(), capabilMsg->GetPackBodySize());

  // Deserialize the transform data
  // If you want to skip CRC check, call Unpack() without argument.
  int c = capabilMsg->Unpack(1);
  
  if (c & igtl::MessageHeader::UNPACK_BODY) // if CRC check is OK
    {
    int nTypes = capabilMsg->GetNumberOfTypes();
    for (int i = 0; i < nTypes; i ++)
      {
      std::cerr << "Typename #" << i << ": " << capabilMsg->GetType(i) << std::endl;
      }
    }

  return 1;
  
}


int ReceiveTrackingData(igtl::Socket * socket, igtl::MessageHeader * header, std::ofstream& ofs)
{
  std::cerr << "Receiving TDATA data type." << std::endl;
  
  //------------------------------------------------------------
  // Allocate TrackingData Message Class

  igtl::TrackingDataMessage::Pointer trackingData;
  trackingData = igtl::TrackingDataMessage::New();
  trackingData->SetMessageHeader(header);
  trackingData->AllocatePack();

  // Receive body from the socket
  socket->Receive(trackingData->GetPackBodyPointer(), trackingData->GetPackBodySize());

  // Deserialize the transform data
  // If you want to skip CRC check, call Unpack() without argument.
  int c = trackingData->Unpack(1);

  if (c & igtl::MessageHeader::UNPACK_BODY) // if CRC check is OK
    {
    int nElements = trackingData->GetNumberOfTrackingDataElements();
    for (int i = 0; i < nElements; i ++)
      {
      igtl::TrackingDataElement::Pointer trackingElement;
      trackingData->GetTrackingDataElement(i, trackingElement);

      igtl::Matrix4x4 matrix;
      trackingElement->GetMatrix(matrix);

      // Format: element#, group name, owner, R, G, B, A, x, y, z, rad, element#, group name, owner, R, G, B, A, x, y, z, rad, ...
      ofs << trackingElement->GetName() << ","
          << "GROUP,"
          << "OWNER,"
          << "0,0,0,0,"
          << std::fixed << matrix[0][3] << "," << matrix[1][3] << "," << matrix[2][3] << ","
          << "1.0,";
      }
    return 1;
    }
  return 0;
}

#endif //OpenIGTLink_PROTOCOL_VERSION >= 2
