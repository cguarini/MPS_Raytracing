//This file contains the code that the master process will execute.

#include <iostream>
#include <mpi.h>
#include <math.h>
#include "RayTrace.h"
#include "slave.h"

void slaveMain(ConfigData* data)
{
    //Depending on the partitioning scheme, different things will happen.
    //You should have a different function for each of the required 
    //schemes that returns some values that you need to handle.
    switch (data->partitioningMode)
    {
        case PART_MODE_NONE:
            //The slave will do nothing since this means sequential operation.
            break;
        case PART_MODE_STATIC_STRIPS_HORIZONTAL:
            renderStaticStripsHorizontal(data);
            break;
        case PART_MODE_STATIC_STRIPS_VERTICAL:
            //The slave will do nothing since this means sequential operation.
            break;
        case PART_MODE_STATIC_BLOCKS:
            slaveStaticBlocks(data);
            break;
        case PART_MODE_STATIC_CYCLES_HORIZONTAL:
            slaveStaticCyclesHorizontal(data);
            break;
        case PART_MODE_STATIC_CYCLES_VERTICAL:
            //The slave will do nothing since this means sequential operation.
            break;
        case PART_MODE_DYNAMIC:
            //The slave will do nothing since this means sequential operation.
            slaveDynamic(data);
            break;
        default:
            std::cout << "This mode (" << data->partitioningMode;
            std::cout << ") is not currently implemented." << std::endl;
            break;
    }
}

void renderStaticStripsHorizontal(ConfigData * data) {
    float startTime, stopTime;
    startTime = MPI_Wtime();
    ///////////////////////////////////////////////////////////////////
    //Index Calculation
    ///////////////////////////////////////////////////////////////////

    //Calculate initial start row for this task
    int startIndex = (data->height / data->mpi_procs) * data->mpi_rank;
    //Distribute leftover strips, since it may not be easily divisible by processors
    int leftoverStrips = data->height % data->mpi_procs;
    if(data->mpi_rank < leftoverStrips){
        startIndex += data->mpi_rank % leftoverStrips;
    }
    else{
        startIndex += leftoverStrips;
    }
    
    //Calculate the index to end for this task
    int endIndex = startIndex + (data->height / data->mpi_procs);
    //Allocate an additional strip if there are leftovers
    if(data->mpi_rank < leftoverStrips){
        endIndex++;
    }

    int rowsToRender = endIndex - startIndex;

    //////////////////////////////////////////////////////////////////////////
    //Rendering
    //////////////////////////////////////////////////////////////////////////

    //Allocate space for the image on the master.
    int renderSize = 3 * (rowsToRender ) * data->width;
    float * pixels = new float[renderSize];
    //Render this processor's portion of the scene.
    for( int row = 0; row < rowsToRender; row++)
    {
        for( int col = 0; col < data->width; col++ )
        {

            //Calculate the index into the array.
            int index = 3 * ( row * data->width + col );

            //Call the function to shade the pixel.
            shadePixel(&(pixels[index]), row + startIndex, col, data);
        }
    }
    stopTime = MPI_Wtime();

    ////////////////////////////////////////////////////////////////////
    //Communication
    ////////////////////////////////////////////////////////////////////

    //Create packet on heap
    //Size of render + timing data + index of rows calculated
    int packetSize = (renderSize) + 3;
    float * packet = new float[packetSize];
    //Communication time
    packet[0] = stopTime - startTime;
    packet[1] = (float)startIndex;
    packet[2] = (float)endIndex;
    //Render
    memcpy(&(packet[3]), pixels, (renderSize * sizeof(float)));

    //Send the packet
    MPI_Send(packet, packetSize, MPI_FLOAT, 0, 'p', MPI_COMM_WORLD);
    ///////////////////////////////////////////////////////////////////
    //Clean-up
    ///////////////////////////////////////////////////////////////////
    delete[] pixels;
    delete[] packet;
}

//Returns true if number is a perfect square
int perfectSquareSlave(double num) {
    double squareRoot = sqrt(num);
    return !(squareRoot - floor(squareRoot));
}

void slaveStaticBlocks(ConfigData * data) {

    double startTime, stopTime, totalCompTime;

    if(!perfectSquareSlave(data->mpi_procs)){
        return;
    }

    startTime = MPI_Wtime();

    ///////////////////////////////////////////////////////////////////
    //Index Calculation
    ///////////////////////////////////////////////////////////////////

    //Calculate initial start row and column (x,y) for this task

    //Should be an int since perfect square
    int procPerRow = (int) sqrt(data->mpi_procs);//How many processors are in a row
    int procPerCol = procPerRow; //How many processors are in a row, perfect square
    int localRowPosition = data->mpi_rank % procPerRow;//Which block of Y values this belongs to
    int localColPosition = data->mpi_rank / procPerCol;//Block of X values

    int startX = ((data->height / procPerRow) * localRowPosition);
    int startY= ((data->width / procPerCol) * localColPosition);
    int endX = ((data->height / procPerRow) * (localRowPosition + 1));
    int endY = ((data->width / procPerCol) * (localColPosition + 1));
    //Useful for packet size calculations
    int maxEndX = endX + (data->width % procPerRow);
    int maxEndY = endY + (data->height % procPerCol);
    int maxRangeX = maxEndX - startX;
    int maxRangeY = maxEndY - startY;
    //Check if height and width are not cleanly divisible, and handle it.
    if(localRowPosition == procPerRow - 1){
        endX = maxEndX;
    }
    if(localColPosition == procPerCol - 1){
        endY = maxEndY;
    }
    int rangeX = endX - startX;
    int rangeY = endY - startY;


    //////////////////////////////////////////////////////////////////////////
    //Rendering
    //////////////////////////////////////////////////////////////////////////

    //Allocate space for the image on the master.
    int renderSize = 3 * rangeX * rangeY;
    float* pixels = new float[renderSize];

    //Render this processor's portion of the scene.
    for( int row = 0; row < rangeY; ++row )
    {
        for( int col = 0; col < rangeX; ++col )
        {

            //Calculate the index into the array.
            int index = 3 * ( row * rangeX + col );

            //Call the function to shade the pixel.
            shadePixel(&(pixels[index]), row + startY, col + startX, data);
        }
    }
    


    stopTime = MPI_Wtime();
    totalCompTime = stopTime - startTime;

    ////////////////////////////////////////////////////////////////////
    //Communication
    ////////////////////////////////////////////////////////////////////

    //Create packet on heap
    //Size of render + timing data
    int packetSize = ( 3 * maxRangeX * maxRangeY) + 5;
    float * packet = new float[packetSize];

    //Communication time
    packet[0] = totalCompTime;
    packet[1] = (float)startX;
    packet[2] = (float)startY;
    packet[3] = (float)endX;
    packet[4] = (float)endY;
    //Render
    memcpy(&(packet[5]), pixels, (renderSize * sizeof(float)));

    //Send the packet
    MPI_Send(packet, packetSize, MPI_FLOAT, 0, 'p', MPI_COMM_WORLD);
    ///////////////////////////////////////////////////////////////////
    //Clean-up
    ///////////////////////////////////////////////////////////////////
    delete[] pixels;
    delete[] packet;
}

void slaveStaticCyclesHorizontal(ConfigData * data) {

    double startTime, stopTime, totalCompTime;

    startTime = MPI_Wtime();

    ///////////////////////////////////////////////////////////////////
    //Index Calculation
    ///////////////////////////////////////////////////////////////////

    

    //Unused for this method

    

    //////////////////////////////////////////////////////////////////////////
    //Rendering
    //////////////////////////////////////////////////////////////////////////

    //Allocate space for the image on the master.
    int renderSize = 3 * data->width * data->height;
    float* pixels = new float[renderSize];


    //Render this processor's portion of the scene.
    for( int row = 0; row < data->height; ++row )
    {
        //Check if this row belongs to this processor
        if( ( (row / data->cycleSize) % data->mpi_procs ) == data->mpi_rank){
            //Row is in a cycle owned by this processor
            for( int col = 0; col < data->width; ++col )
            {

                //Calculate the index into the array.
                int index = 3 * ( row * data->width + col );

                //Call the function to shade the pixel.
                shadePixel(&(pixels[index]), row, col, data);
            }
        }
    }

    stopTime = MPI_Wtime();
    totalCompTime = stopTime - startTime;

    ////////////////////////////////////////////////////////////////////
    //Communication
    ////////////////////////////////////////////////////////////////////
    //Create packet on heap
    //Size of render + timing data
    int packetSize = 1 + ( 3 * data->width * data->height);
    float * packet = new float[packetSize];

    //Communication time
    packet[0] = totalCompTime;

    //Render
    memcpy(&(packet[1]), pixels, (renderSize * sizeof(float)));

    //Send the packet
    MPI_Send(packet, packetSize, MPI_FLOAT, 0, 'p', MPI_COMM_WORLD);
    ///////////////////////////////////////////////////////////////////
    //Clean-up
    ///////////////////////////////////////////////////////////////////
    delete[] pixels;
    delete[] packet;
}

void slaveDynamic(ConfigData * data) {

    double startTime, stopTime, totalCompTime;
    MPI_Status status;


    ///////////////////////////////////////////////////////////////////
    //Packet Allocation
    ///////////////////////////////////////////////////////////////////

    //Allocate space for the image on the master.
    int renderSize = 3 * data->width * data->height;
    float* pixels = new float[renderSize];

    //initialize render to 0
    for(int i = 0; i < renderSize; i++){
        pixels[i] = 0.0;
    }
        
    //Create space for packets on the heap
    //Packet containing rendering instructions
    //0 - finsihedFlag, 1 - startX, 2 - startY, 3 - endX, 4 - endY
    int instructionPacketSize = 5;
    int * instructionPacket = new int[instructionPacketSize];

    //Packet to request new instructions
    int requestPacketSize = 1; //rank
    int * requestPacket = new int[requestPacketSize];
    requestPacket[0] = data->mpi_rank; //never changes

    //Packet containing render data along with computation time
    int renderPacketSize = renderSize + 1;
    float * renderPacket = new float[renderPacketSize];

    

    //////////////////////////////////////////////////////////////////////////
    //Rendering
    //////////////////////////////////////////////////////////////////////////


    startTime = MPI_Wtime();

    int finishedFlag = 0;

    //Send initial request
    MPI_Send(requestPacket, requestPacketSize, MPI_INT, 0, 'r', MPI_COMM_WORLD);

    while(!finishedFlag){
        //Recieve a block from the queue
        MPI_Recv(instructionPacket, instructionPacketSize, MPI_INT, 0, 'i', MPI_COMM_WORLD, &status);

        //Extract instructions
        finishedFlag = instructionPacket[0];
        int startX = instructionPacket[1];
        int startY = instructionPacket[2];
        int endX = instructionPacket[3];
        int endY = instructionPacket[4];

        //Render this processor's portion of the scene.
        for( int row = startY; row < endY; ++row )
        {
            //Only render what is assigned to this processor
            for( int col = startX; col < endX; ++col )
            {

                //Calculate the index into the array.
                int index = 3 * ( row * data->width + col );

                //Call the function to shade the pixel.
                shadePixel(&(pixels[index]), row, col, data);
            }
        }

        //Render is complete, request a new packet
        MPI_Send(requestPacket, requestPacketSize, MPI_INT, 0, 'r', MPI_COMM_WORLD);

    }

    stopTime = MPI_Wtime();
    totalCompTime = stopTime - startTime;


    //Communication time
    renderPacket[0] = totalCompTime;

    //Render
    memcpy(&(renderPacket[1]), pixels, (renderSize * sizeof(float)));

    //Send the packet
    MPI_Send(renderPacket, renderPacketSize, MPI_FLOAT, 0, 'p', MPI_COMM_WORLD);
    ///////////////////////////////////////////////////////////////////
    //Clean-up
    ///////////////////////////////////////////////////////////////////
    delete[] pixels;
    delete[] instructionPacket;
    delete[] requestPacket;
    delete[] renderPacket;
}