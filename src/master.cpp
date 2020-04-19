//This file contains the code that the master process will execute.

#include <iostream>
#include <mpi.h>
#include <math.h>

#include "RayTrace.h"
#include "master.h"

void masterMain(ConfigData* data)
{
    //Depending on the partitioning scheme, different things will happen.
    //You should have a different function for each of the required 
    //schemes that returns some values that you need to handle.
    
    //Allocate space for the image on the master.
    float* pixels = new float[3 * data->width * data->height];
    
    //Execution time will be defined as how long it takes
    //for the given function to execute based on partitioning
    //type.
    double renderTime = 0.0, startTime, stopTime;
    double compTime;

	//Add the required partitioning methods here in the case statement.
	//You do not need to handle all cases; the default will catch any
	//statements that are not specified. This switch/case statement is the
	//only place that you should be adding code in this function. Make sure
	//that you update the header files with the new functions that will be
	//called.
	//It is suggested that you use the same parameters to your functions as shown
	//in the sequential example below.
    switch (data->partitioningMode)
    {
        case PART_MODE_NONE:
            //Call the function that will handle this.
            startTime = MPI_Wtime();
            masterSequential(data, pixels);
            stopTime = MPI_Wtime();
            break;
        case PART_MODE_STATIC_STRIPS_HORIZONTAL:
            //Call the function that will handle this.
            startTime = MPI_Wtime();
            compTime = masterStaticStripsHorizontal(data, pixels);
            stopTime = MPI_Wtime();
            break;
        case PART_MODE_STATIC_STRIPS_VERTICAL:
            //Call the function that will handle this.
            startTime = MPI_Wtime();
            masterSequential(data, pixels);
            stopTime = MPI_Wtime();
            break;
        case PART_MODE_STATIC_BLOCKS:
            //Call the function that will handle this.
            startTime = MPI_Wtime();
            compTime = masterStaticBlocks(data, pixels);
            stopTime = MPI_Wtime();
            break;
        case PART_MODE_STATIC_CYCLES_HORIZONTAL:
            //Call the function that will handle this.
            startTime = MPI_Wtime();
            compTime = masterStaticCyclesHorizontal(data, pixels);
            stopTime = MPI_Wtime();
            break;
        case PART_MODE_STATIC_CYCLES_VERTICAL:
            //Call the function that will handle this.
            startTime = MPI_Wtime();
            masterSequential(data, pixels);
            stopTime = MPI_Wtime();
            break;
        case PART_MODE_DYNAMIC:
            //Call the function that will handle this.
            startTime = MPI_Wtime();
            compTime = masterDynamic(data, pixels);
            stopTime = MPI_Wtime();
            break;
        default:
            std::cout << "This mode (" << data->partitioningMode;
            std::cout << ") is not currently implemented." << std::endl;
            break;
    }

    /////////////////////////////////////////////////////////////
    //Report Results
    /////////////////////////////////////////////////////////////
    renderTime = stopTime - startTime;
    double commTime = renderTime - compTime;
    double commCompRatio = commTime / compTime;
    std::cout << "Execution Time: " << renderTime << " seconds" << std::endl;
    std::cout << "Total Communication Time: " << commTime << " seconds" <<std::endl;
    std::cout << "Total Computation Time: " << compTime << " seconds" << std::endl;
    std::cout << "Communication to Computation Ratio: " <<  commCompRatio << " seconds" <<std::endl << std::endl;

    //After this gets done, save the image.
    std::cout << "Image will be save to: ";
    std::string file = generateFileName(data);
    std::cout << file << std::endl;
    savePixels(file, pixels, data);

    //Delete the pixel data.
    delete[] pixels; 
}

void masterSequential(ConfigData* data, float* pixels)
{
    //Start the computation time timer.
    double computationStart = MPI_Wtime();

    //Render the scene.
    for( int i = 0; i < data->height; ++i )
    {
        for( int j = 0; j < data->width; ++j )
        {
            int row = i;
            int column = j;

            //Calculate the index into the array.
            int baseIndex = 3 * ( row * data->width + column );

            //Call the function to shade the pixel.
            shadePixel(&(pixels[baseIndex]),row,j,data);
        }
    }

    //Stop the comp. timer
    double computationStop = MPI_Wtime();
    double computationTime = computationStop - computationStart;

    //After receiving from all processes, the communication time will
    //be obtained.
    double communicationTime = 0.0;

    //Print the times and the c-to-c ratio
	//This section of printing, IN THIS ORDER, needs to be included in all of the
	//functions that you write at the end of the function.
    std::cout << "Total Computation Time: " << computationTime << " seconds" << std::endl;
    std::cout << "Total Communication Time: " << communicationTime << " seconds" << std::endl;
    double c2cRatio = communicationTime / computationTime;
    std::cout << "C-to-C Ratio: " << c2cRatio << std::endl;
}



double masterStaticStripsHorizontal(ConfigData * data, float * pixels) {

    double startTime, stopTime, totalCompTime;
    MPI_Status status;

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
    totalCompTime = stopTime - startTime;

    ////////////////////////////////////////////////////////////////////
    //Communication
    ////////////////////////////////////////////////////////////////////

    //Create packet on heap
    //Size of render + timing data
    int packetSize = (renderSize) + 3;
    float * packet = new float[packetSize];

    int i = 1;
    //Recieve portions of render from slaves
    for(i = 1; i < leftoverStrips; i++){
        
        //Recieve packet
        MPI_Recv(packet, packetSize, MPI_FLOAT, i, 'p', MPI_COMM_WORLD, &status);

        //Extract results
        //Comp time is the time of the longest render
        if(packet[0] > totalCompTime){
            totalCompTime = packet[0];
        }
        int slaveStart = (int)packet[1];
        int totalSize = 3 * data->width * ((data->height / data->mpi_procs) + 1);
        int slaveIndex = 3 * slaveStart * data->width;
        memcpy(&(pixels[slaveIndex]), &packet[3], totalSize * sizeof(float));

    }
        //Recieve portions of render from slaves
    for(; i < data->mpi_procs; i++){
        //Recieve packet
        MPI_Recv(packet, packetSize, MPI_FLOAT, i, 'p', MPI_COMM_WORLD, &status);
        //Extract results
        //Comp time is the time of the longest render
        if(packet[0] > totalCompTime){
            totalCompTime = packet[0];
        }
        int slaveStart = (int)packet[1];
        int totalSize = 3 * data->width * ((data->height / data->mpi_procs));
        int slaveIndex = 3 * slaveStart * data->width;
        memcpy(&(pixels[slaveIndex]), &packet[3], totalSize * sizeof(float));

    }

    ///////////////////////////////////////////////////////////////////
    //Clean-up
    ///////////////////////////////////////////////////////////////////
    delete[] packet;

    return totalCompTime;
}

//Returns true if number is a perfect square
int perfectSquare(double num) {
    double squareRoot = sqrt(num);
    return !(squareRoot - floor(squareRoot));
}

double masterStaticBlocks(ConfigData * data, float * pixels) {

    double startTime, stopTime, totalCompTime;
    MPI_Status status;

    if(!perfectSquare(data->mpi_procs)){
        std::cout << "Number of processors must be a perfect square!" << std::endl;
        return 0.0;
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

    //Recieve portions of render from slaves
    for(int i = 1; i < data->mpi_procs; i++){
        
        //Recieve packet
        MPI_Recv(packet, packetSize, MPI_FLOAT, i, 'p', MPI_COMM_WORLD, &status);

        //Extract results
        //Comp time is the time of the longest render
        if(packet[0] > totalCompTime){
            totalCompTime = packet[0];
        }
        int slaveStartX = (int)packet[1];
        int slaveStartY = (int)packet[2];
        int slaveEndX = (int)packet[3];
        int slaveEndY = (int)packet[4];
        int slaveRangeX = slaveEndX - slaveStartX;
        int slaveRangeY = slaveEndY - slaveStartY;

        //copy over row by row
        for(int row = 0; row < slaveRangeY; row++){
            for(int col = 0; col < slaveRangeX; col++){

                int pixelsIndex = 3 * (((slaveStartY + row) * data->width) + slaveStartX + col);
                int packetIndex = 5 + (3 * ((row * slaveRangeX)+ col));
                pixels[pixelsIndex] = packet[packetIndex];
                pixels[pixelsIndex + 1] = packet[packetIndex + 1];
                pixels[pixelsIndex + 2] = packet[packetIndex + 2];
            }
        }

    }

    ///////////////////////////////////////////////////////////////////
    //Clean-up
    ///////////////////////////////////////////////////////////////////
    delete[] packet;

    return totalCompTime;
}

double masterStaticCyclesHorizontal(ConfigData * data, float * pixels) {

    double startTime, stopTime, totalCompTime;
    MPI_Status status;

    startTime = MPI_Wtime();

    ///////////////////////////////////////////////////////////////////
    //Index Calculation
    ///////////////////////////////////////////////////////////////////

    

    //Unused for this method

    

    //////////////////////////////////////////////////////////////////////////
    //Rendering
    //////////////////////////////////////////////////////////////////////////
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
    //Packets are the entire picture
    int packetSize = ( 3 * data->width * data->height) + 1;
    float * packet = new float[packetSize];

    //Recieve portions of render from slaves
    for(int i = 1; i < data->mpi_procs; i++){
        
        //Recieve packet
        MPI_Recv(packet, packetSize, MPI_FLOAT, i, 'p', MPI_COMM_WORLD, &status);

        //Extract results
        //Comp time is the time of the longest render
        if(packet[0] > totalCompTime){
            totalCompTime = packet[0];
        }

        //copy over row by row
        for(int row = 0; row < data->height; row++){
            if( ( (row / data->cycleSize) % data->mpi_procs ) == i){
                for(int col = 0; col < data->width; col++){

                    int pixelsIndex = 3 * ((row * data->width) + col);
                    int packetIndex = 1 + pixelsIndex;
                    pixels[pixelsIndex] = packet[packetIndex];
                    pixels[pixelsIndex + 1] = packet[packetIndex + 1];
                    pixels[pixelsIndex + 2] = packet[packetIndex + 2];
                }
            }
        }

    }
    ///////////////////////////////////////////////////////////////////
    //Clean-up
    ///////////////////////////////////////////////////////////////////
    delete[] packet;

    return totalCompTime;
}

double masterDynamic(ConfigData * data, float * pixels) {
    double totalCompTime = 0.0;
    MPI_Status status;


    ///////////////////////////////////////////////////////////////////
    //Index Calculation
    ///////////////////////////////////////////////////////////////////

    //Calculate number of blocks in rows and columns
    int numRows = data->height / data->dynamicBlockHeight;
    if(data->height % data->dynamicBlockHeight){
        numRows++;
    }
    int numCols = data->width / data->dynamicBlockWidth;
    if(data->width % data->dynamicBlockWidth){
        numCols++;
    }

    ///////////////////////////////////////////////////////////////////
    //Packet Allocation
    ///////////////////////////////////////////////////////////////////

    //Allocate space for the image on the master.
    int renderSize = 3 * data->width * data->height;

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
    float * renderPacket = new float[renderSize];

    ////////////////////////////////////////////////////////////////////
    //Communication
    ////////////////////////////////////////////////////////////////////


    //Recieve portions of render from slaves
    for(int row = 0; row < numRows; row++){

        for(int col = 0; col < numCols; col++){

            //Recieve a packet from a slave that is done computation
            MPI_Recv(requestPacket, requestPacketSize, MPI_INT, MPI_ANY_SOURCE, 'r', MPI_COMM_WORLD, &status);
            
            //Parse packet
            int requesterRank = requestPacket[0];
            std::cout << "Received request from " << requesterRank << std::endl;
            //Calculate and send new coordinates
            int startX = col * data->dynamicBlockWidth;
            int startY = row * data->dynamicBlockHeight;
            int endX = (col + 1) * data->dynamicBlockWidth;
            if( endX > data->width){
                endX = data->width;
            }
            int endY = (row + 1) * data->dynamicBlockHeight;
            if( endY > data->height){
                endY = data->height;
            }
            instructionPacket[0] = 0;//Finished flag, if true slave should end
            instructionPacket[1] = startX;
            instructionPacket[2] = startY;
            instructionPacket[3] = endX;
            instructionPacket[4] = endY;
            
            std::cout << "Sending " << requesterRank <<" block " << startX <<"," <<startY << std::endl;
            
            MPI_Send(instructionPacket, instructionPacketSize, MPI_INT, requesterRank, 'i', MPI_COMM_WORLD);
        }

    }

    std::cout << "===================================" << std::endl;

    //Tell slaves computation is finished
    for(int i = 1; i < data->mpi_procs; i++){

        //Recieve a packet from a slave that is done computation
        MPI_Recv(requestPacket, requestPacketSize, MPI_INT, MPI_ANY_SOURCE, 'r', MPI_COMM_WORLD, &status);
        
        std::cout << "Received request from " << requesterRank << std::endl;

        //Parse packet
        int requesterRank = requestPacket[0];

        instructionPacket[0] = 1;//Finished flag, if true slave should end
        instructionPacket[1] = 0;
        instructionPacket[2] = 0;
        instructionPacket[3] = 0;
        instructionPacket[4] = 0;

        std::cout << "Sending " << requesterRank <<" block " << startX <<"," <<startY << std::endl;

        MPI_Send(instructionPacket, instructionPacketSize, MPI_INT, requesterRank, 'i', MPI_COMM_WORLD);
    }

    std::cout << "===================================" << std::endl;

    //Recieve and save final results
    for(int i = 1; i < data->mpi_procs; i++){

        //Recieve a packet from a slave that is done computation
        MPI_Recv(renderPacket, renderPacketSize, MPI_INT, MPI_ANY_SOURCE, 'p', MPI_COMM_WORLD, &status);
        std::cout << "Recieved render " << i << std::endl;

        //Find total computation time
        if(renderPacket[0] > totalCompTime){
            totalCompTime = renderPacket[0];
        }

        //Save pixels
        for(int i = 0; i < renderSize; i++){
            pixels[i] += renderPacket[i + 1];
        }

    }
    ///////////////////////////////////////////////////////////////////
    //Clean-up
    ///////////////////////////////////////////////////////////////////
    delete[] instructionPacket;
    delete[] requestPacket;
    delete[] renderPacket;

    return totalCompTime;
}
