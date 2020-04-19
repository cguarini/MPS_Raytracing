#ifndef __SLAVE_PROCESS_H__
#define __SLAVE_PROCESS_H__

#include "RayTrace.h"

void slaveMain( ConfigData *data );
void renderStaticStripsHorizontal(ConfigData * data);
void slaveStaticBlocks(ConfigData * data);
void slaveStaticCyclesHorizontal(ConfigData * data);
void slaveDynamic(ConfigData * data);

#endif
