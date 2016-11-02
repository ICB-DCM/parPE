#ifndef LOCAL_OPTIMIZATION_H
#define LOCAL_OPTIMIZATION_H

#include <IpStdCInterface.h>

#include "dataprovider.h"
#include <include/udata.h>
#include <include/edata.h>
#include <include/rdata.h>
#include "dataprovider.h"

#include "logger.h"

#define IPTOPT_LOG_FILE "/home/dweindl/src/CanPathProSSH/dw/ipopt.log"

void getLocalOptimum(datapath dataPath);

#endif
