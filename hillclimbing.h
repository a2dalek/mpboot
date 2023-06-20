#ifndef HILLCLIMBING_H
#define HILLCLIMBING_H

#include "iqtree.h"
#include "pllrepo/src/pll.h"

enum class OptimizeMethod {
    TBR,
    SPR,
};

void pllOptimizeParsimony(pllInstance *tr, partitionList *pr, int mintrav,
                            int maxtrav, IQTree *_iqtree, OptimizeMethod _method);

#endif