#pragma once

#ifndef SIMULATEDANNEALING_H
#define SIMULATEDANNEALING_H

#include "iqtree.h"
#include "pllrepo/src/pll.h"

void InitSA(pllInstance *tr, partitionList *pr, unsigned int &bestIterationScoreHits, int maxtrav, bool usingTBR = false);

bool checkContinueSA(pllInstance *tr);

void applyBestTree(pllInstance *tr);

bool checkAcceptTree(pllInstance *tr, partitionList *pr, int perSiteScores, bool& haveChange, unsigned int mp, unsigned long& bestTreeScoreHits);

#endif