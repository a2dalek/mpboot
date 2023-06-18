#include "hillclimbing.h"
#include "sprparsimony.h"
#include "tbrparsimony.h"

void hillClimbing(pllInstance *tr, partitionList *pr, int mintrav,
                  int maxtrav, IQTree *_iqtree, OptimizeMethod method,
                  int &perSiteScores, unsigned int &randomMP,
                  unsigned int &startMP, bool &is_first_loop, unsigned int &bestIterationScoreHits) {

    if (method == OptimizeMethod::TBR) {

    } else if (method == OptimizeMethod::SPR) {
        if (tr->pureSA) tr->usingSA = true;
        pllTreeToNewick(tr->tree_string, tr, pr,
                      tr->start->back, PLL_FALSE, PLL_TRUE, 0, 0, 0,
                      PLL_SUMMARIZE_LH, 0, 0);
        switch (tr->coolingSchedule)
        {
            case LINEAR_ADDITIVE_COOLING_PLL:
                tr->coolingAmount = (tr->startTemp - tr->finalTemp) / tr->maxCoolingTimes;
            break;

            case LINEAR_MULTIPLICATIVE_COOLING_PLL:
                tr->coolingAmount = ((tr->startTemp / tr->finalTemp) - 1.0) / tr->maxCoolingTimes;
            break;


            case EXPONENTIAL_MULTIPLICATIVE_COOLING_PLL:
                tr->coolingAmount = pow(tr->finalTemp / tr->startTemp, 1.0 / tr->maxCoolingTimes);
            break;
        }

        bestIterationScoreHits = 1;
        tr->coolingTimes = 0;
        tr->temperature = tr->startTemp;
        tr->stepCount = 0;
        tr->limitTrees = (tr->mxtips + tr->mxtips - 2) * 5 * (1<<(maxtrav-1)) / tr->maxCoolingTimes;
        tr->numOfDelta = tr->bestParsimony/200;
        tr->sumOfDelta = 1;

        if (tr->usingSA) {
            while (tr->coolingTimes <= tr->maxCoolingTimes) {
                pllOptimizeSprParsimony(tr, pr, mintrav, maxtrav, _iqtree, perSiteScores, randomMP, startMP, is_first_loop, bestIterationScoreHits);
            }
        } else {
            pllOptimizeSprParsimony(tr, pr, mintrav, maxtrav, _iqtree, perSiteScores, randomMP, startMP, is_first_loop, bestIterationScoreHits);
        }

        pllCheckIterImprove(tr);

        char* bestTreeString = tr->tree_string;
        pllNewickTree *tmpTree = pllNewickParseString(bestTreeString);
        pllTreeInitTopologyNewick(tr, tmpTree, PLL_TRUE);
        pllNewickParseDestroy(&tmpTree);
    }
}

void pllOptimizeParsimony(pllInstance *tr, partitionList *pr, int mintrav,
                            int maxtrav, IQTree *_iqtree, OptimizeMethod method) {

    int perSiteScores;
    unsigned int randomMP, startMP, bestIterationScoreHits;

    bool is_first_loop = true;

    if (tr->plusSA) {

    } else if (tr->pureSA) {
        hillClimbing(tr, pr, mintrav, maxtrav, _iqtree, method, perSiteScores, randomMP, startMP, is_first_loop, bestIterationScoreHits);
    } else {
        hillClimbing(tr, pr, mintrav, maxtrav, _iqtree, method, perSiteScores, randomMP, startMP, is_first_loop, bestIterationScoreHits);
    }
}