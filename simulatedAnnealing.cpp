#include "simulatedAnnealing.h"

void InitSA(pllInstance *tr, partitionList *pr, unsigned int &bestIterationScoreHits, int maxtrav, bool usingTBR) {
    tr->usingSA = true;
    pllTreeToNewick(tr->best_tree_string_sa, tr, pr,
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

    
    if (tr->pureSA) {
        tr->numOfDelta = 1; 
        tr->sumOfDelta = tr->bestParsimony/200;
    }

    bestIterationScoreHits = 1;
    tr->coolingTimes = 0;
    tr->temperature = tr->startTemp;
    tr->stepCount = 0;
    // tr->limitTrees = (tr->mxtips + tr->mxtips - 2) * 5 * (1<<(maxtrav-1)) / tr->maxCoolingTimes;
    tr->limitTrees = tr->maxNumTree / 10 / tr->maxCoolingTimes + ((tr->maxNumTree / 10 % tr->maxCoolingTimes) > 0);
    // printf("aaaa %d\n", (tr->mxtips + tr->mxtips - 2) * 5 * (1<<(maxtrav-1)) / tr->maxCoolingTimes);
    // printf("bbbb %d\n", tr->limitTrees);

    if (usingTBR) tr->limitTrees *= maxtrav-1;
}

bool checkContinueSA(pllInstance *tr) {
    return tr->coolingTimes < tr->maxCoolingTimes;
}

void applyBestTree(pllInstance *tr) {
    char* bestTreeString = tr->best_tree_string_sa;
    pllNewickTree *tmpTree = pllNewickParseString(bestTreeString);
    pllTreeInitTopologyNewick(tr, tmpTree, PLL_TRUE);
    pllNewickParseDestroy(&tmpTree);
}

void checkDecreaseTemp(pllInstance *tr) {
    if (tr->stepCount == tr->limitTrees)  {
        switch (tr->coolingSchedule)
        {
            case LINEAR_ADDITIVE_COOLING_PLL:
                tr->temperature -= tr->coolingAmount;
            break;

            case LINEAR_MULTIPLICATIVE_COOLING_PLL:
                tr->temperature = tr->startTemp / (1 + tr->coolingAmount * tr->coolingTimes);
            break;

            case EXPONENTIAL_ADDITIVE_COOLING_PLL:
            {
                double deltaTemp = tr->startTemp - tr->finalTemp;
                tr->temperature = tr->finalTemp + deltaTemp/(1.0 + exp(2.0*log(deltaTemp)/tr->maxCoolingTimes*(tr->coolingTimes - 0.5*tr->maxCoolingTimes))); 
            }
            break;

            case EXPONENTIAL_MULTIPLICATIVE_COOLING_PLL:
                tr->temperature *= tr->coolingAmount;
            break;
        }

        tr->coolingTimes++;
        tr->stepCount = 0;
    }
}

void checkAcceptTree(pllInstance *tr, partitionList *pr, int perSiteScores, bool& haveChange, unsigned int mp, unsigned long& bestTreeScoreHits) {
    if (mp < tr->bestParsimony) {
        pllTreeToNewick(tr->best_tree_string_sa, tr, pr,
            tr->start->back, PLL_FALSE, PLL_TRUE, 0, 0, 0,
            PLL_SUMMARIZE_LH, 0, 0);
        tr->sumOfDelta += (tr->bestParsimony - mp) * tr->numOfDelta;
        tr->bestParsimony = mp;
    } else {
        tr->sumOfDelta += mp - tr->bestParsimony;
        tr->numOfDelta++;
    }

    if (mp < tr->currentParsimony) {
        tr->currentParsimony = mp;
        bestTreeScoreHits = 1;
        haveChange = true;
    } else if (mp == tr->currentParsimony) {
        bestTreeScoreHits++;
        if (random_double() < 1.0 / bestTreeScoreHits) {
            haveChange = true;
        }
    } else {
        int delta = tr->bestParsimony - mp;
        double tmp = double(delta)/double(double(tr->sumOfDelta) / tr->numOfDelta * tr->temperature);
        double probability = exp(tmp);
        if (random_double() <= probability) {
            haveChange = true;
        } 
    }

    tr->stepCount++; 
    checkDecreaseTemp(tr);   
}