#include "hillclimbing.h"
#include "sprparsimony.h"
#include "tbrparsimony.h"
#include "simulatedAnnealing.h"
#include <iostream>

OptimizeMethod method;

bool checkContinue(pllInstance *tr, unsigned int &randomMP, unsigned int &startMP) {
    if (tr->usingSA) {
        return checkContinueSA(tr);
    } else {
        return randomMP < startMP;
    }
}

void doHillClimbing(pllInstance *tr, partitionList *pr, int mintrav,
                  int maxtrav, IQTree *_iqtree, int &perSiteScores, unsigned int &randomMP,
                  unsigned int &startMP, unsigned int &bestIterationScoreHits) {
    do {
        if (method == OptimizeMethod::TBR) {
            pllOptimizeTbrParsimony(tr, pr, mintrav, maxtrav, _iqtree, perSiteScores, randomMP, startMP, bestIterationScoreHits);
        } else {
            pllOptimizeSprParsimony(tr, pr, mintrav, maxtrav, _iqtree, perSiteScores, randomMP, startMP, bestIterationScoreHits);
        }
    } while (checkContinue(tr, randomMP, startMP));
}

void InitPllOptimizeParsimony(pllInstance *tr, partitionList *pr, int mintrav,
                  int maxtrav, IQTree *_iqtree, int &perSiteScores, unsigned int &randomMP,
                  unsigned int &startMP, unsigned int &bestIterationScoreHits) {

    if (method == OptimizeMethod::TBR) {
        InitPllOptimizeTbrParsimony(tr, pr, mintrav, maxtrav, _iqtree, perSiteScores, randomMP, startMP, bestIterationScoreHits);
    } else if (method == OptimizeMethod::SPR) {
        InitPllOptimizeSprParsimony(tr, pr, mintrav, maxtrav, _iqtree, perSiteScores, randomMP, startMP, bestIterationScoreHits);
    }
}

void pllCheckIterImprove(pllInstance *tr) {
    if (method == OptimizeMethod::TBR) {
        pllCheckIterImproveTBR(tr);
    } else if (method == OptimizeMethod::SPR) {
        pllCheckIterImproveSPR(tr);
    }
}

void pllOptimizeParsimony(pllInstance *tr, partitionList *pr, int mintrav,
                            int maxtrav, IQTree *_iqtree, OptimizeMethod _method) {
                        
    int perSiteScores;
    unsigned int randomMP, startMP, bestIterationScoreHits;

    method = _method;

    tr->normalIter++;
    InitPllOptimizeParsimony(tr, pr, mintrav, maxtrav, _iqtree, perSiteScores, randomMP, startMP, bestIterationScoreHits);
    if (tr->plusSA && tr->normalIter > 10) {
        doHillClimbing(tr, pr, mintrav, maxtrav, _iqtree, perSiteScores, randomMP, startMP, bestIterationScoreHits);
        InitSA(tr, pr, bestIterationScoreHits, maxtrav, (method == OptimizeMethod::TBR));
        doHillClimbing(tr, pr, mintrav, maxtrav, _iqtree, perSiteScores, randomMP, startMP, bestIterationScoreHits);

        // pllCheckIterImprove(tr);
        applyBestTree(tr);
    } else if (tr->pureSA && tr->normalIter > 10) {
        InitSA(tr, pr, bestIterationScoreHits, maxtrav, (method == OptimizeMethod::TBR));
        doHillClimbing(tr, pr, mintrav, maxtrav, _iqtree, perSiteScores, randomMP, startMP, bestIterationScoreHits);

        // pllCheckIterImprove(tr);
        applyBestTree(tr);
        // std::cout<<tr->cnt1 << " " << tr->cnt2 << "\n";
    } else {
        doHillClimbing(tr, pr, mintrav, maxtrav, _iqtree, perSiteScores, randomMP, startMP, bestIterationScoreHits);
    }
}