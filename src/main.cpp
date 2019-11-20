#include "readData.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <chrono>
#include <limits>

using namespace std;

enum NL{N1, N2, N3, N4, N5};

double ** costMatrix;
int dimension;

typedef struct{
    double cost;
    int nodeInserted;
    int edgeRemoved;
} InsertionInfo;


void printCostMatrix() {
    std::cout << "dimension: " << dimension << endl;
    for (size_t i = 1; i <= dimension; i++) {
        for (size_t j = 1; j <= dimension; j++) {
            std::cout << costMatrix[i][j] << " ";
        }
        std::cout << endl;
    }
    cout << endl << endl;
}


bool compareCosts(const InsertionInfo &a, const InsertionInfo &b) {
    return (a.cost < b.cost);
}

vector<int> construction(double alpha) {

    vector<int> solution = {1};

    vector<int> candidateList;

    const int subtourSize = 3;

    //Fills candidates list according to dimension
    for(int i = 2; i < dimension + 1; i++){ //dimension -1(node 1 already in) +2(starts at 2)
        candidateList.push_back(i);
    }

    //Initial subtour
    for(int i = 0; i < subtourSize; i++){
        int j = rand() % candidateList.size();
        solution.push_back(candidateList[j]);
        candidateList.erase(candidateList.begin() + j);
    }
    solution.push_back(1);


    while(!candidateList.empty())
    {
        vector<InsertionInfo> insertionCosts((solution.size() - 1) * candidateList.size());

        //Calculates and sorts insertion costs
        for(int i = 0, j = 1, l = 0; i < solution.size() - 1; i++, j++){
            for(auto k : candidateList){
                insertionCosts[l].cost = costMatrix[solution[i]][k] + costMatrix[solution[j]][k] - costMatrix[solution[j]][solution[i]];
                insertionCosts[l].nodeInserted = k;
                insertionCosts[l].edgeRemoved = j;
                l++;
            }
        }

        sort(insertionCosts.begin(), insertionCosts.end(), compareCosts);

        //Inserts according to cost
        int listFraction = alpha * insertionCosts.size();
        if(listFraction < 1)
            listFraction = 1;

        int x = rand() % listFraction;
        solution.insert(solution.begin() + insertionCosts[x].edgeRemoved, insertionCosts[x].nodeInserted);

        //Removes candidate from list
        for(int i = 0; i < candidateList.size(); i++) {
            if(candidateList[i] == insertionCosts[x].nodeInserted) {
                candidateList.erase(candidateList.begin() + i);
            }
        }
    }

    return solution;
}


double solutionCost (vector <int> solution){
    double cost = 0;
    for(int i = 0; i < solution.size() - 1; i++){
        cost += costMatrix[solution[i]][solution[i+1]];
    }
    return cost;
}


vector<int> swap (vector<int> s, double *bestDelta){
    int best_i = 0, best_j = 0;
    *bestDelta = 0; //0 = no change
    double delta = 0, semiDelta = 0;

    for(int j = 3; j < s.size() - 2; j++){

        semiDelta =  -costMatrix[s[j]][s[j-1]] -costMatrix[s[j]][s[j+1]]; //relies only on j

        for(int i = 1; i < j - 1; i++){

            delta = semiDelta + costMatrix[s[i]][s[j-1]] + costMatrix[s[i]][s[j+1]]
            + costMatrix[s[j]][s[i-1]] + costMatrix[s[j]][s[i+1]]
            - costMatrix[s[i]][s[i-1]] - costMatrix[s[i]][s[i+1]];

            if(delta < *bestDelta){
                *bestDelta = delta;

                best_i = i;
                best_j = j;
            }
        }
    }

    if(*bestDelta < 0){
        if(best_i < best_j){
            s.insert(s.begin() + best_i + 1, s[best_j]);
            s.insert(s.begin() + best_j + 2, s[best_i]);
            s.erase(s.begin() + best_i);
            s.erase(s.begin() + best_j);
        }else if(best_j < best_i){
            s.insert(s.begin() + best_j + 1, s[best_i]);
            s.insert(s.begin() + best_i + 2, s[best_j]);
            s.erase(s.begin() + best_j);
            s.erase(s.begin() + best_i);
        }
    }else{
        *bestDelta = 0;
    }

    return s;
}

vector<int> flip (vector<int> s, double *bestDelta){
    int best_i = 0, best_j = 0;
    *bestDelta = 0; //0 = no change
    double delta = 0, semiDelta = 0;

    for(int j = 3; j < s.size() - 1; j++){
        
        semiDelta = -costMatrix[s[j-1]][s[j]]; //relies only on j

        for(int i = 1; i < j - 1; i++){

            delta = semiDelta + costMatrix[s[i]][s[j]] +costMatrix[s[j-1]][s[i-1]] -costMatrix[s[i-1]][s[i]];

            if(delta < *bestDelta){
                *bestDelta = delta;

                best_i = i;
                best_j = j;
            }
        }
    }

    if(*bestDelta < 0){
        for(int k = best_i, l = 0; k < best_j; k++, l++){
            s.insert(s.begin() + best_i, s[k+l]);
        }
        for(int k = best_i; k < best_j; k++){
            s.erase(s.begin() + best_j);
        }
    }else{
        *bestDelta = 0;
    }

    return s;
}

vector<int> reinsertion (vector<int> s, double *bestDelta){
    int best_i = 0, best_j = 0;
    *bestDelta = 0; //0 = no change
    double delta = 0, semiDelta = 0;

    for(int j = 1; j < s.size(); j++){

        semiDelta = -costMatrix[s[j]][s[j-1]]; //relies only on j

        for(int i = 1; i < s.size() - 1; i++){
            if(i == j || j == i+1 || i == j+1) continue;

            delta = semiDelta + costMatrix[s[i]][s[j]] + costMatrix[s[i]][s[j-1]] + costMatrix[s[i-1]][s[i+1]]
            - costMatrix[s[i]][s[i-1]] - costMatrix[s[i]][s[i+1]];

            if(delta < *bestDelta){
                *bestDelta = delta;

                best_i = i;
                best_j = j;
            }

        }
    }

    if(*bestDelta < 0){
        s.insert(s.begin() + best_j, s[best_i]);
        if(best_j < best_i) best_i++;
        s.erase(s.begin() + best_i);
    }else{
        *bestDelta = 0;
    }

    return s;
}

vector<int> oropt2 (vector<int> s, double *bestDelta){
    int best_i = 0, best_j = 0;
    *bestDelta = 0; //0 = no change
    double delta = 0, semiDelta = 0;

    for(int j = 1; j < s.size(); j++){
        semiDelta = -costMatrix[s[j]][s[j-1]]; //relies only on j

        for(int i = 1; i < s.size() - 2; i++){

            if (i == j || j == i+1 || j == i+2) continue;

            delta = semiDelta + costMatrix[s[i]][s[j-1]] + costMatrix[s[i+1]][s[j]] + costMatrix[s[i-1]][s[i+2]]
            - costMatrix[s[i]][s[i-1]] - costMatrix[s[i+1]][s[i+2]];

            if(delta < *bestDelta){
                *bestDelta = delta;

                best_i = i;
                best_j = j;
            }

        }
    }

    if(*bestDelta < 0){
        if(best_i < best_j){
            s.insert(s.begin() + best_j, s[best_i+1]);
            s.insert(s.begin() + best_j, s[best_i]);
            s.erase(s.begin() + best_i);
            s.erase(s.begin() + best_i);
        }
        else if(best_j < best_i){
            s.insert(s.begin() + best_j, s[best_i+1]);
            s.insert(s.begin() + best_j, s[best_i+1]);
            s.erase(s.begin() + best_i+2);
            s.erase(s.begin() + best_i+2);
        }
    }else{
        *bestDelta = 0;
    }

    return s;
}

vector<int> oropt3 (vector<int> s, double *bestDelta){
    int best_i = 0, best_j = 0;
    *bestDelta = 0; //0 = no change
    double delta = 0, semiDelta = 0;

    for(int j = 1; j < s.size(); j++){
        semiDelta = -costMatrix[s[j]][s[j-1]]; //relies only on j

        for(int i = 1; i < s.size() - 3; i++){

            if (i == j || j == i+1 || j == i+2 || j == i+3) continue;

            delta = semiDelta + costMatrix[s[i]][s[j-1]] + costMatrix[s[i+2]][s[j]] + costMatrix[s[i-1]][s[i+3]]
            - costMatrix[s[i]][s[i-1]] - costMatrix[s[i+2]][s[i+3]];

            if(delta < *bestDelta){
                *bestDelta = delta;

                best_i = i;
                best_j = j;
            }
        }
    }

    if(*bestDelta < 0){
        if(best_i < best_j){
            s.insert(s.begin() + best_j, s[best_i+2]);
            s.insert(s.begin() + best_j, s[best_i+1]);
            s.insert(s.begin() + best_j, s[best_i]);
            s.erase(s.begin() + best_i);
            s.erase(s.begin() + best_i);
            s.erase(s.begin() + best_i);
        }
        else if(best_j < best_i){
            s.insert(s.begin() + best_j, s[best_i+2]);
            s.insert(s.begin() + best_j, s[best_i+2]);
            s.insert(s.begin() + best_j, s[best_i+2]);
            s.erase(s.begin() + best_i+3);
            s.erase(s.begin() + best_i+3);
            s.erase(s.begin() + best_i+3);
        }
    }else{
        *bestDelta = 0;
    }

    return s;
}

vector<int> RVND (vector<int> s, double *mainCost){
    vector<int> ngbhList = {N1, N2, N3};
    int ngbh_n;

    vector<int> neighbour_s = s;
    double neighbourCost = numeric_limits<double>::max();
    double delta = 0; //Receives delta from neighbouthood movements

    while(!ngbhList.empty())
    {
        ngbh_n = ngbhList[rand() % ngbhList.size()];

        switch(ngbh_n){
            case N1:
                neighbour_s = swap(s, &delta);
                break;
            case N2:
                neighbour_s = flip(s, &delta);
                break;
            case N3:
                neighbour_s = reinsertion(s, &delta);
                break;
            case N4:
                neighbour_s = oropt2(s, &delta);
                break;
            case N5:
                neighbour_s = oropt3(s, &delta);
                break;
        }
        neighbourCost = *mainCost + delta;

        if(neighbourCost < *mainCost){
            s = neighbour_s;
            *mainCost = neighbourCost;

            ngbhList = {N1, N2, N3};
        }else{
            ngbhList.erase(std::remove(ngbhList.begin(), ngbhList.end(), ngbh_n), ngbhList.end());
        }
    }

    return s;
}

vector<int> perturb (vector<int> s, double *cost){
    int alphaSize = ( rand() % (dimension/10) ) + 2;
    int alphaStart = (rand() % (dimension-2) ) + 1;
    int alphaEnd = alphaStart + alphaSize;

    //Gets a new alpha subsegment when there's no space for beta
    while(alphaEnd >= (dimension - 4))
    {
        alphaSize = ( rand() % (dimension/10) ) + 2;
        alphaStart = (rand() % (dimension-2) ) + 1;
        alphaEnd = alphaStart + alphaSize;
    }

    int betaSize = ( rand() % (dimension/10) ) + 2;
    int betaStart = (rand() % (dimension-2-alphaEnd) ) + alphaEnd + 1;
    int betaEnd = betaStart + betaSize;

    //Gets a new beta subsegment when beta covers the extremes
    while(betaEnd >= dimension)
    {
        betaSize = ( rand() % (dimension/10) ) + 2;
        betaStart = (rand() % (dimension-2-alphaEnd) ) + alphaEnd + 1;
        betaEnd = betaStart + betaSize;
    }

    //Gets cost difference after movement is made
    double delta = costMatrix[s[betaEnd-1]][s[alphaStart-1]] + costMatrix[s[betaStart]][s[alphaEnd]]
    + costMatrix[s[alphaEnd-1]][s[betaStart-1]] + costMatrix[s[alphaStart]][s[betaEnd]]
    - costMatrix[s[alphaStart]][s[alphaStart-1]] - costMatrix[s[alphaEnd]][s[alphaEnd-1]]
    - costMatrix[s[betaStart]][s[betaStart-1]] - costMatrix[s[betaEnd]][s[betaEnd-1]];

    vector<int> s_copy = s;


    //Invert subsegment elements
    for(int i = 0; i < alphaSize; i++){
        s.erase(s.begin() + alphaStart);
    }
    for(int i = 0; i < alphaSize; i++){
        s.insert(s.begin() + alphaStart, s_copy[alphaStart + i]);
    }

    for(int i = 0; i < betaSize; i++){
        s.erase(s.begin() + betaStart);
    }
    for(int i = 0; i < betaSize; i++){
        s.insert(s.begin() + betaStart, s_copy[betaStart + i]);
    }

    //Invert subsegment order
    s_copy = s;

    for(int i = betaEnd-1; i >= betaStart; i--){
        s.insert(s.begin() + alphaStart, s_copy[i]);
    }
    for(int i = 0; i < alphaSize; i++){
        s.erase(s.begin() + alphaStart + betaSize);
    }

    /*At this point, the solution will look like [... beta ... beta ...], and we need to
    insert alpha at the start of the second beta. However, betaStart no longer represents
    where the second beta starts, as this variable considers alpha at the first spot where
    currently there is a beta, and alpha and beta can have different sizes. As such, to find out
    where the new beta starts, we need to offset the betaStart variable by the difference
    between betaSize and alphaSize.*/
    int offset = betaSize - alphaSize;

    for(int i = alphaEnd-1; i >= alphaStart; i--){
        s.insert(s.begin() + betaStart + offset, s_copy[i]);
    }
    for(int i = 0; i < betaSize; i++){
        s.erase(s.begin() + betaStart + offset + alphaSize);
    }

    *cost += delta;

    return s;
}

int main(int argc, char** argv) {
    srand(time(NULL));

    readData(argc, argv, &dimension, &costMatrix);
    //std::cout << "\tCOST MATRIX: \n";
    //printCostMatrix();

    auto timerStart = chrono::system_clock::now();

    const int I_MAX = 50;
    const int I_ILS = (dimension >= 150) ? (dimension/2) : (dimension);

    vector<int> solutionAlpha, solutionBeta, solutionOmega;
    double costAlpha, costBeta, costOmega = numeric_limits<double>::max();

    for(int i = 0; i < I_MAX; i++)
    {
        double alpha = (rand() % 100)/100;
        solutionAlpha = construction(alpha);
        solutionBeta = solutionAlpha;

        costAlpha = solutionCost(solutionAlpha);
        costBeta = costAlpha;

        for(int iterILS = 0; iterILS < I_ILS; iterILS++){
            solutionAlpha = RVND(solutionAlpha, &costAlpha);

            if(costAlpha < costBeta){
                solutionBeta = solutionAlpha;
                costBeta = costAlpha;
                iterILS = 0;
            }

            solutionAlpha = perturb(solutionBeta, &costAlpha);
        }

        if(costBeta < costOmega){
            solutionOmega = solutionBeta;
            costOmega = costBeta;
        }
    }

    auto timerEnd = chrono::system_clock::now();
    chrono::duration<double> elapsedSeconds = timerEnd - timerStart;

    //PRINT COST AND SOLUTION
    std::cout << "SOLUTION:\n";
    for(auto k : solutionOmega){
        std::cout << k << ' ';
    }
    std::cout << "\n\n" << "COST: " << costOmega << "\n";
    std::cout << "TIME: " << elapsedSeconds.count() << "\n";

    return 0;
}
