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
    cout << "dimension: " << dimension << endl;
    for (size_t i = 1; i <= dimension; i++) {
        for (size_t j = 1; j <= dimension; j++) {
            cout << costMatrix[i][j] << " ";
        }
        cout << endl;
    }
}


bool compareCosts(const InsertionInfo &a, const InsertionInfo &b) {
    return (a.cost < b.cost);
}

vector<int> construction(double alpha) {

    vector<int> solution = {1, 1};

    vector<int> candidateList;

    int subtourSize = 3;
    int candidateListSize = dimension - 1; //Node 1 is already in the solution.

    //Fills candidates list according to dimension
    for(int i = 0, j = 2; i < candidateListSize; i++, j++){
        candidateList.insert(candidateList.begin() + i, j);
    }

    //Initial subtour
    for(int i = 0; i < subtourSize; i++){
        int j = rand() % candidateList.size();
        solution.insert(solution.begin() + 1 , candidateList[j]);
        candidateList.erase(candidateList.begin() + j);
    }


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


vector<int> swap (vector<int> s){
    int best_i = 0, best_j = 0;
    double bestDelta = 0; //0 = no change
    double delta;

    for(int j = 3; j < s.size() - 2; j++){
        for(int i = 1; i < j - 1; i++){

            delta = costMatrix[s[i]][s[j-1]] + costMatrix[s[i]][s[j+1]] + costMatrix[s[j]][s[i-1]] + costMatrix[s[j]][s[i+1]]
            - costMatrix[s[i]][s[i-1]] - costMatrix[s[i]][s[i+1]] - costMatrix[s[j]][s[j-1]] - costMatrix[s[j]][s[j+1]];

            if(delta < bestDelta){
                bestDelta = delta;

                best_i = i;
                best_j = j;
            }
        }
    }

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

    return s;
}

vector<int> flip (vector<int> s){
    int best_i = 0, best_j = 0;
    double bestDelta = 0; //0 = no change
    double delta;

    for(int j = 3; j < s.size() - 1; j++){
        for(int i = 1; i < j - 1; i++){

            delta = costMatrix[s[i]][s[j]] + costMatrix[s[j-1]][s[i-1]] - costMatrix[s[i-1]][s[i]] - costMatrix[s[j-1]][s[j]];

            if(delta < bestDelta){
                bestDelta = delta;

                best_i = i;
                best_j = j;                
            }
        }
    }

    for(int k = best_i, l = 0; k < best_j; k++, l++){
        s.insert(s.begin() + best_i, s[k+l]);
    }
    for(int k = best_i; k < best_j; k++){
        s.erase(s.begin() + best_j);
    }

    return s;
}

vector<int> reinsertion (vector<int> s){
    int best_i = 0, best_j = 0;
    double bestDelta = 0; //0 = no change
    double delta;

    for(int j = 1; j < s.size(); j++){
        for(int i = 1; i < s.size() - 1; i++){

            delta = costMatrix[s[i]][s[j]] + costMatrix[s[i]][s[j-1]] + costMatrix[s[i-1]][s[i+1]]
            - costMatrix[s[i]][s[i-1]] - costMatrix[s[i]][s[i+1]] - costMatrix[s[j]][s[j-1]];

            if(delta < bestDelta){
                bestDelta = delta;
                
                best_i = i;
                best_j = j;
            }
            
        }
    }

    s.insert(s.begin() + best_j, s[best_i]);
    if(best_j < best_i) best_i++;
    s.erase(s.begin() + best_i);

    return s;
}

vector<int> oropt2 (vector<int> s){
    int best_i = 0, best_j = 0;
    double bestDelta = 0; //0 = no change
    double delta;

    for(int j = 1; j < s.size(); j++){
        for(int i = 1; i < s.size() - 2; i++){

            if (j == i+1) continue;

            delta = costMatrix[s[i]][s[j-1]] + costMatrix[s[i+1]][s[j]] + costMatrix[s[i-1]][s[i+2]]
            - costMatrix[s[i]][s[i-1]] - costMatrix[s[i+1]][s[i+2]] - costMatrix[s[j]][s[j-1]];

            if(delta < bestDelta){
                bestDelta = delta;

                best_i = i;
                best_j = j;
            }
            
        }
    }

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

    return s;
}

vector<int> oropt3 (vector<int> s){
    int best_i = 0, best_j = 0;
    double bestDelta = 0; //0 = no change
    double delta;

    for(int j = 1; j < s.size(); j++){
        for(int i = 1; i < s.size() - 3; i++){

            if (j == i+1 || j == i+2) continue;
            
            delta = costMatrix[s[i]][s[j-1]] + costMatrix[s[i+2]][s[j]] + costMatrix[s[i-1]][s[i+3]] 
            - costMatrix[s[i]][s[i-1]] - costMatrix[s[i+2]][s[i+3]] - costMatrix[s[j]][s[j-1]];

            if(delta < bestDelta){
                bestDelta = delta;

                best_i = i;
                best_j = j;
            }
        }
    }

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

    return s;
}

vector<int> RVND (vector<int> s){
    vector<int> ngbhList = {N1, N2, N3, N4, N5};
    int ngbh_n;

    vector<int> neighbour_s;
    double s_cost = solutionCost(s);
    double neighbour_s_cost;

    while(!ngbhList.empty())
    {
        ngbh_n = ngbhList[rand() % ngbhList.size()];

        switch(ngbh_n){
            case N1:
                neighbour_s = swap(s);
                break;
            case N2:
                neighbour_s = flip(s);
                break;
            case N3:
                neighbour_s = reinsertion(s);
                break;
            case N4:
                neighbour_s = oropt2(s);
                break;
            case N5:
                neighbour_s = oropt3(s);
                break;
        }

        neighbour_s_cost = solutionCost(neighbour_s);
        if(neighbour_s_cost < s_cost){
            s = neighbour_s;
            s_cost = neighbour_s_cost;

            ngbhList = {N1, N2, N3, N4, N5};
        }else{
            ngbhList.erase(std::remove(ngbhList.begin(), ngbhList.end(), ngbh_n), ngbhList.end());
        }
    }

    return s;
}

vector<int> perturb (vector<int> s){
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
        betaStart = (rand() % (dimension-2-alphaEnd) ) + alphaEnd;
        betaEnd = betaStart + betaSize;
    }
       
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


    return s;
}

int main(int argc, char** argv) {
    auto timerStart = chrono::system_clock::now();

    srand(time(NULL));

    readData(argc, argv, &dimension, &costMatrix);
    //cout << "\tCOST MATRIX: \n";
    //printCostMatrix();

    int I_MAX = 50, I_ILS;
    if(dimension >= 150){
        I_ILS = dimension/2;
    }else{
        I_ILS = dimension;
    }

    vector<int> solutionAlpha, solutionBeta, solutionOmega;
    double costBeta, costOmega = numeric_limits<double>::max();
    
    for(int i = 0; i < I_MAX; i++)
    {
        double alpha = (rand() % 100)/100;
        solutionAlpha = construction(alpha);
        solutionBeta = solutionAlpha;

        for(int iterILS = 0; iterILS < I_ILS; iterILS++){
            solutionAlpha = RVND(solutionAlpha);

            if(solutionCost(solutionAlpha) < solutionCost(solutionBeta)){
                solutionBeta = solutionAlpha;
                iterILS = 0;
            }

            solutionAlpha = perturb(solutionBeta);
        }

        costBeta = solutionCost(solutionBeta);
        if(costBeta < costOmega){
            solutionOmega = solutionBeta;
            costOmega = costBeta;
        }
    }

    auto timerEnd = chrono::system_clock::now();
    chrono::duration<double> elapsedSeconds = timerEnd - timerStart;

    //PRINT COST AND SOLUTION
    cout << "\n\n\n" << "\tSOLUTION:\n";
    for(auto k : solutionOmega){
        cout << k << ' ';
    }
    cout << "\n\n\n\n" << "\tCOST: " << costOmega << "\n\n";
    cout << "\tTIME: " << elapsedSeconds.count() << "s\n\n\n";   

    return 0;
}