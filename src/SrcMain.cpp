#include "SrcMain.h"
#include <iostream>
#include <random>
#include "TSP.h"
#include <fstream>
#include <algorithm>


std::vector<Location>* readLocations(const std::string fileDic){
    std::ifstream inFile(fileDic);
    std::vector<Location> *locVector = new std::vector<Location>;
    std::string tempLine;
    unsigned int apoInd1;
    unsigned int apoInd2;
    std::string tempName;
    double tempLong;
    double tempLat;
    if(inFile.is_open()) {
        while(!inFile.eof()){
            //capture input line
            std::getline(inFile,tempLine);
            //detect different variables in line
            apoInd1 = tempLine.find(",");
            apoInd2 = tempLine.find(",",apoInd1+1);
            //save data from line
            tempName = tempLine.substr(0,apoInd1);
            tempLat = std::stod(tempLine.substr(apoInd1+1,apoInd2));
            tempLong = std::stod(tempLine.substr(apoInd2+1));
            //declare and intialize location struct
            Location vecElement;
            vecElement.mLongitude = tempLong;
            vecElement.mLatitude = tempLat;
            vecElement.mName = tempName;
            //add struct to final vector
            locVector -> push_back(vecElement);
        }
    }
    inFile.close();
    return locVector;
}


//for generate function in generate population

std::vector<std::vector<int>>* generatePopulation(int locSize, int popsize, std::mt19937 &randGener){
    
    std::vector<std::vector<int>>* population = new  std::vector<std::vector<int>>;
    for(int i = 0; i < popsize; i++) {
    std::vector<int> initV(locSize);
    
    int j = 0;
    std::iota(initV.begin(),initV.end(),j++);
    
    //std::uniform_int_distribution<int> randomLoc(0,locSize);
    std::shuffle(initV.begin() + 1, initV.end(), randGener);
    population -> push_back(initV);
    }
    return population;
}


std::vector<std::pair<int,double>>* calcFitness(std::vector<std::vector<int>>* population, int popsize,std::vector<Location> &locVec){
    //output vector to be filled
    std::vector<std::pair<int, double>>* fitnessVec= new std::vector<std::pair<int,double>>;
    
    //convert index vector to location vector
    std::vector<std::vector<Location>> popLocations;//population in the location translation
    std::for_each(population->begin(), population->end(), [&popLocations, locVec](std::vector<int> member){
        std::vector<Location> tempVec;
        std::for_each(member.begin(), member.end(),[&locVec,&tempVec](int ind){
            Location tempLoc = locVec[ind];
            tempVec.push_back(tempLoc);
        } );
        popLocations.push_back(tempVec);
    });
    
    
    //calculate differences
    const double Pi = 0.0174533;
    std::vector<std::vector<double>> calcResults;
    std::for_each(popLocations.begin(), popLocations.end(), [Pi,&calcResults](std::vector<Location> locMem){
        std::vector<Location> calcResultTemp;
        std::adjacent_difference(locMem.begin(), locMem.end(), std::back_inserter(calcResultTemp),  [Pi](Location first, Location second){
            
            double dlon = second.mLongitude - first.mLongitude;
            double dlat = second.mLatitude - first.mLatitude;
            double a = pow(sin(Pi *(dlat/2)),2) + cos(Pi*first.mLatitude) * cos(Pi*second.mLatitude) * pow(sin(Pi *(dlon/2)),2);
            double c = 2 * atan2(sqrt(a),sqrt(1-a));
            double distance = 3961 * c;
            Location tempLoc;
            tempLoc.distanceDiff = distance;
            
            return tempLoc;
        });
        //first element correction
        Location beg = locMem[0];
        Location end = locMem[locMem.size() - 1];
        double dlon = beg.mLongitude - end.mLongitude;
        double dlat = beg.mLatitude - end.mLatitude;
        double a = pow(sin(Pi *(dlat/2)),2) + (cos(Pi*end.mLatitude) * cos(Pi*beg.mLatitude) * pow(sin(Pi *(dlon/2)),2));
        double c = 2 * atan2(sqrt(a),sqrt(1-a));
        calcResultTemp[0].distanceDiff = 3961 * c;
        std::vector<double> doubleTemp;
        //transform location to double vector
        std::transform(calcResultTemp.begin(), calcResultTemp.end(),std::back_inserter(doubleTemp), [](Location loc){return loc.distanceDiff;});
        calcResults.push_back(doubleTemp);
    });
    
    int i = 0;
    std::for_each(calcResults.begin(), calcResults.end(), [fitnessVec,i](std::vector<double> popResult) mutable{
        std::pair<int,double> tempPair;
        tempPair.second = std::accumulate(popResult.begin(), popResult.end(),0.0, [](double a, double b){
            return a + b;
        });
        tempPair.first = i;
        i++;
        fitnessVec->push_back(tempPair);
    });
    
    return fitnessVec;
}


bool fitnessComparator(std::pair<int,double> A, std::pair<int,double> B){
    if(A.second < B.second)
        return true;
    else
        return false;
}
std::vector<std::pair<int,int>>* selectionPairs(std::vector<std::pair<int,double>> fitnessVec, std::mt19937 &randGen, int popsize){
    //sort fitness vector
    std::sort(fitnessVec.begin(),fitnessVec.end(),fitnessComparator);
    std::vector<double> probabilityVec(popsize);
    //generate initial probalility vector
    std::generate(probabilityVec.begin(), probabilityVec.end(), [&popsize](){return 1.0/popsize;});
    //give precedence to top two
    int topTwo1 = fitnessVec[0].first;
    int topTwo2 = fitnessVec[1].first;
    probabilityVec[topTwo1] *= 6.0;
    probabilityVec[topTwo2] *= 6.0;
    //give precendence to remainder of top half
    int tempInd;
    for(int i = 2; i < (popsize/2); i++){
        tempInd = fitnessVec[i].first;
        probabilityVec[tempInd] *= 3.0;
    }
    //renormailize prob vector
    double probSum = std::accumulate(probabilityVec.begin(), probabilityVec.end(),0.0, [](double a, double b){
        return a + b;
    });
    //finalProbV is normalized vector
    std::vector<double> finalProbV;
    std::transform(probabilityVec.begin(), probabilityVec.end(), std::back_inserter(finalProbV), [&probSum](double A){return (A/probSum);});
    //picking pairs
    std::uniform_real_distribution<double> realDis(0,1);
    std::vector<std::pair<int, int>> *pairsVector = new std::vector<std::pair<int,int>>;
    double randProb;
    double pSum = 0;
    int firstParent;
    int secondParent;
    for(int i = 0; i < popsize; i++){
        //first random parent
        randProb = realDis(randGen);
        for(int j = 0; j < popsize; j++){
            pSum += finalProbV[j];
            if(pSum >= randProb){
                firstParent = j;
                pSum = 0;
                break;
            }
        }
        //second random parent
        randProb = realDis(randGen);
        for(int k = 0; k < popsize; k++){
             pSum += finalProbV[k];
            if(pSum >= randProb){
                secondParent = k;
                pSum = 0;
                break;
            }
        }
        std::pair<int,int> parentsPair;
        parentsPair.first = firstParent;
        parentsPair.second = secondParent;
        pairsVector->push_back(parentsPair);
    }
    return pairsVector;
}

std::vector<std::vector<int>>* crossover(std::vector<std::vector<int>> &population, int initialSize, std::mt19937 &randGen, std::vector<std::pair<int,int>> &parents, int popsize, double mutationchance){
    
    
    std::uniform_real_distribution<double> doubleDis(0,1);
    std::uniform_int_distribution<int> mutateDis(1,initialSize-1);

    std::vector<std::vector<int>>* crossoverV = new std::vector<std::vector<int>>;
    std::vector<int> firstParentV;
    std::vector<int> secondParentV;
    for(int i = 0; i < popsize; i++){
    /**std::for_each(parents.begin(), parents.end(), [&crossoverV,&firstParentV,&secondParentV,&intDis,&intDis2,&population,&doubleDis,&randGen,&mutationchance,&initialSize](std::pair<int,int> parentPair){**/
        
        std::uniform_int_distribution<int> intDis(1,initialSize - 2);
        int crossoverIndex = intDis(randGen);
        std::uniform_int_distribution<int> intDis2(0,1);
        int parentLead = intDis2(randGen);
        parentLead = 0;
        std::vector<int> newChild;
        if(parentLead == 0){
            firstParentV = population[parents[i].second];
            secondParentV = population[parents[i].first];
        } else {
            firstParentV = population[parents[i].first];
            secondParentV = population[parents[i].second];
        }
        std::copy_n(firstParentV.begin(), crossoverIndex, std::back_inserter(newChild));
        std::copy_if(secondParentV.begin(), secondParentV.end(), std::back_inserter(newChild), [&newChild](int A){
            auto iter = std::find(newChild.begin(), newChild.end(), A);
            if(iter == newChild.end())
                return true;
            else
                return false;
        });
        //Mutation
        double mutationProb = doubleDis(randGen);
        if(mutationProb <= mutationchance){
            int firstInd = mutateDis(randGen);
            int secondInd = mutateDis(randGen);
            std::swap(newChild[firstInd], newChild[secondInd]);
        }
        crossoverV->push_back(newChild);
        //});
    
    
   }
    return crossoverV;
}

void logOutput(int popsize, std::mt19937 &randGen,int generations ,int initialLocationSize,int mutationchance, std::vector<Location> &initialLoc){
    //print initial log for testing
    std::ofstream oFile("log.txt");
    if(oFile.is_open()){
        //print population
        oFile << "INITIAL POPULATION:" << std::endl;
        std::vector<std::vector<int>>* pop = generatePopulation(initialLocationSize, popsize, randGen);

        std::for_each(pop->begin(), pop->end(), [&oFile](std::vector<int> v){
            std::for_each(v.begin(), v.end(), [&oFile](int elem){
                oFile << elem << ",";
                
            });
            oFile << std::endl;
        });
        //print fitness
        oFile << "FITNESS:" << std::endl;
        std::vector<std::pair<int,double>>* fitness = calcFitness(pop, popsize, initialLoc);

        std::for_each(fitness->begin(), fitness->end(), [&oFile](std::pair<int,double> fitMem){
            oFile << fitMem.first << ":" << fitMem.second << std::endl;
        });
        oFile << "SELECTED PAIRS:" << std::endl;
        std::vector<std::pair<int,int>>* selectionP = selectionPairs(*fitness, randGen, popsize);
        std::for_each(selectionP->begin(), selectionP->end(), [&oFile](std::pair<int,int> selMem){
            oFile << "(" << selMem.first << "," << selMem.second << ")" << std::endl;
        });
        //print Generations
        std::vector<std::vector<int>>* testGeneration = pop;
        std::vector<std::pair<int,double>> finalFitness;
        for(int i = 1; i <= generations; i++){
            oFile << "GENERATION: " << i << std::endl;
            testGeneration = crossover(*testGeneration, initialLocationSize, randGen, *selectionP, popsize, (double)(mutationchance/10.0));
            std::for_each(testGeneration->begin(), testGeneration->end(), [&oFile](std::vector<int> v){
            std::for_each(v.begin(), v.end(), [&oFile](int elem){
                oFile << elem << ",";
                });
            oFile << std::endl;
            });
            //print fitness
            oFile << "FITNESS:" << std::endl;
            std::vector<std::pair<int,double>>* fitness = calcFitness(testGeneration, popsize, initialLoc);
            
            std::for_each(fitness->begin(), fitness->end(), [&oFile](std::pair<int,double> fitMem){
                oFile << fitMem.first << ":" << fitMem.second << std::endl;
            });
            std::sort(fitness->begin(),fitness->end(),fitnessComparator);
            finalFitness = *fitness;
            oFile << "SELECTED PAIRS:" << std::endl;
            std::vector<std::pair<int,int>>* selectionP = selectionPairs(*fitness, randGen, popsize);
            std::for_each(selectionP->begin(), selectionP->end(), [&oFile](std::pair<int,int> selMem){
                oFile << "(" << selMem.first << "," << selMem.second << ")" << std::endl;
            });
        }
        oFile << "SOLUTION:" << std::endl;
        int finalInd = finalFitness[0].first;
        std::vector<std::vector<int>> finalGener = *testGeneration;
        std::for_each(finalGener[finalInd].begin(), finalGener[finalInd].end(), [&initialLoc,&oFile](int locInd){
            oFile << initialLoc[locInd].mName << std::endl;
        });
        oFile << initialLoc[finalGener[finalInd][0]].mName << std::endl;
        oFile << "DISTANCE:" << finalFitness[0].second << "miles";
    }
    oFile.close();
}

void ProcessCommandArgs(int argc, const char* argv[])
{
	// TODO
    std::string inputFile = argv[1];
    int popsize = std::stoi(argv[2]);
    int generations = std::stoi(argv[3]);
    int mutationchance = std::stoi(argv[4]);
    int seed = std::stoi(argv[5]);
    std::mt19937 randGen(seed);
    std::vector<Location> *initialLocations = readLocations(inputFile);
    
    logOutput(popsize,randGen,generations,initialLocations->size(),mutationchance,*initialLocations);
    
}
