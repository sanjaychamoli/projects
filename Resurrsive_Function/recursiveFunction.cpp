#include <iostream>
#include <vector>
#include <fstream>
#include <time.h>
#include <sys/time.h>
#include <map>

using namespace std;

typedef uint64_t size64_t;
uint vectorSize = 200;

/**
* @brief Function to measure the wall clock
* @return The current time
*/
double get_wall_time();

/**
* @brief Function to execute the algorithm
* @param number The input for which the value is to be computer using the algorithm
* @param vectorStore A vector with already computer values
* @param mapStore A map to store computed values
* @return The computed value from the alogrithm
*/
size64_t algorithmFunction(size64_t number, vector<int>* vectorStore, std::map<size64_t, size64_t>*  mapStore);

/**
* @brief Stores a few values into a vector.
* @param number The input for which the value is to be computer using the algorithm
* @return The computed value from the alogrithm
*/
size64_t vectorCreator(uint number);

int main(){
    size64_t number = 123456789012345678;

    // Start of timer
    double wall0 = get_wall_time();
    std::vector<int> vectorStore(vectorSize);
    for (uint i = vectorStore.size()-1; i > 0; --i)
    {
        vectorStore[i] = vectorCreator(i);
    }
    std::map<size64_t, size64_t>  mapStore;
    size64_t output = algorithmFunction(number, &vectorStore, &mapStore);
    double wall1 = get_wall_time();
    // End of timer
    cout <<"\nThe output for "<< number << " is: "<< output;
    cout << "\nWall Time Total= " << wall1 - wall0 << "\n";

    return 0;
}


double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL))
    {
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}


size64_t algorithmFunction(size64_t number, vector<int>* vectorStore, std::map<size64_t, size64_t>* mapStore){

    // Checking if the key is already stored in the map
    auto search = mapStore->find(number);
    if ( search != mapStore->end() )
    {
        return search->second;
    }

    // If not stored in map or vector then solve
    if (number > (vectorSize - 1))
    {
        // If number is even
        if (!(number % 2))
        {
            // If number is a power of 2
            if((number & (number - 1)) == 0)
            {
                mapStore-> operator[](number) = 1;
                return 1;
            }
            // Find solve of a number and add to map
            else
            {
                size64_t newNumber = number/2;
                size64_t temp = algorithmFunction(newNumber, vectorStore, mapStore);
                mapStore-> operator[](number) = temp;
                return temp;
            }
        }
        // Odd number. Solve and vectorStore in map
        else
        {
            size64_t newNumber = number/2;

            size64_t sum1 = algorithmFunction(newNumber, vectorStore, mapStore);
            size64_t sum2 = algorithmFunction(newNumber -1, vectorStore, mapStore);

            sum1 += sum2;
            mapStore->operator[](number) = sum1;
            return  sum1;
        }
    }
    return vectorStore->operator [](number);
}


size64_t vectorCreator(uint number){
    if(number > 2)
    {
        //even
        if (!(number % 2))
        {
            if((number & (number - 1)) == 0)
                return 1;
            else
                return vectorCreator(number/2);
        }
        else
        {
            uint newNumber = number/2;
            uint sum1 = vectorCreator(newNumber);
            uint sum2 = vectorCreator(newNumber - 1);
            return  size64_t(sum1 + sum2);
        }
    }
    else
        return 1;
}
