#include<iostream>
#include<fstream>
#include<ctime>

#ifdef BLAZE_DIAG
#include"blaze_diag.h"
#endif

#ifdef BLAZE_BLOCKDIAG
#include"blaze_blockdiag.h"
#endif

#ifdef BLAZE_SPARSEDIAG
#include"blaze_sparsediag.h"
#endif

int main() {
    const int N = 2; // Number of times to repeat the operation

    // measure:
    clock_t startTime = clock();
    for (int i = 0; i < N; ++i) // Perform the operation N times
    {
        operation();
    }
    clock_t endTime = clock();
    double totalTime = double(endTime - startTime) / CLOCKS_PER_SEC;
    double averageTime = totalTime / N;

    // print:
    std::cout << "Benchmark: " << benchmark_name << "\n";
    std::cout << "Total time taken: " << totalTime << " seconds\n";
    std::cout << "Average time per operation: " << averageTime << " seconds\n";

    // store:
    std::ofstream outFile("results.txt", std::ios::app);
    if (outFile.is_open()) 
    {
        outFile << benchmark_name << " | tot = " << totalTime << " s | av = " << averageTime << " s\n";
        outFile.close();
    } 
    else 
    {
        std::cerr << "Unable to open file for writing\n";
        return 1;
    }
}