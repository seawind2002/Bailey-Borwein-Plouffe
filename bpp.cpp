#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <omp.h>
#include <cstdio>

using namespace std;

// Calculate the term in the BBP formula
double CalculateTerm(unsigned long long k) {
    return 1.0 / pow(16, k) * (
        4.0 / (8.0 * k + 1.0) -
        2.0 / (8.0 * k + 4.0) -
        1.0 / (8.0 * k + 5.0) -
        1.0 / (8.0 * k + 6.0)
    );
}

// Calculate the value of Pi using Bailey–Borwein–Plouffe formula (sequential)
double CalculatePiBBP_Sequential(unsigned long long num_iterations) {
    double pi = 0.0;
    for (unsigned long long k = 0; k < num_iterations; ++k) {
        pi += CalculateTerm(k);
    }
    return pi;
}

// Calculate the value of Pi using Bailey–Borwein–Plouffe formula (parallel)
double CalculatePiBBP_Parallel(unsigned long long num_iterations) {
    double pi = 0.0;
#pragma omp parallel for reduction(+:pi)
    for (unsigned long long k = 0; k < num_iterations; ++k) {
        pi += CalculateTerm(k);
    }
    return pi;
}

int main() {
    int num_thread;
    unsigned long long int j = 0;
    unsigned long long int k;
    int num_runs;

    cout << "Enter the value of k: ";
    cin >> k;

    cout << "Enter the number of threads: ";
    cin >> num_thread;
    omp_set_num_threads(num_thread);

    cout << "Enter the number of runs: ";
    cin >> num_runs;

    // Arrays to store sums for each iteration
    double* sum_duration_sequential = new double[k - j + 1]();
    double* sum_duration_parallel = new double[k - j + 1]();
    double* sum_speedup = new double[k - j + 1]();
    double* sum_efficiency = new double[k - j + 1]();

    // Outer loop for multiple runs
    for (int run = 0; run < num_runs; ++run) {
        // Inner loop for different values of k
        for (unsigned long long i = j; i <= k; ++i) {
            unsigned long long num_iterations = pow(2, i);

            // Sequential
            auto start_time_sequential = chrono::steady_clock::now();
            double pi_sequential = CalculatePiBBP_Sequential(num_iterations);
            auto end_time_sequential = chrono::steady_clock::now();
            auto elapsed_time_sequential = chrono::duration_cast<chrono::duration<double>>(end_time_sequential - start_time_sequential);
            sum_duration_sequential[i - j] += elapsed_time_sequential.count();

            // Parallel
            auto start_time_parallel = chrono::steady_clock::now();
            double pi_parallel = CalculatePiBBP_Parallel(num_iterations);
            auto end_time_parallel = chrono::steady_clock::now();
            auto elapsed_time_parallel = chrono::duration_cast<chrono::duration<double>>(end_time_parallel - start_time_parallel);
            sum_duration_parallel[i - j] += elapsed_time_parallel.count();

            // Calculate and accumulate average speedup and efficiency
            sum_speedup[i - j] += elapsed_time_sequential.count() / elapsed_time_parallel.count();
            sum_efficiency[i - j] += (elapsed_time_sequential.count() / elapsed_time_parallel.count()) / num_thread;

            // Printing results to terminal for each run
            cout << "Run: " << run + 1 << " Iteration: " << i << endl;
            cout << setprecision(15) << "Approximation of Pi (Sequential): " << pi_sequential << endl;
            cout << "Time taken to calculate Pi (Sequential): " << elapsed_time_sequential.count() << " seconds" << endl;

            cout << setprecision(15) << "Approximation of Pi (Parallel): " << pi_parallel << endl;
            cout << "Time taken to calculate Pi (Parallel): " << elapsed_time_parallel.count() << " seconds" << endl;

            cout << "Performance Difference (Speedup): " << elapsed_time_sequential.count() / elapsed_time_parallel.count() << endl;
            cout << "Efficiency: " << (elapsed_time_sequential.count() / elapsed_time_parallel.count()) / num_thread << endl;

            cout << "------------------------------------" << endl;
        }
    }

    // Calculate averages
    for (unsigned long long i = j; i <= k; ++i) {
        sum_duration_sequential[i - j] /= num_runs; // Calculate average duration for sequential
        sum_duration_parallel[i - j] /= num_runs;   // Calculate average duration for parallel
        sum_speedup[i - j] /= num_runs;
        sum_efficiency[i - j] /= num_runs;
    }

    // Remove old CSV file if it exists
    remove("performance_results.csv");

    // Write results to CSV file
    ofstream csvFile("performance_results.csv");
    csvFile << "Iteration,AvgDurationSequential,AvgDurationParallel,AvgSpeedup,AvgEfficiency" << endl;

    for (unsigned long long i = j; i <= k; ++i) {
        csvFile << i << "," << setprecision(15) << sum_duration_sequential[i - j] << ","
                << sum_duration_parallel[i - j] << "," << sum_speedup[i - j] << "," << sum_efficiency[i - j] << endl;
    }

    delete[] sum_duration_sequential;
    delete[] sum_duration_parallel;
    delete[] sum_speedup;
    delete[] sum_efficiency;

    return 0;
}
