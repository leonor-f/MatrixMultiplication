#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <cstdlib>
#include <papi.h>
#include <omp.h>
#include <fstream>

using namespace std;

#define SYSTEMTIME clock_t

double OnMult(int m_ar, int m_br)
{
    SYSTEMTIME Time1, Time2;
    char st[100];
    int i, j, k;
    double *pha, *phb, *phc, temp, time;

    pha = (double *)malloc((m_ar * m_ar) * sizeof(double));
    phb = (double *)malloc((m_ar * m_ar) * sizeof(double));
    phc = (double *)malloc((m_ar * m_ar) * sizeof(double));

    for (i = 0; i < m_ar; i++)
        for (j = 0; j < m_ar; j++)
            pha[i * m_ar + j] = (double)1.0;

    for (i = 0; i < m_br; i++)
        for (j = 0; j < m_br; j++)
            phb[i * m_br + j] = (double)(i + 1);

    Time1 = clock();

    for (i = 0; i < m_ar; i++)
    {
        for (j = 0; j < m_br; j++)
        {
            temp = 0;
            for (k = 0; k < m_ar; k++)
                temp += pha[i * m_ar + k] * phb[k * m_br + j];
            phc[i * m_ar + j] = temp;
        }
    }

    Time2 = clock();
    time = (double)(Time2 - Time1) / CLOCKS_PER_SEC;
    cout << endl;
    sprintf(st, "Time: %3.3f seconds\n", time);
    cout << st;

    // Display first 10 elements of the result matrix to verify correctness
    cout << "Result matrix: " << endl;
    for (i = 0; i < 1; i++)
    {
        for (j = 0; j < min(10, m_br); j++)
            cout << phc[j] << " ";
    }
    cout << endl;

    free(pha);
    free(phb);
    free(phc);

    return time;
}

double OnMultElemLine(int m_ar, int m_br, bool flag)
{
    SYSTEMTIME Time1, Time2;
    char st[100];
    int i, j, k;
    double *pha, *phb, *phc, time;

    pha = (double *)malloc((m_ar * m_ar) * sizeof(double));
    phb = (double *)malloc((m_ar * m_ar) * sizeof(double));
    phc = (double *)malloc((m_ar * m_ar) * sizeof(double));

    // Initialize matrices
    for (i = 0; i < m_ar; i++)
        for (j = 0; j < m_ar; j++)
            pha[i * m_ar + j] = (double)1.0;

    for (i = 0; i < m_br; i++)
        for (j = 0; j < m_br; j++)
            phb[i * m_br + j] = (double)(i + 1);

    for (i = 0; i < m_ar; i++)
        for (j = 0; j < m_ar; j++)
            phc[i * m_ar + j] = (double)0.0;

    // Start timing
    Time1 = clock();

    // Element-by-line multiplication
    for (i = 0; i < m_ar; i++)
    {
        for (j = 0; j < m_br; j++)
        {
            for (k = 0; k < m_ar; k++)
                phc[i * m_ar + k] += pha[i * m_ar + j] * phb[j * m_br + k];
        }
    }

    // End timing
    Time2 = clock();
    time = (double)(Time2 - Time1) / CLOCKS_PER_SEC;
    if (flag)
    {
        cout << endl;
        sprintf(st, "Time: %3.3f seconds\n", time);
        cout << st;

        // Display first 10 elements of the result matrix
        cout << "Result matrix: " << endl;
        for (i = 0; i < 1; i++)
            for (j = 0; j < min(10, m_br); j++)
                cout << phc[j] << " ";
        cout << endl;
    }

    free(pha);
    free(phb);
    free(phc);

    return time;
}

double OnMultBlock(int m_ar, int m_br, int bkSize)
{
    SYSTEMTIME Time1, Time2;
    char st[100];
    int i, j, k, l, m, n;
    double *pha, *phb, *phc, time;

    pha = (double *)malloc((m_ar * m_ar) * sizeof(double));
    phb = (double *)malloc((m_ar * m_ar) * sizeof(double));
    phc = (double *)malloc((m_ar * m_ar) * sizeof(double));

    for (i = 0; i < m_ar; i++)
        for (j = 0; j < m_ar; j++)
            pha[i * m_ar + j] = (double)1.0;

    for (i = 0; i < m_br; i++)
        for (j = 0; j < m_br; j++)
            phb[i * m_br + j] = (double)(i + 1);

    for (i = 0; i < m_ar; i++)
        for (j = 0; j < m_ar; j++)
            phc[i * m_ar + j] = (double)0.0;

    Time1 = clock();

    for (l = 0; l < m_ar; l += bkSize)
    {
        for (m = 0; m < m_br; m += bkSize)
        {
            for (n = 0; n < m_ar; n += bkSize)
            {
                // Perform element-line multiplication for each block
                for (i = l; i < min(l + bkSize, m_ar); i++)
                {
                    for (k = n; k < min(n + bkSize, m_ar); k++)
                    {
                        for (j = m; j < min(m + bkSize, m_br); j++)
                        {
                            phc[i * m_ar + k] += pha[i * m_ar + j] * phb[j * m_br + k];
                        }
                    }
                }
            }
        }
    }

    Time2 = clock();
    time = (double)(Time2 - Time1) / CLOCKS_PER_SEC;
    cout << endl;
    sprintf(st, "Time: %3.3f seconds\n", time);
    cout << st;

    // Display first 10 elements of the result matrix
    cout << "Result matrix: " << endl;
    for (i = 0; i < 1; i++)
        for (j = 0; j < min(10, m_br); j++)
            cout << phc[j] << " ";
    cout << endl;

    free(pha);
    free(phb);
    free(phc);

    return time;
}

double OnMultLineParallel1(int m_ar, int m_br)
{
    char st[100];
    int i, j, k;
    double *pha, *phb, *phc, time, start_time, end_time;

    pha = (double *)malloc((m_ar * m_ar) * sizeof(double));
    phb = (double *)malloc((m_ar * m_ar) * sizeof(double));
    phc = (double *)malloc((m_ar * m_ar) * sizeof(double));

    // Initialize matrices
    for (i = 0; i < m_ar; i++)
        for (j = 0; j < m_ar; j++)
            pha[i * m_ar + j] = (double)1.0;

    for (i = 0; i < m_br; i++)
        for (j = 0; j < m_br; j++)
            phb[i * m_br + j] = (double)(i + 1);

    for (i = 0; i < m_ar; i++)
        for (j = 0; j < m_ar; j++)
            phc[i * m_ar + j] = (double)0.0;

    // Start timing
    start_time = omp_get_wtime();

// Element-line multiplication with OpenMP parallelization (outer for loop version)
#pragma omp parallel for private(i, j, k) shared(pha, phb, phc)
    for (i = 0; i < m_ar; i++)
        for (j = 0; j < m_br; j++)
            for (k = 0; k < m_ar; k++)
                phc[i * m_ar + k] += pha[i * m_ar + j] * phb[j * m_br + k];

    // End timing
    end_time = omp_get_wtime();
    time = end_time - start_time;
    cout << endl;
    sprintf(st, "Time: %3.3f seconds\n", time);
    cout << st;

    // Display first 10 elements of the result matrix
    cout << "Result matrix: " << endl;
    for (i = 0; i < 1; i++)
        for (j = 0; j < min(10, m_br); j++)
            cout << phc[j] << " ";
    cout << endl;

    free(pha);
    free(phb);
    free(phc);

    return time;
}

double OnMultLineParallel2(int m_ar, int m_br)
{
    char st[100];
    int i, j, k;
    double *pha, *phb, *phc, time, start_time, end_time;

    pha = (double *)malloc((m_ar * m_ar) * sizeof(double));
    phb = (double *)malloc((m_ar * m_ar) * sizeof(double));
    phc = (double *)malloc((m_ar * m_ar) * sizeof(double));

    // Initialize matrices
    for (i = 0; i < m_ar; i++)
        for (j = 0; j < m_ar; j++)
            pha[i * m_ar + j] = (double)1.0;

    for (i = 0; i < m_br; i++)
        for (j = 0; j < m_br; j++)
            phb[i * m_br + j] = (double)(i + 1);

    for (i = 0; i < m_ar; i++)
        for (j = 0; j < m_ar; j++)
            phc[i * m_ar + j] = (double)0.0;

    // Start timing
    start_time = omp_get_wtime();

// Element-line multiplication with OpenMP parallelization (inner for loop version)
#pragma omp parallel private(i, j, k) shared(pha, phb, phc)
    {
        for (i = 0; i < m_ar; i++)
        {
            for (k = 0; k < m_ar; k++)
            {
#pragma omp for
                for (j = 0; j < m_br; j++)
                {
                    phc[i * m_ar + j] += pha[i * m_ar + k] * phb[k * m_br + j];
                }
            }
        }
    }

    // End timing
    end_time = omp_get_wtime();
    time = end_time - start_time;
    cout << endl;
    sprintf(st, "Time: %3.3f seconds\n", time);
    cout << st;

    // Display first 10 elements of the result matrix
    cout << "Result matrix: " << endl;
    for (i = 0; i < 1; i++)
        for (j = 0; j < min(10, m_br); j++)
            cout << phc[j] << " ";
    cout << endl;

    free(pha);
    free(phb);
    free(phc);

    return time;
}

void measure_performance(double (*func)(int, int), int m_ar, int m_br)
{
    double t_parallel, mflops, speedup, t_sequential, efficiency;
    int num_threads;
    char st[100];

    t_parallel = func(m_ar, m_br);
    t_sequential = OnMultElemLine(m_ar, m_br, false);

    mflops = (2.0 * m_ar * m_br * m_ar) / (t_parallel * 1e6);
    num_threads = omp_get_max_threads();
    speedup = t_sequential / t_parallel;
    efficiency = speedup / num_threads;

    sprintf(st, "Mflops: %3.3f\n", mflops);
    cout << st;
    sprintf(st, "Speedup: %3.3f\n", speedup);
    cout << st;
    sprintf(st, "Efficiency: %3.3f\n", efficiency);
    cout << st;
}

void handle_error(int retval)
{
    printf("PAPI error %d: %s\n", retval, PAPI_strerror(retval));
    exit(1);
}

void init_papi()
{
    int retval = PAPI_library_init(PAPI_VER_CURRENT);
    if (retval != PAPI_VER_CURRENT && retval < 0)
    {
        printf("PAPI library version mismatch!\n");
        exit(1);
    }
    if (retval < 0)
        handle_error(retval);

    std::cout << "PAPI Version Number: MAJOR: " << PAPI_VERSION_MAJOR(retval)
              << " MINOR: " << PAPI_VERSION_MINOR(retval)
              << " REVISION: " << PAPI_VERSION_REVISION(retval) << "\n";
}

void write_result_mult(int eventSet, long long values[2], double (*func)(int, int), ofstream &file, const string &name, int min, int max, int step, int blockSize)
{
    double time, total_memory_accesses;
    int ret;

    for (int matrixSize = min; matrixSize <= max; matrixSize += step)
    {
        ret = PAPI_start(eventSet);
        if (ret != PAPI_OK)
        {
            cout << "ERROR: Start PAPI" << endl;
            handle_error(ret);
        }
        time = func(matrixSize, matrixSize);
        ret = PAPI_stop(eventSet, values);
        total_memory_accesses = 3 * (long long)matrixSize * matrixSize * matrixSize;
        file << name << ","
             << to_string(matrixSize) << ","
             << (blockSize == 1 ? "-," : "")
             << time << ","
             << values[0] << ","
             << values[1] << ","
             << total_memory_accesses - values[0] << ","
             << total_memory_accesses - values[1] << endl;
        ret = PAPI_reset(eventSet);
    }
}

void write_result_line(int eventSet, long long values[2], double (*func)(int, int, bool), ofstream &file, const string &name, int min, int max, int step, int blockSize)
{
    double time, total_memory_accesses;
    int ret;

    for (int matrixSize = min; matrixSize <= max; matrixSize += step)
    {
        ret = PAPI_start(eventSet);
        if (ret != PAPI_OK)
        {
            cout << "ERROR: Start PAPI" << endl;
            handle_error(ret);
        }
        time = func(matrixSize, matrixSize, true);
        ret = PAPI_stop(eventSet, values);
        total_memory_accesses = 3 * (long long)matrixSize * matrixSize * matrixSize;
        file << name << ","
             << to_string(matrixSize) << ","
             << (blockSize == 1 ? "-," : "")
             << time << ","
             << values[0] << ","
             << values[1] << ","
             << total_memory_accesses - values[0] << ","
             << total_memory_accesses - values[1] << endl;
        ret = PAPI_reset(eventSet);
    }
}

void write_result_block(int eventSet, long long values[2], ofstream &file, int blockSize, const string &name, int min, int max, int step)
{
    double time, total_memory_accesses;
    int ret;

    for (int matrixSize = min; matrixSize <= max; matrixSize += step)
    {
        ret = PAPI_start(eventSet);
        if (ret != PAPI_OK)
        {
            cout << "ERROR: Start PAPI" << endl;
            handle_error(ret);
        }
        time = OnMultBlock(matrixSize, matrixSize, blockSize);
        ret = PAPI_stop(eventSet, values);
        total_memory_accesses = 3 * (long long)matrixSize * matrixSize * matrixSize;
        file << name << ","
             << to_string(matrixSize) << ","
             << to_string(blockSize) << ","
             << time << ","
             << values[0] << ","
             << values[1] << ","
             << total_memory_accesses - values[0] << ","
             << total_memory_accesses - values[1] << endl;
        ret = PAPI_reset(eventSet);
    }
}

void write_all_results_mult(int eventSet, long long values[2])
{
    cout << "Creating file result_multiplication.csv with results..." << endl;

    ofstream file("result_multiplication.csv");
    file << "Name,Matrix Size,Block Size,Time,L1_DCM,L2_DCM,L1_Hits,L2_Hits\n";

    int block_sizes[] = {128, 256, 512};

    // Standard multiplication
    write_result_mult(eventSet, values, OnMult, file, "Standard", 600, 3000, 400, 1);

    // Element-line multiplication
    write_result_line(eventSet, values, OnMultElemLine, file, "Element-Line", 600, 3000, 400, 1);
    write_result_line(eventSet, values, OnMultElemLine, file, "Element-Line", 4096, 10240, 2048, 1);

    // Block multiplication
    for (int block : block_sizes)
    {
        write_result_block(eventSet, values, file, block, "Block", 600, 3000, 400);
        write_result_block(eventSet, values, file, block, "Block", 4096, 10240, 2048);
    }

    cout << endl
         << "Finished writing results into result_multiplication.csv" << endl;
}

void write_result_parallel(int eventSet, long long values[2], double (*func)(int, int), ofstream &file, const string &name, int min, int max, int step)
{
    double t_parallel, t_sequential, mflops, speedup, efficiency, total_memory_accesses;
    int ret, num_threads;

    for (int matrixSize = min; matrixSize <= max; matrixSize += step)
    {
        ret = PAPI_start(eventSet);
        if (ret != PAPI_OK)
        {
            cout << "ERROR: Start PAPI" << endl;
            handle_error(ret);
        }
        t_parallel = func(matrixSize, matrixSize);
        ret = PAPI_stop(eventSet, values);
        t_sequential = OnMultElemLine(matrixSize, matrixSize, false);

        mflops = (2.0 * matrixSize * matrixSize * matrixSize) / (t_parallel * 1e6);
        num_threads = omp_get_max_threads();
        speedup = t_sequential / t_parallel;
        efficiency = speedup / num_threads;
        total_memory_accesses = 3 * (long long)matrixSize * matrixSize * matrixSize;
        
        file << name << ","
             << to_string(matrixSize) << ","
             << t_parallel << ","
             << values[0] << ","
             << values[1] << ","
             << total_memory_accesses - values[0] << ","
             << total_memory_accesses - values[1]<< ","
             << mflops << ","
             << speedup << ","
             << efficiency << endl;
        ret = PAPI_reset(eventSet);
    }
}

void write_all_results_parallel(int eventSet, long long values[2])
{
    cout << endl
         << "Creating file result_parallelization.csv with results..." << endl;

    ofstream file("result_parallelization.csv");
    file << "Name,Matrix Size,Time,L1_DCM,L2_DCM,L1_Hits,L2_Hits,MFlops,Speedup,Efficiency\n";

    // Parallel 1
    write_result_parallel(eventSet, values, OnMultLineParallel1, file, "Parallel 1", 600, 3000, 400);
    write_result_parallel(eventSet, values, OnMultLineParallel1, file, "Parallel 1", 4096, 10240, 2048);

    // Parallel 2
    write_result_parallel(eventSet, values, OnMultLineParallel2, file, "Parallel 2", 600, 3000, 400);
    write_result_parallel(eventSet, values, OnMultLineParallel2, file, "Parallel 2", 4096, 10240, 2048);

    cout << endl
         << "Finished writing results into result_parallelization.csv" << endl;
}

void callPythonScript(const string &scriptName, int matrixSize, double &time)
{
    string command = "python3 " + scriptName + " " + to_string(matrixSize);
    FILE *pipe = popen(command.c_str(), "r");
    if (!pipe)
    {
        cerr << "Failed to run Python script" << endl;
        exit(1);
    }

    char buffer[128];
    string result = "";
    while (fgets(buffer, sizeof(buffer), pipe) != NULL)
        result += buffer;
    pclose(pipe);

    if (result.empty())
    {
        cerr << "Python script did not return a valid result" << endl;
        exit(1);
    }

    // Parse the output
    time = stod(result);
}

void write_all_results_python()
{
    cout << endl
         << "Creating file result_python.csv with results..." << endl;

    ofstream file("result_python.csv");
    file << "Name,Matrix Size,Time\n";

    double time;

    // Standard multiplication
    for (int matrixSize = 600; matrixSize <= 3000; matrixSize += 400)
    {
        callPythonScript("matrix_mult.py", matrixSize, time);
        file << "Standard,"
             << matrixSize << ","
             << time << endl;
    }

    // Element-line multiplication
    for (int matrixSize = 600; matrixSize <= 3000; matrixSize += 400)
    {
        callPythonScript("line_mult.py", matrixSize, time);
        file << "Element-Line,"
             << matrixSize << ","
             << time << endl;
    }

    cout << endl
         << "Finished writing results into result_python.csv" << endl;
}

void generate_results(const string &argv, int eventSet, long long values[2])
{
    if (argv == "1")
        write_all_results_python();
    else if (argv == "2")
        write_all_results_mult(eventSet, values);
    else if (argv == "3")
        write_all_results_parallel(eventSet, values);
    else
        cout << "Invalid argument" << endl;
}

int main(int argc, char *argv[])
{
    char c;
    int lin, col, blockSize;
    int op, sub_op;
    long long values[2], total_memory_accesses, l1_cache_hits, l2_cache_hits;
    int EventSet = PAPI_NULL;
    int ret;

    ret = PAPI_library_init(PAPI_VER_CURRENT);
    if (ret != PAPI_VER_CURRENT)
        std::cout << "FAIL" << endl;

    ret = PAPI_create_eventset(&EventSet);
    if (ret != PAPI_OK)
        cout << "ERROR: create eventset" << endl;

    ret = PAPI_add_event(EventSet, PAPI_L1_DCM);
    if (ret != PAPI_OK)
        cout << "ERROR: PAPI_L1_DCM" << endl;

    ret = PAPI_add_event(EventSet, PAPI_L2_DCM);
    if (ret != PAPI_OK)
        cout << "ERROR: PAPI_L2_DCM" << endl;

    if (argc > 1)
        generate_results(argv[1], EventSet, values);

    op = -1;
    sub_op = -1;
    do
    {
        op = -1;
        sub_op = -1;
        cout << endl
             << "1. Multiplication" << endl;
        cout << "2. Element-Line Multiplication" << endl;
        cout << "3. Block Multiplication" << endl;
        cout << "0. Exit" << endl;
        cout << "Selection?: ";
        cin >> op;

        if (op == 0)
            break;

        if (op == 1)
        {
            cout << endl;
            cout << "1. Python" << endl;
            cout << "2. C++" << endl;
            cout << "Selection?: ";
            cin >> sub_op;

            if (sub_op == 1)
            {
                printf("Dimensions: lins=cols ? ");
                cin >> lin;
                col = lin;
                string command = "python3 matrix_mult.py " + to_string(lin);
                int sys_ret = system(command.c_str());
                (void)sys_ret;
                continue;
            }
            sub_op = -1;
        }

        if (op == 2)
        {
            cout << endl;
            cout << "1. Python" << endl;
            cout << "2. C++" << endl;
            cout << "3. Parallel Implementation (Solution 1)" << endl;
            cout << "4. Parallel Implementation (Solution 2)" << endl;
            cout << "Selection?: ";
            cin >> sub_op;
            op = -1;

            if (sub_op == 1)
            {
                printf("Dimensions: lins=cols ? ");
                cin >> lin;
                col = lin;
                string command = "python3 line_mult.py " + to_string(lin);
                int sys_ret = system(command.c_str());
                (void)sys_ret;
                continue;
            }
        }

        printf("Dimensions: lins=cols ? ");
        cin >> lin;
        col = lin;

        // Start counting
        ret = PAPI_start(EventSet);
        if (ret != PAPI_OK)
            cout << "ERROR: Start PAPI" << endl;

        if (op == 1)
        {
            (void)OnMult(lin, col);
        }
        else if (op == 3)
        {
            cout << "Block Size? ";
            cin >> blockSize;
            (void)OnMultBlock(lin, col, blockSize);
        }
        else if (sub_op == 2)
        {
            (void)OnMultElemLine(lin, col, true);
        }
        else if (sub_op == 3)
        {
            measure_performance(OnMultLineParallel1, lin, col);
        }
        else if (sub_op == 4)
        {
            measure_performance(OnMultLineParallel2, lin, col);
        }

        ret = PAPI_stop(EventSet, values);
        if (ret != PAPI_OK)
            cout << "ERROR: Stop PAPI" << endl;
        printf("L1 DCM: %lld \n", values[0]);
        printf("L2 DCM: %lld \n", values[1]);

        total_memory_accesses = 3 * (long long)lin * lin * lin;
        l1_cache_hits = total_memory_accesses - values[0];
        l2_cache_hits = total_memory_accesses - values[1];
        printf("L1 Hits: %lld \n", l1_cache_hits);
        printf("L2 Hits: %lld \n", l2_cache_hits);

        ret = PAPI_reset(EventSet);
        if (ret != PAPI_OK)
            std::cout << "FAIL reset" << endl;

    } while (op != 0);

    ret = PAPI_remove_event(EventSet, PAPI_L1_DCM);
    if (ret != PAPI_OK)
        std::cout << "FAIL remove event" << endl;

    ret = PAPI_remove_event(EventSet, PAPI_L2_DCM);
    if (ret != PAPI_OK)
        std::cout << "FAIL remove event" << endl;

    ret = PAPI_destroy_eventset(&EventSet);
    if (ret != PAPI_OK)
        std::cout << "FAIL destroy" << endl;

    return 0;
}