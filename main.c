#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#define inv_sqrt_2xPI 0.39894228040143270286

#define tests_size  (int )sizeof(data_init) / sizeof(TestsData)
#define num_threads 16

typedef struct TestsData_ {
    float spotPrice; // spot price
    float strike; // strike price
    float rate; // risk-free interest rate
    float divq; // dividend rate
    float volatility; // volatility
    float otime; // time to maturity or option expiration in years
    //     (1yr = 1.0, 6mos = 0.5, 3mos = 0.25, ..., etc)
    const char *otype; // Option type.  "P"=PUT, "C"=CALL
    float divs; // dividend vals (not used in this test)
    float DGrefval; // DerivaGem Reference Value
} TestsData;

TestsData data_init[] = {
#include "testData.txt"
#include "testData.txt"
#include "testData.txt"
#include "testData.txt"
#include "testData.txt"
#include "testData.txt"
#include "testData.txt"
#include "testData.txt"
#include "testData.txt"
#include "testData.txt"
#include "testData.txt"
};

typedef struct singleThread {
    int id;
    int data_size;
    TestsData data[tests_size];
    float result[tests_size];
} singleThread;

float CNDF(float InputX) {
    int sign;

    float OutputX;
    float xInput;
    float xNPrimeofX;
    float expValues;
    float xK2;
    float xK2_2, xK2_3;
    float xK2_4, xK2_5;
    float xLocal, xLocal_1;
    float xLocal_2, xLocal_3;

    // Check for negative value of InputX
    if (InputX < 0.0) {
        InputX = -InputX;
        sign = 1;
    } else
        sign = 0;

    xInput = InputX;

    // Compute NPrimeX term common to both four & six decimal accuracy calcs
    expValues = exp(-0.5f * InputX * InputX);
    xNPrimeofX = expValues;
    xNPrimeofX = xNPrimeofX * inv_sqrt_2xPI;

    xK2 = 0.2316419 * xInput;
    xK2 = 1.0 + xK2;
    xK2 = 1.0 / xK2;
    xK2_2 = xK2 * xK2;
    xK2_3 = xK2_2 * xK2;
    xK2_4 = xK2_3 * xK2;
    xK2_5 = xK2_4 * xK2;

    xLocal_1 = xK2 * 0.319381530;
    xLocal_2 = xK2_2 * (-0.356563782);
    xLocal_3 = xK2_3 * 1.781477937;
    xLocal_2 = xLocal_2 + xLocal_3;
    xLocal_3 = xK2_4 * (-1.821255978);
    xLocal_2 = xLocal_2 + xLocal_3;
    xLocal_3 = xK2_5 * 1.330274429;
    xLocal_2 = xLocal_2 + xLocal_3;

    xLocal_1 = xLocal_2 + xLocal_1;
    xLocal = xLocal_1 * xNPrimeofX;
    xLocal = 1.0 - xLocal;

    OutputX = xLocal;

    if (sign) {
        OutputX = 1.0 - OutputX;
    }

    return OutputX;
}


float blackScholes(float sptprice, float strike, float rate, float volatility,
                   float otime, int otype, float timet) {
    float OptionPrice;

    // local private working variables for the calculation
    float xStockPrice;
    float xStrikePrice;
    float xRiskFreeRate;
    float xVolatility;
    float xTime;
    float xSqrtTime;

    float logValues;
    float xLogTerm;
    float xD1;
    float xD2;
    float xPowerTerm;
    float xDen;
    float d1;
    float d2;
    float FutureValueX;
    float NofXd1;
    float NofXd2;
    float NegNofXd1;
    float NegNofXd2;

    xStockPrice = sptprice;
    xStrikePrice = strike;
    xRiskFreeRate = rate;
    xVolatility = volatility;

    xTime = otime;
    xSqrtTime = sqrt(xTime);

    logValues = log(sptprice / strike);

    xLogTerm = logValues;


    xPowerTerm = xVolatility * xVolatility;
    xPowerTerm = xPowerTerm * 0.5;

    xD1 = xRiskFreeRate + xPowerTerm;
    xD1 = xD1 * xTime;
    xD1 = xD1 + xLogTerm;

    xDen = xVolatility * xSqrtTime;
    xD1 = xD1 / xDen;
    xD2 = xD1 - xDen;

    d1 = xD1;
    d2 = xD2;

    NofXd1 = CNDF(d1);
    NofXd2 = CNDF(d2);

    FutureValueX = strike * (exp(-(rate) * (otime)));
    if (otype == 0) {
        OptionPrice = (sptprice * NofXd1) - (FutureValueX * NofXd2);
    } else {
        NegNofXd1 = (1.0 - NofXd1);
        NegNofXd2 = (1.0 - NofXd2);
        OptionPrice = (FutureValueX * NegNofXd2) - (sptprice * NegNofXd1);
    }

    return OptionPrice;
}

int oTypeCharToInt(const char *c) {
    if (strcmp(c, "P") == 0)
        return 1;
    return 0;
}

void *blackScholesCaller(void *args) {
    singleThread *thread_args = (singleThread *) args;
    int id = thread_args->id;
    int data_size = thread_args->data_size;
    TestsData *data = thread_args->data;
    float *result = thread_args->result;
    for (int i = 0; i < data_size; i++) {
        result[i] = blackScholes(data[i].spotPrice, data[i].strike, data[i].rate,
                                 data[i].volatility, data[i].otime, oTypeCharToInt(data[i].otype),
                                 0);
    }
}

int main(int args, char *argv[]) {
    pthread_t thread[num_threads];
    singleThread thread_args[num_threads];
    float sisd_result[tests_size];
    float mimd_result[tests_size];
    bool is_valid = true;
    int normalizer = 100;

    // Time keeping
    struct timespec b_start_t, a_start_t, b_end_t, a_end_t;
    unsigned long long b_time = 0, a_time = 0;


    //loop normlizer time to get the sum time for all the tests for normalization
    for (int k = 0; k < normalizer; k++) {
        // Warm up
        clock_gettime(CLOCK_MONOTONIC, &b_start_t);
        clock_gettime(CLOCK_MONOTONIC, &b_end_t);

        // Perform SISD calculation on test data
        clock_gettime(CLOCK_MONOTONIC, &b_start_t);
        for (int i = 0; i < tests_size; i++) {
            sisd_result[i] = blackScholes(data_init[i].spotPrice, data_init[i].strike, data_init[i].rate,
                                          data_init[i].volatility, data_init[i].otime,
                                          oTypeCharToInt(data_init[i].otype),
                                          0);
        }
        clock_gettime(CLOCK_MONOTONIC, &b_end_t);


        clock_gettime(CLOCK_MONOTONIC, &a_start_t);

        // Perform MIMD calculation on test data that is divided into num_threads threads
        for (int i = 0; i < num_threads; i++) {
            // Initialize the argument struture
            thread_args[i].id = i;
            thread_args[i].data_size = tests_size / num_threads;
            // printf("data size %d\n", thread_args[i].data_size);
            for (int j = 0; j < thread_args[i].data_size; j++) {
                // printf("for thread %d, data %d\n", i, j);
                int originIndex = i * thread_args[i].data_size;
                thread_args[i].data[j].spotPrice = data_init[originIndex + j].spotPrice;
                thread_args[i].data[j].strike = data_init[originIndex + j].strike;
                thread_args[i].data[j].rate = data_init[originIndex + j].rate;
                thread_args[i].data[j].volatility = data_init[originIndex + j].volatility;
                thread_args[i].data[j].otime = data_init[originIndex + j].otime;
                thread_args[i].data[j].otype = data_init[originIndex + j].otype;
            }

            pthread_create(&thread[i], NULL, blackScholesCaller, (void *) &thread_args[i]);
        }

        //calculate the rest of the data
        for (int i = 0; i < tests_size % num_threads; i++) {
            int originIndex = tests_size - tests_size % num_threads;
            mimd_result[originIndex + i] = blackScholes(data_init[originIndex + i].spotPrice,
                                                        data_init[originIndex + i].strike,
                                                        data_init[originIndex + i].rate,
                                                        data_init[originIndex + i].volatility,
                                                        data_init[originIndex + i].otime,
                                                        oTypeCharToInt(data_init[originIndex + i].otype),
                                                        0);
        }

        for (int i = 0; i < num_threads; i++) {
            pthread_join(thread[i], NULL);
            for (int j = 0; j < thread_args[i].data_size; j++) {
                mimd_result[i * thread_args[i].data_size + j] = thread_args[i].result[j];
            }
        }
        clock_gettime(CLOCK_MONOTONIC, &a_end_t);

        //Check if the results are maching or not
        for (int i = 0; i < tests_size; i++) {
            if (sisd_result[i] != mimd_result[i]) {
                printf("sisd_result[%d] = %.5f\n", i, sisd_result[i]);
                printf("mimd_result[%d] = %.5f\n", i, mimd_result[i]);
                is_valid = false;
                break;
            }
        }

        //sum the time
        b_time = b_time + ((b_end_t.tv_sec - b_start_t.tv_sec) * 1000000000) + (
                     b_end_t.tv_nsec - b_start_t.tv_nsec);
        a_time = a_time + ((a_end_t.tv_sec - a_start_t.tv_sec) * 1000000000) + (
                     a_end_t.tv_nsec - a_start_t.tv_nsec);
    }

    //normalize the time
    b_time = b_time / normalizer;
    a_time = a_time / normalizer;

    // Display information
    printf("Result validation:\n");
    if (is_valid)
        printf("      -> Results are matching\n");
    else
        printf("      -> Results are not matching\n");


    printf("*** SISD  Implementation: %lld ns\n", b_time);
    printf("*** MIMD Implementation: %lld ns\n", a_time);
    printf("\n");
    printf("\n");
    printf("      -> Speedup = %.2fx\n", ((b_time * 1.0f) / a_time));
    return 0;
}
