#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#define inv_sqrt_2xPI 0.39894228040143270286


typedef struct TestsData_ {
    float spotPrice;          // spot price
    float strike;     // strike price
    float rate;          // risk-free interest rate
    float divq;       // dividend rate
    float volatility;          // volatility
    float otime;          // time to maturity or option expiration in years
    //     (1yr = 1.0, 6mos = 0.5, 3mos = 0.25, ..., etc)
    const char *otype;  // Option type.  "P"=PUT, "C"=CALL
    float divs;       // dividend vals (not used in this test)
    float DGrefval;   // DerivaGem Reference Value
} TestsData;

TestsData data_init[] = {
#include "TestData.txt"
};

float CNDF (float InputX)
{
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
    xLocal   = xLocal_1 * xNPrimeofX;
    xLocal   = 1.0 - xLocal;

    OutputX  = xLocal;

    if (sign) {
        OutputX = 1.0 - OutputX;
    }

    return OutputX;
}


float blackScholes(float sptprice, float strike, float rate, float volatility,
                   float otime, int otype, float timet)
{
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

    FutureValueX = strike * (exp(-(rate)*(otime)));
    if (otype == 0) {
        OptionPrice = (sptprice * NofXd1) - (FutureValueX * NofXd2);
    } else {
        NegNofXd1 = (1.0 - NofXd1);
        NegNofXd2 = (1.0 - NofXd2);
        OptionPrice = (FutureValueX * NegNofXd2) - (sptprice * NegNofXd1);
    }

    return OptionPrice;
}

int main(int args, char *argv[]) {
    printf("hello world\n");
    int tests_size = sizeof(data_init)/sizeof(TestsData);
    int num_threads = 16;

    // Time keeping
    struct timespec b_start_t, a_start_t, b_end_t, a_end_t;


    // Warm up
    clock_gettime(CLOCK_MONOTONIC, &b_start_t);
    clock_gettime(CLOCK_MONOTONIC, &b_end_t);

    // Perform SISD calculation on test data
    clock_gettime(CLOCK_MONOTONIC, &b_start_t);
    for(int i=0; i<tests_size; i++) {
        float oprice = blackScholes(data_init[i].spotPrice, data_init[i].strike, data_init[i].rate, data_init[i].volatility, data_init[i].otime, data_init[i].otype, 0);
        // printf("price of %d = %.5f\n", i+1, oprice);
    }
    clock_gettime(CLOCK_MONOTONIC, &b_end_t);

    unsigned long long b_time = ((b_end_t.tv_sec - b_start_t.tv_sec) * 1000000000) + (b_end_t.tv_nsec - b_start_t.tv_nsec);


    // Display information
    printf("Time for running:\n");
    printf("*** SISD  Implementation: %lld ns\n", b_time);
    // printf("*** MIMD Implementation: %lld ns\n", a_time);
    printf("\n");
    printf("\n");
    // printf("      -> Speedup = %.2fx\n", ((b_time * 1.0f) / a_time));
    return 0;
}
