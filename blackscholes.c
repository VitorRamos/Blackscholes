// Copyright (c) 2007 Intel Corp.

// Black-Scholes
// Analytical method for calculating European Options
//
//
// Reference Source: Options, Futures, and Other Derivatives, 3rd Edition, Prentice
// Hall, John C. Hull,


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <string.h>

//Precision to use for calculations
#define fptype float

#define NUM_RUNS 100


typedef struct ReadData_
{
    fptype s;          // spot price
    fptype strike;     // strike price
    fptype r;          // risk-free interest rate
    fptype divq;       // dividend rate (not used)
    fptype v;          // volatility
    fptype t;          // time to maturity or option expiration in years
    //     (1yr = 1.0, 6mos = 0.5, 3mos = 0.25, ..., etc)
    char OptionType;   // Option type.  "P"=PUT, "C"=CALL
    fptype divs;       // dividend vals (not used in this test)
    fptype DGrefval;   // DerivaGem Reference Value
} ReadData;


typedef struct OptionData_
{
    fptype s;          // spot price
    fptype strike;     // strike price
    fptype r;          // risk-free interest rate
//    fptype divq;       // dividend rate (not used)
    fptype v;          // volatility
    fptype t;          // time to maturity or option expiration in years
    //     (1yr = 1.0, 6mos = 0.5, 3mos = 0.25, ..., etc)
    int OptionType;   // Option type.  "P"=PUT, "C"=CALL
//    fptype divs;       // dividend vals (not used in this test)
//    fptype DGrefval;   // DerivaGem Reference Value
} OptionData;


OptionData *local_data;
fptype *local_prices;
int local_numOptions;

/*
int    * otype;
fptype * sptprice;
fptype * strike;
fptype * rate;
fptype * volatility;
fptype * otime;
*/
int numError = 0;
int nThreads;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Cumulative Normal Distribution Function
// See Hull, Section 11.8, P.243-244
#define inv_sqrt_2xPI 0.39894228040143270286

fptype CNDF ( fptype InputX )
{
    int sign;

    fptype OutputX;
    fptype xInput;
    fptype xNPrimeofX;
    fptype expValues;
    fptype xK2;
    fptype xK2_2, xK2_3;
    fptype xK2_4, xK2_5;
    fptype xLocal, xLocal_1;
    fptype xLocal_2, xLocal_3;

    // Check for negative value of InputX
    if (InputX < 0.0)
    {
        InputX = -InputX;
        sign = 1;
    }
    else
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

    if (sign)
        OutputX = 1.0 - OutputX;

    return OutputX;
}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
fptype BlkSchlsEqEuroNoDiv( fptype sptprice,
                            fptype strike, fptype rate, fptype volatility,
                            fptype time, int otype, float timet )
{
    fptype OptionPrice;

    // local private working variables for the calculation
    fptype xStockPrice;
    fptype xStrikePrice;
    fptype xRiskFreeRate;
    fptype xVolatility;
    fptype xTime;
    fptype xSqrtTime;

    fptype logValues;
    fptype xLogTerm;
    fptype xD1;
    fptype xD2;
    fptype xPowerTerm;
    fptype xDen;
    fptype d1;
    fptype d2;
    fptype FutureValueX;
    fptype NofXd1;
    fptype NofXd2;
    fptype NegNofXd1;
    fptype NegNofXd2;

    xStockPrice = sptprice;
    xStrikePrice = strike;
    xRiskFreeRate = rate;
    xVolatility = volatility;

    xTime = time;
    xSqrtTime = sqrt(xTime);

    logValues = log( sptprice / strike );

    xLogTerm = logValues;


    xPowerTerm = xVolatility * xVolatility;
    xPowerTerm = xPowerTerm * 0.5;

    xD1 = xRiskFreeRate + xPowerTerm;
    xD1 = xD1 * xTime;
    xD1 = xD1 + xLogTerm;

    xDen = xVolatility * xSqrtTime;
    xD1 = xD1 / xDen;
    xD2 = xD1 -  xDen;

    d1 = xD1;
    d2 = xD2;

    NofXd1 = CNDF( d1 );
    NofXd2 = CNDF( d2 );

    FutureValueX = strike * ( exp( -(rate)*(time) ) );
    if (otype == 0)
    {
        OptionPrice = (sptprice * NofXd1) - (FutureValueX * NofXd2);
    }
    else
    {
        NegNofXd1 = (1.0 - NofXd1);
        NegNofXd2 = (1.0 - NofXd2);
        OptionPrice = (FutureValueX * NegNofXd2) - (sptprice * NegNofXd1);
    }

    return OptionPrice;
}


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////


int bs_thread(void *tid_ptr)
{
    int i, j;
    fptype price;
    fptype priceDelta;
    int tid = *(int *)tid_ptr;

    for (j=0; j<NUM_RUNS; j++)
    {
        #pragma omp parallel for private(i, price, priceDelta)
        for (i=0; i<local_numOptions; i++)
        {
            /* Calling main function to calculate option value based on
             * Black & Scholes's equation.
             */
            price = BlkSchlsEqEuroNoDiv( local_data[i].s, local_data[i].strike,
                                         local_data[i].r, local_data[i].v, local_data[i].t,
                                         local_data[i].OptionType, 0);
            local_prices[i] = price;
            #ifdef ERR_CHK
                priceDelta = data[i].DGrefval - price;
                if( fabs(priceDelta) >= 1e-4 )
                {
                    printf("Error on %d. Computed=%.5f, Ref=%.5f, Delta=%.5f\n",
                           i, price, data[i].DGrefval, priceDelta);
                    numError ++;
                }
            #endif
        }
    }

    return 0;
}


int main (int argc, char **argv)
{
    int rank, comm_sz;
    int loopnum, rv;

    int numOptions;
    char *inputFile, *outputFile;
    MPI_File in_file, out_file;
    MPI_Offset offset;

    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    double ttotal= MPI_Wtime();

    if(argc < 4)
    {
        MPI_Finalize();
    }

    nThreads= atoi(argv[1]);
    inputFile= argv[2];
    outputFile= argv[3];

    if(rank == 0)
    {
        printf("PARSEC Benchmark Suite\n");
    }
    MPI_File_open(MPI_COMM_WORLD, inputFile, MPI_MODE_RDONLY, MPI_INFO_NULL, &in_file);
    MPI_File_read(in_file, &numOptions, 1, MPI_INT, MPI_STATUS_IGNORE);

    int r= numOptions%comm_sz;
    local_numOptions= (numOptions/comm_sz);
    local_numOptions+= (r-rank)>0?1:0;

    // alloc spaces for the option data
    local_data = (OptionData*)malloc(local_numOptions*sizeof(OptionData));
    local_prices = (fptype*)malloc(local_numOptions*sizeof(fptype));

    ReadData A, *local_read;
    MPI_Datatype mpi_read_data, mpi_read_data_;
    int blocklen[9]= {1, 1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype types[9]= {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,
                            MPI_FLOAT, MPI_FLOAT, MPI_CHAR, MPI_FLOAT, MPI_FLOAT};
    MPI_Aint base_addr, addr, disp[9];
    MPI_Get_address(&A.s, &base_addr);
    disp[0]= base_addr-base_addr;
    MPI_Get_address(&A.strike, &addr);
    disp[1]= addr-base_addr;
    MPI_Get_address(&A.r, &addr);
    disp[2]= addr-base_addr;
    MPI_Get_address(&A.divq, &addr);
    disp[3]= addr-base_addr;
    MPI_Get_address(&A.v, &addr);
    disp[4]= addr-base_addr;
    MPI_Get_address(&A.t, &addr);
    disp[5]= addr-base_addr;
    MPI_Get_address(&A.OptionType, &addr);
    disp[6]= addr-base_addr;
    MPI_Get_address(&A.divs, &addr);
    disp[7]= addr-base_addr;
    MPI_Get_address(&A.DGrefval, &addr);
    disp[8]= addr-base_addr;
    MPI_Type_create_struct(9, blocklen, disp, types, &mpi_read_data_);
    MPI_Type_create_resized(mpi_read_data_, 0, 1*sizeof(ReadData), &mpi_read_data);
    MPI_Type_commit(&mpi_read_data);

    int data_len= 8*sizeof(float)+sizeof(char);
    local_read= (ReadData*)malloc(local_numOptions*sizeof(ReadData));
    offset= rank*local_numOptions*(data_len)+sizeof(int);

    MPI_File_seek(in_file,offset,MPI_SEEK_SET);
    MPI_File_read(in_file,local_read,local_numOptions,mpi_read_data,MPI_STATUS_IGNORE);

    fptype divq, divs, DGrefval;
    char OptionType;
    for (loopnum= 0; loopnum<local_numOptions; ++loopnum)
    {
        local_data[loopnum].s= local_read[loopnum].s;
        local_data[loopnum].strike= local_read[loopnum].strike;
        local_data[loopnum].r= local_read[loopnum].r;
        local_data[loopnum].v= local_read[loopnum].v;
        local_data[loopnum].t= local_read[loopnum].t;
        local_data[loopnum].OptionType= (local_read[loopnum].OptionType == 'P') ? 1 : 0;
    }
    MPI_File_close(&in_file);
    free(local_read);

    if(rank == 0)
    {
        printf("Input file %s\n", inputFile);
        printf("Num of Options: %d\n", numOptions);
        printf("Num of Local Options: %d\n", local_numOptions);
        printf("Num of Nos: %d\n", comm_sz);
        printf("Num of Threads: %d\n", nThreads);
        printf("Num of Runs: %d\n", NUM_RUNS);
        printf("Size of data: %ld\n", numOptions * (sizeof(OptionData) + sizeof(int)));
    }

    double t1= MPI_Wtime(), time;
    {
        int tid=0;
        omp_set_num_threads(nThreads);
        bs_thread(&tid);
    }
    t1= MPI_Wtime()-t1;
    MPI_Reduce(&t1,&time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    if(rank == 0)
    {
        printf("Tempo dos calculos : %lf\n", time);
    }

    char* string_out= malloc(13*(local_numOptions+1)*sizeof(char));
    offset=rank*local_numOptions*sizeof(char)*(13);
    MPI_File_open(MPI_COMM_WORLD,outputFile,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL, &out_file);
    for(loopnum= 0; loopnum<local_numOptions; ++loopnum)
    {
        if(local_prices[loopnum]>=0)
            sprintf(string_out+(loopnum*13),"+%.5e\n",local_prices[loopnum]);
        else
            sprintf(string_out+(loopnum*13),"%.5e\n",local_prices[loopnum]);
    }
    MPI_File_write_ordered(out_file,string_out,13*local_numOptions,MPI_CHAR,MPI_STATUS_IGNORE);
    MPI_File_close(&out_file);

#ifdef ERR_CHK
    printf("Num Errors: %d\n", numError);
#endif
    //free(string_out);
    free(local_data);
    free(local_prices);

    ttotal= MPI_Wtime()-ttotal;
    MPI_Reduce(&ttotal,&time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    if(rank == 0)
        printf("Tempo Total : %lf\n", time);

    MPI_Finalize();
    return 0;
}
