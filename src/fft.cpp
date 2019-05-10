#include <windows.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MatPlot.h"
using namespace MatPlot;

#define round(x) (int)(((x)>=0 ? (x)+0.5 : (x)-0.5))
struct cmplx{
    double real;
    double imag;
};
void r2dif_nr(int N, struct cmplx *wTb, struct cmplx *b);
void r2dit_nr(int N, struct cmplx *wTb, struct cmplx *b);
void r4dif_nr(int N, struct cmplx *wTb, struct cmplx *b);
int bit_rev(int a, int radix, int len);
struct cmplx cmplx_add(struct cmplx a, struct cmplx b);
struct cmplx cmplx_minus(struct cmplx a, struct cmplx b);
struct cmplx cmplx_mul(struct cmplx a, struct cmplx b);
struct cmplx wTb4[4] = {/*{{{*/
                            { 1.0000,  0.0000},
                            { 0.0000, -1.0000},
                            {-1.0000,  0.0000},
                            { 0.0000,  1.0000},
};/*}}}*/
struct cmplx wTb8[8] = {/*{{{*/
                            { 1.0000,  0.0000},
                            { 0.7071,  0.7071},
                            { 0.0000,  1.0000},
                            {-0.7071,  0.7071},
                            {-1.0000,  0.0000},
                            {-0.7071, -0.7071},
                            { 0.0000, -1.0000},
                            { 0.7071, -0.7071}
};/*}}}*/
struct cmplx wTb16[16] = {/*{{{*/
                            { 1.0000,  0.0000},
                            { 0.9239,  0.3827},
                            { 0.7071,  0.7071},
                            { 0.3827,  0.9239},
                            { 0.0000,  1.0000},
                            {-0.3827,  0.9239},
                            {-0.7071,  0.7071},
                            {-0.9239,  0.3827},
                            {-1.0000,  0.0000},
                            {-0.9239, -0.3827},
                            {-0.7071, -0.7071},
                            {-0.3827, -0.9239},
                            {-0.0000, -1.0000},
                            { 0.3827, -0.9239},
                            { 0.7071, -0.7071},
                            { 0.9239, -0.3827}
};/*}}}*/
struct cmplx wTb32[32] = {/*{{{*/
                            { 1.0000,  0.0000},
                            { 0.9808,  0.1951},
                            { 0.9239,  0.3827},
                            { 0.8315,  0.5556},
                            { 0.7071,  0.7071},
                            { 0.5556,  0.8315},
                            { 0.3827,  0.9239},
                            { 0.1951,  0.9808},
                            { 0.0000,  1.0000},
                            {-0.1951,  0.9808},
                            {-0.3827,  0.9239},
                            {-0.5556,  0.8315},
                            {-0.7071,  0.7071},
                            {-0.8315,  0.5556},
                            {-0.9239,  0.3827},
                            {-0.9808,  0.1951},
                            {-1.0000,  0.0000},
                            {-0.9808, -0.1951},
                            {-0.9239, -0.3827},
                            {-0.8315, -0.5556},
                            {-0.7071, -0.7071},
                            {-0.5556, -0.8315},
                            {-0.3827, -0.9239},
                            {-0.1951, -0.9808},
                            {-0.0000, -1.0000},
                            { 0.1951, -0.9808},
                            { 0.3827, -0.9239},
                            { 0.5556, -0.8315},
                            { 0.7071, -0.7071},
                            { 0.8315, -0.5556},
                            { 0.9239, -0.3827},
                            { 0.9808, -0.1951}
};/*}}}*/
struct cmplx wTb64[64] = {/*{{{*/
                            { 1.0000,  0.0000},
                            { 0.9952, -0.0980},
                            { 0.9808, -0.1951},
                            { 0.9569, -0.2903},
                            { 0.9239, -0.3827},
                            { 0.8819, -0.4714},
                            { 0.8315, -0.5556},
                            { 0.7730, -0.6344},
                            { 0.7071, -0.7071},
                            { 0.6344, -0.7730},
                            { 0.5556, -0.8315},
                            { 0.4714, -0.8819},
                            { 0.3827, -0.9239},
                            { 0.2903, -0.9569},
                            { 0.1951, -0.9808},
                            { 0.0980, -0.9952},
                            { 0.0000, -1.0000},
                            {-0.0980, -0.9952},
                            {-0.1951, -0.9808},
                            {-0.2903, -0.9569},
                            {-0.3827, -0.9239},
                            {-0.4714, -0.8819},
                            {-0.5556, -0.8315},
                            {-0.6344, -0.7730},
                            {-0.7071, -0.7071},
                            {-0.7730, -0.6344},
                            {-0.8315, -0.5556},
                            {-0.8819, -0.4714},
                            {-0.9239, -0.3827},
                            {-0.9569, -0.2903},
                            {-0.9808, -0.1951},
                            {-0.9952, -0.0980},
                            {-1.0000,  0.0000},
                            {-0.9952,  0.0980},
                            {-0.9808,  0.1951},
                            {-0.9569,  0.2903},
                            {-0.9239,  0.3827},
                            {-0.8819,  0.4714},
                            {-0.8315,  0.5556},
                            {-0.7730,  0.6344},
                            {-0.7071,  0.7071},
                            {-0.6344,  0.7730},
                            {-0.5556,  0.8315},
                            {-0.4714,  0.8819},
                            {-0.3827,  0.9239},
                            {-0.2903,  0.9569},
                            {-0.1951,  0.9808},
                            {-0.0980,  0.9952},
                            { 0.0000,  1.0000},
                            { 0.0980,  0.9952},
                            { 0.1951,  0.9808},
                            { 0.2903,  0.9569},
                            { 0.3827,  0.9239},
                            { 0.4714,  0.8819},
                            { 0.5556,  0.8315},
                            { 0.6344,  0.7730},
                            { 0.7071,  0.7071},
                            { 0.7730,  0.6344},
                            { 0.8315,  0.5556},
                            { 0.8819,  0.4714},
                            { 0.9239,  0.3827},
                            { 0.9569,  0.2903},
                            { 0.9808,  0.1951},
                            { 0.9952,  0.0980}
};/*}}}*/
#if 1
struct cmplx data_a[64] = {/*{{{*/
                          { 0.0000,  0.0000},
                          { 1.5801,  0.0000},
                          {-1.4266,  0.0000},
                          { 2.0566,  0.0000},
                          {-1.4695,  0.0000},
                          {-0.5000,  0.0000},
                          { 0.8817,  0.0000},
                          {-1.7476,  0.0000},
                          { 2.3776,  0.0000},
                          {-0.7711,  0.0000},
                          { 0.0000,  0.0000},
                          { 0.7711,  0.0000},
                          {-2.3776,  0.0000},
                          { 1.7476,  0.0000},
                          {-0.8817,  0.0000},
                          { 0.5000,  0.0000},
                          { 1.4695,  0.0000},
                          {-2.0566,  0.0000},
                          { 1.4266,  0.0000},
                          {-1.5801,  0.0000},
                          { 0.0000,  0.0000},
                          { 1.5801,  0.0000},
                          {-1.4266,  0.0000},
                          { 2.0566,  0.0000},
                          {-1.4695,  0.0000},
                          {-0.5000,  0.0000},
                          { 0.8817,  0.0000},
                          {-1.7476,  0.0000},
                          { 2.3776,  0.0000},
                          {-0.7711,  0.0000},
                          { 0.0000,  0.0000},
                          { 0.7711,  0.0000},
                          {-2.3776,  0.0000},
                          { 1.7476,  0.0000},
                          {-0.8817,  0.0000},
                          { 0.5000,  0.0000},
                          { 1.4695,  0.0000},
                          {-2.0566,  0.0000},
                          { 1.4266,  0.0000},
                          {-1.5801,  0.0000},
                          { 0.0000,  0.0000},
                          { 1.5801,  0.0000},
                          {-1.4266,  0.0000},
                          { 2.0566,  0.0000},
                          {-1.4695,  0.0000},
                          {-0.5000,  0.0000},
                          { 0.8817,  0.0000},
                          {-1.7476,  0.0000},
                          { 2.3776,  0.0000},
                          {-0.7711,  0.0000},
                          { 0.0000,  0.0000},
                          { 0.7711,  0.0000},
                          {-2.3776,  0.0000},
                          { 1.7476,  0.0000},
                          {-0.8817,  0.0000},
                          { 0.5000,  0.0000},
                          { 1.4695,  0.0000},
                          {-2.0566,  0.0000},
                          { 1.4266,  0.0000},
                          {-1.5801,  0.0000},
                          { 0.0000,  0.0000},
                          { 1.5801,  0.0000},
                          {-1.4266,  0.0000},
                          { 2.0566,  0.0000}
};/*}}}*/
struct cmplx data_b[64] = {/*{{{*/
                          { 0.0000,  0.0000},
                          { 1.5801,  0.0000},
                          {-1.4266,  0.0000},
                          { 2.0566,  0.0000},
                          {-1.4695,  0.0000},
                          {-0.5000,  0.0000},
                          { 0.8817,  0.0000},
                          {-1.7476,  0.0000},
                          { 2.3776,  0.0000},
                          {-0.7711,  0.0000},
                          { 0.0000,  0.0000},
                          { 0.7711,  0.0000},
                          {-2.3776,  0.0000},
                          { 1.7476,  0.0000},
                          {-0.8817,  0.0000},
                          { 0.5000,  0.0000},
                          { 1.4695,  0.0000},
                          {-2.0566,  0.0000},
                          { 1.4266,  0.0000},
                          {-1.5801,  0.0000},
                          { 0.0000,  0.0000},
                          { 1.5801,  0.0000},
                          {-1.4266,  0.0000},
                          { 2.0566,  0.0000},
                          {-1.4695,  0.0000},
                          {-0.5000,  0.0000},
                          { 0.8817,  0.0000},
                          {-1.7476,  0.0000},
                          { 2.3776,  0.0000},
                          {-0.7711,  0.0000},
                          { 0.0000,  0.0000},
                          { 0.7711,  0.0000},
                          {-2.3776,  0.0000},
                          { 1.7476,  0.0000},
                          {-0.8817,  0.0000},
                          { 0.5000,  0.0000},
                          { 1.4695,  0.0000},
                          {-2.0566,  0.0000},
                          { 1.4266,  0.0000},
                          {-1.5801,  0.0000},
                          { 0.0000,  0.0000},
                          { 1.5801,  0.0000},
                          {-1.4266,  0.0000},
                          { 2.0566,  0.0000},
                          {-1.4695,  0.0000},
                          {-0.5000,  0.0000},
                          { 0.8817,  0.0000},
                          {-1.7476,  0.0000},
                          { 2.3776,  0.0000},
                          {-0.7711,  0.0000},
                          { 0.0000,  0.0000},
                          { 0.7711,  0.0000},
                          {-2.3776,  0.0000},
                          { 1.7476,  0.0000},
                          {-0.8817,  0.0000},
                          { 0.5000,  0.0000},
                          { 1.4695,  0.0000},
                          {-2.0566,  0.0000},
                          { 1.4266,  0.0000},
                          {-1.5801,  0.0000},
                          { 0.0000,  0.0000},
                          { 1.5801,  0.0000},
                          {-1.4266,  0.0000},
                          { 2.0566,  0.0000}
};/*}}}*/
#endif
#if 0
struct cmplx data_a[64] = {/*{{{*/
                          { 1.5000, -1.0000},
                          {-2.5000,  2.5000},
                          { 4.5000,  3.5000},
                          {-3.5000, -2.5000}
};/*}}}*/
struct cmplx data_b[64] = {/*{{{*/
                          { 1.5000, -1.0000},
                          {-2.5000,  2.5000},
                          { 4.5000,  3.5000},
                          {-3.5000, -2.5000}
};/*}}}*/
#endif
void main()
#if 1
{
    MatPlotInit();

    int N=64;
    double x[64], y[64];

    r2dif_nr(N, wTb64, data_a);
    for(int i=0; i<N; i++){
        x[i] = i;
        y[i] = sqrt(pow(data_a[bit_rev(i, 2, round(log(N)/log(2)))].real, 2)+
                    pow(data_a[bit_rev(i, 2, round(log(N)/log(2)))].imag, 2));
    }
    plot(x, y, N);
    printf("////////////////////////////////////////////////////////////////////////////////\n");
    r4dif_nr(N, wTb64, data_b);
    for(int i=0; i<N; i++){
        x[i] = i;
        y[i] = sqrt(pow(data_b[bit_rev(i, 4, round(log(N)/log(4)))].real, 2)+
                    pow(data_b[bit_rev(i, 4, round(log(N)/log(4)))].imag, 2));
    }
    figure(1);
    plot(x, y, N);

    getchar();
    MatPlotClose();
}
#endif
#if 0
{
    double x[5] = {0,1,2,3,4};
    double y[5] = {1,3,2,4,3};
    int N=5;
    MatPlotInit();
    plot(x, y, N);
    getchar();
    MatPlotClose();
}
#endif
void r4dif_nr(int N, struct cmplx *wTb, struct cmplx *b)/*{{{*/
{
    int numOfProblems;
    int problemSize;
    int jFirst, jLast;
    int jTwiddle1, jTwiddle2, jTwiddle3, jTwiddle4;
    int quarterSize;
    int j, k, nr;
    int power=round(log(N)/log(4));
    struct cmplx temp1, temp2, temp3, temp4;
    struct cmplx j1 = {0.0,  1.0};
    for(problemSize=N, numOfProblems=1; problemSize>1; problemSize=quarterSize, numOfProblems*=4){
        quarterSize = problemSize/4;
        for(k=0; k<numOfProblems; k++){
            jFirst = k * problemSize;
            jLast  = jFirst + quarterSize;
            jTwiddle1 = 0;
            jTwiddle2 = 0;
            jTwiddle3 = 0;
            jTwiddle4 = 0;
            for(j=jFirst; j<jLast; j++){
                temp1 = b[j];
                temp2 = b[j+quarterSize];
                temp3 = b[j+quarterSize*2];
                temp4 = b[j+quarterSize*3];
                b[j]                = cmplx_add(
                                                cmplx_add(
                                                            cmplx_add(temp1, temp2),
                                                temp3),
                                      temp4);
                b[j+quarterSize]    = cmplx_mul(
                                                cmplx_add(
                                                            cmplx_minus(
                                                                        cmplx_minus(
                                                                                    temp1,
                                                                        cmplx_mul(j1, temp2)),
                                                            temp3),
                                                cmplx_mul(j1, temp4)),
                                      wTb[jTwiddle2]);
                b[j+quarterSize*2]  = cmplx_mul(
                                                cmplx_minus(
                                                            cmplx_add(
                                                                        cmplx_minus(
                                                                                    temp1,
                                                                        temp2),
                                                            temp3),
                                                temp4),
                                      wTb[jTwiddle3]);
                b[j+quarterSize*3]  = cmplx_mul(
                                                cmplx_minus(
                                                            cmplx_minus(
                                                                        cmplx_add(
                                                                                    temp1,
                                                                        cmplx_mul(j1, temp2)),
                                                            temp3),
                                                cmplx_mul(j1, temp4)),
                                      wTb[jTwiddle4]);
                jTwiddle1 += numOfProblems;
                jTwiddle2 += numOfProblems*1;
                jTwiddle3 += numOfProblems*2;
                jTwiddle4 += numOfProblems*3;
            }
        }
    }
#if 1
    for(nr=0, j=0; j<N; j++){
        nr=bit_rev(j, 4, power);
        if(b[nr].real>=0)
            printf("X(%2d)= %02.6f", j+1, b[nr].real);
        else
            printf("X(%2d)=%02.6f", j+1, b[nr].real);
        if(b[nr].imag>=0)
            printf("+%02.6fi\n", b[nr].imag);
        else
            printf("%02.6fi\n", b[nr].imag);
    }
#endif
}/*}}}*/
void r2dif_nr(int N, struct cmplx *wTb, struct cmplx *b)/*{{{*/
{
    int numOfProblems;
    int problemSize;
    int jFirst, jLast, jTwiddle;
    int halfSize;
    int j, k;
    int power=round(log(N)/log(2));
    struct cmplx temp1, temp2;
    for(problemSize=N, numOfProblems=1; problemSize>1; problemSize /=2, numOfProblems*=2){
        halfSize = problemSize/2;
        for(k=0; k<numOfProblems; k++){
            jFirst = k * problemSize;
            jLast  = jFirst + halfSize;
            jTwiddle = 0;
            for(j=jFirst; j<jLast; j++, jTwiddle+=numOfProblems){
                temp1 = cmplx_add(b[j], b[j+halfSize]);
                temp2 = cmplx_mul(cmplx_minus(b[j], b[j+halfSize]), wTb[jTwiddle]);
                b[j]          = temp1;
                b[j+halfSize] = temp2;
            }
        }
    }
#if 1
    for(j=0; j<N; j++){
        if(b[bit_rev(j, 2, power)].real>=0)
            printf("X(%2d)= %02.6f", j+1, b[bit_rev(j, 2, power)].real);
        else
            printf("X(%2d)=%02.6f", j+1, b[bit_rev(j, 2, power)].real);
        if(b[bit_rev(j, 2, power)].imag>=0)
            printf("+%02.6fi\n", b[bit_rev(j, 2, power)].imag);
        else
            printf("%02.6fi\n", b[bit_rev(j, 2, power)].imag);
    }
#endif
}/*}}}*/
void r2dit_nr(int N, struct cmplx *wTb, struct cmplx *b)/*{{{*/
{
    int pairsInGroup;
    int numOfGroups;
    int distance;
    int jFirst, jLast, jTwiddle;
    int power=round(log(N)/log(2));
    struct cmplx temp1, temp2, temp3;
    for(pairsInGroup=N/2, numOfGroups=1, distance=N/2; numOfGroups<N; pairsInGroup/=2, numOfGroups*=2, distance/=2){
        for(int k=0; k<numOfGroups; k++){
            jFirst = 2 * k * pairsInGroup;
            jLast  = jFirst + pairsInGroup;
            jTwiddle = k;
            for(int j=jFirst; j<jLast; j++){
                temp1 = cmplx_mul(b[j+distance], wTb[bit_rev(jTwiddle, 2, power-1)]);
                temp2 = cmplx_add(b[j], temp1);
                temp3 = cmplx_minus(b[j], temp1);
                b[j]          = temp2;
                b[j+distance] = temp3;
            }
        }
    }
#if 1
    for(int j=0; j<N; j++){
        printf("X[%d]=%f+j%f\n", j, b[bit_rev(j, 2, power)].real, b[bit_rev(j, 2, power)].imag);
    }
#endif
}/*}}}*/
int bit_rev(int din, int radix, int len)/*{{{*/
{
	int i, rslt;
    int sft=round(log(radix)/log(2));
    for(i=0, rslt=0; i<len; i++){
        rslt <<= sft;
        rslt |= din & ((1<<sft)-1);
        din  >>= sft;
    }
    return rslt;
}/*}}}*/
struct cmplx cmplx_add(struct cmplx a, struct cmplx b)/*{{{*/
{
    struct cmplx rslt;
    rslt.real = a.real + b.real;
    rslt.imag = a.imag + b.imag;
    return rslt;
}/*}}}*/
struct cmplx cmplx_minus(struct cmplx a, struct cmplx b)/*{{{*/
{
    struct cmplx rslt;
    rslt.real = a.real - b.real;
    rslt.imag = a.imag - b.imag;
    return rslt;
}/*}}}*/
struct cmplx cmplx_mul(struct cmplx a, struct cmplx b)/*{{{*/
{
    struct cmplx rslt;
    double temp1, temp2, temp3, temp4, temp5, temp6;
    temp1 = a.real - a.imag;    // A-B
    temp2 = b.real - b.imag;    // C-D
    temp3 = a.real + a.imag;    // A+B

    temp4 = b.real * temp1;     // C(A-B)
    temp5 = a.imag * temp2;     // B(C-D)
    temp6 = b.imag * temp3;     // D(A+B)

    rslt.real = temp4 + temp5;  // real = C(A-B) + B(C-D)
    rslt.imag = temp6 + temp5;  // imag = D(A+B) + B(C-D)
    return rslt;
}/*}}}*/
// vim: fdm=marker
