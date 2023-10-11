
void    rhsrecip(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double *that,int* teid,double* rhs);
void rhsrecip2nd(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double *that,int* teid,double* rhs);
void rhsrecip3rd(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double *that,int* teid,double* rhs);


void    conductioncurr(const  int* te2p,const  int nte,const double* p,const double* sigma, const double* u,  double* rv,  double* jv);
void conductioncurr2nd(const  int* te2p,const  int nte,const double* p,const double* sigma, const double* u,  double* rv,  double* jv);
void conductioncurr3rd(const  int* te2p,const  int nte,const double* p,const double* sigma, const double* u,  double* rv,  double* jv);

void    conductioncurrTensor(const  int* te2p,const  int nte,const double* p,const double* sigma, const double* u,  double* rv,  double* jv);
void conductioncurr2ndTensor(const  int* te2p,const  int nte,const double* p,const double* sigma, const double* u,  double* rv,  double* jv);
void conductioncurr3rdTensor(const  int* te2p,const  int nte,const double* p,const double* sigma, const double* u,  double* rv,  double* jv);

void impressedcurr(const  int* te2p,const  int nte,const double* p , const double* that,const int* teid,int nsource,  double* rv,  double* jv);

void Hprim(int nsource,int ntarget, double *rs,double* js,double* robs, double* Hprimary,double iprec);
void Eprim(int nsource,int ntarget, double *rs,double* js,double* robs, double* Eprimary,double iprec);
void Aprim(int nsource,int ntarget, double *rs,double* js,double* robs, double* Aprimary,double* curlAprimary,double iprec);

   void rhsfunc(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double* rs,double* js,double* rhs,double iprec);
void rhsfunc2nd(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double *rs,double* js,double* rhs,double iprec);
void rhsfunc3rd(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double *rs,double* js,double* rhs,double iprec);

void rhsfuncks(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double* rs,double* js,double* rhs,double iprec);
void rhsfunc2ndks(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double *rs,double* js,double* rhs,double iprec);
void rhsfunc3rdks(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double *rs,double* js,double* rhs,double iprec);
void rhsfunckstensor(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double *rs,double* js,double* rhs,double iprec);


void assemblemat(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int *crvindex,double* val);
void assemblemat2nd(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int *crvindex,double* val);
void assemblemat3rd(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int *crvindex,double* val);


void assemblemattensor(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int* crvindex,double* val);
void assemblemat2ndtensor(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int* crvindex,double* val);
void assemblemat3rdtensor(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int* crvindex,double* val);


void pointLocation_f(const  double* robs,const  int nobs,const  int* te2p,const  int nte,const double* p
  ,const int te_stride,const int* te2te,  int* teid, double* bari);

void secondaryEfield(const  double* robs,const  int nobs,const int* te2te,const  int* te2p,const  int nte,const double* p
  ,const double* u,  double* soln,const double iprec);
void secondaryEfield2nd(const  double* robs,const  int nobs,const int* te2te,const  int* te2p,const  int nte,const double* p
  ,const double* u,  double* soln,const double iprec);
void secondaryEfield3rd(const  double* robs,const  int nobs,const int* te2te,const  int* te2p,const  int nte,const double* p
  ,const double* u,  double* soln,const double iprec);


void secondaryEfieldc(const  double* robs,const  int nobs,const int* te2te,const  int* te2p,const  int nte,const double* p
  ,const double* u,  double* soln,const double iprec);
void secondaryEfield2ndc(const  double* robs,const  int nobs,const int* te2te,const  int* te2p,const  int nte,const double* p
  ,const double* u,  double* soln,const double iprec);
void secondaryEfield3rdc(const  double* robs,const  int nobs,const int* te2te,const  int* te2p,const  int nte,const double* p
  ,const double* u,  double* soln,const double iprec);



void generatemesh3rd(const int nte,const  int* te2p,const double* p,const int np,const int ne,int* te2p2,double* p2);
void generatemesh2nd(const int nte,const  int* te2p,const double* p,const int np,int* te2p2,double* p2);

int generateedge_ct(const int nte,const  int* te2p);
int generateface_ct(const int nte,const  int* te2p);

void generatetetrapaths(const int nte,const  int* te2p,int* te2te);
#include <uchar.h>
/* --------------------------------------------------- */
/* Automatically generated by mwrap                    */
/* --------------------------------------------------- */

/* Code generated by mwrap */
/*
  Copyright statement for mwrap:

  mwrap -- MEX file generation for MATLAB and Octave
  Copyright (c) 2007-2008 David Bindel

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  You may distribute a work that contains part or all of the source code
  generated by mwrap under the terms of your choice.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  THE SOFTWARE.
*/

#include <mex.h>
#include <stdio.h>
#include <string.h>


#ifndef ulong
#  define ulong unsigned long
#endif
#ifndef uint
#  define uint  unsigned int
#endif
#ifndef uchar
#  define uchar unsigned char
#endif


/*
 * Records for call profile.
 */
int* mexprofrecord_= NULL;


/*
 * Support routines for copying data into and out of the MEX stubs
 */

void* mxWrapGetP(const mxArray* a, const char* fmt, const char** e)
{
    void* p = 0;
    mxArray* ap;
    if (mxGetClassID(a) == mxDOUBLE_CLASS && 
        mxGetM(a)*mxGetN(a) == 1 && *mxGetPr(a) == 0)
        return p;
    if (mxIsChar(a)) {
        char pbuf[128];
        mxGetString(a, pbuf, sizeof(pbuf));
        sscanf(pbuf, fmt, &p);
    } 
#ifdef R2008OO
    else if (ap = mxGetProperty(a, 0, "mwptr")) {
        return mxWrapGetP(ap, fmt, e);
    }
#endif
    if (p == 0)
        *e = "Invalid pointer";
    return p;
}

mxArray* mxWrapCreateP(void* p, const char* fmt)
{
    if (p == 0) {
        mxArray* z = mxCreateDoubleMatrix(1,1, mxREAL);
        *mxGetPr(z) = 0;
        return z;
    } else {
        char pbuf[128];
        sprintf(pbuf, fmt, p);
        return mxCreateString(pbuf);
    }
}

mxArray* mxWrapStrncpy(const char* s)
{
    if (s) {
        return mxCreateString(s);
    } else {
        mxArray* z = mxCreateDoubleMatrix(1,1, mxREAL);
        *mxGetPr(z) = 0;
        return z;
    }
}

double mxWrapGetScalar(const mxArray* a, const char** e)
{
    if (!a || mxGetClassID(a) != mxDOUBLE_CLASS || mxGetM(a)*mxGetN(a) != 1) {
        *e = "Invalid scalar argument";
        return 0;
    }
    return *mxGetPr(a);
}

char* mxWrapGetString(const mxArray* a, const char** e)
{
    char* s;
    int slen;
    if (!a || (!mxIsChar(a) && mxGetM(a)*mxGetN(a) > 0)) {
        *e = "Invalid string argument";
        return NULL;
    }
    slen = mxGetM(a)*mxGetN(a) + 1;
    s = (char*) mxMalloc(slen);
    if (mxGetM(a)*mxGetN(a) == 0)
        *s = 0;
    else
        mxGetString(a, s, slen);
    return s;
}


#define mxWrapGetArrayDef(func, T) \
T* func(const mxArray* a, const char** e)     \
{ \
    T* array; \
    int arraylen; \
    int i; \
    T* p; \
    double* q; \
    if (!a || mxGetClassID(a) != mxDOUBLE_CLASS) { \
        *e = "Invalid array argument"; \
        return 0; \
    } \
    arraylen = mxGetM(a)*mxGetN(a); \
    array = (T*) mxMalloc(mxGetM(a)*mxGetN(a) * sizeof(T)); \
    p = array; \
    q = mxGetPr(a); \
    for (i = 0; i < arraylen; ++i) \
        *p++ = (T) (*q++); \
    return array; \
}


#define mxWrapCopyDef(func, T) \
void func(mxArray* a, const T* q, int n) \
{ \
    int i; \
    double* p = mxGetPr(a); \
    for (i = 0; i < n; ++i) \
        *p++ = *q++; \
}


#define mxWrapReturnDef(func, T) \
mxArray* func(const T* q, int m, int n) \
{ \
    int i; \
    double* p; \
    if (!q) { \
        return mxCreateDoubleMatrix(0,0, mxREAL); \
    } else { \
        mxArray* a = mxCreateDoubleMatrix(m,n, mxREAL); \
        p = mxGetPr(a); \
        for (i = 0; i < m*n; ++i) \
            *p++ = *q++; \
        return a; \
    } \
}


#define mxWrapGetScalarZDef(func, T, ZT, setz) \
void func(T* z, const mxArray* a) \
{ \
    double* pr = mxGetPr(a); \
    double* pi = mxGetPi(a); \
    setz(z, (ZT) *pr, (pi ? (ZT) *pi : (ZT) 0)); \
}


#define mxWrapGetArrayZDef(func, T, ZT, setz) \
T* func(const mxArray* a, const char** e) \
{ \
    T* array; \
    int arraylen; \
    int i; \
    T* p; \
    double* qr; \
    double* qi; \
    if (!a || mxGetClassID(a) != mxDOUBLE_CLASS) { \
        *e = "Invalid array argument"; \
        return 0; \
    } \
    arraylen = mxGetM(a)*mxGetN(a); \
    array = (T*) mxMalloc(mxGetM(a)*mxGetN(a) * sizeof(T)); \
    p = array; \
    qr = mxGetPr(a); \
    qi = mxGetPi(a); \
    for (i = 0; i < arraylen; ++i) { \
        ZT val_qr = *qr++; \
        ZT val_qi = (qi ? (ZT) *qi++ : (ZT) 0); \
        setz(p, val_qr, val_qi); \
        ++p; \
    } \
    return array; \
}


#define mxWrapCopyZDef(func, T, real, imag) \
void func(mxArray* a, const T* q, int n) \
{ \
    int i; \
    double* pr = mxGetPr(a); \
    double* pi = mxGetPi(a); \
    for (i = 0; i < n; ++i) { \
        *pr++ = real(*q); \
        *pi++ = imag(*q); \
        ++q; \
    } \
}


#define mxWrapReturnZDef(func, T, real, imag) \
mxArray* func(const T* q, int m, int n) \
{ \
    int i; \
    double* pr; \
    double* pi; \
    if (!q) { \
        return mxCreateDoubleMatrix(0,0, mxCOMPLEX); \
    } else { \
        mxArray* a = mxCreateDoubleMatrix(m,n, mxCOMPLEX); \
        pr = mxGetPr(a); \
        pi = mxGetPi(a); \
        for (i = 0; i < m*n; ++i) { \
            *pr++ = real(*q); \
            *pi++ = imag(*q); \
            ++q; \
        } \
        return a; \
    } \
}

/* Array copier definitions */
mxWrapGetArrayDef(mxWrapGetArray_bool, bool)
mxWrapCopyDef    (mxWrapCopy_bool,     bool)
mxWrapReturnDef  (mxWrapReturn_bool,   bool)
mxWrapGetArrayDef(mxWrapGetArray_char, char)
mxWrapCopyDef    (mxWrapCopy_char,     char)
mxWrapReturnDef  (mxWrapReturn_char,   char)
mxWrapGetArrayDef(mxWrapGetArray_double, double)
mxWrapCopyDef    (mxWrapCopy_double,     double)
mxWrapReturnDef  (mxWrapReturn_double,   double)
mxWrapGetArrayDef(mxWrapGetArray_float, float)
mxWrapCopyDef    (mxWrapCopy_float,     float)
mxWrapReturnDef  (mxWrapReturn_float,   float)
mxWrapGetArrayDef(mxWrapGetArray_int, int)
mxWrapCopyDef    (mxWrapCopy_int,     int)
mxWrapReturnDef  (mxWrapReturn_int,   int)
mxWrapGetArrayDef(mxWrapGetArray_long, long)
mxWrapCopyDef    (mxWrapCopy_long,     long)
mxWrapReturnDef  (mxWrapReturn_long,   long)
mxWrapGetArrayDef(mxWrapGetArray_size_t, size_t)
mxWrapCopyDef    (mxWrapCopy_size_t,     size_t)
mxWrapReturnDef  (mxWrapReturn_size_t,   size_t)
mxWrapGetArrayDef(mxWrapGetArray_uchar, uchar)
mxWrapCopyDef    (mxWrapCopy_uchar,     uchar)
mxWrapReturnDef  (mxWrapReturn_uchar,   uchar)
mxWrapGetArrayDef(mxWrapGetArray_uint, uint)
mxWrapCopyDef    (mxWrapCopy_uint,     uint)
mxWrapReturnDef  (mxWrapReturn_uint,   uint)
mxWrapGetArrayDef(mxWrapGetArray_ulong, ulong)
mxWrapCopyDef    (mxWrapCopy_ulong,     ulong)
mxWrapReturnDef  (mxWrapReturn_ulong,   ulong)

/* ---- demo.mw: 100 ----
 * rhsrecip(int nte, int np, int[4, nte] te2p2, double[3, np] p2, double[9, nte] conductivity, int nsource, inout double[3, nsource] that, int[nsource] teid, inout double[np] rhs);
 */
const char* stubids1_ = "rhsrecip(i int, i int, i int[xx], i double[xx], i double[xx], i int, io double[xx], i int[x], io double[x])";

void mexStub1(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nte        */
    int         in1_;    /* np         */
    int*        in2_ =0; /* te2p2      */
    double*     in3_ =0; /* p2         */
    double*     in4_ =0; /* conductivity */
    int         in5_;    /* nsource    */
    double*     in6_ =0; /* that       */
    int*        in7_ =0; /* teid       */
    double*     in8_ =0; /* rhs        */
    int         dim9_;   /* 4          */
    int         dim10_;   /* nte        */
    int         dim11_;   /* 3          */
    int         dim12_;   /* np         */
    int         dim13_;   /* 9          */
    int         dim14_;   /* nte        */
    int         dim15_;   /* 3          */
    int         dim16_;   /* nsource    */
    int         dim17_;   /* nsource    */
    int         dim18_;   /* np         */

    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);
    dim18_ = (int) mxWrapGetScalar(prhs[18], &mw_err_txt_);

    if (mxGetM(prhs[2]) != dim9_ ||
        mxGetN(prhs[2]) != dim10_) {
        mw_err_txt_ = "Bad argument size: te2p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim11_ ||
        mxGetN(prhs[3]) != dim12_) {
        mw_err_txt_ = "Bad argument size: p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4]) != dim13_ ||
        mxGetN(prhs[4]) != dim14_) {
        mw_err_txt_ = "Bad argument size: conductivity";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6]) != dim15_ ||
        mxGetN(prhs[6]) != dim16_) {
        mw_err_txt_ = "Bad argument size: that";
        goto mw_err_label;
    }

    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != dim17_) {
        mw_err_txt_ = "Bad argument size: teid";        goto mw_err_label;
    }

    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != dim18_) {
        mw_err_txt_ = "Bad argument size: rhs";        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    in5_ = (int) mxWrapGetScalar(prhs[5], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_double(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxWrapGetArray_int(prhs[7], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in7_ = NULL;
    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != 0) {
        in8_ = mxWrapGetArray_double(prhs[8], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in8_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[1]++;
    rhsrecip(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_);
    plhs[0] = mxCreateDoubleMatrix(dim15_, dim16_, mxREAL);
    mxWrapCopy_double(plhs[0], in6_, dim15_*dim16_);
    plhs[1] = mxCreateDoubleMatrix(dim18_, 1, mxREAL);
    mxWrapCopy_double(plhs[1], in8_, dim18_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in6_)  mxFree(in6_);
    if (in7_)  mxFree(in7_);
    if (in8_)  mxFree(in8_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 102 ----
 * rhsrecip2nd(int nte, int np, int[10, nte] te2p2, double[3, np] p2, double[9, nte] conductivity, int nsource, inout double[3, nsource] that, int[nsource] teid, inout double[np] rhs);
 */
const char* stubids2_ = "rhsrecip2nd(i int, i int, i int[xx], i double[xx], i double[xx], i int, io double[xx], i int[x], io double[x])";

void mexStub2(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nte        */
    int         in1_;    /* np         */
    int*        in2_ =0; /* te2p2      */
    double*     in3_ =0; /* p2         */
    double*     in4_ =0; /* conductivity */
    int         in5_;    /* nsource    */
    double*     in6_ =0; /* that       */
    int*        in7_ =0; /* teid       */
    double*     in8_ =0; /* rhs        */
    int         dim9_;   /* 10         */
    int         dim10_;   /* nte        */
    int         dim11_;   /* 3          */
    int         dim12_;   /* np         */
    int         dim13_;   /* 9          */
    int         dim14_;   /* nte        */
    int         dim15_;   /* 3          */
    int         dim16_;   /* nsource    */
    int         dim17_;   /* nsource    */
    int         dim18_;   /* np         */

    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);
    dim18_ = (int) mxWrapGetScalar(prhs[18], &mw_err_txt_);

    if (mxGetM(prhs[2]) != dim9_ ||
        mxGetN(prhs[2]) != dim10_) {
        mw_err_txt_ = "Bad argument size: te2p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim11_ ||
        mxGetN(prhs[3]) != dim12_) {
        mw_err_txt_ = "Bad argument size: p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4]) != dim13_ ||
        mxGetN(prhs[4]) != dim14_) {
        mw_err_txt_ = "Bad argument size: conductivity";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6]) != dim15_ ||
        mxGetN(prhs[6]) != dim16_) {
        mw_err_txt_ = "Bad argument size: that";
        goto mw_err_label;
    }

    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != dim17_) {
        mw_err_txt_ = "Bad argument size: teid";        goto mw_err_label;
    }

    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != dim18_) {
        mw_err_txt_ = "Bad argument size: rhs";        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    in5_ = (int) mxWrapGetScalar(prhs[5], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_double(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxWrapGetArray_int(prhs[7], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in7_ = NULL;
    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != 0) {
        in8_ = mxWrapGetArray_double(prhs[8], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in8_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[2]++;
    rhsrecip2nd(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_);
    plhs[0] = mxCreateDoubleMatrix(dim15_, dim16_, mxREAL);
    mxWrapCopy_double(plhs[0], in6_, dim15_*dim16_);
    plhs[1] = mxCreateDoubleMatrix(dim18_, 1, mxREAL);
    mxWrapCopy_double(plhs[1], in8_, dim18_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in6_)  mxFree(in6_);
    if (in7_)  mxFree(in7_);
    if (in8_)  mxFree(in8_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 104 ----
 * rhsrecip3rd(int nte, int np, int[20, nte] te2p2, double[3, np] p2, double[9, nte] conductivity, int nsource, inout double[3, nsource] that, int[nsource] teid, inout double[np] rhs);
 */
const char* stubids3_ = "rhsrecip3rd(i int, i int, i int[xx], i double[xx], i double[xx], i int, io double[xx], i int[x], io double[x])";

void mexStub3(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nte        */
    int         in1_;    /* np         */
    int*        in2_ =0; /* te2p2      */
    double*     in3_ =0; /* p2         */
    double*     in4_ =0; /* conductivity */
    int         in5_;    /* nsource    */
    double*     in6_ =0; /* that       */
    int*        in7_ =0; /* teid       */
    double*     in8_ =0; /* rhs        */
    int         dim9_;   /* 20         */
    int         dim10_;   /* nte        */
    int         dim11_;   /* 3          */
    int         dim12_;   /* np         */
    int         dim13_;   /* 9          */
    int         dim14_;   /* nte        */
    int         dim15_;   /* 3          */
    int         dim16_;   /* nsource    */
    int         dim17_;   /* nsource    */
    int         dim18_;   /* np         */

    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);
    dim18_ = (int) mxWrapGetScalar(prhs[18], &mw_err_txt_);

    if (mxGetM(prhs[2]) != dim9_ ||
        mxGetN(prhs[2]) != dim10_) {
        mw_err_txt_ = "Bad argument size: te2p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim11_ ||
        mxGetN(prhs[3]) != dim12_) {
        mw_err_txt_ = "Bad argument size: p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4]) != dim13_ ||
        mxGetN(prhs[4]) != dim14_) {
        mw_err_txt_ = "Bad argument size: conductivity";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6]) != dim15_ ||
        mxGetN(prhs[6]) != dim16_) {
        mw_err_txt_ = "Bad argument size: that";
        goto mw_err_label;
    }

    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != dim17_) {
        mw_err_txt_ = "Bad argument size: teid";        goto mw_err_label;
    }

    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != dim18_) {
        mw_err_txt_ = "Bad argument size: rhs";        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    in5_ = (int) mxWrapGetScalar(prhs[5], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_double(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxWrapGetArray_int(prhs[7], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in7_ = NULL;
    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != 0) {
        in8_ = mxWrapGetArray_double(prhs[8], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in8_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[3]++;
    rhsrecip3rd(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_);
    plhs[0] = mxCreateDoubleMatrix(dim15_, dim16_, mxREAL);
    mxWrapCopy_double(plhs[0], in6_, dim15_*dim16_);
    plhs[1] = mxCreateDoubleMatrix(dim18_, 1, mxREAL);
    mxWrapCopy_double(plhs[1], in8_, dim18_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in6_)  mxFree(in6_);
    if (in7_)  mxFree(in7_);
    if (in8_)  mxFree(in8_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 122 ----
 * conductioncurrTensor(int[4, nte] te2p2, int nte, double[3, np] p2, double[9, nte] conductivity, double[np] x, inout double[3, ntens] rv, inout double[3, ntens] jv);
 */
const char* stubids4_ = "conductioncurrTensor(i int[xx], i int, i double[xx], i double[xx], i double[x], io double[xx], io double[xx])";

void mexStub4(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int*        in0_ =0; /* te2p2      */
    int         in1_;    /* nte        */
    double*     in2_ =0; /* p2         */
    double*     in3_ =0; /* conductivity */
    double*     in4_ =0; /* x          */
    double*     in5_ =0; /* rv         */
    double*     in6_ =0; /* jv         */
    int         dim7_;   /* 4          */
    int         dim8_;   /* nte        */
    int         dim9_;   /* 3          */
    int         dim10_;   /* np         */
    int         dim11_;   /* 9          */
    int         dim12_;   /* nte        */
    int         dim13_;   /* np         */
    int         dim14_;   /* 3          */
    int         dim15_;   /* ntens      */
    int         dim16_;   /* 3          */
    int         dim17_;   /* ntens      */

    dim7_ = (int) mxWrapGetScalar(prhs[7], &mw_err_txt_);
    dim8_ = (int) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);

    if (mxGetM(prhs[0]) != dim7_ ||
        mxGetN(prhs[0]) != dim8_) {
        mw_err_txt_ = "Bad argument size: te2p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[2]) != dim9_ ||
        mxGetN(prhs[2]) != dim10_) {
        mw_err_txt_ = "Bad argument size: p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim11_ ||
        mxGetN(prhs[3]) != dim12_) {
        mw_err_txt_ = "Bad argument size: conductivity";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim13_) {
        mw_err_txt_ = "Bad argument size: x";        goto mw_err_label;
    }

    if (mxGetM(prhs[5]) != dim14_ ||
        mxGetN(prhs[5]) != dim15_) {
        mw_err_txt_ = "Bad argument size: rv";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6]) != dim16_ ||
        mxGetN(prhs[6]) != dim17_) {
        mw_err_txt_ = "Bad argument size: jv";
        goto mw_err_label;
    }

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != 0) {
        in0_ = mxWrapGetArray_int(prhs[0], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in0_ = NULL;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxGetPr(prhs[2]);
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxWrapGetArray_double(prhs[5], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_double(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[4]++;
    conductioncurrTensor(in0_, in1_, in2_, in3_, in4_, in5_, in6_);
    plhs[0] = mxCreateDoubleMatrix(dim14_, dim15_, mxREAL);
    mxWrapCopy_double(plhs[0], in5_, dim14_*dim15_);
    plhs[1] = mxCreateDoubleMatrix(dim16_, dim17_, mxREAL);
    mxWrapCopy_double(plhs[1], in6_, dim16_*dim17_);

mw_err_label:
    if (in0_)  mxFree(in0_);
    if (in5_)  mxFree(in5_);
    if (in6_)  mxFree(in6_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 124 ----
 * conductioncurr2ndTensor(int[10, nte] te2p2, int nte, double[3, np] p2, double[9, nte] conductivity, double[np] x, inout double[3, ntens] rv, inout double[3, ntens] jv);
 */
const char* stubids5_ = "conductioncurr2ndTensor(i int[xx], i int, i double[xx], i double[xx], i double[x], io double[xx], io double[xx])";

void mexStub5(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int*        in0_ =0; /* te2p2      */
    int         in1_;    /* nte        */
    double*     in2_ =0; /* p2         */
    double*     in3_ =0; /* conductivity */
    double*     in4_ =0; /* x          */
    double*     in5_ =0; /* rv         */
    double*     in6_ =0; /* jv         */
    int         dim7_;   /* 10         */
    int         dim8_;   /* nte        */
    int         dim9_;   /* 3          */
    int         dim10_;   /* np         */
    int         dim11_;   /* 9          */
    int         dim12_;   /* nte        */
    int         dim13_;   /* np         */
    int         dim14_;   /* 3          */
    int         dim15_;   /* ntens      */
    int         dim16_;   /* 3          */
    int         dim17_;   /* ntens      */

    dim7_ = (int) mxWrapGetScalar(prhs[7], &mw_err_txt_);
    dim8_ = (int) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);

    if (mxGetM(prhs[0]) != dim7_ ||
        mxGetN(prhs[0]) != dim8_) {
        mw_err_txt_ = "Bad argument size: te2p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[2]) != dim9_ ||
        mxGetN(prhs[2]) != dim10_) {
        mw_err_txt_ = "Bad argument size: p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim11_ ||
        mxGetN(prhs[3]) != dim12_) {
        mw_err_txt_ = "Bad argument size: conductivity";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim13_) {
        mw_err_txt_ = "Bad argument size: x";        goto mw_err_label;
    }

    if (mxGetM(prhs[5]) != dim14_ ||
        mxGetN(prhs[5]) != dim15_) {
        mw_err_txt_ = "Bad argument size: rv";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6]) != dim16_ ||
        mxGetN(prhs[6]) != dim17_) {
        mw_err_txt_ = "Bad argument size: jv";
        goto mw_err_label;
    }

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != 0) {
        in0_ = mxWrapGetArray_int(prhs[0], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in0_ = NULL;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxGetPr(prhs[2]);
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxWrapGetArray_double(prhs[5], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_double(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[5]++;
    conductioncurr2ndTensor(in0_, in1_, in2_, in3_, in4_, in5_, in6_);
    plhs[0] = mxCreateDoubleMatrix(dim14_, dim15_, mxREAL);
    mxWrapCopy_double(plhs[0], in5_, dim14_*dim15_);
    plhs[1] = mxCreateDoubleMatrix(dim16_, dim17_, mxREAL);
    mxWrapCopy_double(plhs[1], in6_, dim16_*dim17_);

mw_err_label:
    if (in0_)  mxFree(in0_);
    if (in5_)  mxFree(in5_);
    if (in6_)  mxFree(in6_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 126 ----
 * conductioncurr3rdTensor(int[20, nte] te2p2, int nte, double[3, np] p2, double[9, nte] conductivity, double[np] x, inout double[3, ntens] rv, inout double[3, ntens] jv);
 */
const char* stubids6_ = "conductioncurr3rdTensor(i int[xx], i int, i double[xx], i double[xx], i double[x], io double[xx], io double[xx])";

void mexStub6(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int*        in0_ =0; /* te2p2      */
    int         in1_;    /* nte        */
    double*     in2_ =0; /* p2         */
    double*     in3_ =0; /* conductivity */
    double*     in4_ =0; /* x          */
    double*     in5_ =0; /* rv         */
    double*     in6_ =0; /* jv         */
    int         dim7_;   /* 20         */
    int         dim8_;   /* nte        */
    int         dim9_;   /* 3          */
    int         dim10_;   /* np         */
    int         dim11_;   /* 9          */
    int         dim12_;   /* nte        */
    int         dim13_;   /* np         */
    int         dim14_;   /* 3          */
    int         dim15_;   /* ntens      */
    int         dim16_;   /* 3          */
    int         dim17_;   /* ntens      */

    dim7_ = (int) mxWrapGetScalar(prhs[7], &mw_err_txt_);
    dim8_ = (int) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);

    if (mxGetM(prhs[0]) != dim7_ ||
        mxGetN(prhs[0]) != dim8_) {
        mw_err_txt_ = "Bad argument size: te2p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[2]) != dim9_ ||
        mxGetN(prhs[2]) != dim10_) {
        mw_err_txt_ = "Bad argument size: p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim11_ ||
        mxGetN(prhs[3]) != dim12_) {
        mw_err_txt_ = "Bad argument size: conductivity";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim13_) {
        mw_err_txt_ = "Bad argument size: x";        goto mw_err_label;
    }

    if (mxGetM(prhs[5]) != dim14_ ||
        mxGetN(prhs[5]) != dim15_) {
        mw_err_txt_ = "Bad argument size: rv";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6]) != dim16_ ||
        mxGetN(prhs[6]) != dim17_) {
        mw_err_txt_ = "Bad argument size: jv";
        goto mw_err_label;
    }

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != 0) {
        in0_ = mxWrapGetArray_int(prhs[0], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in0_ = NULL;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxGetPr(prhs[2]);
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxWrapGetArray_double(prhs[5], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_double(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[6]++;
    conductioncurr3rdTensor(in0_, in1_, in2_, in3_, in4_, in5_, in6_);
    plhs[0] = mxCreateDoubleMatrix(dim14_, dim15_, mxREAL);
    mxWrapCopy_double(plhs[0], in5_, dim14_*dim15_);
    plhs[1] = mxCreateDoubleMatrix(dim16_, dim17_, mxREAL);
    mxWrapCopy_double(plhs[1], in6_, dim16_*dim17_);

mw_err_label:
    if (in0_)  mxFree(in0_);
    if (in5_)  mxFree(in5_);
    if (in6_)  mxFree(in6_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 130 ----
 * impressedcurr(int[4, nte] te2p, int nte, double[3, np] p, double[3, nsource] that, int[nsource] teid, int nsource, inout double[3, ntens] rv, inout double[3, ntens] jv);
 * Also at demo.mw: 187
 */
const char* stubids7_ = "impressedcurr(i int[xx], i int, i double[xx], i double[xx], i int[x], i int, io double[xx], io double[xx])";

void mexStub7(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int*        in0_ =0; /* te2p       */
    int         in1_;    /* nte        */
    double*     in2_ =0; /* p          */
    double*     in3_ =0; /* that       */
    int*        in4_ =0; /* teid       */
    int         in5_;    /* nsource    */
    double*     in6_ =0; /* rv         */
    double*     in7_ =0; /* jv         */
    int         dim8_;   /* 4          */
    int         dim9_;   /* nte        */
    int         dim10_;   /* 3          */
    int         dim11_;   /* np         */
    int         dim12_;   /* 3          */
    int         dim13_;   /* nsource    */
    int         dim14_;   /* nsource    */
    int         dim15_;   /* 3          */
    int         dim16_;   /* ntens      */
    int         dim17_;   /* 3          */
    int         dim18_;   /* ntens      */

    dim8_ = (int) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);
    dim18_ = (int) mxWrapGetScalar(prhs[18], &mw_err_txt_);

    if (mxGetM(prhs[0]) != dim8_ ||
        mxGetN(prhs[0]) != dim9_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[2]) != dim10_ ||
        mxGetN(prhs[2]) != dim11_) {
        mw_err_txt_ = "Bad argument size: p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim12_ ||
        mxGetN(prhs[3]) != dim13_) {
        mw_err_txt_ = "Bad argument size: that";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim14_) {
        mw_err_txt_ = "Bad argument size: teid";        goto mw_err_label;
    }

    if (mxGetM(prhs[6]) != dim15_ ||
        mxGetN(prhs[6]) != dim16_) {
        mw_err_txt_ = "Bad argument size: rv";
        goto mw_err_label;
    }

    if (mxGetM(prhs[7]) != dim17_ ||
        mxGetN(prhs[7]) != dim18_) {
        mw_err_txt_ = "Bad argument size: jv";
        goto mw_err_label;
    }

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != 0) {
        in0_ = mxWrapGetArray_int(prhs[0], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in0_ = NULL;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxGetPr(prhs[2]);
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxWrapGetArray_int(prhs[4], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in4_ = NULL;
    in5_ = (int) mxWrapGetScalar(prhs[5], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_double(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxWrapGetArray_double(prhs[7], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in7_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[7]++;
    impressedcurr(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_);
    plhs[0] = mxCreateDoubleMatrix(dim15_, dim16_, mxREAL);
    mxWrapCopy_double(plhs[0], in6_, dim15_*dim16_);
    plhs[1] = mxCreateDoubleMatrix(dim17_, dim18_, mxREAL);
    mxWrapCopy_double(plhs[1], in7_, dim17_*dim18_);

mw_err_label:
    if (in0_)  mxFree(in0_);
    if (in4_)  mxFree(in4_);
    if (in6_)  mxFree(in6_);
    if (in7_)  mxFree(in7_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 157 ----
 * rhsrecip(int nte, int np, int[4, nte] te2p2, double[3, np] p2, double[nte] conductivity, int nsource, inout double[3, nsource] that, int[nsource] teid, inout double[np] rhs);
 */
const char* stubids8_ = "rhsrecip(i int, i int, i int[xx], i double[xx], i double[x], i int, io double[xx], i int[x], io double[x])";

void mexStub8(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nte        */
    int         in1_;    /* np         */
    int*        in2_ =0; /* te2p2      */
    double*     in3_ =0; /* p2         */
    double*     in4_ =0; /* conductivity */
    int         in5_;    /* nsource    */
    double*     in6_ =0; /* that       */
    int*        in7_ =0; /* teid       */
    double*     in8_ =0; /* rhs        */
    int         dim9_;   /* 4          */
    int         dim10_;   /* nte        */
    int         dim11_;   /* 3          */
    int         dim12_;   /* np         */
    int         dim13_;   /* nte        */
    int         dim14_;   /* 3          */
    int         dim15_;   /* nsource    */
    int         dim16_;   /* nsource    */
    int         dim17_;   /* np         */

    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);

    if (mxGetM(prhs[2]) != dim9_ ||
        mxGetN(prhs[2]) != dim10_) {
        mw_err_txt_ = "Bad argument size: te2p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim11_ ||
        mxGetN(prhs[3]) != dim12_) {
        mw_err_txt_ = "Bad argument size: p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim13_) {
        mw_err_txt_ = "Bad argument size: conductivity";        goto mw_err_label;
    }

    if (mxGetM(prhs[6]) != dim14_ ||
        mxGetN(prhs[6]) != dim15_) {
        mw_err_txt_ = "Bad argument size: that";
        goto mw_err_label;
    }

    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != dim16_) {
        mw_err_txt_ = "Bad argument size: teid";        goto mw_err_label;
    }

    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != dim17_) {
        mw_err_txt_ = "Bad argument size: rhs";        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    in5_ = (int) mxWrapGetScalar(prhs[5], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_double(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxWrapGetArray_int(prhs[7], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in7_ = NULL;
    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != 0) {
        in8_ = mxWrapGetArray_double(prhs[8], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in8_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[8]++;
    rhsrecip(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_);
    plhs[0] = mxCreateDoubleMatrix(dim14_, dim15_, mxREAL);
    mxWrapCopy_double(plhs[0], in6_, dim14_*dim15_);
    plhs[1] = mxCreateDoubleMatrix(dim17_, 1, mxREAL);
    mxWrapCopy_double(plhs[1], in8_, dim17_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in6_)  mxFree(in6_);
    if (in7_)  mxFree(in7_);
    if (in8_)  mxFree(in8_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 159 ----
 * rhsrecip2nd(int nte, int np, int[10, nte] te2p2, double[3, np] p2, double[nte] conductivity, int nsource, inout double[3, nsource] that, int[nsource] teid, inout double[np] rhs);
 */
const char* stubids9_ = "rhsrecip2nd(i int, i int, i int[xx], i double[xx], i double[x], i int, io double[xx], i int[x], io double[x])";

void mexStub9(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nte        */
    int         in1_;    /* np         */
    int*        in2_ =0; /* te2p2      */
    double*     in3_ =0; /* p2         */
    double*     in4_ =0; /* conductivity */
    int         in5_;    /* nsource    */
    double*     in6_ =0; /* that       */
    int*        in7_ =0; /* teid       */
    double*     in8_ =0; /* rhs        */
    int         dim9_;   /* 10         */
    int         dim10_;   /* nte        */
    int         dim11_;   /* 3          */
    int         dim12_;   /* np         */
    int         dim13_;   /* nte        */
    int         dim14_;   /* 3          */
    int         dim15_;   /* nsource    */
    int         dim16_;   /* nsource    */
    int         dim17_;   /* np         */

    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);

    if (mxGetM(prhs[2]) != dim9_ ||
        mxGetN(prhs[2]) != dim10_) {
        mw_err_txt_ = "Bad argument size: te2p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim11_ ||
        mxGetN(prhs[3]) != dim12_) {
        mw_err_txt_ = "Bad argument size: p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim13_) {
        mw_err_txt_ = "Bad argument size: conductivity";        goto mw_err_label;
    }

    if (mxGetM(prhs[6]) != dim14_ ||
        mxGetN(prhs[6]) != dim15_) {
        mw_err_txt_ = "Bad argument size: that";
        goto mw_err_label;
    }

    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != dim16_) {
        mw_err_txt_ = "Bad argument size: teid";        goto mw_err_label;
    }

    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != dim17_) {
        mw_err_txt_ = "Bad argument size: rhs";        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    in5_ = (int) mxWrapGetScalar(prhs[5], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_double(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxWrapGetArray_int(prhs[7], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in7_ = NULL;
    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != 0) {
        in8_ = mxWrapGetArray_double(prhs[8], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in8_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[9]++;
    rhsrecip2nd(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_);
    plhs[0] = mxCreateDoubleMatrix(dim14_, dim15_, mxREAL);
    mxWrapCopy_double(plhs[0], in6_, dim14_*dim15_);
    plhs[1] = mxCreateDoubleMatrix(dim17_, 1, mxREAL);
    mxWrapCopy_double(plhs[1], in8_, dim17_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in6_)  mxFree(in6_);
    if (in7_)  mxFree(in7_);
    if (in8_)  mxFree(in8_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 161 ----
 * rhsrecip3rd(int nte, int np, int[20, nte] te2p2, double[3, np] p2, double[nte] conductivity, int nsource, inout double[3, nsource] that, int[nsource] teid, inout double[np] rhs);
 */
const char* stubids10_ = "rhsrecip3rd(i int, i int, i int[xx], i double[xx], i double[x], i int, io double[xx], i int[x], io double[x])";

void mexStub10(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nte        */
    int         in1_;    /* np         */
    int*        in2_ =0; /* te2p2      */
    double*     in3_ =0; /* p2         */
    double*     in4_ =0; /* conductivity */
    int         in5_;    /* nsource    */
    double*     in6_ =0; /* that       */
    int*        in7_ =0; /* teid       */
    double*     in8_ =0; /* rhs        */
    int         dim9_;   /* 20         */
    int         dim10_;   /* nte        */
    int         dim11_;   /* 3          */
    int         dim12_;   /* np         */
    int         dim13_;   /* nte        */
    int         dim14_;   /* 3          */
    int         dim15_;   /* nsource    */
    int         dim16_;   /* nsource    */
    int         dim17_;   /* np         */

    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);

    if (mxGetM(prhs[2]) != dim9_ ||
        mxGetN(prhs[2]) != dim10_) {
        mw_err_txt_ = "Bad argument size: te2p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim11_ ||
        mxGetN(prhs[3]) != dim12_) {
        mw_err_txt_ = "Bad argument size: p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim13_) {
        mw_err_txt_ = "Bad argument size: conductivity";        goto mw_err_label;
    }

    if (mxGetM(prhs[6]) != dim14_ ||
        mxGetN(prhs[6]) != dim15_) {
        mw_err_txt_ = "Bad argument size: that";
        goto mw_err_label;
    }

    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != dim16_) {
        mw_err_txt_ = "Bad argument size: teid";        goto mw_err_label;
    }

    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != dim17_) {
        mw_err_txt_ = "Bad argument size: rhs";        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    in5_ = (int) mxWrapGetScalar(prhs[5], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_double(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxWrapGetArray_int(prhs[7], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in7_ = NULL;
    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != 0) {
        in8_ = mxWrapGetArray_double(prhs[8], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in8_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[10]++;
    rhsrecip3rd(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_);
    plhs[0] = mxCreateDoubleMatrix(dim14_, dim15_, mxREAL);
    mxWrapCopy_double(plhs[0], in6_, dim14_*dim15_);
    plhs[1] = mxCreateDoubleMatrix(dim17_, 1, mxREAL);
    mxWrapCopy_double(plhs[1], in8_, dim17_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in6_)  mxFree(in6_);
    if (in7_)  mxFree(in7_);
    if (in8_)  mxFree(in8_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 179 ----
 * conductioncurr(int[4, nte] te2p2, int nte, double[3, np] p2, double[nte] conductivity, double[np] x, inout double[3, ntens] rv, inout double[3, ntens] jv);
 */
const char* stubids11_ = "conductioncurr(i int[xx], i int, i double[xx], i double[x], i double[x], io double[xx], io double[xx])";

void mexStub11(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int*        in0_ =0; /* te2p2      */
    int         in1_;    /* nte        */
    double*     in2_ =0; /* p2         */
    double*     in3_ =0; /* conductivity */
    double*     in4_ =0; /* x          */
    double*     in5_ =0; /* rv         */
    double*     in6_ =0; /* jv         */
    int         dim7_;   /* 4          */
    int         dim8_;   /* nte        */
    int         dim9_;   /* 3          */
    int         dim10_;   /* np         */
    int         dim11_;   /* nte        */
    int         dim12_;   /* np         */
    int         dim13_;   /* 3          */
    int         dim14_;   /* ntens      */
    int         dim15_;   /* 3          */
    int         dim16_;   /* ntens      */

    dim7_ = (int) mxWrapGetScalar(prhs[7], &mw_err_txt_);
    dim8_ = (int) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);

    if (mxGetM(prhs[0]) != dim7_ ||
        mxGetN(prhs[0]) != dim8_) {
        mw_err_txt_ = "Bad argument size: te2p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[2]) != dim9_ ||
        mxGetN(prhs[2]) != dim10_) {
        mw_err_txt_ = "Bad argument size: p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != dim11_) {
        mw_err_txt_ = "Bad argument size: conductivity";        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim12_) {
        mw_err_txt_ = "Bad argument size: x";        goto mw_err_label;
    }

    if (mxGetM(prhs[5]) != dim13_ ||
        mxGetN(prhs[5]) != dim14_) {
        mw_err_txt_ = "Bad argument size: rv";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6]) != dim15_ ||
        mxGetN(prhs[6]) != dim16_) {
        mw_err_txt_ = "Bad argument size: jv";
        goto mw_err_label;
    }

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != 0) {
        in0_ = mxWrapGetArray_int(prhs[0], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in0_ = NULL;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxGetPr(prhs[2]);
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxWrapGetArray_double(prhs[5], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_double(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[11]++;
    conductioncurr(in0_, in1_, in2_, in3_, in4_, in5_, in6_);
    plhs[0] = mxCreateDoubleMatrix(dim13_, dim14_, mxREAL);
    mxWrapCopy_double(plhs[0], in5_, dim13_*dim14_);
    plhs[1] = mxCreateDoubleMatrix(dim15_, dim16_, mxREAL);
    mxWrapCopy_double(plhs[1], in6_, dim15_*dim16_);

mw_err_label:
    if (in0_)  mxFree(in0_);
    if (in5_)  mxFree(in5_);
    if (in6_)  mxFree(in6_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 181 ----
 * conductioncurr2nd(int[10, nte] te2p2, int nte, double[3, np] p2, double[nte] conductivity, double[np] x, inout double[3, ntens] rv, inout double[3, ntens] jv);
 */
const char* stubids12_ = "conductioncurr2nd(i int[xx], i int, i double[xx], i double[x], i double[x], io double[xx], io double[xx])";

void mexStub12(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int*        in0_ =0; /* te2p2      */
    int         in1_;    /* nte        */
    double*     in2_ =0; /* p2         */
    double*     in3_ =0; /* conductivity */
    double*     in4_ =0; /* x          */
    double*     in5_ =0; /* rv         */
    double*     in6_ =0; /* jv         */
    int         dim7_;   /* 10         */
    int         dim8_;   /* nte        */
    int         dim9_;   /* 3          */
    int         dim10_;   /* np         */
    int         dim11_;   /* nte        */
    int         dim12_;   /* np         */
    int         dim13_;   /* 3          */
    int         dim14_;   /* ntens      */
    int         dim15_;   /* 3          */
    int         dim16_;   /* ntens      */

    dim7_ = (int) mxWrapGetScalar(prhs[7], &mw_err_txt_);
    dim8_ = (int) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);

    if (mxGetM(prhs[0]) != dim7_ ||
        mxGetN(prhs[0]) != dim8_) {
        mw_err_txt_ = "Bad argument size: te2p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[2]) != dim9_ ||
        mxGetN(prhs[2]) != dim10_) {
        mw_err_txt_ = "Bad argument size: p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != dim11_) {
        mw_err_txt_ = "Bad argument size: conductivity";        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim12_) {
        mw_err_txt_ = "Bad argument size: x";        goto mw_err_label;
    }

    if (mxGetM(prhs[5]) != dim13_ ||
        mxGetN(prhs[5]) != dim14_) {
        mw_err_txt_ = "Bad argument size: rv";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6]) != dim15_ ||
        mxGetN(prhs[6]) != dim16_) {
        mw_err_txt_ = "Bad argument size: jv";
        goto mw_err_label;
    }

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != 0) {
        in0_ = mxWrapGetArray_int(prhs[0], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in0_ = NULL;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxGetPr(prhs[2]);
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxWrapGetArray_double(prhs[5], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_double(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[12]++;
    conductioncurr2nd(in0_, in1_, in2_, in3_, in4_, in5_, in6_);
    plhs[0] = mxCreateDoubleMatrix(dim13_, dim14_, mxREAL);
    mxWrapCopy_double(plhs[0], in5_, dim13_*dim14_);
    plhs[1] = mxCreateDoubleMatrix(dim15_, dim16_, mxREAL);
    mxWrapCopy_double(plhs[1], in6_, dim15_*dim16_);

mw_err_label:
    if (in0_)  mxFree(in0_);
    if (in5_)  mxFree(in5_);
    if (in6_)  mxFree(in6_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 183 ----
 * conductioncurr3rd(int[20, nte] te2p2, int nte, double[3, np] p2, double[nte] conductivity, double[np] x, inout double[3, ntens] rv, inout double[3, ntens] jv);
 */
const char* stubids13_ = "conductioncurr3rd(i int[xx], i int, i double[xx], i double[x], i double[x], io double[xx], io double[xx])";

void mexStub13(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int*        in0_ =0; /* te2p2      */
    int         in1_;    /* nte        */
    double*     in2_ =0; /* p2         */
    double*     in3_ =0; /* conductivity */
    double*     in4_ =0; /* x          */
    double*     in5_ =0; /* rv         */
    double*     in6_ =0; /* jv         */
    int         dim7_;   /* 20         */
    int         dim8_;   /* nte        */
    int         dim9_;   /* 3          */
    int         dim10_;   /* np         */
    int         dim11_;   /* nte        */
    int         dim12_;   /* np         */
    int         dim13_;   /* 3          */
    int         dim14_;   /* ntens      */
    int         dim15_;   /* 3          */
    int         dim16_;   /* ntens      */

    dim7_ = (int) mxWrapGetScalar(prhs[7], &mw_err_txt_);
    dim8_ = (int) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);

    if (mxGetM(prhs[0]) != dim7_ ||
        mxGetN(prhs[0]) != dim8_) {
        mw_err_txt_ = "Bad argument size: te2p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[2]) != dim9_ ||
        mxGetN(prhs[2]) != dim10_) {
        mw_err_txt_ = "Bad argument size: p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != dim11_) {
        mw_err_txt_ = "Bad argument size: conductivity";        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim12_) {
        mw_err_txt_ = "Bad argument size: x";        goto mw_err_label;
    }

    if (mxGetM(prhs[5]) != dim13_ ||
        mxGetN(prhs[5]) != dim14_) {
        mw_err_txt_ = "Bad argument size: rv";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6]) != dim15_ ||
        mxGetN(prhs[6]) != dim16_) {
        mw_err_txt_ = "Bad argument size: jv";
        goto mw_err_label;
    }

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != 0) {
        in0_ = mxWrapGetArray_int(prhs[0], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in0_ = NULL;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxGetPr(prhs[2]);
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxWrapGetArray_double(prhs[5], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_double(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[13]++;
    conductioncurr3rd(in0_, in1_, in2_, in3_, in4_, in5_, in6_);
    plhs[0] = mxCreateDoubleMatrix(dim13_, dim14_, mxREAL);
    mxWrapCopy_double(plhs[0], in5_, dim13_*dim14_);
    plhs[1] = mxCreateDoubleMatrix(dim15_, dim16_, mxREAL);
    mxWrapCopy_double(plhs[1], in6_, dim15_*dim16_);

mw_err_label:
    if (in0_)  mxFree(in0_);
    if (in5_)  mxFree(in5_);
    if (in6_)  mxFree(in6_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 337 ----
 * int ne = generateedge_ct(int nte, int[4, nte] te2p);
 * Also at demo.mw: 342
 */
const char* stubids15_ = "o int = generateedge_ct(i int, i int[xx])";

void mexStub15(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nte        */
    int*        in1_ =0; /* te2p       */
    int         out0_;   /* ne         */
    int         dim2_;   /* 4          */
    int         dim3_;   /* nte        */

    dim2_ = (int) mxWrapGetScalar(prhs[2], &mw_err_txt_);
    dim3_ = (int) mxWrapGetScalar(prhs[3], &mw_err_txt_);

    if (mxGetM(prhs[1]) != dim2_ ||
        mxGetN(prhs[1]) != dim3_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[1])*mxGetN(prhs[1]) != 0) {
        in1_ = mxWrapGetArray_int(prhs[1], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in1_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[15]++;
    out0_ = generateedge_ct(in0_, in1_);
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(plhs[0]) = out0_;

mw_err_label:
    if (in1_)  mxFree(in1_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 340 ----
 * generatemesh2nd(int nte, int[4, nte] te2p, double[3, np] p, int np, inout int[10, nte] te2p2, inout double[3, npne] p2);
 */
const char* stubids16_ = "generatemesh2nd(i int, i int[xx], i double[xx], i int, io int[xx], io double[xx])";

void mexStub16(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nte        */
    int*        in1_ =0; /* te2p       */
    double*     in2_ =0; /* p          */
    int         in3_;    /* np         */
    int*        in4_ =0; /* te2p2      */
    double*     in5_ =0; /* p2         */
    int         dim6_;   /* 4          */
    int         dim7_;   /* nte        */
    int         dim8_;   /* 3          */
    int         dim9_;   /* np         */
    int         dim10_;   /* 10         */
    int         dim11_;   /* nte        */
    int         dim12_;   /* 3          */
    int         dim13_;   /* npne       */

    dim6_ = (int) mxWrapGetScalar(prhs[6], &mw_err_txt_);
    dim7_ = (int) mxWrapGetScalar(prhs[7], &mw_err_txt_);
    dim8_ = (int) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);

    if (mxGetM(prhs[1]) != dim6_ ||
        mxGetN(prhs[1]) != dim7_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[2]) != dim8_ ||
        mxGetN(prhs[2]) != dim9_) {
        mw_err_txt_ = "Bad argument size: p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4]) != dim10_ ||
        mxGetN(prhs[4]) != dim11_) {
        mw_err_txt_ = "Bad argument size: te2p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[5]) != dim12_ ||
        mxGetN(prhs[5]) != dim13_) {
        mw_err_txt_ = "Bad argument size: p2";
        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[1])*mxGetN(prhs[1]) != 0) {
        in1_ = mxWrapGetArray_int(prhs[1], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in1_ = NULL;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxGetPr(prhs[2]);
    } else
        in2_ = NULL;
    in3_ = (int) mxWrapGetScalar(prhs[3], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxWrapGetArray_int(prhs[4], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in4_ = NULL;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxWrapGetArray_double(prhs[5], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in5_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[16]++;
    generatemesh2nd(in0_, in1_, in2_, in3_, in4_, in5_);
    plhs[0] = mxCreateDoubleMatrix(dim10_, dim11_, mxREAL);
    mxWrapCopy_int(plhs[0], in4_, dim10_*dim11_);
    plhs[1] = mxCreateDoubleMatrix(dim12_, dim13_, mxREAL);
    mxWrapCopy_double(plhs[1], in5_, dim12_*dim13_);

mw_err_label:
    if (in1_)  mxFree(in1_);
    if (in4_)  mxFree(in4_);
    if (in5_)  mxFree(in5_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 343 ----
 * int nf = generateface_ct(int nte, int[4, nte] te2p);
 */
const char* stubids18_ = "o int = generateface_ct(i int, i int[xx])";

void mexStub18(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nte        */
    int*        in1_ =0; /* te2p       */
    int         out0_;   /* nf         */
    int         dim2_;   /* 4          */
    int         dim3_;   /* nte        */

    dim2_ = (int) mxWrapGetScalar(prhs[2], &mw_err_txt_);
    dim3_ = (int) mxWrapGetScalar(prhs[3], &mw_err_txt_);

    if (mxGetM(prhs[1]) != dim2_ ||
        mxGetN(prhs[1]) != dim3_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[1])*mxGetN(prhs[1]) != 0) {
        in1_ = mxWrapGetArray_int(prhs[1], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in1_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[18]++;
    out0_ = generateface_ct(in0_, in1_);
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(plhs[0]) = out0_;

mw_err_label:
    if (in1_)  mxFree(in1_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 346 ----
 * generatemesh3rd(int nte, int[4, nte] te2p, double[3, np] p, int np, int ne, inout int[20, nte] te2p2, inout double[3, npne] p2);
 */
const char* stubids19_ = "generatemesh3rd(i int, i int[xx], i double[xx], i int, i int, io int[xx], io double[xx])";

void mexStub19(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nte        */
    int*        in1_ =0; /* te2p       */
    double*     in2_ =0; /* p          */
    int         in3_;    /* np         */
    int         in4_;    /* ne         */
    int*        in5_ =0; /* te2p2      */
    double*     in6_ =0; /* p2         */
    int         dim7_;   /* 4          */
    int         dim8_;   /* nte        */
    int         dim9_;   /* 3          */
    int         dim10_;   /* np         */
    int         dim11_;   /* 20         */
    int         dim12_;   /* nte        */
    int         dim13_;   /* 3          */
    int         dim14_;   /* npne       */

    dim7_ = (int) mxWrapGetScalar(prhs[7], &mw_err_txt_);
    dim8_ = (int) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);

    if (mxGetM(prhs[1]) != dim7_ ||
        mxGetN(prhs[1]) != dim8_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[2]) != dim9_ ||
        mxGetN(prhs[2]) != dim10_) {
        mw_err_txt_ = "Bad argument size: p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[5]) != dim11_ ||
        mxGetN(prhs[5]) != dim12_) {
        mw_err_txt_ = "Bad argument size: te2p2";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6]) != dim13_ ||
        mxGetN(prhs[6]) != dim14_) {
        mw_err_txt_ = "Bad argument size: p2";
        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[1])*mxGetN(prhs[1]) != 0) {
        in1_ = mxWrapGetArray_int(prhs[1], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in1_ = NULL;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxGetPr(prhs[2]);
    } else
        in2_ = NULL;
    in3_ = (int) mxWrapGetScalar(prhs[3], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    in4_ = (int) mxWrapGetScalar(prhs[4], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxWrapGetArray_int(prhs[5], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_double(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[19]++;
    generatemesh3rd(in0_, in1_, in2_, in3_, in4_, in5_, in6_);
    plhs[0] = mxCreateDoubleMatrix(dim11_, dim12_, mxREAL);
    mxWrapCopy_int(plhs[0], in5_, dim11_*dim12_);
    plhs[1] = mxCreateDoubleMatrix(dim13_, dim14_, mxREAL);
    mxWrapCopy_double(plhs[1], in6_, dim13_*dim14_);

mw_err_label:
    if (in1_)  mxFree(in1_);
    if (in5_)  mxFree(in5_);
    if (in6_)  mxFree(in6_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 358 ----
 * assemblemattensor(int nte, int np, int[4, nte] te2p, double[3, np] p, double[9, nte] reg, inout int[3, nte16] crvindex, inout double[nte16] val);
 */
const char* stubids20_ = "assemblemattensor(i int, i int, i int[xx], i double[xx], i double[xx], io int[xx], io double[x])";

void mexStub20(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nte        */
    int         in1_;    /* np         */
    int*        in2_ =0; /* te2p       */
    double*     in3_ =0; /* p          */
    double*     in4_ =0; /* reg        */
    int*        in5_ =0; /* crvindex   */
    double*     in6_ =0; /* val        */
    int         dim7_;   /* 4          */
    int         dim8_;   /* nte        */
    int         dim9_;   /* 3          */
    int         dim10_;   /* np         */
    int         dim11_;   /* 9          */
    int         dim12_;   /* nte        */
    int         dim13_;   /* 3          */
    int         dim14_;   /* nte16      */
    int         dim15_;   /* nte16      */

    dim7_ = (int) mxWrapGetScalar(prhs[7], &mw_err_txt_);
    dim8_ = (int) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);

    if (mxGetM(prhs[2]) != dim7_ ||
        mxGetN(prhs[2]) != dim8_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim9_ ||
        mxGetN(prhs[3]) != dim10_) {
        mw_err_txt_ = "Bad argument size: p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4]) != dim11_ ||
        mxGetN(prhs[4]) != dim12_) {
        mw_err_txt_ = "Bad argument size: reg";
        goto mw_err_label;
    }

    if (mxGetM(prhs[5]) != dim13_ ||
        mxGetN(prhs[5]) != dim14_) {
        mw_err_txt_ = "Bad argument size: crvindex";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != dim15_) {
        mw_err_txt_ = "Bad argument size: val";        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxWrapGetArray_int(prhs[5], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_double(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[20]++;
    assemblemattensor(in0_, in1_, in2_, in3_, in4_, in5_, in6_);
    plhs[0] = mxCreateDoubleMatrix(dim13_, dim14_, mxREAL);
    mxWrapCopy_int(plhs[0], in5_, dim13_*dim14_);
    plhs[1] = mxCreateDoubleMatrix(dim15_, 1, mxREAL);
    mxWrapCopy_double(plhs[1], in6_, dim15_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in5_)  mxFree(in5_);
    if (in6_)  mxFree(in6_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 363 ----
 * assemblemat2ndtensor(int nte, int np, int[10, nte] te2p, double[3, np] p, double[9, nte] reg, inout int[3, nte100] crvindex, inout double[nte100] val);
 */
const char* stubids21_ = "assemblemat2ndtensor(i int, i int, i int[xx], i double[xx], i double[xx], io int[xx], io double[x])";

void mexStub21(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nte        */
    int         in1_;    /* np         */
    int*        in2_ =0; /* te2p       */
    double*     in3_ =0; /* p          */
    double*     in4_ =0; /* reg        */
    int*        in5_ =0; /* crvindex   */
    double*     in6_ =0; /* val        */
    int         dim7_;   /* 10         */
    int         dim8_;   /* nte        */
    int         dim9_;   /* 3          */
    int         dim10_;   /* np         */
    int         dim11_;   /* 9          */
    int         dim12_;   /* nte        */
    int         dim13_;   /* 3          */
    int         dim14_;   /* nte100     */
    int         dim15_;   /* nte100     */

    dim7_ = (int) mxWrapGetScalar(prhs[7], &mw_err_txt_);
    dim8_ = (int) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);

    if (mxGetM(prhs[2]) != dim7_ ||
        mxGetN(prhs[2]) != dim8_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim9_ ||
        mxGetN(prhs[3]) != dim10_) {
        mw_err_txt_ = "Bad argument size: p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4]) != dim11_ ||
        mxGetN(prhs[4]) != dim12_) {
        mw_err_txt_ = "Bad argument size: reg";
        goto mw_err_label;
    }

    if (mxGetM(prhs[5]) != dim13_ ||
        mxGetN(prhs[5]) != dim14_) {
        mw_err_txt_ = "Bad argument size: crvindex";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != dim15_) {
        mw_err_txt_ = "Bad argument size: val";        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxWrapGetArray_int(prhs[5], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_double(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[21]++;
    assemblemat2ndtensor(in0_, in1_, in2_, in3_, in4_, in5_, in6_);
    plhs[0] = mxCreateDoubleMatrix(dim13_, dim14_, mxREAL);
    mxWrapCopy_int(plhs[0], in5_, dim13_*dim14_);
    plhs[1] = mxCreateDoubleMatrix(dim15_, 1, mxREAL);
    mxWrapCopy_double(plhs[1], in6_, dim15_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in5_)  mxFree(in5_);
    if (in6_)  mxFree(in6_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 368 ----
 * assemblemat3rdtensor(int nte, int np, int[20, nte] te2p, double[3, np] p, double[9, nte] reg, inout int[3, nte400] crvindex, inout double[nte400] val);
 */
const char* stubids22_ = "assemblemat3rdtensor(i int, i int, i int[xx], i double[xx], i double[xx], io int[xx], io double[x])";

void mexStub22(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nte        */
    int         in1_;    /* np         */
    int*        in2_ =0; /* te2p       */
    double*     in3_ =0; /* p          */
    double*     in4_ =0; /* reg        */
    int*        in5_ =0; /* crvindex   */
    double*     in6_ =0; /* val        */
    int         dim7_;   /* 20         */
    int         dim8_;   /* nte        */
    int         dim9_;   /* 3          */
    int         dim10_;   /* np         */
    int         dim11_;   /* 9          */
    int         dim12_;   /* nte        */
    int         dim13_;   /* 3          */
    int         dim14_;   /* nte400     */
    int         dim15_;   /* nte400     */

    dim7_ = (int) mxWrapGetScalar(prhs[7], &mw_err_txt_);
    dim8_ = (int) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);

    if (mxGetM(prhs[2]) != dim7_ ||
        mxGetN(prhs[2]) != dim8_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim9_ ||
        mxGetN(prhs[3]) != dim10_) {
        mw_err_txt_ = "Bad argument size: p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4]) != dim11_ ||
        mxGetN(prhs[4]) != dim12_) {
        mw_err_txt_ = "Bad argument size: reg";
        goto mw_err_label;
    }

    if (mxGetM(prhs[5]) != dim13_ ||
        mxGetN(prhs[5]) != dim14_) {
        mw_err_txt_ = "Bad argument size: crvindex";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != dim15_) {
        mw_err_txt_ = "Bad argument size: val";        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxWrapGetArray_int(prhs[5], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_double(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[22]++;
    assemblemat3rdtensor(in0_, in1_, in2_, in3_, in4_, in5_, in6_);
    plhs[0] = mxCreateDoubleMatrix(dim13_, dim14_, mxREAL);
    mxWrapCopy_int(plhs[0], in5_, dim13_*dim14_);
    plhs[1] = mxCreateDoubleMatrix(dim15_, 1, mxREAL);
    mxWrapCopy_double(plhs[1], in6_, dim15_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in5_)  mxFree(in5_);
    if (in6_)  mxFree(in6_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 382 ----
 * assemblemat(int nte, int np, int[4, nte] te2p, double[3, np] p, double[nte] reg, inout int[3, nte16] crvindex, inout double[nte16] val);
 */
const char* stubids23_ = "assemblemat(i int, i int, i int[xx], i double[xx], i double[x], io int[xx], io double[x])";

void mexStub23(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nte        */
    int         in1_;    /* np         */
    int*        in2_ =0; /* te2p       */
    double*     in3_ =0; /* p          */
    double*     in4_ =0; /* reg        */
    int*        in5_ =0; /* crvindex   */
    double*     in6_ =0; /* val        */
    int         dim7_;   /* 4          */
    int         dim8_;   /* nte        */
    int         dim9_;   /* 3          */
    int         dim10_;   /* np         */
    int         dim11_;   /* nte        */
    int         dim12_;   /* 3          */
    int         dim13_;   /* nte16      */
    int         dim14_;   /* nte16      */

    dim7_ = (int) mxWrapGetScalar(prhs[7], &mw_err_txt_);
    dim8_ = (int) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);

    if (mxGetM(prhs[2]) != dim7_ ||
        mxGetN(prhs[2]) != dim8_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim9_ ||
        mxGetN(prhs[3]) != dim10_) {
        mw_err_txt_ = "Bad argument size: p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim11_) {
        mw_err_txt_ = "Bad argument size: reg";        goto mw_err_label;
    }

    if (mxGetM(prhs[5]) != dim12_ ||
        mxGetN(prhs[5]) != dim13_) {
        mw_err_txt_ = "Bad argument size: crvindex";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != dim14_) {
        mw_err_txt_ = "Bad argument size: val";        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxWrapGetArray_int(prhs[5], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_double(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[23]++;
    assemblemat(in0_, in1_, in2_, in3_, in4_, in5_, in6_);
    plhs[0] = mxCreateDoubleMatrix(dim12_, dim13_, mxREAL);
    mxWrapCopy_int(plhs[0], in5_, dim12_*dim13_);
    plhs[1] = mxCreateDoubleMatrix(dim14_, 1, mxREAL);
    mxWrapCopy_double(plhs[1], in6_, dim14_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in5_)  mxFree(in5_);
    if (in6_)  mxFree(in6_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 387 ----
 * assemblemat2nd(int nte, int np, int[10, nte] te2p, double[3, np] p, double[nte] reg, inout int[3, nte100] crvindex, inout double[nte100] val);
 */
const char* stubids24_ = "assemblemat2nd(i int, i int, i int[xx], i double[xx], i double[x], io int[xx], io double[x])";

void mexStub24(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nte        */
    int         in1_;    /* np         */
    int*        in2_ =0; /* te2p       */
    double*     in3_ =0; /* p          */
    double*     in4_ =0; /* reg        */
    int*        in5_ =0; /* crvindex   */
    double*     in6_ =0; /* val        */
    int         dim7_;   /* 10         */
    int         dim8_;   /* nte        */
    int         dim9_;   /* 3          */
    int         dim10_;   /* np         */
    int         dim11_;   /* nte        */
    int         dim12_;   /* 3          */
    int         dim13_;   /* nte100     */
    int         dim14_;   /* nte100     */

    dim7_ = (int) mxWrapGetScalar(prhs[7], &mw_err_txt_);
    dim8_ = (int) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);

    if (mxGetM(prhs[2]) != dim7_ ||
        mxGetN(prhs[2]) != dim8_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim9_ ||
        mxGetN(prhs[3]) != dim10_) {
        mw_err_txt_ = "Bad argument size: p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim11_) {
        mw_err_txt_ = "Bad argument size: reg";        goto mw_err_label;
    }

    if (mxGetM(prhs[5]) != dim12_ ||
        mxGetN(prhs[5]) != dim13_) {
        mw_err_txt_ = "Bad argument size: crvindex";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != dim14_) {
        mw_err_txt_ = "Bad argument size: val";        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxWrapGetArray_int(prhs[5], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_double(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[24]++;
    assemblemat2nd(in0_, in1_, in2_, in3_, in4_, in5_, in6_);
    plhs[0] = mxCreateDoubleMatrix(dim12_, dim13_, mxREAL);
    mxWrapCopy_int(plhs[0], in5_, dim12_*dim13_);
    plhs[1] = mxCreateDoubleMatrix(dim14_, 1, mxREAL);
    mxWrapCopy_double(plhs[1], in6_, dim14_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in5_)  mxFree(in5_);
    if (in6_)  mxFree(in6_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 392 ----
 * assemblemat3rd(int nte, int np, int[20, nte] te2p, double[3, np] p, double[nte] reg, inout int[3, nte400] crvindex, inout double[nte400] val);
 */
const char* stubids25_ = "assemblemat3rd(i int, i int, i int[xx], i double[xx], i double[x], io int[xx], io double[x])";

void mexStub25(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nte        */
    int         in1_;    /* np         */
    int*        in2_ =0; /* te2p       */
    double*     in3_ =0; /* p          */
    double*     in4_ =0; /* reg        */
    int*        in5_ =0; /* crvindex   */
    double*     in6_ =0; /* val        */
    int         dim7_;   /* 20         */
    int         dim8_;   /* nte        */
    int         dim9_;   /* 3          */
    int         dim10_;   /* np         */
    int         dim11_;   /* nte        */
    int         dim12_;   /* 3          */
    int         dim13_;   /* nte400     */
    int         dim14_;   /* nte400     */

    dim7_ = (int) mxWrapGetScalar(prhs[7], &mw_err_txt_);
    dim8_ = (int) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);

    if (mxGetM(prhs[2]) != dim7_ ||
        mxGetN(prhs[2]) != dim8_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim9_ ||
        mxGetN(prhs[3]) != dim10_) {
        mw_err_txt_ = "Bad argument size: p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim11_) {
        mw_err_txt_ = "Bad argument size: reg";        goto mw_err_label;
    }

    if (mxGetM(prhs[5]) != dim12_ ||
        mxGetN(prhs[5]) != dim13_) {
        mw_err_txt_ = "Bad argument size: crvindex";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != dim14_) {
        mw_err_txt_ = "Bad argument size: val";        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxWrapGetArray_int(prhs[5], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_double(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[25]++;
    assemblemat3rd(in0_, in1_, in2_, in3_, in4_, in5_, in6_);
    plhs[0] = mxCreateDoubleMatrix(dim12_, dim13_, mxREAL);
    mxWrapCopy_int(plhs[0], in5_, dim12_*dim13_);
    plhs[1] = mxCreateDoubleMatrix(dim14_, 1, mxREAL);
    mxWrapCopy_double(plhs[1], in6_, dim14_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in5_)  mxFree(in5_);
    if (in6_)  mxFree(in6_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 406 ----
 * rhsfunc(int nte, int np, int[4, nte] te2p, double[3, np] p, double[nte] reg, int nsource, double[3, nsource] rs, double[3, nsource] js, inout double[np] rhs, double iprec);
 */
const char* stubids26_ = "rhsfunc(i int, i int, i int[xx], i double[xx], i double[x], i int, i double[xx], i double[xx], io double[x], i double)";

void mexStub26(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nte        */
    int         in1_;    /* np         */
    int*        in2_ =0; /* te2p       */
    double*     in3_ =0; /* p          */
    double*     in4_ =0; /* reg        */
    int         in5_;    /* nsource    */
    double*     in6_ =0; /* rs         */
    double*     in7_ =0; /* js         */
    double*     in8_ =0; /* rhs        */
    double      in9_;    /* iprec      */
    int         dim10_;   /* 4          */
    int         dim11_;   /* nte        */
    int         dim12_;   /* 3          */
    int         dim13_;   /* np         */
    int         dim14_;   /* nte        */
    int         dim15_;   /* 3          */
    int         dim16_;   /* nsource    */
    int         dim17_;   /* 3          */
    int         dim18_;   /* nsource    */
    int         dim19_;   /* np         */

    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);
    dim18_ = (int) mxWrapGetScalar(prhs[18], &mw_err_txt_);
    dim19_ = (int) mxWrapGetScalar(prhs[19], &mw_err_txt_);

    if (mxGetM(prhs[2]) != dim10_ ||
        mxGetN(prhs[2]) != dim11_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim12_ ||
        mxGetN(prhs[3]) != dim13_) {
        mw_err_txt_ = "Bad argument size: p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim14_) {
        mw_err_txt_ = "Bad argument size: reg";        goto mw_err_label;
    }

    if (mxGetM(prhs[6]) != dim15_ ||
        mxGetN(prhs[6]) != dim16_) {
        mw_err_txt_ = "Bad argument size: rs";
        goto mw_err_label;
    }

    if (mxGetM(prhs[7]) != dim17_ ||
        mxGetN(prhs[7]) != dim18_) {
        mw_err_txt_ = "Bad argument size: js";
        goto mw_err_label;
    }

    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != dim19_) {
        mw_err_txt_ = "Bad argument size: rhs";        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    in5_ = (int) mxWrapGetScalar(prhs[5], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxGetPr(prhs[6]);
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxGetPr(prhs[7]);
    } else
        in7_ = NULL;
    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != 0) {
        in8_ = mxWrapGetArray_double(prhs[8], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in8_ = NULL;
    in9_ = (double) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mexprofrecord_)
        mexprofrecord_[26]++;
    rhsfunc(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_, in9_);
    plhs[0] = mxCreateDoubleMatrix(dim19_, 1, mxREAL);
    mxWrapCopy_double(plhs[0], in8_, dim19_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in8_)  mxFree(in8_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 408 ----
 * rhsfunc2nd(int nte, int np, int[10, nte] te2p, double[3, np] p, double[nte] reg, int nsource, double[3, nsource] rs, double[3, nsource] js, inout double[np] rhs, double iprec);
 */
const char* stubids27_ = "rhsfunc2nd(i int, i int, i int[xx], i double[xx], i double[x], i int, i double[xx], i double[xx], io double[x], i double)";

void mexStub27(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nte        */
    int         in1_;    /* np         */
    int*        in2_ =0; /* te2p       */
    double*     in3_ =0; /* p          */
    double*     in4_ =0; /* reg        */
    int         in5_;    /* nsource    */
    double*     in6_ =0; /* rs         */
    double*     in7_ =0; /* js         */
    double*     in8_ =0; /* rhs        */
    double      in9_;    /* iprec      */
    int         dim10_;   /* 10         */
    int         dim11_;   /* nte        */
    int         dim12_;   /* 3          */
    int         dim13_;   /* np         */
    int         dim14_;   /* nte        */
    int         dim15_;   /* 3          */
    int         dim16_;   /* nsource    */
    int         dim17_;   /* 3          */
    int         dim18_;   /* nsource    */
    int         dim19_;   /* np         */

    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);
    dim18_ = (int) mxWrapGetScalar(prhs[18], &mw_err_txt_);
    dim19_ = (int) mxWrapGetScalar(prhs[19], &mw_err_txt_);

    if (mxGetM(prhs[2]) != dim10_ ||
        mxGetN(prhs[2]) != dim11_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim12_ ||
        mxGetN(prhs[3]) != dim13_) {
        mw_err_txt_ = "Bad argument size: p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim14_) {
        mw_err_txt_ = "Bad argument size: reg";        goto mw_err_label;
    }

    if (mxGetM(prhs[6]) != dim15_ ||
        mxGetN(prhs[6]) != dim16_) {
        mw_err_txt_ = "Bad argument size: rs";
        goto mw_err_label;
    }

    if (mxGetM(prhs[7]) != dim17_ ||
        mxGetN(prhs[7]) != dim18_) {
        mw_err_txt_ = "Bad argument size: js";
        goto mw_err_label;
    }

    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != dim19_) {
        mw_err_txt_ = "Bad argument size: rhs";        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    in5_ = (int) mxWrapGetScalar(prhs[5], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxGetPr(prhs[6]);
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxGetPr(prhs[7]);
    } else
        in7_ = NULL;
    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != 0) {
        in8_ = mxWrapGetArray_double(prhs[8], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in8_ = NULL;
    in9_ = (double) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mexprofrecord_)
        mexprofrecord_[27]++;
    rhsfunc2nd(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_, in9_);
    plhs[0] = mxCreateDoubleMatrix(dim19_, 1, mxREAL);
    mxWrapCopy_double(plhs[0], in8_, dim19_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in8_)  mxFree(in8_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 410 ----
 * rhsfunc3rd(int nte, int np, int[20, nte] te2p, double[3, np] p, double[nte] reg, int nsource, double[3, nsource] rs, double[3, nsource] js, inout double[np] rhs, double iprec);
 */
const char* stubids28_ = "rhsfunc3rd(i int, i int, i int[xx], i double[xx], i double[x], i int, i double[xx], i double[xx], io double[x], i double)";

void mexStub28(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nte        */
    int         in1_;    /* np         */
    int*        in2_ =0; /* te2p       */
    double*     in3_ =0; /* p          */
    double*     in4_ =0; /* reg        */
    int         in5_;    /* nsource    */
    double*     in6_ =0; /* rs         */
    double*     in7_ =0; /* js         */
    double*     in8_ =0; /* rhs        */
    double      in9_;    /* iprec      */
    int         dim10_;   /* 20         */
    int         dim11_;   /* nte        */
    int         dim12_;   /* 3          */
    int         dim13_;   /* np         */
    int         dim14_;   /* nte        */
    int         dim15_;   /* 3          */
    int         dim16_;   /* nsource    */
    int         dim17_;   /* 3          */
    int         dim18_;   /* nsource    */
    int         dim19_;   /* np         */

    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);
    dim18_ = (int) mxWrapGetScalar(prhs[18], &mw_err_txt_);
    dim19_ = (int) mxWrapGetScalar(prhs[19], &mw_err_txt_);

    if (mxGetM(prhs[2]) != dim10_ ||
        mxGetN(prhs[2]) != dim11_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim12_ ||
        mxGetN(prhs[3]) != dim13_) {
        mw_err_txt_ = "Bad argument size: p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim14_) {
        mw_err_txt_ = "Bad argument size: reg";        goto mw_err_label;
    }

    if (mxGetM(prhs[6]) != dim15_ ||
        mxGetN(prhs[6]) != dim16_) {
        mw_err_txt_ = "Bad argument size: rs";
        goto mw_err_label;
    }

    if (mxGetM(prhs[7]) != dim17_ ||
        mxGetN(prhs[7]) != dim18_) {
        mw_err_txt_ = "Bad argument size: js";
        goto mw_err_label;
    }

    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != dim19_) {
        mw_err_txt_ = "Bad argument size: rhs";        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    in5_ = (int) mxWrapGetScalar(prhs[5], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxGetPr(prhs[6]);
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxGetPr(prhs[7]);
    } else
        in7_ = NULL;
    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != 0) {
        in8_ = mxWrapGetArray_double(prhs[8], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in8_ = NULL;
    in9_ = (double) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mexprofrecord_)
        mexprofrecord_[28]++;
    rhsfunc3rd(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_, in9_);
    plhs[0] = mxCreateDoubleMatrix(dim19_, 1, mxREAL);
    mxWrapCopy_double(plhs[0], in8_, dim19_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in8_)  mxFree(in8_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 424 ----
 * rhsfuncks(int nte, int np, int[4, nte] te2p, double[3, np] p, double[nte] reg, int nsource, double[3, nsource] rs, double[3, nsource] js, inout double[np] rhs, double iprec);
 */
const char* stubids29_ = "rhsfuncks(i int, i int, i int[xx], i double[xx], i double[x], i int, i double[xx], i double[xx], io double[x], i double)";

void mexStub29(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nte        */
    int         in1_;    /* np         */
    int*        in2_ =0; /* te2p       */
    double*     in3_ =0; /* p          */
    double*     in4_ =0; /* reg        */
    int         in5_;    /* nsource    */
    double*     in6_ =0; /* rs         */
    double*     in7_ =0; /* js         */
    double*     in8_ =0; /* rhs        */
    double      in9_;    /* iprec      */
    int         dim10_;   /* 4          */
    int         dim11_;   /* nte        */
    int         dim12_;   /* 3          */
    int         dim13_;   /* np         */
    int         dim14_;   /* nte        */
    int         dim15_;   /* 3          */
    int         dim16_;   /* nsource    */
    int         dim17_;   /* 3          */
    int         dim18_;   /* nsource    */
    int         dim19_;   /* np         */

    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);
    dim18_ = (int) mxWrapGetScalar(prhs[18], &mw_err_txt_);
    dim19_ = (int) mxWrapGetScalar(prhs[19], &mw_err_txt_);

    if (mxGetM(prhs[2]) != dim10_ ||
        mxGetN(prhs[2]) != dim11_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim12_ ||
        mxGetN(prhs[3]) != dim13_) {
        mw_err_txt_ = "Bad argument size: p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim14_) {
        mw_err_txt_ = "Bad argument size: reg";        goto mw_err_label;
    }

    if (mxGetM(prhs[6]) != dim15_ ||
        mxGetN(prhs[6]) != dim16_) {
        mw_err_txt_ = "Bad argument size: rs";
        goto mw_err_label;
    }

    if (mxGetM(prhs[7]) != dim17_ ||
        mxGetN(prhs[7]) != dim18_) {
        mw_err_txt_ = "Bad argument size: js";
        goto mw_err_label;
    }

    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != dim19_) {
        mw_err_txt_ = "Bad argument size: rhs";        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    in5_ = (int) mxWrapGetScalar(prhs[5], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxGetPr(prhs[6]);
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxGetPr(prhs[7]);
    } else
        in7_ = NULL;
    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != 0) {
        in8_ = mxWrapGetArray_double(prhs[8], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in8_ = NULL;
    in9_ = (double) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mexprofrecord_)
        mexprofrecord_[29]++;
    rhsfuncks(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_, in9_);
    plhs[0] = mxCreateDoubleMatrix(dim19_, 1, mxREAL);
    mxWrapCopy_double(plhs[0], in8_, dim19_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in8_)  mxFree(in8_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 426 ----
 * rhsfunc2ndks(int nte, int np, int[10, nte] te2p, double[3, np] p, double[nte] reg, int nsource, double[3, nsource] rs, double[3, nsource] js, inout double[np] rhs, double iprec);
 * Also at demo.mw: 448
 */
const char* stubids30_ = "rhsfunc2ndks(i int, i int, i int[xx], i double[xx], i double[x], i int, i double[xx], i double[xx], io double[x], i double)";

void mexStub30(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nte        */
    int         in1_;    /* np         */
    int*        in2_ =0; /* te2p       */
    double*     in3_ =0; /* p          */
    double*     in4_ =0; /* reg        */
    int         in5_;    /* nsource    */
    double*     in6_ =0; /* rs         */
    double*     in7_ =0; /* js         */
    double*     in8_ =0; /* rhs        */
    double      in9_;    /* iprec      */
    int         dim10_;   /* 10         */
    int         dim11_;   /* nte        */
    int         dim12_;   /* 3          */
    int         dim13_;   /* np         */
    int         dim14_;   /* nte        */
    int         dim15_;   /* 3          */
    int         dim16_;   /* nsource    */
    int         dim17_;   /* 3          */
    int         dim18_;   /* nsource    */
    int         dim19_;   /* np         */

    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);
    dim18_ = (int) mxWrapGetScalar(prhs[18], &mw_err_txt_);
    dim19_ = (int) mxWrapGetScalar(prhs[19], &mw_err_txt_);

    if (mxGetM(prhs[2]) != dim10_ ||
        mxGetN(prhs[2]) != dim11_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim12_ ||
        mxGetN(prhs[3]) != dim13_) {
        mw_err_txt_ = "Bad argument size: p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim14_) {
        mw_err_txt_ = "Bad argument size: reg";        goto mw_err_label;
    }

    if (mxGetM(prhs[6]) != dim15_ ||
        mxGetN(prhs[6]) != dim16_) {
        mw_err_txt_ = "Bad argument size: rs";
        goto mw_err_label;
    }

    if (mxGetM(prhs[7]) != dim17_ ||
        mxGetN(prhs[7]) != dim18_) {
        mw_err_txt_ = "Bad argument size: js";
        goto mw_err_label;
    }

    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != dim19_) {
        mw_err_txt_ = "Bad argument size: rhs";        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    in5_ = (int) mxWrapGetScalar(prhs[5], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxGetPr(prhs[6]);
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxGetPr(prhs[7]);
    } else
        in7_ = NULL;
    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != 0) {
        in8_ = mxWrapGetArray_double(prhs[8], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in8_ = NULL;
    in9_ = (double) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mexprofrecord_)
        mexprofrecord_[30]++;
    rhsfunc2ndks(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_, in9_);
    plhs[0] = mxCreateDoubleMatrix(dim19_, 1, mxREAL);
    mxWrapCopy_double(plhs[0], in8_, dim19_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in8_)  mxFree(in8_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 428 ----
 * rhsfunc3rdks(int nte, int np, int[20, nte] te2p, double[3, np] p, double[nte] reg, int nsource, double[3, nsource] rs, double[3, nsource] js, inout double[np] rhs, double iprec);
 * Also at demo.mw: 450
 */
const char* stubids31_ = "rhsfunc3rdks(i int, i int, i int[xx], i double[xx], i double[x], i int, i double[xx], i double[xx], io double[x], i double)";

void mexStub31(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nte        */
    int         in1_;    /* np         */
    int*        in2_ =0; /* te2p       */
    double*     in3_ =0; /* p          */
    double*     in4_ =0; /* reg        */
    int         in5_;    /* nsource    */
    double*     in6_ =0; /* rs         */
    double*     in7_ =0; /* js         */
    double*     in8_ =0; /* rhs        */
    double      in9_;    /* iprec      */
    int         dim10_;   /* 20         */
    int         dim11_;   /* nte        */
    int         dim12_;   /* 3          */
    int         dim13_;   /* np         */
    int         dim14_;   /* nte        */
    int         dim15_;   /* 3          */
    int         dim16_;   /* nsource    */
    int         dim17_;   /* 3          */
    int         dim18_;   /* nsource    */
    int         dim19_;   /* np         */

    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);
    dim18_ = (int) mxWrapGetScalar(prhs[18], &mw_err_txt_);
    dim19_ = (int) mxWrapGetScalar(prhs[19], &mw_err_txt_);

    if (mxGetM(prhs[2]) != dim10_ ||
        mxGetN(prhs[2]) != dim11_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim12_ ||
        mxGetN(prhs[3]) != dim13_) {
        mw_err_txt_ = "Bad argument size: p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim14_) {
        mw_err_txt_ = "Bad argument size: reg";        goto mw_err_label;
    }

    if (mxGetM(prhs[6]) != dim15_ ||
        mxGetN(prhs[6]) != dim16_) {
        mw_err_txt_ = "Bad argument size: rs";
        goto mw_err_label;
    }

    if (mxGetM(prhs[7]) != dim17_ ||
        mxGetN(prhs[7]) != dim18_) {
        mw_err_txt_ = "Bad argument size: js";
        goto mw_err_label;
    }

    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != dim19_) {
        mw_err_txt_ = "Bad argument size: rhs";        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    in5_ = (int) mxWrapGetScalar(prhs[5], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxGetPr(prhs[6]);
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxGetPr(prhs[7]);
    } else
        in7_ = NULL;
    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != 0) {
        in8_ = mxWrapGetArray_double(prhs[8], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in8_ = NULL;
    in9_ = (double) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mexprofrecord_)
        mexprofrecord_[31]++;
    rhsfunc3rdks(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_, in9_);
    plhs[0] = mxCreateDoubleMatrix(dim19_, 1, mxREAL);
    mxWrapCopy_double(plhs[0], in8_, dim19_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in8_)  mxFree(in8_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 446 ----
 * rhsfunckstensor(int nte, int np, int[4, nte] te2p, double[3, np] p, double[9, nte] reg, int nsource, double[3, nsource] rs, double[3, nsource] js, inout double[np] rhs, double iprec);
 */
const char* stubids32_ = "rhsfunckstensor(i int, i int, i int[xx], i double[xx], i double[xx], i int, i double[xx], i double[xx], io double[x], i double)";

void mexStub32(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nte        */
    int         in1_;    /* np         */
    int*        in2_ =0; /* te2p       */
    double*     in3_ =0; /* p          */
    double*     in4_ =0; /* reg        */
    int         in5_;    /* nsource    */
    double*     in6_ =0; /* rs         */
    double*     in7_ =0; /* js         */
    double*     in8_ =0; /* rhs        */
    double      in9_;    /* iprec      */
    int         dim10_;   /* 4          */
    int         dim11_;   /* nte        */
    int         dim12_;   /* 3          */
    int         dim13_;   /* np         */
    int         dim14_;   /* 9          */
    int         dim15_;   /* nte        */
    int         dim16_;   /* 3          */
    int         dim17_;   /* nsource    */
    int         dim18_;   /* 3          */
    int         dim19_;   /* nsource    */
    int         dim20_;   /* np         */

    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);
    dim18_ = (int) mxWrapGetScalar(prhs[18], &mw_err_txt_);
    dim19_ = (int) mxWrapGetScalar(prhs[19], &mw_err_txt_);
    dim20_ = (int) mxWrapGetScalar(prhs[20], &mw_err_txt_);

    if (mxGetM(prhs[2]) != dim10_ ||
        mxGetN(prhs[2]) != dim11_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim12_ ||
        mxGetN(prhs[3]) != dim13_) {
        mw_err_txt_ = "Bad argument size: p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4]) != dim14_ ||
        mxGetN(prhs[4]) != dim15_) {
        mw_err_txt_ = "Bad argument size: reg";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6]) != dim16_ ||
        mxGetN(prhs[6]) != dim17_) {
        mw_err_txt_ = "Bad argument size: rs";
        goto mw_err_label;
    }

    if (mxGetM(prhs[7]) != dim18_ ||
        mxGetN(prhs[7]) != dim19_) {
        mw_err_txt_ = "Bad argument size: js";
        goto mw_err_label;
    }

    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != dim20_) {
        mw_err_txt_ = "Bad argument size: rhs";        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    in5_ = (int) mxWrapGetScalar(prhs[5], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxGetPr(prhs[6]);
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxGetPr(prhs[7]);
    } else
        in7_ = NULL;
    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != 0) {
        in8_ = mxWrapGetArray_double(prhs[8], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in8_ = NULL;
    in9_ = (double) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mexprofrecord_)
        mexprofrecord_[32]++;
    rhsfunckstensor(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_, in9_);
    plhs[0] = mxCreateDoubleMatrix(dim20_, 1, mxREAL);
    mxWrapCopy_double(plhs[0], in8_, dim20_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in8_)  mxFree(in8_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 457 ----
 * Eprim(int nsource, int ntarget, double[3, nsource] rs, double[3, nsource] js, double[3, ntarget] robs, inout double[3, ntarget] Eout, double iprecEp);
 * Also at demo.mw: 482
 */
const char* stubids35_ = "Eprim(i int, i int, i double[xx], i double[xx], i double[xx], io double[xx], i double)";

void mexStub35(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nsource    */
    int         in1_;    /* ntarget    */
    double*     in2_ =0; /* rs         */
    double*     in3_ =0; /* js         */
    double*     in4_ =0; /* robs       */
    double*     in5_ =0; /* Eout       */
    double      in6_;    /* iprecEp    */
    int         dim7_;   /* 3          */
    int         dim8_;   /* nsource    */
    int         dim9_;   /* 3          */
    int         dim10_;   /* nsource    */
    int         dim11_;   /* 3          */
    int         dim12_;   /* ntarget    */
    int         dim13_;   /* 3          */
    int         dim14_;   /* ntarget    */

    dim7_ = (int) mxWrapGetScalar(prhs[7], &mw_err_txt_);
    dim8_ = (int) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);

    if (mxGetM(prhs[2]) != dim7_ ||
        mxGetN(prhs[2]) != dim8_) {
        mw_err_txt_ = "Bad argument size: rs";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim9_ ||
        mxGetN(prhs[3]) != dim10_) {
        mw_err_txt_ = "Bad argument size: js";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4]) != dim11_ ||
        mxGetN(prhs[4]) != dim12_) {
        mw_err_txt_ = "Bad argument size: robs";
        goto mw_err_label;
    }

    if (mxGetM(prhs[5]) != dim13_ ||
        mxGetN(prhs[5]) != dim14_) {
        mw_err_txt_ = "Bad argument size: Eout";
        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxGetPr(prhs[2]);
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxWrapGetArray_double(prhs[5], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in5_ = NULL;
    in6_ = (double) mxWrapGetScalar(prhs[6], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mexprofrecord_)
        mexprofrecord_[35]++;
    Eprim(in0_, in1_, in2_, in3_, in4_, in5_, in6_);
    plhs[0] = mxCreateDoubleMatrix(dim13_, dim14_, mxREAL);
    mxWrapCopy_double(plhs[0], in5_, dim13_*dim14_);

mw_err_label:
    if (in5_)  mxFree(in5_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 464 ----
 * Hprim(int nsource, int ntarget, double[3, nsource] rs, double[3, nsource] js, double[3, ntarget] robs, inout double[3, ntarget] Hout, double iprecEp);
 * Also at demo.mw: 503
 */
const char* stubids36_ = "Hprim(i int, i int, i double[xx], i double[xx], i double[xx], io double[xx], i double)";

void mexStub36(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nsource    */
    int         in1_;    /* ntarget    */
    double*     in2_ =0; /* rs         */
    double*     in3_ =0; /* js         */
    double*     in4_ =0; /* robs       */
    double*     in5_ =0; /* Hout       */
    double      in6_;    /* iprecEp    */
    int         dim7_;   /* 3          */
    int         dim8_;   /* nsource    */
    int         dim9_;   /* 3          */
    int         dim10_;   /* nsource    */
    int         dim11_;   /* 3          */
    int         dim12_;   /* ntarget    */
    int         dim13_;   /* 3          */
    int         dim14_;   /* ntarget    */

    dim7_ = (int) mxWrapGetScalar(prhs[7], &mw_err_txt_);
    dim8_ = (int) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);

    if (mxGetM(prhs[2]) != dim7_ ||
        mxGetN(prhs[2]) != dim8_) {
        mw_err_txt_ = "Bad argument size: rs";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim9_ ||
        mxGetN(prhs[3]) != dim10_) {
        mw_err_txt_ = "Bad argument size: js";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4]) != dim11_ ||
        mxGetN(prhs[4]) != dim12_) {
        mw_err_txt_ = "Bad argument size: robs";
        goto mw_err_label;
    }

    if (mxGetM(prhs[5]) != dim13_ ||
        mxGetN(prhs[5]) != dim14_) {
        mw_err_txt_ = "Bad argument size: Hout";
        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxGetPr(prhs[2]);
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxWrapGetArray_double(prhs[5], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in5_ = NULL;
    in6_ = (double) mxWrapGetScalar(prhs[6], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mexprofrecord_)
        mexprofrecord_[36]++;
    Hprim(in0_, in1_, in2_, in3_, in4_, in5_, in6_);
    plhs[0] = mxCreateDoubleMatrix(dim13_, dim14_, mxREAL);
    mxWrapCopy_double(plhs[0], in5_, dim13_*dim14_);

mw_err_label:
    if (in5_)  mxFree(in5_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 473 ----
 * Aprim(int nsource, int ntarget, double[3, nsource] rs, double[3, nsource] js, double[3, ntarget] robs, inout double[3, ntarget] Aout, inout double[3, ntarget] curlAout, double iprecEp);
 */
const char* stubids37_ = "Aprim(i int, i int, i double[xx], i double[xx], i double[xx], io double[xx], io double[xx], i double)";

void mexStub37(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nsource    */
    int         in1_;    /* ntarget    */
    double*     in2_ =0; /* rs         */
    double*     in3_ =0; /* js         */
    double*     in4_ =0; /* robs       */
    double*     in5_ =0; /* Aout       */
    double*     in6_ =0; /* curlAout   */
    double      in7_;    /* iprecEp    */
    int         dim8_;   /* 3          */
    int         dim9_;   /* nsource    */
    int         dim10_;   /* 3          */
    int         dim11_;   /* nsource    */
    int         dim12_;   /* 3          */
    int         dim13_;   /* ntarget    */
    int         dim14_;   /* 3          */
    int         dim15_;   /* ntarget    */
    int         dim16_;   /* 3          */
    int         dim17_;   /* ntarget    */

    dim8_ = (int) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);

    if (mxGetM(prhs[2]) != dim8_ ||
        mxGetN(prhs[2]) != dim9_) {
        mw_err_txt_ = "Bad argument size: rs";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim10_ ||
        mxGetN(prhs[3]) != dim11_) {
        mw_err_txt_ = "Bad argument size: js";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4]) != dim12_ ||
        mxGetN(prhs[4]) != dim13_) {
        mw_err_txt_ = "Bad argument size: robs";
        goto mw_err_label;
    }

    if (mxGetM(prhs[5]) != dim14_ ||
        mxGetN(prhs[5]) != dim15_) {
        mw_err_txt_ = "Bad argument size: Aout";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6]) != dim16_ ||
        mxGetN(prhs[6]) != dim17_) {
        mw_err_txt_ = "Bad argument size: curlAout";
        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxGetPr(prhs[2]);
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxGetPr(prhs[3]);
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxWrapGetArray_double(prhs[5], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_double(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    in7_ = (double) mxWrapGetScalar(prhs[7], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mexprofrecord_)
        mexprofrecord_[37]++;
    Aprim(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_);
    plhs[0] = mxCreateDoubleMatrix(dim14_, dim15_, mxREAL);
    mxWrapCopy_double(plhs[0], in5_, dim14_*dim15_);
    plhs[1] = mxCreateDoubleMatrix(dim16_, dim17_, mxREAL);
    mxWrapCopy_double(plhs[1], in6_, dim16_*dim17_);

mw_err_label:
    if (in5_)  mxFree(in5_);
    if (in6_)  mxFree(in6_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 486 ----
 * secondaryEfield(double[3, nobs] robs, int nobs, int[4, nte] te2te, int[4, nte] te2p, int nte, double[3, np] p, double[np] xval, inout double[3, nobs] Eout, double iprecEp);
 * Also at demo.mw: 508
 */
const char* stubids39_ = "secondaryEfield(i double[xx], i int, i int[xx], i int[xx], i int, i double[xx], i double[x], io double[xx], i double)";

void mexStub39(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    double*     in0_ =0; /* robs       */
    int         in1_;    /* nobs       */
    int*        in2_ =0; /* te2te      */
    int*        in3_ =0; /* te2p       */
    int         in4_;    /* nte        */
    double*     in5_ =0; /* p          */
    double*     in6_ =0; /* xval       */
    double*     in7_ =0; /* Eout       */
    double      in8_;    /* iprecEp    */
    int         dim9_;   /* 3          */
    int         dim10_;   /* nobs       */
    int         dim11_;   /* 4          */
    int         dim12_;   /* nte        */
    int         dim13_;   /* 4          */
    int         dim14_;   /* nte        */
    int         dim15_;   /* 3          */
    int         dim16_;   /* np         */
    int         dim17_;   /* np         */
    int         dim18_;   /* 3          */
    int         dim19_;   /* nobs       */

    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);
    dim18_ = (int) mxWrapGetScalar(prhs[18], &mw_err_txt_);
    dim19_ = (int) mxWrapGetScalar(prhs[19], &mw_err_txt_);

    if (mxGetM(prhs[0]) != dim9_ ||
        mxGetN(prhs[0]) != dim10_) {
        mw_err_txt_ = "Bad argument size: robs";
        goto mw_err_label;
    }

    if (mxGetM(prhs[2]) != dim11_ ||
        mxGetN(prhs[2]) != dim12_) {
        mw_err_txt_ = "Bad argument size: te2te";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim13_ ||
        mxGetN(prhs[3]) != dim14_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[5]) != dim15_ ||
        mxGetN(prhs[5]) != dim16_) {
        mw_err_txt_ = "Bad argument size: p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != dim17_) {
        mw_err_txt_ = "Bad argument size: xval";        goto mw_err_label;
    }

    if (mxGetM(prhs[7]) != dim18_ ||
        mxGetN(prhs[7]) != dim19_) {
        mw_err_txt_ = "Bad argument size: Eout";
        goto mw_err_label;
    }

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != 0) {
        in0_ = mxGetPr(prhs[0]);
    } else
        in0_ = NULL;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxWrapGetArray_int(prhs[3], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in3_ = NULL;
    in4_ = (int) mxWrapGetScalar(prhs[4], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxGetPr(prhs[5]);
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxGetPr(prhs[6]);
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxWrapGetArray_double(prhs[7], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in7_ = NULL;
    in8_ = (double) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mexprofrecord_)
        mexprofrecord_[39]++;
    secondaryEfield(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_);
    plhs[0] = mxCreateDoubleMatrix(dim18_, dim19_, mxREAL);
    mxWrapCopy_double(plhs[0], in7_, dim18_*dim19_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in3_)  mxFree(in3_);
    if (in7_)  mxFree(in7_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 489 ----
 * secondaryEfield2nd(double[3, nobs] robs, int nobs, int[4, nte] te2te, int[10, nte] te2p, int nte, double[3, np] p, double[np] xval, inout double[3, nobs] Eout, double iprecEp);
 * Also at demo.mw: 511
 */
const char* stubids40_ = "secondaryEfield2nd(i double[xx], i int, i int[xx], i int[xx], i int, i double[xx], i double[x], io double[xx], i double)";

void mexStub40(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    double*     in0_ =0; /* robs       */
    int         in1_;    /* nobs       */
    int*        in2_ =0; /* te2te      */
    int*        in3_ =0; /* te2p       */
    int         in4_;    /* nte        */
    double*     in5_ =0; /* p          */
    double*     in6_ =0; /* xval       */
    double*     in7_ =0; /* Eout       */
    double      in8_;    /* iprecEp    */
    int         dim9_;   /* 3          */
    int         dim10_;   /* nobs       */
    int         dim11_;   /* 4          */
    int         dim12_;   /* nte        */
    int         dim13_;   /* 10         */
    int         dim14_;   /* nte        */
    int         dim15_;   /* 3          */
    int         dim16_;   /* np         */
    int         dim17_;   /* np         */
    int         dim18_;   /* 3          */
    int         dim19_;   /* nobs       */

    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);
    dim18_ = (int) mxWrapGetScalar(prhs[18], &mw_err_txt_);
    dim19_ = (int) mxWrapGetScalar(prhs[19], &mw_err_txt_);

    if (mxGetM(prhs[0]) != dim9_ ||
        mxGetN(prhs[0]) != dim10_) {
        mw_err_txt_ = "Bad argument size: robs";
        goto mw_err_label;
    }

    if (mxGetM(prhs[2]) != dim11_ ||
        mxGetN(prhs[2]) != dim12_) {
        mw_err_txt_ = "Bad argument size: te2te";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim13_ ||
        mxGetN(prhs[3]) != dim14_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[5]) != dim15_ ||
        mxGetN(prhs[5]) != dim16_) {
        mw_err_txt_ = "Bad argument size: p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != dim17_) {
        mw_err_txt_ = "Bad argument size: xval";        goto mw_err_label;
    }

    if (mxGetM(prhs[7]) != dim18_ ||
        mxGetN(prhs[7]) != dim19_) {
        mw_err_txt_ = "Bad argument size: Eout";
        goto mw_err_label;
    }

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != 0) {
        in0_ = mxGetPr(prhs[0]);
    } else
        in0_ = NULL;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxWrapGetArray_int(prhs[3], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in3_ = NULL;
    in4_ = (int) mxWrapGetScalar(prhs[4], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxGetPr(prhs[5]);
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxGetPr(prhs[6]);
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxWrapGetArray_double(prhs[7], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in7_ = NULL;
    in8_ = (double) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mexprofrecord_)
        mexprofrecord_[40]++;
    secondaryEfield2nd(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_);
    plhs[0] = mxCreateDoubleMatrix(dim18_, dim19_, mxREAL);
    mxWrapCopy_double(plhs[0], in7_, dim18_*dim19_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in3_)  mxFree(in3_);
    if (in7_)  mxFree(in7_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 492 ----
 * secondaryEfield3rd(double[3, nobs] robs, int nobs, int[4, nte] te2te, int[20, nte] te2p, int nte, double[3, np] p, double[np] xval, inout double[3, nobs] Eout, double iprecEp);
 * Also at demo.mw: 514
 */
const char* stubids41_ = "secondaryEfield3rd(i double[xx], i int, i int[xx], i int[xx], i int, i double[xx], i double[x], io double[xx], i double)";

void mexStub41(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    double*     in0_ =0; /* robs       */
    int         in1_;    /* nobs       */
    int*        in2_ =0; /* te2te      */
    int*        in3_ =0; /* te2p       */
    int         in4_;    /* nte        */
    double*     in5_ =0; /* p          */
    double*     in6_ =0; /* xval       */
    double*     in7_ =0; /* Eout       */
    double      in8_;    /* iprecEp    */
    int         dim9_;   /* 3          */
    int         dim10_;   /* nobs       */
    int         dim11_;   /* 4          */
    int         dim12_;   /* nte        */
    int         dim13_;   /* 20         */
    int         dim14_;   /* nte        */
    int         dim15_;   /* 3          */
    int         dim16_;   /* np         */
    int         dim17_;   /* np         */
    int         dim18_;   /* 3          */
    int         dim19_;   /* nobs       */

    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);
    dim18_ = (int) mxWrapGetScalar(prhs[18], &mw_err_txt_);
    dim19_ = (int) mxWrapGetScalar(prhs[19], &mw_err_txt_);

    if (mxGetM(prhs[0]) != dim9_ ||
        mxGetN(prhs[0]) != dim10_) {
        mw_err_txt_ = "Bad argument size: robs";
        goto mw_err_label;
    }

    if (mxGetM(prhs[2]) != dim11_ ||
        mxGetN(prhs[2]) != dim12_) {
        mw_err_txt_ = "Bad argument size: te2te";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim13_ ||
        mxGetN(prhs[3]) != dim14_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[5]) != dim15_ ||
        mxGetN(prhs[5]) != dim16_) {
        mw_err_txt_ = "Bad argument size: p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != dim17_) {
        mw_err_txt_ = "Bad argument size: xval";        goto mw_err_label;
    }

    if (mxGetM(prhs[7]) != dim18_ ||
        mxGetN(prhs[7]) != dim19_) {
        mw_err_txt_ = "Bad argument size: Eout";
        goto mw_err_label;
    }

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != 0) {
        in0_ = mxGetPr(prhs[0]);
    } else
        in0_ = NULL;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxWrapGetArray_int(prhs[3], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in3_ = NULL;
    in4_ = (int) mxWrapGetScalar(prhs[4], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxGetPr(prhs[5]);
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxGetPr(prhs[6]);
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxWrapGetArray_double(prhs[7], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in7_ = NULL;
    in8_ = (double) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mexprofrecord_)
        mexprofrecord_[41]++;
    secondaryEfield3rd(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_);
    plhs[0] = mxCreateDoubleMatrix(dim18_, dim19_, mxREAL);
    mxWrapCopy_double(plhs[0], in7_, dim18_*dim19_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in3_)  mxFree(in3_);
    if (in7_)  mxFree(in7_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 530 ----
 * secondaryEfieldc(double[3, nobs] robs, int nobs, int[4, nte] te2te, int[4, nte] te2p, int nte, double[3, np] p, double[np] xval, inout double[3, nobs] Eout, double iprecEp);
 * Also at demo.mw: 552
 */
const char* stubids47_ = "secondaryEfieldc(i double[xx], i int, i int[xx], i int[xx], i int, i double[xx], i double[x], io double[xx], i double)";

void mexStub47(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    double*     in0_ =0; /* robs       */
    int         in1_;    /* nobs       */
    int*        in2_ =0; /* te2te      */
    int*        in3_ =0; /* te2p       */
    int         in4_;    /* nte        */
    double*     in5_ =0; /* p          */
    double*     in6_ =0; /* xval       */
    double*     in7_ =0; /* Eout       */
    double      in8_;    /* iprecEp    */
    int         dim9_;   /* 3          */
    int         dim10_;   /* nobs       */
    int         dim11_;   /* 4          */
    int         dim12_;   /* nte        */
    int         dim13_;   /* 4          */
    int         dim14_;   /* nte        */
    int         dim15_;   /* 3          */
    int         dim16_;   /* np         */
    int         dim17_;   /* np         */
    int         dim18_;   /* 3          */
    int         dim19_;   /* nobs       */

    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);
    dim18_ = (int) mxWrapGetScalar(prhs[18], &mw_err_txt_);
    dim19_ = (int) mxWrapGetScalar(prhs[19], &mw_err_txt_);

    if (mxGetM(prhs[0]) != dim9_ ||
        mxGetN(prhs[0]) != dim10_) {
        mw_err_txt_ = "Bad argument size: robs";
        goto mw_err_label;
    }

    if (mxGetM(prhs[2]) != dim11_ ||
        mxGetN(prhs[2]) != dim12_) {
        mw_err_txt_ = "Bad argument size: te2te";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim13_ ||
        mxGetN(prhs[3]) != dim14_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[5]) != dim15_ ||
        mxGetN(prhs[5]) != dim16_) {
        mw_err_txt_ = "Bad argument size: p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != dim17_) {
        mw_err_txt_ = "Bad argument size: xval";        goto mw_err_label;
    }

    if (mxGetM(prhs[7]) != dim18_ ||
        mxGetN(prhs[7]) != dim19_) {
        mw_err_txt_ = "Bad argument size: Eout";
        goto mw_err_label;
    }

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != 0) {
        in0_ = mxGetPr(prhs[0]);
    } else
        in0_ = NULL;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxWrapGetArray_int(prhs[3], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in3_ = NULL;
    in4_ = (int) mxWrapGetScalar(prhs[4], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxGetPr(prhs[5]);
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxGetPr(prhs[6]);
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxWrapGetArray_double(prhs[7], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in7_ = NULL;
    in8_ = (double) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mexprofrecord_)
        mexprofrecord_[47]++;
    secondaryEfieldc(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_);
    plhs[0] = mxCreateDoubleMatrix(dim18_, dim19_, mxREAL);
    mxWrapCopy_double(plhs[0], in7_, dim18_*dim19_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in3_)  mxFree(in3_);
    if (in7_)  mxFree(in7_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 533 ----
 * secondaryEfield2ndc(double[3, nobs] robs, int nobs, int[4, nte] te2te, int[10, nte] te2p, int nte, double[3, np] p, double[np] xval, inout double[3, nobs] Eout, double iprecEp);
 * Also at demo.mw: 555
 */
const char* stubids48_ = "secondaryEfield2ndc(i double[xx], i int, i int[xx], i int[xx], i int, i double[xx], i double[x], io double[xx], i double)";

void mexStub48(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    double*     in0_ =0; /* robs       */
    int         in1_;    /* nobs       */
    int*        in2_ =0; /* te2te      */
    int*        in3_ =0; /* te2p       */
    int         in4_;    /* nte        */
    double*     in5_ =0; /* p          */
    double*     in6_ =0; /* xval       */
    double*     in7_ =0; /* Eout       */
    double      in8_;    /* iprecEp    */
    int         dim9_;   /* 3          */
    int         dim10_;   /* nobs       */
    int         dim11_;   /* 4          */
    int         dim12_;   /* nte        */
    int         dim13_;   /* 10         */
    int         dim14_;   /* nte        */
    int         dim15_;   /* 3          */
    int         dim16_;   /* np         */
    int         dim17_;   /* np         */
    int         dim18_;   /* 3          */
    int         dim19_;   /* nobs       */

    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);
    dim18_ = (int) mxWrapGetScalar(prhs[18], &mw_err_txt_);
    dim19_ = (int) mxWrapGetScalar(prhs[19], &mw_err_txt_);

    if (mxGetM(prhs[0]) != dim9_ ||
        mxGetN(prhs[0]) != dim10_) {
        mw_err_txt_ = "Bad argument size: robs";
        goto mw_err_label;
    }

    if (mxGetM(prhs[2]) != dim11_ ||
        mxGetN(prhs[2]) != dim12_) {
        mw_err_txt_ = "Bad argument size: te2te";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim13_ ||
        mxGetN(prhs[3]) != dim14_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[5]) != dim15_ ||
        mxGetN(prhs[5]) != dim16_) {
        mw_err_txt_ = "Bad argument size: p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != dim17_) {
        mw_err_txt_ = "Bad argument size: xval";        goto mw_err_label;
    }

    if (mxGetM(prhs[7]) != dim18_ ||
        mxGetN(prhs[7]) != dim19_) {
        mw_err_txt_ = "Bad argument size: Eout";
        goto mw_err_label;
    }

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != 0) {
        in0_ = mxGetPr(prhs[0]);
    } else
        in0_ = NULL;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxWrapGetArray_int(prhs[3], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in3_ = NULL;
    in4_ = (int) mxWrapGetScalar(prhs[4], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxGetPr(prhs[5]);
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxGetPr(prhs[6]);
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxWrapGetArray_double(prhs[7], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in7_ = NULL;
    in8_ = (double) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mexprofrecord_)
        mexprofrecord_[48]++;
    secondaryEfield2ndc(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_);
    plhs[0] = mxCreateDoubleMatrix(dim18_, dim19_, mxREAL);
    mxWrapCopy_double(plhs[0], in7_, dim18_*dim19_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in3_)  mxFree(in3_);
    if (in7_)  mxFree(in7_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 536 ----
 * secondaryEfield3rdc(double[3, nobs] robs, int nobs, int[4, nte] te2te, int[20, nte] te2p, int nte, double[3, np] p, double[np] xval, inout double[3, nobs] Eout, double iprecEp);
 * Also at demo.mw: 558
 */
const char* stubids49_ = "secondaryEfield3rdc(i double[xx], i int, i int[xx], i int[xx], i int, i double[xx], i double[x], io double[xx], i double)";

void mexStub49(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    double*     in0_ =0; /* robs       */
    int         in1_;    /* nobs       */
    int*        in2_ =0; /* te2te      */
    int*        in3_ =0; /* te2p       */
    int         in4_;    /* nte        */
    double*     in5_ =0; /* p          */
    double*     in6_ =0; /* xval       */
    double*     in7_ =0; /* Eout       */
    double      in8_;    /* iprecEp    */
    int         dim9_;   /* 3          */
    int         dim10_;   /* nobs       */
    int         dim11_;   /* 4          */
    int         dim12_;   /* nte        */
    int         dim13_;   /* 20         */
    int         dim14_;   /* nte        */
    int         dim15_;   /* 3          */
    int         dim16_;   /* np         */
    int         dim17_;   /* np         */
    int         dim18_;   /* 3          */
    int         dim19_;   /* nobs       */

    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);
    dim18_ = (int) mxWrapGetScalar(prhs[18], &mw_err_txt_);
    dim19_ = (int) mxWrapGetScalar(prhs[19], &mw_err_txt_);

    if (mxGetM(prhs[0]) != dim9_ ||
        mxGetN(prhs[0]) != dim10_) {
        mw_err_txt_ = "Bad argument size: robs";
        goto mw_err_label;
    }

    if (mxGetM(prhs[2]) != dim11_ ||
        mxGetN(prhs[2]) != dim12_) {
        mw_err_txt_ = "Bad argument size: te2te";
        goto mw_err_label;
    }

    if (mxGetM(prhs[3]) != dim13_ ||
        mxGetN(prhs[3]) != dim14_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[5]) != dim15_ ||
        mxGetN(prhs[5]) != dim16_) {
        mw_err_txt_ = "Bad argument size: p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != dim17_) {
        mw_err_txt_ = "Bad argument size: xval";        goto mw_err_label;
    }

    if (mxGetM(prhs[7]) != dim18_ ||
        mxGetN(prhs[7]) != dim19_) {
        mw_err_txt_ = "Bad argument size: Eout";
        goto mw_err_label;
    }

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != 0) {
        in0_ = mxGetPr(prhs[0]);
    } else
        in0_ = NULL;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxWrapGetArray_int(prhs[3], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in3_ = NULL;
    in4_ = (int) mxWrapGetScalar(prhs[4], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxGetPr(prhs[5]);
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxGetPr(prhs[6]);
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxWrapGetArray_double(prhs[7], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in7_ = NULL;
    in8_ = (double) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mexprofrecord_)
        mexprofrecord_[49]++;
    secondaryEfield3rdc(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_);
    plhs[0] = mxCreateDoubleMatrix(dim18_, dim19_, mxREAL);
    mxWrapCopy_double(plhs[0], in7_, dim18_*dim19_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in3_)  mxFree(in3_);
    if (in7_)  mxFree(in7_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 570 ----
 * pointLocation_f(double[3, nobs] robs, int nobs, int[te_stride, nte] te2p, int nte, double[3, np] p, int te_stride, int[4, nte] te2te, inout int[nobs] teid, inout double[4, nobs] bari);
 */
const char* stubids54_ = "pointLocation_f(i double[xx], i int, i int[xx], i int, i double[xx], i int, i int[xx], io int[x], io double[xx])";

void mexStub54(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    double*     in0_ =0; /* robs       */
    int         in1_;    /* nobs       */
    int*        in2_ =0; /* te2p       */
    int         in3_;    /* nte        */
    double*     in4_ =0; /* p          */
    int         in5_;    /* te_stride  */
    int*        in6_ =0; /* te2te      */
    int*        in7_ =0; /* teid       */
    double*     in8_ =0; /* bari       */
    int         dim9_;   /* 3          */
    int         dim10_;   /* nobs       */
    int         dim11_;   /* te_stride  */
    int         dim12_;   /* nte        */
    int         dim13_;   /* 3          */
    int         dim14_;   /* np         */
    int         dim15_;   /* 4          */
    int         dim16_;   /* nte        */
    int         dim17_;   /* nobs       */
    int         dim18_;   /* 4          */
    int         dim19_;   /* nobs       */

    dim9_ = (int) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (int) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (int) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (int) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (int) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (int) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (int) mxWrapGetScalar(prhs[15], &mw_err_txt_);
    dim16_ = (int) mxWrapGetScalar(prhs[16], &mw_err_txt_);
    dim17_ = (int) mxWrapGetScalar(prhs[17], &mw_err_txt_);
    dim18_ = (int) mxWrapGetScalar(prhs[18], &mw_err_txt_);
    dim19_ = (int) mxWrapGetScalar(prhs[19], &mw_err_txt_);

    if (mxGetM(prhs[0]) != dim9_ ||
        mxGetN(prhs[0]) != dim10_) {
        mw_err_txt_ = "Bad argument size: robs";
        goto mw_err_label;
    }

    if (mxGetM(prhs[2]) != dim11_ ||
        mxGetN(prhs[2]) != dim12_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[4]) != dim13_ ||
        mxGetN(prhs[4]) != dim14_) {
        mw_err_txt_ = "Bad argument size: p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[6]) != dim15_ ||
        mxGetN(prhs[6]) != dim16_) {
        mw_err_txt_ = "Bad argument size: te2te";
        goto mw_err_label;
    }

    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != dim17_) {
        mw_err_txt_ = "Bad argument size: teid";        goto mw_err_label;
    }

    if (mxGetM(prhs[8]) != dim18_ ||
        mxGetN(prhs[8]) != dim19_) {
        mw_err_txt_ = "Bad argument size: bari";
        goto mw_err_label;
    }

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != 0) {
        in0_ = mxGetPr(prhs[0]);
    } else
        in0_ = NULL;
    in1_ = (int) mxWrapGetScalar(prhs[1], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    in3_ = (int) mxWrapGetScalar(prhs[3], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxGetPr(prhs[4]);
    } else
        in4_ = NULL;
    in5_ = (int) mxWrapGetScalar(prhs[5], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_int(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxWrapGetArray_int(prhs[7], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in7_ = NULL;
    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != 0) {
        in8_ = mxWrapGetArray_double(prhs[8], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in8_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[54]++;
    pointLocation_f(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_);
    plhs[0] = mxCreateDoubleMatrix(dim17_, 1, mxREAL);
    mxWrapCopy_int(plhs[0], in7_, dim17_);
    plhs[1] = mxCreateDoubleMatrix(dim18_, dim19_, mxREAL);
    mxWrapCopy_double(plhs[1], in8_, dim18_*dim19_);

mw_err_label:
    if (in2_)  mxFree(in2_);
    if (in6_)  mxFree(in6_);
    if (in7_)  mxFree(in7_);
    if (in8_)  mxFree(in8_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- demo.mw: 578 ----
 * generatetetrapaths(int nte, int[4, nte] te2p, inout int[4, nte] te2te);
 */
const char* stubids55_ = "generatetetrapaths(i int, i int[xx], io int[xx])";

void mexStub55(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int         in0_;    /* nte        */
    int*        in1_ =0; /* te2p       */
    int*        in2_ =0; /* te2te      */
    int         dim3_;   /* 4          */
    int         dim4_;   /* nte        */
    int         dim5_;   /* 4          */
    int         dim6_;   /* nte        */

    dim3_ = (int) mxWrapGetScalar(prhs[3], &mw_err_txt_);
    dim4_ = (int) mxWrapGetScalar(prhs[4], &mw_err_txt_);
    dim5_ = (int) mxWrapGetScalar(prhs[5], &mw_err_txt_);
    dim6_ = (int) mxWrapGetScalar(prhs[6], &mw_err_txt_);

    if (mxGetM(prhs[1]) != dim3_ ||
        mxGetN(prhs[1]) != dim4_) {
        mw_err_txt_ = "Bad argument size: te2p";
        goto mw_err_label;
    }

    if (mxGetM(prhs[2]) != dim5_ ||
        mxGetN(prhs[2]) != dim6_) {
        mw_err_txt_ = "Bad argument size: te2te";
        goto mw_err_label;
    }

    in0_ = (int) mxWrapGetScalar(prhs[0], &mw_err_txt_);
    if (mw_err_txt_)
        goto mw_err_label;
    if (mxGetM(prhs[1])*mxGetN(prhs[1]) != 0) {
        in1_ = mxWrapGetArray_int(prhs[1], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in1_ = NULL;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxWrapGetArray_int(prhs[2], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in2_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[55]++;
    generatetetrapaths(in0_, in1_, in2_);
    plhs[0] = mxCreateDoubleMatrix(dim5_, dim6_, mxREAL);
    mxWrapCopy_int(plhs[0], in2_, dim5_*dim6_);

mw_err_label:
    if (in1_)  mxFree(in1_);
    if (in2_)  mxFree(in2_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ----
 */
void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
    char id[512];
    if (nrhs == 0) {
        mexPrintf("Mex function installed\n");
        return;
    }

    if (mxGetString(prhs[0], id, sizeof(id)) != 0)
        mexErrMsgTxt("Identifier should be a string");
    else if (strcmp(id, stubids1_) == 0)
        mexStub1(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids2_) == 0)
        mexStub2(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids3_) == 0)
        mexStub3(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids4_) == 0)
        mexStub4(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids5_) == 0)
        mexStub5(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids6_) == 0)
        mexStub6(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids7_) == 0)
        mexStub7(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids8_) == 0)
        mexStub8(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids9_) == 0)
        mexStub9(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids10_) == 0)
        mexStub10(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids11_) == 0)
        mexStub11(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids12_) == 0)
        mexStub12(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids13_) == 0)
        mexStub13(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids15_) == 0)
        mexStub15(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids16_) == 0)
        mexStub16(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids18_) == 0)
        mexStub18(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids19_) == 0)
        mexStub19(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids20_) == 0)
        mexStub20(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids21_) == 0)
        mexStub21(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids22_) == 0)
        mexStub22(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids23_) == 0)
        mexStub23(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids24_) == 0)
        mexStub24(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids25_) == 0)
        mexStub25(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids26_) == 0)
        mexStub26(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids27_) == 0)
        mexStub27(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids28_) == 0)
        mexStub28(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids29_) == 0)
        mexStub29(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids30_) == 0)
        mexStub30(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids31_) == 0)
        mexStub31(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids32_) == 0)
        mexStub32(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids35_) == 0)
        mexStub35(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids36_) == 0)
        mexStub36(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids37_) == 0)
        mexStub37(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids39_) == 0)
        mexStub39(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids40_) == 0)
        mexStub40(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids41_) == 0)
        mexStub41(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids47_) == 0)
        mexStub47(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids48_) == 0)
        mexStub48(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids49_) == 0)
        mexStub49(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids54_) == 0)
        mexStub54(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids55_) == 0)
        mexStub55(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, "*profile on*") == 0) {
        if (!mexprofrecord_) {
            mexprofrecord_ = (int*) malloc(56 * sizeof(int));
            mexLock();
        }
        memset(mexprofrecord_, 0, 56 * sizeof(int));
    } else if (strcmp(id, "*profile off*") == 0) {
        if (mexprofrecord_) {
            free(mexprofrecord_);
            mexUnlock();
        }
        mexprofrecord_ = NULL;
    } else if (strcmp(id, "*profile report*") == 0) {
        if (!mexprofrecord_)
            mexPrintf("Profiler inactive\n");
        mexPrintf("%d calls to demo.mw:100\n", mexprofrecord_[1]);
        mexPrintf("%d calls to demo.mw:102\n", mexprofrecord_[2]);
        mexPrintf("%d calls to demo.mw:104\n", mexprofrecord_[3]);
        mexPrintf("%d calls to demo.mw:122\n", mexprofrecord_[4]);
        mexPrintf("%d calls to demo.mw:124\n", mexprofrecord_[5]);
        mexPrintf("%d calls to demo.mw:126\n", mexprofrecord_[6]);
        mexPrintf("%d calls to demo.mw:130 (demo.mw:187)\n", mexprofrecord_[7]);
        mexPrintf("%d calls to demo.mw:157\n", mexprofrecord_[8]);
        mexPrintf("%d calls to demo.mw:159\n", mexprofrecord_[9]);
        mexPrintf("%d calls to demo.mw:161\n", mexprofrecord_[10]);
        mexPrintf("%d calls to demo.mw:179\n", mexprofrecord_[11]);
        mexPrintf("%d calls to demo.mw:181\n", mexprofrecord_[12]);
        mexPrintf("%d calls to demo.mw:183\n", mexprofrecord_[13]);
        mexPrintf("%d calls to demo.mw:337 (demo.mw:342)\n", mexprofrecord_[15]);
        mexPrintf("%d calls to demo.mw:340\n", mexprofrecord_[16]);
        mexPrintf("%d calls to demo.mw:343\n", mexprofrecord_[18]);
        mexPrintf("%d calls to demo.mw:346\n", mexprofrecord_[19]);
        mexPrintf("%d calls to demo.mw:358\n", mexprofrecord_[20]);
        mexPrintf("%d calls to demo.mw:363\n", mexprofrecord_[21]);
        mexPrintf("%d calls to demo.mw:368\n", mexprofrecord_[22]);
        mexPrintf("%d calls to demo.mw:382\n", mexprofrecord_[23]);
        mexPrintf("%d calls to demo.mw:387\n", mexprofrecord_[24]);
        mexPrintf("%d calls to demo.mw:392\n", mexprofrecord_[25]);
        mexPrintf("%d calls to demo.mw:406\n", mexprofrecord_[26]);
        mexPrintf("%d calls to demo.mw:408\n", mexprofrecord_[27]);
        mexPrintf("%d calls to demo.mw:410\n", mexprofrecord_[28]);
        mexPrintf("%d calls to demo.mw:424\n", mexprofrecord_[29]);
        mexPrintf("%d calls to demo.mw:426 (demo.mw:448)\n", mexprofrecord_[30]);
        mexPrintf("%d calls to demo.mw:428 (demo.mw:450)\n", mexprofrecord_[31]);
        mexPrintf("%d calls to demo.mw:446\n", mexprofrecord_[32]);
        mexPrintf("%d calls to demo.mw:457 (demo.mw:482)\n", mexprofrecord_[35]);
        mexPrintf("%d calls to demo.mw:464 (demo.mw:503)\n", mexprofrecord_[36]);
        mexPrintf("%d calls to demo.mw:473\n", mexprofrecord_[37]);
        mexPrintf("%d calls to demo.mw:486 (demo.mw:508)\n", mexprofrecord_[39]);
        mexPrintf("%d calls to demo.mw:489 (demo.mw:511)\n", mexprofrecord_[40]);
        mexPrintf("%d calls to demo.mw:492 (demo.mw:514)\n", mexprofrecord_[41]);
        mexPrintf("%d calls to demo.mw:530 (demo.mw:552)\n", mexprofrecord_[47]);
        mexPrintf("%d calls to demo.mw:533 (demo.mw:555)\n", mexprofrecord_[48]);
        mexPrintf("%d calls to demo.mw:536 (demo.mw:558)\n", mexprofrecord_[49]);
        mexPrintf("%d calls to demo.mw:570\n", mexprofrecord_[54]);
        mexPrintf("%d calls to demo.mw:578\n", mexprofrecord_[55]);
    } else if (strcmp(id, "*profile log*") == 0) {
        FILE* logfp;
        if (nrhs != 2 || mxGetString(prhs[1], id, sizeof(id)) != 0)
            mexErrMsgTxt("Must have two string arguments");
        logfp = fopen(id, "w+");
        if (!logfp)
            mexErrMsgTxt("Cannot open log for output");
        if (!mexprofrecord_)
            fprintf(logfp, "Profiler inactive\n");
        fprintf(logfp, "%d calls to demo.mw:100\n", mexprofrecord_[1]);
        fprintf(logfp, "%d calls to demo.mw:102\n", mexprofrecord_[2]);
        fprintf(logfp, "%d calls to demo.mw:104\n", mexprofrecord_[3]);
        fprintf(logfp, "%d calls to demo.mw:122\n", mexprofrecord_[4]);
        fprintf(logfp, "%d calls to demo.mw:124\n", mexprofrecord_[5]);
        fprintf(logfp, "%d calls to demo.mw:126\n", mexprofrecord_[6]);
        fprintf(logfp, "%d calls to demo.mw:130 (demo.mw:187)\n", mexprofrecord_[7]);
        fprintf(logfp, "%d calls to demo.mw:157\n", mexprofrecord_[8]);
        fprintf(logfp, "%d calls to demo.mw:159\n", mexprofrecord_[9]);
        fprintf(logfp, "%d calls to demo.mw:161\n", mexprofrecord_[10]);
        fprintf(logfp, "%d calls to demo.mw:179\n", mexprofrecord_[11]);
        fprintf(logfp, "%d calls to demo.mw:181\n", mexprofrecord_[12]);
        fprintf(logfp, "%d calls to demo.mw:183\n", mexprofrecord_[13]);
        fprintf(logfp, "%d calls to demo.mw:337 (demo.mw:342)\n", mexprofrecord_[15]);
        fprintf(logfp, "%d calls to demo.mw:340\n", mexprofrecord_[16]);
        fprintf(logfp, "%d calls to demo.mw:343\n", mexprofrecord_[18]);
        fprintf(logfp, "%d calls to demo.mw:346\n", mexprofrecord_[19]);
        fprintf(logfp, "%d calls to demo.mw:358\n", mexprofrecord_[20]);
        fprintf(logfp, "%d calls to demo.mw:363\n", mexprofrecord_[21]);
        fprintf(logfp, "%d calls to demo.mw:368\n", mexprofrecord_[22]);
        fprintf(logfp, "%d calls to demo.mw:382\n", mexprofrecord_[23]);
        fprintf(logfp, "%d calls to demo.mw:387\n", mexprofrecord_[24]);
        fprintf(logfp, "%d calls to demo.mw:392\n", mexprofrecord_[25]);
        fprintf(logfp, "%d calls to demo.mw:406\n", mexprofrecord_[26]);
        fprintf(logfp, "%d calls to demo.mw:408\n", mexprofrecord_[27]);
        fprintf(logfp, "%d calls to demo.mw:410\n", mexprofrecord_[28]);
        fprintf(logfp, "%d calls to demo.mw:424\n", mexprofrecord_[29]);
        fprintf(logfp, "%d calls to demo.mw:426 (demo.mw:448)\n", mexprofrecord_[30]);
        fprintf(logfp, "%d calls to demo.mw:428 (demo.mw:450)\n", mexprofrecord_[31]);
        fprintf(logfp, "%d calls to demo.mw:446\n", mexprofrecord_[32]);
        fprintf(logfp, "%d calls to demo.mw:457 (demo.mw:482)\n", mexprofrecord_[35]);
        fprintf(logfp, "%d calls to demo.mw:464 (demo.mw:503)\n", mexprofrecord_[36]);
        fprintf(logfp, "%d calls to demo.mw:473\n", mexprofrecord_[37]);
        fprintf(logfp, "%d calls to demo.mw:486 (demo.mw:508)\n", mexprofrecord_[39]);
        fprintf(logfp, "%d calls to demo.mw:489 (demo.mw:511)\n", mexprofrecord_[40]);
        fprintf(logfp, "%d calls to demo.mw:492 (demo.mw:514)\n", mexprofrecord_[41]);
        fprintf(logfp, "%d calls to demo.mw:530 (demo.mw:552)\n", mexprofrecord_[47]);
        fprintf(logfp, "%d calls to demo.mw:533 (demo.mw:555)\n", mexprofrecord_[48]);
        fprintf(logfp, "%d calls to demo.mw:536 (demo.mw:558)\n", mexprofrecord_[49]);
        fprintf(logfp, "%d calls to demo.mw:570\n", mexprofrecord_[54]);
        fprintf(logfp, "%d calls to demo.mw:578\n", mexprofrecord_[55]);
        fclose(logfp);
    } else
        mexErrMsgTxt("Unknown identifier");
}

