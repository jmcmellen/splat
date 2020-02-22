/********************************************************************************
 * Irregular Terrain With Obstructions Model (ITWOM)                             *
 *                                                                               *
 * Version 3.0b, January 2, 2019  File: itwom3.0.c                               *
 *                                                                               *
 * Provenance:   Further test version of itwom2.0m re adj to Hrzn range factors  *
 * 1. This file is based on a thorough debugging, completion, and update of the  *
 * ITM, based on an original, public domain version of this file obtained from:  *
 * ftp://flattop.its.bldrdoc.gov/itm/ITMDLL.cpp prior to May, 2007. C++ routines *
 * for this program are taken from a translation of the FORTRAN code written by  *
 * U.S. Department of Commerce NTIA/ITS Institute for Telecommunication Sciences    *
 * Irregular Terrain Model (ITM) (Longley-Rice).                                 *
 * 2. The Linux version of this file incorporates improvements suggested by a    *
 * study of changes made to file itm.cpp by J. D. McDonald to remove Microsoft   *
 * Windows dll-isms and to debug an ambguity in overloaded calls.                *
 * 3. The Linux version of this file also incorporates improvements suggested by *
 * a study of further modifications made to itm.cpp by John A. Magliacane to     *
 * remove unused variables, unneeded #includes, and to replace pow() statements     *
 * with explicit multiplications to improve execution speed and accuracy.        *
 * 4. On August 19, 2007 this file was modified by Sid Shumate to include        *
 * changes and updates included in version 7.0 of ITMDLL.cpp, which was released *
 * by the NTIA/ITS on June 26, 2007. With correction set SS1 and SS2: itm71.cpp.    *
 * 5. On Feb. 5, 2008 this file became v.1.0 of the ITWOM with the addition, by     *
 * Sid Shumate, of multiple corrections, the replacement of subroutines lrprop   *
 * and alos with lrprop2 and alos2, and the addition of subroutine saalos to     *
 * incorporate Radiative Transfer Engine (RTE) computations in the line of sight *
 * range.                                                                        *
 * Update 8 Jun 2010 to modify alos to match 2010 series of IEEE-BTS             *
 * newsletter articles                                                           *
 * Update June 12, 2010 to z version to change test outputs                      *
 * Offshoot start date June 23, 2010 to start itwom2.0 dual version for FCC.     *
 * Update to 2.0b July 25 to correct if statement errors in adiff2 re two peak   *
 * calculations starting at line 525                                             *
 * Development to 2.0c 8 Aug 2010 after modifying saalos and adiff for full      *
 * addition of saalos treatment to post obstruction calculations and debugging.  *
 * Modified to make 1st obs loss=5.8 only, no clutter loss considered            *
 * Updated to 3.1 after much cleanup, conversion to pure C, and added           *
 * documentation by Michel Hoche-Mong.                                           *
 *                                                                               *
 * Commented out unused variables and calculations to eliminate gcc warnings     *
 *    (-Wunused-but-set-variable)  -- John A. Magliacane -- July 25, 2013        *
 ********************************************************************************/

#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#ifdef ITM_ELEV_DOUBLE
#define elev_t double
#else
#define elev_t float
#endif

/* 
 * prop_type
 *
 * Propagation Type
 *
 * This is used as the main state container throughout all these calls. It
 * is initialized in the main call to point_to_point() or area() and used
 * until those calls return.
 *
 * Heights and distances are in meters, except where noted.
 * Notes/Questions:
 *     What's the difference between rch and he? What about tgh?
 */
typedef struct prop_type
{
    double aref;     /* reference attenuation                                  */
    double dist;     /* propagation path distance (in km?) */
    double hg[2];    /* hg[0] = transmitter site tx height, AGL                */
                     /* hg[1] = receive site tx height, AGL                    */
    double rch[2];   /* rch[0] = effective height of tx antenna                */
                     /* rch[1] = effective height of rx antenna                */
    double wn;       /* wave number: freq/47.7MHz*meter; units are 1/meters    */
    double dh;       /* delta H, terrain irregularity factor                   */
    double dhd;      /* override for dh, passed in when calling into library   */
                     /* in point_to_point() or area(). If >1, this will be     */
                     /* used instead of the calculated dh.                     */
    double ens;      /* refractivity of the atmosphere at sea level            */
    double encc;     /* average clutter canopy refractivity constant. Default  */
                     /* setting is 1000 (parts per million).                   */
    double cch;      /* average height of clutter, AGL                         */
    double cd;       /* direct signal clutter exposure distance                */
    double gme;      /* effective earth curvature (actual+refraction): 1/radius */
    double zgndreal; /* resistive component of earth impedance                 */
    double zgndimag; /* reactive component of earth impedance                  */
    double he[2];    /* he[0] = effective height of tx antenna                 */
                     /* he[1] = effective height ot rx antenna                 */
    double dl[2];    /* dl[0] = tx antenna horizon or highest obstacle dist    */
                     /* dl[1] = rx antenna horizon or highest obstacle dist    */
    double the[2];   /* the[0] = take-off angle from tx ant to dl[0]           */
                     /* the[1] = take-off angle from rx ant to dl[1]           */
    double tiw;      /* incremental distance between pts in pfl[] elev array   */
    double ght;      /* tx antenna height above MSL                            */
    double ghr;      /* rx antenna height above MSL                            */
    double rph;      /* height of 2-ray reflection point                       */
    double hht;      /* height of horizon or highest obstacle visible to tx    */
    double hhr;      /* height of horizon or highest obstacle visible to rx    */
    double tgh;      /* tx antenna height, AGL                                 */
    double tsgh;     /*  the transmitter site ground height (AMSL) of either   */
                     /* the main transmitter in a call from alos2, or of an    */
                     /* obstruction peak in a call from adiff2.                */
    double thera;    /* receive approach angle, radians. used in alos()        */
    double thenr;    /* theta near receiver, slope of line from last elev pt   */
                     /* to receiver. used in alos().                           */
    double xae;      /* relationship of frequency to earth's curvature         */
    int rpl;         /* loc of 2-ray reflection point in pfl[] elev array      */
    int kwx;         /* error indicator                                        */
    int mdp;         /* propagation model:                                     */
                     /*      -1   point-to-point                               */
                     /*       0   area prediction mode                         */
                     /*       1   initializing area prediction mode            */
    int ptx;         /* polarity of transmitted signal, passed in main call    */
                     /* to point_to_point() or area().                         */
                     /*       0   horizontal                                   */
                     /*       1   vertical                                     */
                     /*       2   circular                                     */
    bool los;        /* true if line-of-sight mode                             */
    bool wlos;       /* true if line-if-sight coefficients have been calculated */
    bool wscat;      /* true if troposcatter coefficients have been calculated */
} prop_type;

/*
 * propv_type
 *
 * Propagation Variability
 *
 */
typedef struct propv_type
{
    double sgc;      /* stddev of situation variability (confidence)           */
    int mdvar;       /* mode of variability. Range 0-23. Preset to 12.         */
    int klim;        /* radio climate, set by the user from a preset list:     */
                     /* 1. Equatorial; (Africa, along the equator)             */
                     /* 2. Continental Subtropical; (Sudan region)             */
                     /* 3. Subtropical (a.k.a. Maritime Subtropical (West      */
                     /*    Coast of Africa);                                   */
                     /* 4. Desert (Death Valley, NV; Sahara);                  */
                     /* 5. Continental Temperate (usual general U.S. default); */
                     /* 6. Maritime Temperate Over Land (California to State   */
                     /*    of Washington; West Coast of Europe including U.K.) */
                     /* 7. Maritime Temperate, Over Sea.                       */
} propv_type;

/*
 * propa_type
 *
 * Propagation ...Attenuation?
 *
 */
typedef struct propa_type
{
    double dlsa;   /* the sum of the smooth-earth horizon distance           */
    double dx;     /* distance where diffraction mode gives way to scatter mode */
    double ael;
    double ak1;
    double ak2;
    double aed;
    double emd;
    double aes;
    double ems;
    double dls[2]; /* dls[0] = smooth earth horizon distance for tx antenna   */
                   /* dls[1] = smooth earth horizon distance for rx antenna   */
    double dla;    /* sum of the two horizon distances                        */
    double tha;    /* total bending angle in radians                          */
} propa_type;

#ifndef min
#define min(i, j) ( i < j ? i : j)
#endif
#ifndef max
#define max(i, j) ( i > j ? i : j)
#endif

#define FORTRAN_DIM(x, y) ( (x) > (y) ? ((x)-(y)) : (0.0) )

#define THIRD (1.0/3.0)

/* mini complex number library */
typedef struct tcomplex {
    double _real;
    double _imag;
} tcomplex;

double tcreal(tcomplex tc)
{
    return tc._real;
}

tcomplex tcadd(float r, tcomplex tc) {
    tc._real = r + tc._real;
    return tc;
}

tcomplex tcsubtract(float r, tcomplex tc) {
    tc._real = r - tc._real;
    return tc;
}
double tcimag(tcomplex tc)
{
    return tc._imag;
}

double tcarg(tcomplex tc){
    return atan2(tc._imag, tc._real);
}

double tcabs(tcomplex tc)
{
    return sqrt(tcreal(tc)*tcreal(tc) + tcimag(tc)*tcimag(tc));
}

tcomplex tcmult(tcomplex a, tcomplex b){
    tcomplex ret = { a._real*b._real - a._imag*b._imag,
        a._real*b._imag + a._imag*b._real };
    return ret;
}

tcomplex tcdiv(tcomplex a, tcomplex b){
    tcomplex ret = { (a._real*b._real + a._imag*b._imag)/(b._real*b._real + b._imag*b._imag),
        (a._imag*b._real - a._real*b._imag)/(b._real*b._real + b._imag*b._imag) };
    return ret;
}

tcomplex tcsqrt(tcomplex tc){
    tcomplex ret = { sqrt(tcabs(tc)) * cos(tcarg(tc)/2),  sqrt(tcabs(tc)) * sin(tcarg(tc)/2)};
    return ret;
}


/* 
 * Attenuation from Knife Edge diffraction
 *
 * v2: the square of the value v, where v is the "internal wedge phase angle" of an
 *     obstruction
 *
 * Returns the attenuation factor
 *
 * See ITWOM-SUB-ROUTINES.pdf, p10-11 for the definition of v2, and p405 for the analysis of
 * the function.
 *
 */
double aknfe(const double v2)
{
    double a;

    if (v2<5.76)
        a=6.02+9.11*sqrt(v2)-1.27*v2;
    else
        a=12.953+10*log10(v2);
    return a;
}

/*
 * Function Height-Gain for Three-Radii
 *
 * Calculates the height-gain over a smooth spherical earth for use when
 * calculating the rounded earth attenuation. This is described in 
 * Hufford (itm_alg.pdf, p17, formula 6.4) as an approximation to an
 * Airy function used as part of L.E. Vogler's 1964 formulation of the
 * solution to the smooth, spherical earth problem (itm_alg.pdf, p7-8)
 * and TN101, Sec 8. This used in the three-radii smooth earth analysis.
 *
 *  x:
 * pk:
 *
 * See ITWOM-SUB-ROUTINES.pdf, pp12-15 for a general description
 * See ITWOM-SUB-ROUTINES.pdf p143, for a function description
 * ...which refers to itm_alg.pdf, p15, specifically section 6.4.
 *
 */
double fht(const double x, const double pk)
{
    double w, fhtv;

    if (x<200.0)
    {
        w=-log(pk);

        if (pk<1.0e-5 || x*w*w*w > 5495.0)
        {
            fhtv=-117.0;

            if (x>1.0)
                fhtv=40.0*log10(x)+fhtv;
        }
        else
            fhtv=2.5e-5*x*x/pk-8.686*w-15.0;
    }
    else
    {
        fhtv=0.05751*x-10.0*log10(x);

        if (x<2000.0)
        {
            w=0.0134*x*exp(-0.005*x);
            fhtv=(1.0-w)*fhtv+w*(40.0*log10(x)-117.0);
        }
    }
    return fhtv;
}

/*
 * H0 Frequency Gain for Scatter Fields
 *
 *  r: twice the angular distance
 * et: scatter efficiency factor
 *
 * See ITWOM-SUB-ROUTINES.pdf, p 146
 */
double h0f(double r, double et)
{
    double a[5]={25.0, 80.0, 177.0, 395.0, 705.0};
    double b[5]={24.0, 45.0,  68.0,  80.0, 105.0};
    double q, x;
    double h0fv, temp; 
    int it;

    it=(int)et;

    if (it<=0)
    {
        it=1;
        q=0.0;
    }

    else if (it>=5)
    {
        it=5;
        q=0.0;
    }

    else
        q=et-it;

    /* x=pow(1.0/r,2.0); */

    temp=1.0/r;
    x=temp*temp;

    h0fv=4.343*log((a[it-1]*x+b[it-1])*x+1.0);

    if (q!=0.0)
        h0fv=(1.0-q)*h0fv+q*4.343*log((a[it]*x+b[it])*x+1.0);

    return h0fv;
}

/*
 * Approximate theta-d function for scatter fields
 *
 * td:
 *
 * See ITWOM-SUB_ROUTINES.pdf, p 50
 */
double ahd(double td)
{
    int i;
    double a[3]={   133.4,    104.6,     71.8};
    double b[3]={0.332e-3, 0.212e-3, 0.157e-3};
    double c[3]={  -4.343,   -1.086,    2.171};

    if (td<=10e3)
        i=0;

    else if (td<=70e3)
        i=1;

    else
        i=2;

    return a[i]+b[i]*td+c[i]*log(td);
}

double abq_alos(tcomplex r)
{
    return tcreal(r)*tcreal(r) + tcimag(r)*tcimag(r);
}

/* saalos()
 *
 * Shumate's Approximations of Attenuation in the Line Of Sight prior to
 * an obstruction.
 *
 * See ITWOM-SUB-ROUTINES.pdf, p 275
 */
double saalos(double dist, prop_type *prop)
{
    double ensa, encca, q, dp, dx, tde, hc, ucrpc, ctip, tip, tic, stic, ctic, sta;
    double ttc, cttc, crpc, ssnps, d1a, rsp, tsp, arte, zi, pd, pdk, hone, tvsr;
    double saalosv=0.0;

    q=0.0;

    if (dist==0.0)
    {
        tsp=1.0;
        rsp=0.0;
        d1a=50.0;
        saalosv=0.0;
    }
    else if(prop->hg[1] > prop->cch)
    {
        saalosv=0.0;
    }
    else
    {
        pd=dist;
        pdk=pd/1000.0;
        tsp=1.0;
        rsp=0.0;
        d1a=pd;
        /* at first, hone is transmitter antenna height 
           relative to receive site ground level. */
        hone=prop->tgh+prop->tsgh-(prop->rch[1]-prop->hg[1]);  

        if(prop->tgh>prop->cch)  /* for TX ant above all clutter height*/
        {
            ensa=1+prop->ens*0.000001;
            encca=1+prop->encc*0.000001;
            dp=pd;

            for (int j=0; j<5; ++j)
            {
                tde=dp/6378137.0;
                hc=(prop->cch+6378137.0)*(1-cos(tde));
                dx=(prop->cch+6378137.0)*sin(tde);
                ucrpc=sqrt((hone-prop->cch+hc)*(hone-prop->cch+hc)+(dx*dx));
                ctip=(hone-prop->cch+hc)/ucrpc;
                tip=acos(ctip);
                tic=tip+tde;
                tic=max(0.0,tic);
                stic=sin(tic);
                sta=(ensa/encca)*stic;
                ttc=asin(sta);
                cttc=sqrt(1-(sin(ttc))*(sin(ttc)));
                crpc=(prop->cch-prop->hg[1])/cttc;
                if(crpc>=dp)
                {
                    crpc=dp-1/dp;
                }

                ssnps=(3.1415926535897/2)-tic;
                d1a=(crpc*sin(ttc))/(1-1/6378137.0);
                dp=pd-d1a;

            }

            ctic=cos(tic);

            /* if the ucrpc path touches the canopy before reaching the
               end of the ucrpc, the entry point moves toward the
               transmitter, extending the crpc and d1a. Estimating the d1a: */

            if(ssnps<=0.0)
            {
                d1a=min(0.1*pd,600.0);
                crpc=d1a;
                /* hone must be redefined as being barely above
                   the canopy height with respect to the receiver
                   canopy height, which despite the earth curvature
                   is at or above the transmitter antenna height. */
                hone=prop->cch+1;           
                rsp=.997;
                tsp=1-rsp;
            }
            else
            {

                if (prop->ptx>=1)  /* polarity ptx is vertical or circular */
                {
                    q=((ensa*cttc-encca*ctic)/(ensa*cttc+encca*ctic));
                    rsp=q*q;
                    tsp=1-rsp;

                    if (prop->ptx==2)  /* polarity is circular - new */
                    {
                        q=((ensa*ctic-encca*cttc)/(ensa*ctic+encca*cttc));
                        rsp=((ensa*cttc-encca*ctic)/(ensa*cttc+encca*ctic));
                        rsp=(q*q+rsp*rsp)/2;
                        tsp=1-rsp;
                    }
                }
                else    /* ptx is 0, horizontal, or undefined */
                {
                    q=((ensa*ctic-encca*cttc)/(ensa*ctic+encca*cttc));
                    rsp=q*q;
                    tsp=1-rsp;
                }
            }
            /* tvsr is defined as tx ant height above receiver ant height */
            tvsr= max(0.0,prop->tgh+prop->tsgh-prop->rch[1]);  

            if (d1a<50.0)
            {
                arte=0.0195*crpc-20*log10(tsp);
            }

            else
            {
                if (d1a<225.0)
                {

                    if (tvsr>1000.0)
                    {
                        q=d1a*(0.03*exp(-0.14*pdk));
                    }
                    else 
                    {
                        q=d1a*(0.07*exp(-0.17*pdk));
                    }

                    arte=q+(0.7*pdk-max(0.01,log10(prop->wn*47.7)-2))*(prop->hg[1]/hone);
                }

                else
                {
                    q=0.00055*(pdk)+log10(pdk)*(0.041-0.0017*sqrt(hone)+0.019);

                    arte=d1a*q-(18*log10(rsp))/(exp(hone/37.5));

                    zi=1.5*sqrt(hone-prop->cch);

                    if(pdk>zi)
                    {
                        q=(pdk-zi)*10.2*((sqrt(max(0.01,log10(prop->wn*47.7)-2.0)))/(100-zi));
                    }
                    else
                    {
                        q=((zi-pdk)/zi)*(-20.0*max(0.01,log10(prop->wn*47.7)-2.0))/sqrt(hone);
                    }
                    arte=arte+q;

                }
            }    
        }
        else  /* for TX at or below clutter height */
        {
            q=(prop->cch-prop->tgh)*(2.06943-1.56184*exp(1/prop->cch-prop->tgh));
            q=q+(17.98-0.84224*(prop->cch-prop->tgh))*exp(-0.00000061*pd);
            arte=q+1.34795*20*log10(pd+1.0);
            arte=arte-(max(0.01,log10(prop->wn*47.7)-2))*(prop->hg[1]/prop->tgh);
        }
        saalosv=arte;
    }
    return saalosv;
}

typedef struct adiff_state
{
    double wd1;     /* 1/wd, the inverse of the weighting factor wd */
    double xd1;     /* dla added to a curvature adjustment tha/gme  */
    double afo;     /* attenuation from absorption/scattering from oxygen, water vapor, and terrain clutter */
    double qk;
    double aht;
    double xht;
} adiff_state;

/* Attenuation from Diffraction, Initialization
 *
 * This is intended to be called before adiff() is called. Replaces the old
 * setup that was done by calling adiff() with a distance of 0.
 *
 * See ITWOM-SUB-ROUTINES.pdf, p2
 */
void adiff_init(adiff_state *state, prop_type *prop, propa_type *propa)
{
    tcomplex prop_zgnd = {prop->zgndreal, prop->zgndimag};
    double a, q, pk, wa;

    q=prop->hg[0]*prop->hg[1];
    state->qk=prop->he[0]*prop->he[1]-q;

    if (prop->mdp<0.0)
        q+=10.0;

    state->wd1=sqrt(1.0+state->qk/q);
    state->xd1=propa->dla+propa->tha/prop->gme;
    q=(1.0-0.8*exp(-propa->dlsa/50e3))*prop->dh;
    q*=0.78*exp(-pow(q/16.0,0.25));
    state->afo=min(15.0,2.171*log(1.0+4.77e-4*prop->hg[0]*prop->hg[1]*prop->wn*q));
    state->qk=1.0/tcabs(prop_zgnd);
    state->aht=20.0;
    state->xht=0.0;

    for (int j=0; j<2; ++j)
    {
        /* a=0.5*pow(prop->dl[j],2.0)/prop->he[j]; */
        a=0.5*(prop->dl[j]*prop->dl[j])/prop->he[j];
        wa=pow(a*prop->wn,THIRD);
        pk=state->qk/wa;
        q=(1.607-pk)*151.0*wa*prop->dl[j]/a;
        state->xht+=q;
        state->aht+=fht(q,pk);
    }
}

/* Attenuation from Diffraction
 *
 * Finds the diffraction attenuation for the path distance d.
 *
 * Note that before this can be called, adiff_init() should be called to
 * set up adiff_state.
 *
 * See ITWOM-SUB-ROUTINES.pdf, p2
 */
double adiff(adiff_state *state, double dist, prop_type *prop, propa_type *propa)
{
    double a, q, pk, ds, th, wa, ar, wd, adiffv;

    th=propa->tha+dist*prop->gme;
    ds=dist-propa->dla;
    /* q=0.0795775*prop->wn*ds*pow(th,2.0); */
    q=0.0795775*prop->wn*ds*th*th;
    adiffv=aknfe(q*prop->dl[0]/(ds+prop->dl[0]))+aknfe(q*prop->dl[1]/(ds+prop->dl[1]));
    a=ds/th;
    wa=pow(a*prop->wn,THIRD);
    pk=state->qk/wa;
    q=(1.607-pk)*151.0*wa*th+state->xht;
    ar=0.05751*q-4.343*log(q)-state->aht;
    q=(state->wd1+state->xd1/dist)*min(((1.0-0.8*exp(-dist/50e3))*prop->dh*prop->wn),6283.2);
    wd=25.1/(25.1+sqrt(q));
    adiffv=ar*wd+(1.0-wd)*adiffv+state->afo;

    return adiffv;
}

typedef struct adiff2_state {
    double wd1;
    double xd1;
    double qk;
    double aht;
    double xht;
} adiff2_state;

/* Attenuation from Diffraction Version 2, Initialization
 *
 * This is intended to be called before adiff2() is called. Replaces the old
 * setup that was done by calling adiff2() with a distance of 0.
 *
 * See ITWOM-SUB-ROUTINES.pdf, p10
 */
void adiff2_init(adiff2_state *state, prop_type *prop, propa_type *propa)
{
    tcomplex prop_zgnd = { prop->zgndreal, prop->zgndimag };
    double a, q, pk, wa;

    q=prop->hg[0]*prop->hg[1];
    state->qk=prop->he[0]*prop->he[1]-q;
    /* dhec=2.73; */

    if (prop->mdp<0.0)
        q+=10.0;

    /* coefficients for a standard four radii, rounded earth computation are prepared */
    state->wd1=sqrt(1.0+state->qk/q);
    state->xd1=propa->dla+propa->tha/prop->gme;
    q=(1.0-0.8*exp(-propa->dlsa/50e3))*prop->dh;
    q*=0.78*exp(-pow(q/16.0,0.25));
    state->qk=1.0/tcabs(prop_zgnd);
    state->aht=20.0;
    state->xht=0.0;
    a=0.5*(prop->dl[0]*prop->dl[0])/prop->he[0];
    wa=pow(a*prop->wn,THIRD);
    pk=state->qk/wa;
    q=(1.607-pk)*151.0*wa*prop->dl[0]/a;
    state->xht=q;
    state->aht+=fht(q,pk);


    if ((trunc(prop->dl[1])==0.0) || (prop->the[1]>0.2))
    { 
        state->xht+=state->xht;
        state->aht+=(state->aht-20.0);
    }
    else
    {
        a=0.5*(prop->dl[1]*prop->dl[1])/prop->he[1];
        wa=pow(a*prop->wn,THIRD);
        pk=state->qk/wa;
        q=(1.607-pk)*151.0*wa*prop->dl[1]/a;
        state->xht+=q;
        state->aht+=fht(q,pk);
    }
}

/* Attenuation from Diffraction, Version 2
 *
 * Finds the diffraction attenuation for the path distance d.
 *
 * Note that before this can be called, adiff2_init() should be called to
 * set up adiff2_state.
 *
 * See ITWOM-SUB-ROUTINES.pdf, p10
 */
double adiff2(adiff2_state *state, double dist, prop_type *prop, propa_type *propa)
{
    double toh, toho, roh, roho;
    double dtr, dtro, dto, dto1, dtof, dto1f;
    double dro, dro2, drof, dro2f;
    double drto;
    double dhh1, dhh2;
    double q, rd, ds, dsl, sf2, vv;
    double kedr=0.0, arp=0.0, sdr=0.0, pd=0.0, srp=0.0, kem=0.0, csd=0.0, sdl=0.0, adiffv2=0.0, closs=0.0;

    /* sf1=1.0; */ /* average empirical hilltop foliage scatter factor for 1 obstruction  */
    sf2=1.0;  /* average empirical hilltop foliage scatter factor for 2 obstructions */

    dsl=max(dist-propa->dla,0.0);
    ds=dist-propa->dla;
    toh=prop->hht-(prop->rch[0]-prop->dl[0]*((prop->rch[1]-prop->rch[0])/prop->dist));
    roh=prop->hhr-(prop->rch[0]-(prop->dist-prop->dl[1])*((prop->rch[1]-prop->rch[0])/prop->dist));
    toho=prop->hht-(prop->rch[0]-(prop->dl[0]+dsl)*((prop->hhr-prop->rch[0])/(prop->dist-prop->dl[1])));
    roho=prop->hhr-(prop->hht-dsl*((prop->rch[1]-prop->hht)/dsl));
    dto=sqrt(prop->dl[0]*prop->dl[0]+toh*toh);
    dto+=prop->gme*prop->dl[0];
    dto1=sqrt(prop->dl[0]*prop->dl[0]+toho*toho);
    dto1+=prop->gme*prop->dl[0];
    dtro=sqrt((prop->dl[0]+dsl)*(prop->dl[0]+dsl)+prop->hhr*prop->hhr);
    dtro+=prop->gme*(prop->dl[0]+dsl);
    drto=sqrt((prop->dl[1]+dsl)*(prop->dl[1]+dsl)+prop->hht*prop->hht);
    drto+=prop->gme*(prop->dl[1]+dsl);
    dro=sqrt(prop->dl[1]*prop->dl[1]+roh*roh);
    dro+=prop->gme*(prop->dl[1]);
    dro2=sqrt(prop->dl[1]*prop->dl[1]+roho*roho);
    dro2+=prop->gme*(prop->dl[1]);
    dtr=sqrt(prop->dist*prop->dist+(prop->rch[0]-prop->rch[1])*(prop->rch[0]-prop->rch[1]));
    dtr+=prop->gme*prop->dist;
    dhh1=sqrt((prop->dist-propa->dla)*(prop->dist-propa->dla)+toho*toho);
    dhh1+=prop->gme*(prop->dist-propa->dla);
    dhh2=sqrt((prop->dist-propa->dla)*(prop->dist-propa->dla)+roho*roho);
    dhh2+=prop->gme*(prop->dist-propa->dla);

    /* for 1 obst tree base path */
    dtof=sqrt(prop->dl[0]*prop->dl[0]+(toh-prop->cch)*(toh-prop->cch));
    dtof+=prop->gme*prop->dl[0];
    dto1f=sqrt(prop->dl[0]*prop->dl[0]+(toho-prop->cch)*(toho-prop->cch));
    dto1f+=prop->gme*prop->dl[0];
    drof=sqrt(prop->dl[1]*prop->dl[1]+(roh-prop->cch)*(roh-prop->cch));
    drof+=prop->gme*(prop->dl[1]);
    dro2f=sqrt(prop->dl[1]*prop->dl[1]+(roho-prop->cch)*(roho-prop->cch));
    dro2f+=prop->gme*(prop->dl[1]);

    /* saalos coefficients preset for post-obstacle receive path */
    prop->tgh=prop->cch+1.0;
    prop->tsgh=prop->hhr;
    rd=prop->dl[1];

    /* two obstacle diffraction calculation */ 
    if (trunc(ds)>0)   /* there are 2 obstacles */ 
    {
        if(trunc(prop->dl[1])>0.0) /* receive site past 2nd peak */
        {
            /* knife edge vs round weighting */
            q=(1.0-0.8*exp(-dist/50e3))*prop->dh;
            q=(state->wd1+state->xd1/dist)*min((q*prop->wn),6283.2);
            /* wd=25.1/(25.1+sqrt(q)); */

            q=0.6365*prop->wn;

            if(prop->the[1]<0.2)  /* receive grazing angle below 0.2 rad */
            {
                /* knife edge attenuation for two obstructions */

                if(prop->hht < 3400)  /* if below tree line, foliage top loss */
                {
                    vv=q*fabs(dto1+dhh1-dtro);
                    adiffv2=-18.0+sf2*aknfe(vv);
                }
                else
                {
                    vv=q*fabs(dto1+dhh1-dtro);
                    adiffv2=aknfe(vv);
                }

                if(prop->hhr < 3400)
                {
                    vv=q*fabs(dro2+dhh2-drto);
                    adiffv2+=(-18.0+sf2*aknfe(vv));
                }
                else
                {
                    vv=q*fabs(dro2+dhh2-drto);
                    adiffv2+=aknfe(vv);
                }
                /* finally, add clutter loss */
                closs=saalos(rd, prop);  
                adiffv2+=min(22.0,closs); 

            }
            else     /* rcvr site too close to 2nd obs */
            {
                /* knife edge attenuation for 1st obs */

                if(prop->hht < 3400)
                {
                    vv=q*fabs(dto1+dhh1-dtro);
                    adiffv2=-18.0+sf2*aknfe(vv);
                }
                else
                {
                    vv=q*fabs(dto1+dhh1-dtro);
                    adiffv2=aknfe(vv);
                }

                /* weighted calc. of knife vs rounded edge 
                   adiffv2=ar*wd+(1.0-wd)*adiffv2; */

                /* clutter path loss past 2nd peak */
                if(prop->the[1]<1.22)
                {
                    rd=prop->dl[1];

                    if(prop->the[1]>0.6) /* through foliage downhill */
                    {
                        prop->tgh=prop->cch;
                    }
                    else    /* close to foliage, rcvr in foliage downslope */
                    {
                        vv=0.6365*prop->wn*fabs(dro2+dhh2-drto);
                    }
                    adiffv2+=aknfe(vv);
                    closs=saalos(rd, prop);
                    adiffv2+=min(closs,22.0); 
                }
                else    /* rcvr very close to bare cliff or skyscraper */
                {
                    adiffv2=5.8+25.0;
                }
            }
        }
        else /* receive site is atop a 2nd peak */
        {
            vv=0.6365*prop->wn*fabs(dto+dro-dtr);
            adiffv2=5.8 + aknfe(vv);
        }
    }
    else /* for single obstacle */
    {

        if(trunc(prop->dl[1])>0.0) /* receive site past 1st peak */
        {

            if(prop->the[1]<0.2) /* receive grazing angle less than .2 radians */
            {
                vv=0.6365*prop->wn*fabs(dto+dro-dtr);

                if(prop->hht < 3400)
                {
                    sdl=18.0;
                    sdl=pow(10,(-sdl/20));
                    /* ke phase difference with respect to direct t-r line */
                    kedr=0.159155*prop->wn*fabs(dto+dro-dtr);
                    arp=fabs(kedr-(trunc(kedr)));
                    kem=aknfe(vv);
                    kem= pow(10,(-kem/20));
                    /* scatter path phase with respect to direct t-r line */
                    sdr=0.5+0.159155*prop->wn*fabs(dtof+drof-dtr);
                    srp=fabs(sdr-(trunc(sdr)));
                    /* difference between scatter and ke phase in radians */
                    pd=6.283185307*fabs(srp-arp);
                    /* report pd prior to restriction 
                       keep pd between 0 and pi radians and adjust for 3&4 quadrant */ 
                    if(pd>=3.141592654)
                    {
                        pd=6.283185307-pd;
                        /* csd=abq_alos(complex<double>(sdl,0)+complex<double>(kem*-cos(pd), kem*-sin(pd)));  */
                        tcomplex val = { sdl + kem*-cos(pd), kem*-sin(pd) };
                        csd=abq_alos(val);
                    } 
                    else
                    {
                        /* csd=abq_alos(complex<double>(sdl,0)+complex<double>(kem*cos(pd), kem*sin(pd))); */
                        tcomplex val = { sdl + kem*cos(pd), kem*sin(pd) };
                        csd=abq_alos(val);
                    }
                    /*csd=max(csd,0.0009); limits maximum loss value to 30.45 db */
                    adiffv2=-3.71-10*log10(csd);
                }
                else
                {
                    adiffv2=aknfe(vv);
                }
                /* finally, add clutter loss */
                closs=saalos(rd, prop);
                adiffv2+=min(closs,22.0);   
            }
            else    /* receive grazing angle too high */
            {
                if(prop->the[1]<1.22)
                {
                    rd=prop->dl[1];

                    if(prop->the[1]>0.6)  /* through foliage downhill */
                    {
                        prop->tgh=prop->cch;
                    }
                    else    /* downhill slope just above foliage  */
                    {
                        vv=0.6365*prop->wn*fabs(dto+dro-dtr);
                        adiffv2=aknfe(vv);
                    }
                    closs=saalos(rd, prop);
                    adiffv2+=min(22.0,closs); 
                }
                else    /* receiver very close to bare cliff or skyscraper */
                {
                    adiffv2=5.8+25.0;
                }
            }
        }
        else    /* if occurs, receive site atop first peak  */
        {
            adiffv2=5.8;
        }
    }

    return adiffv2;
}

typedef struct ascat_state {
    double ad;
    double rr;
    double etq;
    double h0s;
} ascat_state;

/* Attenuation From Scatter
 *
 * This is intended to be called before ascat() is called. Replaces
 * the old setup that was done by calling ascat() with a distance of 0.
 *
 * See ITWOM-SUB-ROUTINES.pdf, p 81
 */
void ascat_init(ascat_state *state, prop_type *prop)
{
    state->ad=prop->dl[0]-prop->dl[1];
    state->rr=prop->he[1]/prop->rch[0];

    if (state->ad<0.0)
    {
        state->ad=-state->ad;
        state->rr=1.0/state->rr;
    }

    state->etq=(5.67e-6*prop->ens-2.32e-3)*prop->ens+0.031;
    state->h0s=-15.0;
}

/* Attenuation From Scatter
 *
 * Finds the "scatter attenuation" for the path distance d.
 *
 * Note that before this can be called, ascat_init() should be called to set
 * up ascat_state.
 *
 * See ITWOM-SUB-ROUTINES.pdf, p81
 */
double ascat(ascat_state *state, double dist, prop_type *prop, propa_type *propa)
{
    double h0, r1, r2, z0, ss, et, ett, th, q;
    double ascatv, temp;

    if (state->h0s>15.0)
        h0=state->h0s;
    else
    {
        th=prop->the[0]+prop->the[1]+dist*prop->gme;
        r2=2.0*prop->wn*th;
        r1=r2*prop->he[0];
        r2*=prop->he[1];

        if (r1<0.2 && r2<0.2)
            return 1001.0;

        ss=(dist-(state->ad))/(dist+(state->ad));
        q=state->rr/ss;
        ss=max(0.1,ss);
        q=min(max(0.1,q),10.0);
        z0=(dist-(state->ad))*(dist+(state->ad))*th*0.25/dist;
        /* et=(etq*exp(-pow(min(1.7,z0/8.0e3),6.0))+1.0)*z0/1.7556e3; */

        temp=min(1.7,z0/8.0e3);
        temp=temp*temp*temp*temp*temp*temp;
        et=(state->etq*exp(-temp)+1.0)*z0/1.7556e3;

        ett=max(et,1.0);
        h0=(h0f(r1,ett)+h0f(r2,ett))*0.5;
        h0+=min(h0,(1.38-log(ett))*log(ss)*log(q)*0.49);
        h0=FORTRAN_DIM(h0,0.0);

        if (et<1.0)
        {
            /* h0=et*h0+(1.0-et)*4.343*log(pow((1.0+1.4142/r1)*(1.0+1.4142/r2),2.0)*(r1+r2)/(r1+r2+2.8284)); */

            temp=((1.0+1.4142/r1)*(1.0+1.4142/r2));
            h0=et*h0+(1.0-et)*4.343*log((temp*temp)*(r1+r2)/(r1+r2+2.8284));
        }

        if (h0>15.0 && state->h0s>=0.0)
            h0=state->h0s;
    }

    state->h0s=h0;
    th=propa->tha+dist*prop->gme;
    /* ascatv=ahd(th*d)+4.343*log(47.7*prop->wn*pow(th,4.0))-0.1*(prop->ens-301.0)*exp(-th*d/40e3)+h0; */
    ascatv=ahd(th*dist)+4.343*log(47.7*prop->wn*(th*th*th*th))-0.1*(prop->ens-301.0)*exp(-th*dist/40e3)+h0;

    return ascatv;
}

/*
 * qerfi()
 *
 * Inverse of qerf(), the standard normal complementary probability.
 *
 * See ITWOM-SUB-ROUTINES, p466
 */
double qerfi(double q)
{
    double x, t, v;
    double c0=2.515516698;
    double c1=0.802853;
    double c2=0.010328;
    double d1=1.432788;
    double d2=0.189269;
    double d3=0.001308;

    x=0.5-q;
    t=max(0.5-fabs(x),0.000001);
    t=sqrt(-2.0*log(t));
    v=t-((c2*t+c1)*t+c0)/(((d3*t+d2)*t+d1)*t+1.0);

    if (x<0.0)
        v=-v;

    return v;
}

/*
 * Quick Longley-Rice Profile Setup
 *
 * Set up function for qlrpfl and qlrpfl2.
 *
 * IN:
 * fmhz: frequency at which to do the study, in mhz. Passed in to point_to_point() or area()
 * zsys: average elevation of the study, usually calculated by averaging the elevation array
 *  en0: atmospheric refractivity at surface level, passed in to point_to_point()
 * ipol: polarity (0 or 1). passed in to point_to_point()
 *  eps: polarization constant. passed in to point_to_point()
 *  sgm: ground constant. passed in to point_to_point()
 *
 * IN/OUT:
 * prop: filled in by various calculations on these values.
 * 
 * Note that nothing here is particularly compute-intensive to calculate. However the values
 * usually stay the same over a study. A (minor) performance enhancement would be to calculate
 * them once and save them for all studies of that tx point.
 *
 * The exception is zsys, which depends on the path.
 *
 *
 * See ITWOM-SUB-ROUTINES.pdf, p262
 */
void qlrps(double fmhz, double zsys, double en0, int ipol, double eps, double sgm, prop_type *prop)
{
    double gma=157e-9;

    prop->wn=fmhz/47.7;
    prop->ens=en0;

    if (zsys!=0.0)
        prop->ens*=exp(-zsys/9460.0);

    prop->gme=gma*(1.0-0.04665*exp(prop->ens/179.3));

    /*complex<double> zq, prop_zgnd(prop->zgndreal,prop->zgndimag); */
    /*zq=complex<double> (eps,376.62*sgm/prop->wn); */
    /*prop_zgnd=sqrt(zq-1.0); */
    tcomplex prop_zgnd = { prop->zgndreal, prop->zgndimag };
    tcomplex zq = { eps, 376.62*sgm/prop->wn };
    prop_zgnd=tcsqrt(tcadd(-1.0, zq));

    if (ipol!=0.0) {
        /* prop_zgnd=prop_zgnd/zq; */
        prop_zgnd=tcdiv(prop_zgnd, zq);
    }

    prop->zgndreal=tcreal(prop_zgnd);
    prop->zgndimag=tcimag(prop_zgnd);
}


typedef struct alos_state {
    double wls;
} alos_state;

/*
 * Attenuation, Line of Sight, Initialization
 *
 * Set up for alos_state. Call this before calling alos()
 */
void alos_init(alos_state *state, prop_type *prop, propa_type *propa)
{
    state->wls=0.021/(0.021+prop->wn*prop->dh/max(10e3,propa->dlsa));
}

/*
 * Attenuation, Line Of Sight
 *
 * See ITWOM-SUB-ROUTINES, p54
 */
double alos(alos_state *state, double dist, prop_type *prop, propa_type *propa)
{
    tcomplex prop_zgnd = { prop->zgndreal, prop->zgndimag };
    tcomplex r = { 0.0, 0.0 };
    double s, sps, q;
    double alosv;

    q=(1.0-0.8*exp(-dist/50e3))*prop->dh;
    s=0.78*q*exp(-pow(q/16.0,0.25));
    q=prop->he[0]+prop->he[1];
    sps=q/sqrt(dist*dist+q*q);
    /* r=(sps-prop_zgnd)/(sps+prop_zgnd)*exp(-min(10.0,prop->wn*s*sps)); */
    double numerator = sps-tcreal(prop_zgnd);
    double denominator = sps+tcreal(prop_zgnd);
    r._real = (numerator / denominator) * exp(-min(10.0,prop->wn*s*sps));
    q=abq_alos(r);

    if (q<0.25 || q<sps) {
        /* r=r*sqrt(sps/q); */
        r._real = r._real * sqrt(sps/q);
    }

    alosv=propa->emd*dist+propa->aed;
    q=prop->wn*prop->he[0]*prop->he[1]*2.0/dist;

    if (q>1.57)
        q=3.14-2.4649/q;

    /*alosv=(-4.343*log(abq_alos(complex<double>(cos(q),-sin(q))+r))-alosv)*wls+alosv; */
    tcomplex val = { cos(q) , -sin(q) };
    val._real += r._real;
    val._imag += r._imag;
    alosv=(-4.343*log(abq_alos(val))-alosv)*(state->wls)+alosv;

    return alosv;
}

/*
 * Attenuation, Line Of Sight, Version 2
 *
 * See ITWOM-SUB-ROUTINES.pdf, p58
 */
double alos2(prop_type *prop, propa_type *propa)
{
    tcomplex prop_zgnd = { prop->zgndreal, prop->zgndimag };
    tcomplex r = { 0.0, 0.0 };
    double cd, cr, dr, hr, hrg, ht, htg, hrp, re, s, sps, q, pd, drh;
    /* int rp; */
    double alosv;

    cd=0.0;
    cr=0.0;
    htg=prop->hg[0];
    hrg=prop->hg[1];
    ht=prop->ght;
    hr=prop->ghr;
    /* rp=prop->rpl; */
    hrp=prop->rph;
    pd=prop->dist;

    q=prop->he[0]+prop->he[1];
    sps=q/sqrt(pd*pd+q*q);
    q=(1.0-0.8*exp(-pd/50e3))*prop->dh;

    if (prop->mdp<0)
    {
        dr=pd/(1+hrg/htg);

        if (dr<(0.5*pd))
        {
            drh=6378137.0-sqrt(-(0.5*pd)*(0.5*pd)+6378137.0*6378137.0+(0.5*pd-dr)*(0.5*pd-dr));  
        }
        else
        {
            drh=6378137.0-sqrt(-(0.5*pd)*(0.5*pd)+6378137.0*6378137.0+(dr-0.5*pd)*(dr-0.5*pd));  
        }

        if ((sps<0.05) && (prop->cch>hrg) && (prop->dist< prop->dl[0])) /* if far from transmitter and receiver below canopy */  
        {
            cd=max(0.01,pd*(prop->cch-hrg)/(htg-hrg));
            cr=max(0.01,pd-dr+dr*(prop->cch-drh)/htg);
            q=((1.0-0.8*exp(-pd/50e3))*prop->dh*(min(-20*log10(cd/cr),1.0)));
        }
    }

    s=0.78*q*exp(-pow(q/16.0,0.25));
    q=exp(-min(10.0,prop->wn*s*sps));
    /*r=q*(sps-prop_zgnd)/(sps+prop_zgnd); */
    r._real=q*(sps-prop_zgnd._real)/(sps+prop_zgnd._real);
    q=abq_alos(r);
    q=min(q,1.0);

    if (q<0.25 || q<sps)
    {
        /*r=r*sqrt(sps/q); */
        r._real=r._real*sqrt(sps/q);
    }
    q=prop->wn*prop->he[0]*prop->he[1]/(pd*3.1415926535897);

    if (prop->mdp<0)
    {
        q=prop->wn*((ht-hrp)*(hr-hrp))/(pd*3.1415926535897);
    }
    q-=floor(q);

    if (q<0.5)
    {
        q*=3.1415926535897;
    }

    else
    {
        q=(1-q)*3.1415926535897;
    }
    /* no longer valid complex conjugate removed 
       by removing minus sign from in front of sin function */
    /* re=abq_alos(complex<double>(cos(q),sin(q))+r); */
    tcomplex val = { cos(q), sin(q) };
    val._real += r._real;
    val._imag += r._imag;
    re=abq_alos(val);
    alosv=-10*log10(re);
    prop->tgh=prop->hg[0];  /*tx above gnd hgt set to antenna height AGL */
    prop->tsgh=prop->rch[0]-prop->hg[0]; /* tsgh set to tx site gl AMSL */

    if ((prop->hg[1]<prop->cch) && (prop->thera<0.785) && (prop->thenr<0.785))
    {
        if (sps<0.05)
        {
            alosv=alosv+saalos(pd, prop);
        }
        else
        {
            alosv=saalos(pd, prop);
        }
    }

    alosv = min(22.0,alosv);
    return alosv;
}

/* Quick Longley-Rice Area, Initialization
 *
 * Parameter preparation subroutine for area mode
 *
 * See ITWOM-SUB-ROUTINES.pdf p468
 *
 * Used only with ITM 1.2.2
 */
void qlra(int kst[], int klimx, int mdvarx, prop_type *prop, propv_type *propv)
{
    double q;

    for (int j=0; j<2; ++j)
    {
        if (kst[j]<=0)
            prop->he[j]=prop->hg[j];
        else
        {
            q=4.0;

            if (kst[j]!=1)
                q=9.0;

            if (prop->hg[j]<5.0)
                q*=sin(0.3141593*prop->hg[j]);

            prop->he[j]=prop->hg[j]+(1.0+q)*exp(-min(20.0,2.0*prop->hg[j]/max(1e-3,prop->dh)));
        }

        q=sqrt(2.0*prop->he[j]/prop->gme);
        prop->dl[j]=q*exp(-0.07*sqrt(prop->dh/max(prop->he[j],5.0)));
        prop->the[j]=(0.65*prop->dh*(q/prop->dl[j]-1.0)-2.0*prop->he[j])/q;
    }

    prop->mdp=1;

    if (mdvarx>=0)
    {
        propv->mdvar=mdvarx;
    }

    if (klimx>0)
    {
        propv->klim=klimx;
    }
}


/*
 * Longley-Rice Propagation
 *
 * This is the PaulM version of lrprop.
 *
 * This was replaced by lrprop2 because no provision was made to accommodate the changes
 * required in the computational procedure that must occur when the line-of-sight and
 * diffraction ranges no longer gently merge, but are, instead, abruptly separated by a major
 * obstruction.
 *
 * See ITWOM-SUB-ROUTINES.pdf, page 185.
 */
void lrprop(double dist, prop_type *prop, propa_type *propa)
{
    /* PaulM_lrprop used for ITM */
    tcomplex prop_zgnd = {prop->zgndreal, prop->zgndimag};
    adiff_state state = {0};
    double a0, a1, a2, a3, a4, a5, a6;
    double d0, d1, d2, d3, d4, d5, d6;
    bool wq;
    double q;
    int j;

    /* if the propagation model is not in area mode we're either initializing
     * area mode or we're running one-shot in point-to-point mode. Either way,
     * set up some coefficients in prop. */
    if (prop->mdp!=0)
    {
        for (j=0; j<2; j++)
            propa->dls[j]=sqrt(2.0*prop->he[j]/prop->gme);

        propa->dlsa=propa->dls[0]+propa->dls[1];
        propa->dla=prop->dl[0]+prop->dl[1];
        propa->tha=max(prop->the[0]+prop->the[1],-propa->dla*prop->gme);
        prop->wlos=false;
        prop->wscat=false;

        /*checking for parameters-in-range, error codes set if not */

        /* Warning: wave number */
        if (prop->wn<0.838 || prop->wn>210.0)
            prop->kwx=max(prop->kwx,1);

        /* Warning: radiation heights, AGL */
        for (j=0; j<2; j++)
            if (prop->hg[j]<1.0 || prop->hg[j]>1000.0)
                prop->kwx=max(prop->kwx,1);

        /* Warning: tx take-off angle, antenna horizon < 1/10 of horizon distance */
        for (j=0; j<2; j++)
            if (fabs(prop->the[j]) >200e-3 || prop->dl[j]<0.1*propa->dls[j] || prop->dl[j]>3.0*propa->dls[j] )
                prop->kwx=max(prop->kwx,3);

        /* Out-of-Range: refractivity, earth curvature, earth impedance/reactance, wave number */
        if (prop->ens < 250.0 || prop->ens > 400.0 || prop->gme < 75e-9 || prop->gme > 250e-9 || (tcreal(prop_zgnd) <= fabs(tcimag(prop_zgnd))) || prop->wn < 0.419 || prop->wn > 420.0)
            prop->kwx=4;

        /* Out-of-Range: radiation heights, AGL */
        for (j=0; j<2; j++)
            if (prop->hg[j]<0.5 || prop->hg[j]>3000.0)
                prop->kwx=4;


        adiff_init(&state, prop, propa);

        /* xae=pow(prop->wn*pow(prop->gme,2.),-THIRD); -- JDM made argument 2 a double */
        prop->xae=pow(prop->wn*(prop->gme*prop->gme),-THIRD);  /* No 2nd pow() */
        d3=max(propa->dlsa,1.3787*prop->xae+propa->dla);
        d4=d3+2.7574*prop->xae;
        a3=adiff(&state, d3,prop,propa);  // use
        a4=adiff(&state, d4,prop,propa);  // use
        propa->emd=(a4-a3)/(d4-d3);
        propa->aed=a3-propa->emd*d3;
    }

    /* if the propagation model is in "initializing area mode", reset it to normal running */
    if (prop->mdp>=0)
    {
        prop->mdp=0;
        prop->dist=dist;
    }

    if (prop->dist>0.0)
    {
        if (prop->dist>1000e3)
            prop->kwx=max(prop->kwx,1);

        double dmin=fabs(prop->he[0]-prop->he[1])/200e-3; /* minimum acceptable path distance length (meters) */
        if (prop->dist<dmin)
            prop->kwx=max(prop->kwx,3);

        if (prop->dist<1e3 || prop->dist>2000e3)
            prop->kwx=4;
    }

    /* if the distance is less than the smooth earth horizon distance, we
     * do some line-of-site setup calculations. Otherwise we skip that
     * and go to scattering. */
    if (prop->dist < propa->dlsa)
    {
        if (!prop->wlos)
        {
            alos_state alosstate = {0};
            alos_init(&alosstate, prop, propa);

            d2=propa->dlsa;
            a2=propa->aed+d2*propa->emd;
            d0=1.908*prop->wn*prop->he[0]*prop->he[1];

            if (propa->aed>=0.0)
            {
                d0=min(d0,0.5*propa->dla);
                d1=d0+0.25*(propa->dla-d0);
            }

            else
                d1=max(-propa->aed/propa->emd,0.25*propa->dla);

            a1=alos(&alosstate, d1,prop,propa);
            wq=false;

            if (d0<d1)
            {
                a0=alos(&alosstate, d0,prop,propa);
                q=log(d2/d0);
                propa->ak2=max(0.0,((d2-d0)*(a1-a0)-(d1-d0)*(a2-a0))/((d2-d0)*log(d1/d0)-(d1-d0)*q));
                wq=propa->aed>=0.0 || propa->ak2>0.0;

                if (wq)
                { 
                    propa->ak1=(a2-a0-propa->ak2*q)/(d2-d0);

                    if (propa->ak1<0.0)
                    {
                        propa->ak1=0.0;
                        propa->ak2=FORTRAN_DIM(a2,a0)/q;

                        if (propa->ak2==0.0)
                            propa->ak1=propa->emd;
                    }
                }

                else
                {
                    propa->ak2=0.0;
                    propa->ak1=(a2-a1)/(d2-d1);

                    if (propa->ak1<=0.0)
                        propa->ak1=propa->emd;
                }
            }

            else
            {
                propa->ak1=(a2-a1)/(d2-d1);
                propa->ak2=0.0;

                if (propa->ak1<=0.0)
                    propa->ak1=propa->emd;
            }

            propa->ael=a2-propa->ak1*d2-propa->ak2*log(d2);
            prop->wlos=true;
        }

        if (prop->dist>0.0)
            prop->aref=propa->ael+propa->ak1*prop->dist+propa->ak2*log(prop->dist);

    }

    if (prop->dist <= 0.0 || prop->dist>=propa->dlsa)
    {
        if (!prop->wscat)
        { 
            ascat_state scatterstate = {0};
            ascat_init(&scatterstate, prop);

            d5=propa->dla+200e3;
            d6=d5+200e3;
            a6=ascat(&scatterstate, d6,prop,propa);
            a5=ascat(&scatterstate, d5,prop,propa);

            if (a5<1000.0)
            {
                propa->ems=(a6-a5)/200e3;
                propa->dx=max(propa->dlsa,max(propa->dla+0.3*prop->xae*log(47.7*prop->wn),(a5-propa->aed-propa->ems*d5)/(propa->emd-propa->ems)));
                propa->aes=(propa->emd-propa->ems)*propa->dx+propa->aed;
            }

            else
            {
                propa->ems=propa->emd;
                propa->aes=propa->aed;
                propa->dx=10.e6;
            }

            prop->wscat=true;
        }

        if (prop->dist>propa->dx)
            prop->aref=propa->aes+propa->ems*prop->dist;
        else
            prop->aref=propa->aed+propa->emd*prop->dist;
    }

    prop->aref=max(prop->aref,0.0);
}


/*
 * Longley-Rice Propagation, the updated version.
 *
 * See ITWOM-SUB-ROUTINES.pdf, page 186.
 */
void lrprop2(double dist, prop_type *prop, propa_type *propa)
{
    /* ITWOM_lrprop2 */
    tcomplex prop_zgnd = {prop->zgndreal, prop->zgndimag };
    adiff2_state state = {0};
    double pd1;
    double a0, a1, a2, a3, a4, a5, a6, iw;
    double d0, d1, d2, d3, d4, d5, d6;
    bool wq;
    double q;
    int j;

    iw=prop->tiw;
    pd1=prop->dist;
    propa->dx=2000000.0;

    /* if the propagation model is not in area mode we're either initializing
     * area mode or we're running one-shot in point-to-point mode. Either way,
     * set up some coefficients in prop. */
    if (prop->mdp!=0)
    {
        for (j=0; j<2; j++)
            propa->dls[j]=sqrt(2.0*prop->he[j]/prop->gme);

        propa->dlsa=propa->dls[0]+propa->dls[1];
        propa->dlsa=min(propa->dlsa,1000000.0);
        propa->dla=prop->dl[0]+prop->dl[1];
        propa->tha=max(prop->the[0]+prop->the[1],-propa->dla*prop->gme);
        prop->wlos=false;
        prop->wscat=false;

        /*checking for parameters-in-range, error codes set if not */

        /* Warning: wave number */
        if (prop->wn<0.838 || prop->wn>210.0)
            prop->kwx=max(prop->kwx,1);

        /* Warning: radiation heights, AGL */
        for (j=0; j<2; j++)
            if (prop->hg[j]<1.0 || prop->hg[j]>1000.0)
                prop->kwx=max(prop->kwx,1);

        /* Warning: tx take-off angle */
        if(fabs(prop->the[0])>200e-3)
            prop->kwx=max(prop->kwx,3);

        /* Warning: rx take-off angle */
        if(fabs(prop->the[1])>1.220)
            prop->kwx=max(prop->kwx,3);

        /*for (j=0; j<2; j++)
          if (prop->dl[j]<0.1*propa->dls[j] || prop->dl[j]>3.0*propa->dls[j])
          prop->kwx=max(prop->kwx,3);  */

        /* Out-of-Range: refractivity, earth curvature, earth impedance/reactance, wave number */
        if (prop->ens<250.0 || prop->ens>400.0 || prop->gme<75e-9 || prop->gme>250e-9 || (tcreal(prop_zgnd)<=fabs(tcimag(prop_zgnd))) || prop->wn<0.419 || prop->wn>420.0)
            prop->kwx=4;

        /* Out-of-Range: radiation heights, AGL */
        for (j=0; j<2; j++)
            if (prop->hg[j]<0.5 || prop->hg[j]>3000.0)
                prop->kwx=4;


        adiff2_init(&state, prop, propa);

        prop->xae=pow(prop->wn*(prop->gme*prop->gme),-THIRD);  
        d3=max(propa->dlsa,1.3787*prop->xae+propa->dla);
        d4=d3+2.7574*prop->xae;
        a3=adiff2(&state, d3,prop,propa);
        a4=adiff2(&state, d4,prop,propa);
        propa->emd=(a4-a3)/(d4-d3);
        propa->aed=a3-propa->emd*d3;
    }

    if (prop->mdp>=0) /* if initializing the area mode */
    {
        prop->mdp=0;   /* area mode is initialized */
        prop->dist=dist;
    }

    if (prop->dist>0.0)
    {
        if (prop->dist>1000e3)   /* prop->dist being in meters, if greater than 1000 km, kwx=1 */
            prop->kwx=max(prop->kwx,1);

        double dmin=fabs(prop->he[0]-prop->he[1])/200e-3; /* minimum acceptable path distance length (meters) */
        if (prop->dist<dmin)
            prop->kwx=max(prop->kwx,3);

        if (prop->dist<1e3 || prop->dist>2000e3)
            prop->kwx=4;
    }

    if (prop->dist<propa->dlsa)
    {

        if (iw<=0.0)   /* if interval width is zero or less, used for area mode */
        {

            if (!prop->wlos)
            {
                d2=propa->dlsa;
                a2=propa->aed+d2*propa->emd;
                d0=1.908*prop->wn*prop->he[0]*prop->he[1];

                if (propa->aed>0.0)
                {
                    prop->aref=propa->aed+propa->emd*prop->dist;
                }
                else
                {
                    if (propa->aed==0.0)
                    {
                        d0=min(d0,0.5*propa->dla);
                        d1=d0+0.25*(propa->dla-d0);
                    }
                    else    /* aed less than zero */
                    {
                        d1=max(-propa->aed/propa->emd,0.25*propa->dla);
                    }
                    a1=alos2(prop,propa);
                    wq=false;

                    if (d0<d1)
                    {
                        a0=alos2(prop,propa);
                        a2=min(a2,alos2(prop,propa));
                        q=log(d2/d0);
                        propa->ak2=max(0.0,((d2-d0)*(a1-a0)-(d1-d0)*(a2-a0))/((d2-d0)*log(d1/d0)-(d1-d0)*q));
                        wq=propa->aed>=0.0 || propa->ak2>0.0;

                        if (wq)
                        { 
                            propa->ak1=(a2-a0-propa->ak2*q)/(d2-d0);

                            if (propa->ak1<0.0)
                            {
                                propa->ak1=0.0;
                                propa->ak2=FORTRAN_DIM(a2,a0)/q;

                                if (propa->ak2==0.0)
                                    propa->ak1=propa->emd;
                            }
                        }
                    }

                    if(!wq)
                    {       
                        propa->ak1=FORTRAN_DIM(a2,a1)/(d2-d1);
                        propa->ak2=0.0;

                        if (propa->ak1==0.0)
                            propa->ak1=propa->emd;

                    }
                    propa->ael=a2-propa->ak1*d2-propa->ak2*log(d2);
                    prop->wlos=true;
                }
            }
        }
        else    /* for ITWOM point-to-point mode */
        {

            if (!prop->wlos)
            {
                q=alos2(prop,propa); /* coefficient setup */
                prop->wlos=true;
            }

            if (prop->los==true)    /* if line of sight */
            {
                prop->aref=alos2(prop,propa);
            }
            else
            {
                if (trunc(prop->dist-prop->dl[0])==0)  /* if at 1st horiz */
                {
                    prop->aref=5.8+alos2(prop,propa);
                }
                else if (trunc(prop->dist-prop->dl[0])>0.0)    /* if past 1st horiz */
                {
                    adiff2_init(&state, prop, propa);
                    prop->aref = adiff2(&state, pd1,prop,propa);
                }
                else
                {
                    prop->aref=1.0;
                }

            }
        }
    }

    /* los and diff. range coefficents done. Starting troposcatter */
    if (prop->dist<=0.0 || prop->dist>=propa->dlsa)
    {
        if (iw==0.0)  /* area mode */
        {
            if(!prop->wscat)
            { 
                ascat_state scatterstate = {0};
                ascat_init(&scatterstate, prop);

                d5=propa->dla+200e3;
                d6=d5+200e3;
                a6=ascat(&scatterstate, d6,prop,propa);
                a5=ascat(&scatterstate, d5,prop,propa);

                if (a5<1000.0)
                {
                    propa->ems=(a6-a5)/200e3;
                    propa->dx=max(propa->dlsa,max(propa->dla+0.3*prop->xae*log(47.7*prop->wn),(a5-propa->aed-propa->ems*d5)/(propa->emd-propa->ems)));

                    propa->aes=(propa->emd-propa->ems)*propa->dx+propa->aed;
                }

                else
                {
                    propa->ems=propa->emd;
                    propa->aes=propa->aed;
                    propa->dx=10000000;
                }
                prop->wscat=true;
            }

            if (prop->dist>propa->dx)
            {
                prop->aref=propa->aes+propa->ems*prop->dist;
            }
            else
            {
                prop->aref=propa->aed+propa->emd*prop->dist;
            }
        }
        else   /* ITWOM mode  q used to preset coefficients with zero input */
        {
            if(!prop->wscat)
            { 
                d5=0.0;
                d6=0.0;
                ascat_state scatterstate = {0};
                ascat_init(&scatterstate, prop);

                a6=ascat(&scatterstate, pd1,prop,propa);

                adiff2_init(&state, prop, propa);
                a5=adiff2(&state, pd1,prop,propa);

                if (a5<=a6)
                {
                    propa->dx=10000000;
                    prop->aref=a5;
                }
                else
                {
                    propa->dx=propa->dlsa;
                    prop->aref=a6;
                }

                prop->wscat=true;
            }
        }
    }
    prop->aref=max(prop->aref,0.0);
}

/*
 * curve()
 *
 * Used in avar
 */
double curve (double const c1, double const c2, double const x1,
              double const x2, double const x3, double const de)
{
    /* return (c1+c2/(1.0+pow((de-x2)/x3,2.0)))*pow(de/x1,2.0)/(1.0+pow(de/x1,2.0)); */
    double temp1, temp2;

    temp1=(de-x2)/x3;
    temp2=de/x1;

    temp1*=temp1;
    temp2*=temp2;

    return (c1+c2/(1.0+temp1))*temp2/(1.0+temp2);
}

/*
 * avar()
 *
 * Analysis of Variants
 *
 * This subroutine calculates the additional loss resulting from statistical analysis
 * of long term variability due to time (reliability), location, and/or situational
 * variability (confidence), depending upon the operating mode determined by the mode of
 * variability code input. The subroutine reports out a single value of loss, avarv, in
 * decibels (dB), which includes aref, the "reference attenuation" summed with the
 * additional statistical loss.
 *
 * See ITWOM-SUB-ROUTINES.pdf, page 444
 *
 * All these const arrays are just lookup tables.
 */
const double bv1[7]={-9.67,-0.62,1.26,-9.21,-0.62,-0.39,3.15};
const double bv2[7]={12.7,9.19,15.5,9.05,9.19,2.86,857.9};
const double xv1[7]={144.9e3,228.9e3,262.6e3,84.1e3,228.9e3,141.7e3,2222.e3};
const double xv2[7]={190.3e3,205.2e3,185.2e3,101.1e3,205.2e3,315.9e3,164.8e3};
const double xv3[7]={133.8e3,143.6e3,99.8e3,98.6e3,143.6e3,167.4e3,116.3e3};
const double bsm1[7]={2.13,2.66,6.11,1.98,2.68,6.86,8.51};
const double bsm2[7]={159.5,7.67,6.65,13.11,7.16,10.38,169.8};
const double xsm1[7]={762.2e3,100.4e3,138.2e3,139.1e3,93.7e3,187.8e3,609.8e3};
const double xsm2[7]={123.6e3,172.5e3,242.2e3,132.7e3,186.8e3,169.6e3,119.9e3};
const double xsm3[7]={94.5e3,136.4e3,178.6e3,193.5e3,133.5e3,108.9e3,106.6e3};
const double bsp1[7]={2.11,6.87,10.08,3.68,4.75,8.58,8.43};
const double bsp2[7]={102.3,15.53,9.60,159.3,8.12,13.97,8.19};
const double xsp1[7]={636.9e3,138.7e3,165.3e3,464.4e3,93.2e3,216.0e3,136.2e3};
const double xsp2[7]={134.8e3,143.7e3,225.7e3,93.1e3,135.9e3,152.0e3,188.5e3};
const double xsp3[7]={95.6e3,98.6e3,129.7e3,94.2e3,113.4e3,122.7e3,122.9e3};
const double bsd1[7]={1.224,0.801,1.380,1.000,1.224,1.518,1.518};
const double bzd1[7]={1.282,2.161,1.282,20.,1.282,1.282,1.282};
const double bfm1[7]={1.0,1.0,1.0,1.0,0.92,1.0,1.0};
const double bfm2[7]={0.0,0.0,0.0,0.0,0.25,0.0,0.0};
const double bfm3[7]={0.0,0.0,0.0,0.0,1.77,0.0,0.0};
const double bfp1[7]={1.0,0.93,1.0,0.93,0.93,1.0,1.0};
const double bfp2[7]={0.0,0.31,0.0,0.19,0.31,0.0,0.0};
const double bfp3[7]={0.0,2.00,0.0,1.79,2.00,0.0,0.0};

double avar(double zzt, double zzl, double zzc, prop_type *prop, propv_type *propv)
{
    int kdv;
    bool ws, w1;
    double dexa, de, vmd, vs0, sgl, sgtm, sgtp, sgtd, tgtd,
           gm, gp, cv1, cv2, yv1, yv2, yv3, csm1, csm2, ysm1, ysm2,
           ysm3, csp1, csp2, ysp1, ysp2, ysp3, csd1, zd, cfm1, cfm2,
           cfm3, cfp1, cfp2, cfp3;

    double rt=7.8, rl=24.0, avarv, q, vs, zt, zl, zc;
    double sgt, yr, temp1, temp2;
    int temp_klim=propv->klim-1;

    if (propv->klim<=0 || propv->klim>7)
    {
        propv->klim=5;
        temp_klim=4;
        prop->kwx=max(prop->kwx,2);
    }

    cv1=bv1[temp_klim];
    cv2=bv2[temp_klim];
    yv1=xv1[temp_klim];
    yv2=xv2[temp_klim];
    yv3=xv3[temp_klim];
    csm1=bsm1[temp_klim];
    csm2=bsm2[temp_klim];
    ysm1=xsm1[temp_klim];
    ysm2=xsm2[temp_klim];
    ysm3=xsm3[temp_klim];
    csp1=bsp1[temp_klim];
    csp2=bsp2[temp_klim];
    ysp1=xsp1[temp_klim];
    ysp2=xsp2[temp_klim];
    ysp3=xsp3[temp_klim];
    csd1=bsd1[temp_klim];
    zd=bzd1[temp_klim];
    cfm1=bfm1[temp_klim];
    cfm2=bfm2[temp_klim];
    cfm3=bfm3[temp_klim];
    cfp1=bfp1[temp_klim];
    cfp2=bfp2[temp_klim];
    cfp3=bfp3[temp_klim];

    kdv=propv->mdvar;
    ws=kdv>=20;

    if (ws)
        kdv-=20;

    w1=kdv>=10;

    if (w1)
        kdv-=10;

    if (kdv<0 || kdv>3)
    {
        kdv=0;
        prop->kwx=max(prop->kwx,2);
    }

    q=log(0.133*prop->wn);

    /* gm=cfm1+cfm2/(pow(cfm3*q,2.0)+1.0); */
    /* gp=cfp1+cfp2/(pow(cfp3*q,2.0)+1.0); */

    gm=cfm1+cfm2/((cfm3*q*cfm3*q)+1.0);
    gp=cfp1+cfp2/((cfp3*q*cfp3*q)+1.0);

    dexa=sqrt(18e6*prop->he[0])+sqrt(18e6*prop->he[1])+pow((575.7e12/prop->wn),THIRD);

    if (prop->dist<dexa)
        de=130e3*prop->dist/dexa;
    else
        de=130e3+prop->dist-dexa;

    vmd=curve(cv1,cv2,yv1,yv2,yv3,de);
    sgtm=curve(csm1,csm2,ysm1,ysm2,ysm3,de)*gm;
    sgtp=curve(csp1,csp2,ysp1,ysp2,ysp3,de)*gp;
    sgtd=sgtp*csd1;
    tgtd=(sgtp-sgtd)*zd;

    if (w1)
        sgl=0.0;
    else
    {
        q=(1.0-0.8*exp(-prop->dist/50e3))*prop->dh*prop->wn;
        sgl=10.0*q/(q+13.0);
    }

    if (ws)
        vs0=0.0;
    else
    {
        /* vs0=pow(5.0+3.0*exp(-de/100e3),2.0); */
        temp1=(5.0+3.0*exp(-de/100e3));
        vs0=temp1*temp1;

    }

    zt=zzt;
    zl=zzl;
    zc=zzc;

    switch (kdv)
    {
        case 0:
            zt=zc;
            zl=zc;
            break;

        case 1:
            zl=zc;
            break;

        case 2:
            zl=zt;
    }

    if (fabs(zt)>3.1 || fabs(zl)>3.1 || fabs(zc)>3.1)
        prop->kwx=max(prop->kwx,1);

    if (zt<0.0)
        sgt=sgtm;

    else if (zt<=zd)
        sgt=sgtp;

    else
        sgt=sgtd+tgtd/zt;

    /* vs=vs0+pow(sgt*zt,2.0)/(rt+zc*zc)+pow(sgl*zl,2.0)/(rl+zc*zc); */

    temp1=sgt*zt;
    temp2=sgl*zl;

    vs=vs0+(temp1*temp1)/(rt+zc*zc)+(temp2*temp2)/(rl+zc*zc);

    if (kdv==0)
    {
        yr=0.0;
        propv->sgc=sqrt(sgt*sgt+sgl*sgl+vs);
    }

    else if (kdv==1)
    {
        yr=sgt*zt;
        propv->sgc=sqrt(sgl*sgl+vs);
    }

    else if (kdv==2)
    {
        yr=sqrt(sgt*sgt+sgl*sgl)*zt;
        propv->sgc=sqrt(vs);
    }

    else
    {
        yr=sgt*zt+sgl*zl;
        propv->sgc=sqrt(vs);
    }

    avarv=prop->aref-vmd-yr-propv->sgc*zc;

    if (avarv<0.0)
        avarv=avarv*(29.0-avarv)/(29.0-10.0*avarv);

    return avarv;
}

/*
 * Horizon calculations.
 *
 * Fills out a bunch of elements in prop:
 *    prop->the[0]: take-off angle at start
 *    prop->the[1]: take-off angle at end
 *
 *    prop->dl[0]: start horizon or highest obstacle distance
 *    prop->dl[1]: end horizon or highest obstacle distance
 *
 * This version is only used only with ITM 1.2.2. An updated version with improved obstacle
 * handling and that does caching for later calculations is below (hzns2).
 *
 * See ITWOM-SUB-ROUTINES.pdf, page 150
 */
void hzns(const elev_t pfl[], prop_type *prop)
{
    bool wq;
    int np;
    double xi, za, zb, qc, q, sb, sa;

    np=(int)pfl[0];                     /* number of elements in pfl array                      */
    xi=pfl[1];                          /* x increment, the distance between elements           */

    za=pfl[2]+prop->hg[0];              /* transmitter location height plus ant height AGL      */
    zb=pfl[np+2]+prop->hg[1];           /* receiver location height plus ant height AGL         */

    qc=0.5*prop->gme;                   /* half the earth's curvature, or 1/(2r)                */
    q=qc*prop->dist;                    /* q is now dist/2r (i.e. radians)                      */

    prop->the[1]=(zb-za)/prop->dist;    /* temporarily set slope from tx to rx.                 */
    /* XXX BUG. This really should be an arctan. See        */
    /* ITWOM pp156-157.                                     */

    prop->the[0]=prop->the[1]-q;        /* set takeoff angle of receiver                        */
    prop->the[1]=-prop->the[1]-q;       /* and set takeoff angle of transmitter to the opposite */

    prop->dl[0]=prop->dist;             /* optimistically set tx horizon to longest path        */
    prop->dl[1]=prop->dist;             /* ditto for rx horizon                                 */

    if (np<2)                           /* short path? all done.                                */
        return;

    sa=0.0;                             /* sa is at tx loc                                      */
    sb=prop->dist;                      /* sb is at rx loc                                      */
    wq=true;                            /* line-of-sight study ok?                              */

    for (int i=1; i<np; i++)
    {
        sa+=xi;                                         /* move sa towards rx loc               */
        sb-=xi;                                         /* move sb towards tx loc               */

        q=pfl[i+2] - ((qc*sa+prop->the[0])*sa) - za;    /* height difference from pfl pt to sa, */
        /* adjusted for curvature               */

        if (q>0.0)                                      /* tx obstruction found?                */
        {
            prop->the[0]+=q/sa;
            prop->dl[0]=sa;
            wq=false;
        }

        if (!wq)                                       /* no tx obstructions, so examine rx end */
        {
            q=pfl[i+2] - ((qc*sb+prop->the[1])*sb) - zb; /* height diff from pfl pt to sb,      */
            /* adjusted for curvature              */

            if (q>0.0)                                 /* rx obstruction found?                 */
            {
                prop->the[1]+=q/sb;
                prop->dl[1]=sb;
            }
        }
    }
}

/*
 * Horizon calculations.
 *
 * Fills out a bunch of elements in prop:
 *    prop->the[0]: take-off angle at start
 *    prop->the[1]: take-off angle at end
 *
 *    prop->dl[0]: distance to start horizon or highest obstacle
 *    prop->dl[1]: distance to end horizon or highest obstacle
 *
 *
 * This version has improved obstacle handling and does some caching for later calculations.
 *
 * See ITWOM-SUB-ROUTINES.pdf, page 150
 */
void hzns2(const elev_t pfl[], prop_type *prop)
{
    bool wq;
    int np, rp, i, j;
    double xi, za, zb, qc, q, sb, sa, dr, dshh;

    np=(int)pfl[0];                     /* number of elements in pfl array                      */
    xi=pfl[1];                          /* x increment, the distance between elements           */

    za=pfl[2]+prop->hg[0];              /* transmitter location height plus ant height AGL      */
    zb=pfl[np+2]+prop->hg[1];           /* receiver location height plus ant height AGL         */

    prop->tiw=xi;                       /* cache some values for later                          */
    prop->ght=za;
    prop->ghr=zb;

    qc=0.5*prop->gme;                   /* gme is 1/r, where r is the earth's radius            */
    q=qc*prop->dist;                    /* q is now dist/2r (i.e. radians)                      */

    prop->the[1]=atan((zb-za)/prop->dist);  /* temporarily set the slope from tx to rx          */

    prop->the[0]=prop->the[1]-q;        /* set takeoff angle of receiver                        */
    prop->the[1]=-prop->the[1]-q;       /* and set takeoff angle of transmitter to the opposite */

    prop->dl[0]=prop->dist;             /* optimistically set tx horizon to longest path        */
    prop->dl[1]=prop->dist;             /* ditto for rx horizon                                 */

    prop->hht=0.0;                      /* set tx obstacle height to 0                          */
    prop->hhr=0.0;                      /* set rx obstacle height to 0                          */
    prop->los=true;                     /* line-of-sight study ok?                              */

    if(np>=2)
    {
        sa=0.0;                         /* sa is at tx loc                                     */
        sb=prop->dist;                  /* sb is at rx loc                                     */
        wq=true;                        /* assume entire path is Line-Of-Sight                 */

        for(j=1; j<np; j++)
        {
            sa+=xi;                                       /* add to the distance from tx       */
            q=pfl[j+2] - ((qc*sa+prop->the[0])*sa) - za;  /* height diff from pfl pt to sa,    */
            /* middle argument adjusts for earth's curvature   */

            if(q>0.0)                                  /* tx obstruction found?                */
            {
                prop->los=false;
                prop->the[0]+=q/sa;
                prop->dl[0]=sa;
                prop->the[0]=min(prop->the[0],1.569); /* limit the max radians of takeoff angle */
                prop->hht=pfl[j+2];         /* Note: this can be 0, causing dr (below) to blow up! */
                wq=false;
            }
        }

        if(!wq)                               /* we know it's not LOS so examine the rx end   */
        {
            for(i=1; i<np; i++)
            {
                sb-=xi;                          /* move sb towards tx                        */

                /* and calc height diff                       */
                q=pfl[np+2-i] - ((qc*(prop->dist-sb)+prop->the[1]) * (prop->dist-sb)) - zb;

                if(q>0.0)
                {
                    prop->the[1]+=q/(prop->dist-sb);

                    prop->the[1]=min(prop->the[1],1.57); /* limit the rx take-off angle range */
                    prop->the[1]=max(prop->the[1],-1.568);

                    prop->hhr=pfl[np+2-i]; /* NOTE: this can be 0, causing dr (below) to blow up! */
                    prop->dl[1]=max(0.0,prop->dist-sb);
                }
            }

            /* recalculate using atan for more accuracy */
            prop->the[0]=atan((prop->hht-za)/prop->dl[0])-0.5*prop->gme*prop->dl[0];
            prop->the[1]=atan((prop->hhr-zb)/prop->dl[1])-0.5*prop->gme*prop->dl[1];
        }
    }

    /* if the receiver horizon path is shorter than the total path, then there's
     * some obstacle, so let's calculate the 2-ray reflection point.
     *
     * dr is the distance in meters from the last obstruction peak to the reflection
     *    point.
     */
    dr = 0.0;
    if((prop->dl[1])<(prop->dist))
    {
        dshh=prop->dist-prop->dl[0]-prop->dl[1]; /* distance between obstacle peaks */

        /* see ITWOM p63 and p162 */
        if(trunc(dshh)==0) /* one obstacle */
        {
            if (prop->hht > 0.0) {
                dr=prop->dl[1]/(1+zb/prop->hht);
            }
        }
        else    /* two obstacles */
        {
            if (prop->hhr > 0.0) {   /* if the antenna has 0 height, there is no 2-ray reflection */
                dr=prop->dl[1]/(1+zb/prop->hhr);
            }
        }
    }
    else    /* line of sight  */
    {
        if (za > 0.0) {
            dr=(prop->dist)/(1+zb/za);   /* if the antenna has 0 height, there is no 2-ray reflection */
        }
    }

    rp=2+(int)(floor(0.5+dr/xi));

    prop->rpl=rp;
    prop->rph=pfl[rp];
}

/* Linear least-squares fit
 *
 * This is trying to create a line fitting the equation y=a+bx using
 * a least squares fit (where b is derived by doing b=sum(x*y)/sum(x*x)
 *
 * THIS HAS A MAJOR FLAW IN IT AND SHOULD NOT BE USED. See the comments
 * towards the end of the function.
 *
 * z[]: input array of data points, where
 *      z[0]: number of points in the array
 *      z[1]: length between the points
 *  x1: length from the start of the array to the start of the section of the
 *      path we want to fit.
 *  x2: length from the end of the array to the end of the section of the
 *      path we want to fit.
 *
 *  i.e. if the path is 20 units long and we want fit the line 2 units in
 *    to 15 units in, x1 will be 2 and x2 will be 5.
 *
 *  Returns:
 *  z0: calculated z-axis value (height) at x1
 *  zn: calculated z-axis value (height) at x2
 *
 * See ITWOM-SUB-ROUTINES.pdf p289
 *
 * Used only with ITM 1.2.2
 */
void z1sq1 (const elev_t z[], const double x1, const double x2, double *z0, double *zn)
{
    double xn, xa, xb, x, a, b;
    int n, ja, jb;

    xn=z[0];                               /* does this really need to be a double? */
    xa=trunc(FORTRAN_DIM(x1/z[1],0.0));    /* set our start interval count */
    xb=xn-trunc(FORTRAN_DIM(xn,x2/z[1]));  /* set our end interval count  */

    if (xb<=xa)     /* if we have a very short path, xb can be less than xa, so... */
    {
        xa=FORTRAN_DIM(xa,1.0);         /* move xa back by one (or set it to 0) */
        xb=xn-FORTRAN_DIM(xn,xb+1.0);   /* move xb up by one (or set to the end) */
    }

    ja=(int)xa;
    jb=(int)xb;

    n=jb-ja;                    /* total number of points as an integer.          */
    xa=xb-xa;                    /* also the total number of points, but as a float. */
    x=-0.5*xa;                  /* number of intervals below and above the midpoint */
    /* crossing. We made it negative because we'll be */
    /* from -x to +x, past the middle at 0.              */

    xb+=x;                 /* xb is  now be the number of intervals to the middle */

    a=0.5*(z[ja+2]+z[jb+2]);   /* avg of the start and end heights, i.e. initial intercept, a */
    b=0.5*(z[ja+2]-z[jb+2])*x; /* avg of the start/end height difference times the run, i.e initial slope */

    for (int i=2; i<=n; ++i)
    {
        ++ja;
        x+=1.0;                /* we're moving from -x to +x, one by one      */
        a+=z[ja+2];               /* a=sum(y)                                 */
        b+=z[ja+2]*x;           /* b=sum(x*y)                                 */
    }

    a/=xa;                        /* a = sum(y)/numpoints                     */

    /* XXX The following is a MAJOR bug. The derived slope is nonsensical!    */
    b=b*12.0/((xa*xa+2.0)*xa);  /* b =  (huh, wha?)                           */

    /* ok, now we have y=a+bx, so let's use it:         */
    *z0=a-b*xb;             /* set height of z0                               */
    *zn=a+b*(xn-xb);        /* set height of zn                               */
}

/* Linear least-squares fit
 *
 * This is trying to create a line fitting the equation y=a+bx using
 * a least squares fit (where b is derived by doing b=sum(x*y)/sum(x*x).
 * This is using a simple linear regression without the intercept term - i.e.
 * we're forcing the line to go through the origin at the midpoint of the
 * path.
 *
 * z[]: input array of data points, where
 *      z[0]: number of points in the array
 *      z[1]: length between the points
 *  x1: length from the start of the array to the start of the section of the
 *      path we want to fit.
 *  x2: length from the end of the array to the end of the section of the
 *      path we want to fit.
 *
 *  i.e. if the path is 20 units long and we want fit the line 2 units in
 *    to 15 units in, x1 will be 2 and x2 will be 5.
 *
 *  Returns:
 *  z0: calculated z-axis value (height) at x1
 *  zn: calculated z-axis value (height) at x2
 *
 *  See ITWOM-SUB-ROUTINES.pdf p298
 */
void z1sq2(const elev_t z[], const double x1, const double x2, double *z0, double *zn)
{
    /* corrected for use with ITWOM */
    double xn, xa, xb, x, a, b, bn;
    int n, ja, jb;

    xn=z[0];                               /* does this really need to be a double? */
    xa=trunc(FORTRAN_DIM(x1/z[1],0.0));    /* set our start interval count */
    xb=xn-trunc(FORTRAN_DIM(xn,x2/z[1]));  /* set our end interval count  */

    if (xb<=xa)     /* if we have a very short path, xb can be less than xa, so... */
    {
        xa=FORTRAN_DIM(xa,1.0);         /* move xa back by one (or set it to 0) */
        xb=xn-FORTRAN_DIM(xn,xb+1.0);   /* move xb up by one (or set to the end) */
    }

    ja=(int)xa;
    jb=(int)xb;

    xa=(2*trunc((xb-xa)/2))-1;    /* make xa an odd integer representing the number */
    /* of intervals to consider. It has to be odd so  */
    /* we have a midpoint zero crossing.              */

    x=-0.5*(xa+1);                /* number of intervals below and above the midpoint*/
    /* crossing. We made it negative because we'll be */
    /* from -x to +x, past the middle at 0.              */

    xb+=x;                        /* xb is the number of intervals to the middle    */
    ja=jb-1-(int)xa;            /* fix up the start point because we insisted on  */
    /* an odd number of intervals. We might need to   */
    /* move it a bit to accomodate.                   */

    n=jb-ja;                   /* total number of intervals, which should be 2*xb */

    a=(z[ja+2]+z[jb+2]);        /* sum of the start and end heights */
    b=(z[ja+2]-z[jb+2])*x;        /* initial slope */

    bn=2*(x*x);                /* initial value for our slope denominator, sum(x*x)   */

    for (int i=2; i<=n; ++i)
    {
        ++ja;
        x+=1.0;            /* we're moving from -x to +x, one by one            */
        bn+=(x*x);         /* bn = sum(x*x)                                     */
        a+=z[ja+2];           /* a = sum(y)                                        */
        b+=z[ja+2]*x;      /* b = sum(x*y)                                      */
    }

    a/=(xa+2);                /* a = sum(y)/numpoints                                */
    b=b/bn;                    /* b = sum(x*y)/sum(x*x)                            */

    /* ok, now we have y=a+bx, so let's use it:         */
    *z0=a-(b*xb);            /* set height of z0                                    */
    *zn=a+(b*(xn-xb));      /* set height of zn                                 */
}


/* Used to find a quantile.  It reorders the array a so that all the elements
 * before ir are greater than or equal to all the elements after ir. In particular,
 * a(ir) will have the same value it would have if a were completely sorted in
 * descending order.
 *
 * For example, if you have an array of 700 numbers (nn=700), and you ask for
 * the 630th quantile (ir=630) (which happens to be the 90% quantile), the
 * function will sort the array just enough to figure out that it should return
 * the value at which 90% the of the things are less. In a truly random array,
 * this should return a value of about 630.
 *
 *  nn: number of elements in the array a[]
 * a[]: data array, of size nn
 *  ir: the quantile desired.
 *
 * returns the value of a(ir)
 *
 * This is essentially a quickselect with Hoare's partition scheme. See
 * https://en.wikipedia.org/wiki/Quickselect for info on the algorithm.
 * 
 * Also see ITWOM-SUB-ROUTINES.pdf p265
 *
 */
double qtile (const int nn, elev_t a[], const int ir)
{
    double q=0.0, r; /* q initialization -- KD2BD */
    int m, n, i, j, j1=0, i0=0, k;  /* more initializations -- KD2BD */
    bool done=false;
    bool goto10=true;

    m=0;
    n=nn;
    k=min(max(0,ir),n-1);

    while (!done)
    {
        if (goto10)
        {
            q=a[k];
            i0=m;
            j1=n;
        }

        i=i0;

        while (i<=n && a[i]>=q)
            i++;

        if (i>n)
            i=n;

        j=j1;

        while (j>=m && a[j]<=q)
            j--;

        if (j<m)
            j=m;

        if (i<j)
        {
            r=a[i];
            a[i]=a[j];
            a[j]=r;
            i0=i+1;
            j1=j-1;
            goto10=false;
        }

        else if (i<k)
        {
            a[k]=a[i];
            a[i]=q;
            m=i+1;
            goto10=true;
        }

        else if (j>k)
        {
            a[k]=a[j];
            a[j]=q;
            n=j-1;
            goto10=true;
        }

        else
            done=true;
    }

    return q;
}

/*
 * qerf()
 *
 * Standard normal complementary probability.
 *
 * qerfi() is the inverse of this.
 *
 * Used only with ITM 1.2.2
 */
double qerf(const double z)
{
    double b1=0.319381530, b2=-0.356563782, b3=1.781477937;
    double b4=-1.821255987, b5=1.330274429;
    double rp=4.317008, rrt2pi=0.398942280;
    double t, x, qerfv;

    x=z;
    t=fabs(x);

    if (t>=10.0)
        qerfv=0.0;
    else
    {
        t=rp/(t+rp);
        qerfv=exp(-0.5*x*x)*rrt2pi*((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
    }

    if (x<0.0)
        qerfv=1.0-qerfv;

    return qerfv;
}

/*
 * Delta h, experimental
 *
 * Used to find delta h, the terrain irregularity factor. delta h is defined
 * as the interdecile range of elevations between two points. The interdecile
 * range is defined as the difference between the 10th and 90th percentiles.
 *
 * pfl[]: profile elevation array, with:
 *          pfl[0]: number of points in pfl
 *          pfl[1]: distance between points in pfl (in meters)
 *    x1: start point along the pfl path, in meters
 *    x2: end point along the pfl path, in meters
 *
 * Returns the delta h value.
 *
 * See notes for d1thx2(), the newer version, below. This version is only used
 * with ITM 1.2.2 and is limited to using a maximum of 245 pfl array elements.
 * It is thus faster but slightly less accurate.
 */
double d1thx(const elev_t pfl[], const double x1, const double x2)
{
    int np, ka, kb, n, k, j;
    double d1thxv, sn, xa, xb;
    elev_t *s;

    np=(int)pfl[0];
    xa=x1/pfl[1];                  /* start point's array element */
    xb=x2/pfl[1];                  /* end point's array element */
    d1thxv=0.0;

    if (xb-xa<2.0) 
        return d1thxv;

    ka=(int)(0.1*(xb-xa+8.0));     /* from 32 to 242 points */
    ka=min(max(4,ka),25);          /* ka can range from 4-25. */
    n=10*ka-5;                     /* n can range from 35-245 */
    kb=n-ka+1;                     /* kb can range from 32-221 */
    sn=n-1;                        /* index of last path element to consider */

    s = (elev_t*)calloc(n+2, sizeof(elev_t));
    assert(s!=NULL);
    s[0]=sn;
    s[1]=1.0;

    xb=(xb-xa)/sn;                 /* rescale xb to the new distance */

    k=(int)(xa+1.0);
    xa-=(double)k;                 /* xa now ranges from -1.0 to near 0 */

    /* Load up a temporary height array, smoothed */
    for (j=0; j<n; j++)
    {
        while (xa>0.0 && k<np)
        {
            xa-=1.0;
            ++k;
        }

        s[j+2]=pfl[k+2]+(pfl[k+2]-pfl[k+1])*xa;
        xa=xa+xb;
    }

    /* Do a least-squares fit of the terrain from point 0.0 to sn meters out
       returning the heights of the endpoints in xa and xb */
    z1sq1(s,0.0,sn,&xa,&xb);

    xb=(xb-xa)/sn;      /* xa and xb are now heights, so xb is now a slope */

    for (j=0; j<n; j++)
    {
        s[j+2]-=xa; /* difference between the original height and new slope */
        xa=xa+xb;
    }

    /* Now find the difference between the 90% and 10% heights.
     * These two qtiles calls take a (relatively) enormous amount of
     * processing time. See notes in d1thx2() below. Here they are limited to
     * a max of 245 elements so it's not so bad. */
    d1thxv=qtile(n-1,s+2,ka-1)-qtile(n-1,s+2,kb-1);

    /* apply empirical data matching magic scaling. See ITWOM p124. */
    d1thxv/=1.0-0.8*exp(-(x2-x1)/50.0e3);

    free(s);

    return d1thxv;
}

/*
 * Delta h, experimental
 *
 * Used to find delta h, the terrain irregularity factor. delta h is defined
 * as the interdecile range of elevations between two points. The interdecile
 * range is defined as the difference between the 10th and 90th percentiles.
 * 
 * pfl[]: profile elevation array, with:
 *          pfl[0]: number of points in pfl
 *          pfl[1]: distance between points in pfl (in meters)
 *    x1: start point along the pfl path, in meters
 *    x2: end point along the pfl path, in meters
 *
 * Returns the delta h value.
 *
 * See ITWOM-SUB-ROUTINES.pdf p125. This version has been modified by Sid
 * Shumate to use the entire range of the pfl array if needed.
 */
double d1thx2(const elev_t pfl[], const double x1, const double x2)
{
    int np, ka, kb, n, k, kmx, j;
    double d1thx2v, sn, xa, xb, xc;
    elev_t *s;

    np=(int)pfl[0];
    xa=x1/pfl[1];                  /* start point's array element */
    xb=x2/pfl[1];                  /* end point's array element */
    d1thx2v=0.0;

    if (xb-xa<2.0)
        return d1thx2v;

    /* Compare this section to d1thx() above. */
    ka=(int)(0.1*(xb-xa+8.0));
    kmx=max(25,(int)(83350/(pfl[1]))); /* Sid's fix: dynamically figure out the max elements */
    ka=min(max(4,ka),kmx);         /* range from 4 to kmax */
    n=10*ka-5;
    kb=n-ka+1;
    sn=n-1;

    s = (elev_t*)calloc(n+2, sizeof(elev_t));
    assert(s!=NULL);
    s[0]=sn;
    s[1]=1.0;

    xb=(xb-xa)/sn;                  /* rescale xb to the new distance */

    k=(int)(trunc(xa+1.0));
    xc=xa-((double)k);              /* xa now ranges from -1.0 to near 0                   */
    /* This new variable, xc, was added by Sid because     */
    /* xa's value is later passed to z1sq2 and he wanted   */
    /* to preserve the original xa value. However, z1sq2() */
    /* never uses the value and just sets it so this is    */
    /* pointless. There's no harm in it though.            */

    /* Load up a temporary height array, smoothed */
    for (j=0; j<n; j++)
    {
        while (xc>0.0 && k<np)
        {
            xc-=1.0;
            ++k;
        }

        s[j+2]=pfl[k+2]+(pfl[k+2]-pfl[k+1])*xc; 
        xc=xc+xb;
    }

    /* Do a least-squares fit of the terrain from point 0.0 to sn meters out
       returning the heights of the endpoints in xa and xb */
    z1sq2(s,0.0,sn,&xa,&xb);

    xb=(xb-xa)/sn;      /* xa and xb are now heights, so xb is now a slope */

    for (j=0; j<n; j++)
    {
        s[j+2]-=xa; /* difference between the original height and new slope */
        xa=xa+xb;
    }

    /* Now find the difference between the 90% and 10% heights.
     * These two qtiles calls take a (relatively) enormous amount of
     * processing time. */
    d1thx2v=qtile(n-1,s+2,ka-1)-qtile(n-1,s+2,kb-1);

    /* apply empirical data matching magic scaling. See ITWOM p124. */
    d1thx2v/=1.0-0.8*exp(-(x2-x1)/50.0e3);

    free(s);

    return d1thx2v;
}

/*
 * qlrpfl()
 *
 * Quick Longley-Rice Profile
 *
 * pfl[]: profile elevation array, with:
 *          pfl[0]: number of points in pfl
 *          pfl[1]: distance between points in pfl (in meters)
 * klimx: climate code
 * mdvarx: mode of variability. Usually 12, can be set to 1 for FCC mode.
 * prop:  prop_type state variable
 * propa: propa_type state variable
 * propv: propv_type state variable
 *
 *
 *
 * See ITWOM-SUB-ROUTINES.pdf p233
 */
void qlrpfl(const elev_t pfl[], int klimx, int mdvarx, prop_type *prop, propa_type *propa, propv_type *propv)
{
    int np, j;
    double xl[2], q, za, zb, temp;

    prop->dist=pfl[0]*pfl[1]; /* total distance of the pfl array */
    np=(int)pfl[0];           /* number of points in the pfl array */

    hzns(pfl,prop); /* analyse pfl and store horizon/obstruction info in prop */

    for (j=0; j<2; j++)                               /* for both tx and rx... */
        xl[j]=min(15.0*prop->hg[j], 0.1*prop->dl[j]); /* ...set xl to min of 15x ant height or 1/10 horizon dist */

    xl[1]=prop->dist-xl[1]; /* adjust the rx distance to be from the far end */

    prop->dh=d1thx(pfl,xl[0],xl[1]);  /* calculate the terrain irregularity factor */

    if (prop->dl[0]+prop->dl[1] > 1.5*prop->dist)
    {
        /* the horizon (or obstruction) is far away... */

        z1sq1(pfl,xl[0],xl[1],&za,&zb);                 /* do a linear least-squares fit */
        /* za has height of line at xl[0] */
        /* zb has height of line at xl[1] */

        /* set effective heights to endpoint heights plus an offset. See ITWOM p236 for discussion */
        prop->he[0]=prop->hg[0]+FORTRAN_DIM(pfl[2],za);
        prop->he[1]=prop->hg[1]+FORTRAN_DIM(pfl[np+2],zb);

        /* arcana to (re)determine dl values that we initially got from hzns. See ITWOM. */
        for (j=0; j<2; j++)
            prop->dl[j]=sqrt(2.0*prop->he[j]/prop->gme)*exp(-0.07*sqrt(prop->dh/max(prop->he[j],5.0)));

        q=prop->dl[0]+prop->dl[1];

        if (q<=prop->dist) /* if there is a rounded horizon, or two obstructions, in the path */
        {
            /* q=pow(prop->dist/q,2.0); */
            temp=prop->dist/q;
            q=temp*temp;

            for (j=0; j<2; j++)
            {
                prop->he[j]*=q; /* tx effective height set to be path dist/distance between obstacles */
                prop->dl[j]=sqrt(2.0*prop->he[j]/prop->gme)*exp(-0.07*sqrt(prop->dh/max(prop->he[j],5.0)));
            }
        }

        for (j=0; j<2; j++) /* original empirical adjustment?  uses delta-h to adjust grazing angles */
        {
            q=sqrt(2.0*prop->he[j]/prop->gme);
            prop->the[j]=(0.65*prop->dh*(q/prop->dl[j]-1.0)-2.0*prop->he[j])/q;
        }
    }
    else
    {
        /* the horizon (or obstruction) is nearish... */

        z1sq1(pfl,xl[0],0.9*prop->dl[0],&za,&q);
        z1sq1(pfl,prop->dist-0.9*prop->dl[1],xl[1],&q,&zb);

        prop->he[0]=prop->hg[0]+FORTRAN_DIM(pfl[2],za);
        prop->he[1]=prop->hg[1]+FORTRAN_DIM(pfl[np+2],zb);
    }

    prop->mdp=-1;

    if (mdvarx>=0)
    {
        propv->mdvar=mdvarx;
    }

    if (klimx>0)
    {
        propv->klim=klimx;
    }

    lrprop(0.0,prop,propa);
}


/*
 * qlrpfl2()
 *
 * Quick Longley-Rice Profile
 *
 * pfl[]: profile elevation array, with:
 *          pfl[0]: number of points in pfl
 *          pfl[1]: distance between points in pfl (in meters)
 * klimx: climate code
 * mdvarx: mode of variability. Usually 12, can be set to 1 for FCC mode.
 * prop:  prop_type state variable
 * propa: propa_type state variable
 * propv: propv_type state variable
 *
 * See ITWOM-SUB-ROUTINES.pdf p247
 */
void qlrpfl2(const elev_t pfl[], int klimx, int mdvarx, prop_type *prop, propa_type *propa, propv_type *propv)
{
    int np, j;
    double xl[2], dlb, za, zb, temp, rad, rae1, rae2;
    double q = 1;

    prop->dist=pfl[0]*pfl[1];  /* total distance of the pfl array */
    np=(int)pfl[0];            /* number of points in the pfl array */

    hzns2(pfl,prop); /* analyse pfl and store horizon/obstruction info in prop */

    dlb=prop->dl[0]+prop->dl[1];
    prop->rch[0]=prop->hg[0]+pfl[2];
    prop->rch[1]=prop->hg[1]+pfl[np+2];

    for (j=0; j<2; j++)                              /* for both tx and rx... */
        xl[j]=min(15.0*prop->hg[j],0.1*prop->dl[j]); /* ...set xl to min of 15x ant height or 1/10 horizon dist */

    xl[1]=prop->dist-xl[1];           /* adjust the rx distance to be from the far end */
    prop->dh=d1thx2(pfl,xl[0],xl[1]); /* calculate the terrain irregularity factor */

    /* The first branch of this if statement is for when there are no points in array (e.g. called in the now-deprecated
     * "area" mode, or when the distance between points is very large.
     *
     * In other words, it's not really ever used in modern code.
     */
    if ((np<1) || (pfl[1]>150.0))
    {
        /* TRANSHORIZON; diffraction over a mutual horizon, or for one or more obstructions */
        if (dlb<1.5*prop->dist)  
        {
            z1sq2(pfl,xl[0],0.9*prop->dl[0],&za,&q);       /* fit line to terrain from start to 90% of start horizon/obstacle */
            z1sq2(pfl,prop->dist-0.9*prop->dl[1],xl[1],&q,&zb); /* ditto, but from the other end */

            prop->he[0]=prop->hg[0]+FORTRAN_DIM(pfl[2],za);    /* set effective height at start (ASL) */
            prop->he[1]=prop->hg[1]+FORTRAN_DIM(pfl[np+2],zb); /* ...and at end                 (ASL) */
        }

        /* Line-of-Sight path */
        else
        {
            z1sq2(pfl,xl[0],xl[1],&za,&zb);                           /* fit a line */
            prop->he[0]=prop->hg[0]+FORTRAN_DIM(pfl[2],za);        /* set effective height at start (ASL) */
            prop->he[1]=prop->hg[1]+FORTRAN_DIM(pfl[np+2],zb);     /* set effective height at end (ASL) */

            /* this is essentially calculating the horizon distances, using the earth's curvature as
             * an approximation. The calculation can be done much more accurately using the hzns call
             * with terrain data, rendering this whole section obsolete. */
            for (j=0; j<2; j++)
                prop->dl[j]=sqrt(2.0*prop->he[j]/prop->gme)*exp(-0.07*sqrt(prop->dh/max(prop->he[j],5.0)));

            /* for one or more obstructions only
             * This was inoperative in the ITM FORTRAN and DLL. It was reactivated and tested and found
             * to be unstable. It was therefore reburied (made nonfunctional).
             * See the discussion in ITWOM p250.
             */
            if ((prop->dl[0]+prop->dl[1])<=prop->dist)
            {
                /* q=pow(prop->dist/(dl[0]+dl[1])),2.0); */
                temp=prop->dist/(prop->dl[0]+prop->dl[1]);
                q=temp*temp;
            }

            /* re-recalculate the horizon distances, using the earth's curvature as an approximation. */
            /* Obsolete. */
            for (j=0; j<2; j++)
            {
                prop->he[j]*=q;
                prop->dl[j]=sqrt(2.0*prop->he[j]/prop->gme)*exp(-0.07*sqrt(prop->dh/max(prop->he[j],5.0)));
            }

            /* this sets (or resets) prop->the, and is not in The Guide FORTRAN QLRPFL */
            for (j=0; j<2; j++)
            {
                q=sqrt(2.0*prop->he[j]/prop->gme);
                prop->the[j]=(0.65*prop->dh*(q/prop->dl[j]-1.0)-2.0*prop->he[j])/q;
            }
        }
    }
    else    /* for ITWOM ,computes he for tx, rcvr, and the receiver approach angles for use in saalos */ 
    {

        prop->he[0]=prop->hg[0]+(pfl[2]);     /* set the effective height at start (ASL)  */
        prop->he[1]=prop->hg[1]+(pfl[np+2]);  /* set the effective height at end (ASL)    */

        rad=(prop->dist-500.0);    /* receiver approach path start location. We consider the last half-km only */

        if (prop->dist>550.0)
        {
            z1sq2(pfl,rad,prop->dist,&rae1,&rae2); /* do a least-squares fit on that last half-km.*/
        }
        else
        {
            rae1=0.0;
            rae2=0.0;
        }

        /* theta, receive angle, is the slope of the least-squares line. This could be returned by
         * z1sq2, but recovering it isn't that hard.
         * It's only used in alos(). */
        prop->thera=atan(fabs(rae2-rae1)/prop->dist);

        if (rae2<rae1)
        {
            prop->thera=-prop->thera;
        }

        /* theta, near receiver, is the slope of the terrain right next to the receiver. It's
         * only used in alos(). */
        prop->thenr=atan(max(0.0,(pfl[np+2]-pfl[np+1]))/pfl[1]);

    }

    prop->mdp=-1;

    if (mdvarx>=0)
    {
        propv->mdvar=mdvarx;
    }

    if (klimx>0)
    {
        propv->klim=klimx;
    }

    lrprop2(0.0,prop,propa);
}


/***************************************************************************************
 * Point-To-Point Mode Calculations 
 ***************************************************************************************/


/******************************************************************************
  point_to_point_ITM()

  Note that point_to_point has become point_to_point_ITM for use as the old ITM 

pol:
0-Horizontal, 1-Vertical

radio_climate:
1-Equatorial, 2-Continental Subtropical,
3-Maritime Tropical, 4-Desert, 5-Continental Temperate,
6-Maritime Temperate, Over Land, 7-Maritime Temperate,
Over Sea

conf, rel: .01 to .99

elev[]: [num points - 1], [delta dist(meters)],
[height(meters) point 1], ..., [height(meters) point n]

errnum: 0- No Error.
1- Warning: Some parameters are nearly out of range.
Results should be used with caution.
2- Note: Default parameters have been substituted for
impossible ones.
3- Warning: A combination of parameters is out of range.
Results are probably invalid.
Other-  Warning: Some parameters are out of range.
Results are probably invalid.

 *****************************************************************************/
void point_to_point_ITM(const elev_t elev[], double tht_m, double rht_m, double eps_dielect, double sgm_conductivity, double eno_ns_surfref, double frq_mhz, int radio_climate, int pol, double conf, double rel, double *dbloss, char *strmode, int *errnum)

{
    prop_type   prop = {0};
    propv_type  propv = {0};
    propa_type  propa = {0};
    double zsys=0;
    double zc, zr;
    double eno, enso, q;
    long ja, jb, i, np;
    /* double dkm, xkm; */
    double fs;

    prop.hg[0]=tht_m;
    prop.hg[1]=rht_m;
    propv.klim=radio_climate;
    prop.kwx=0;
    prop.mdp=-1;
    zc=qerfi(conf);
    zr=qerfi(rel);
    np=(long)elev[0];
    eno=eno_ns_surfref;
    enso=0.0;
    q=enso;

    if (q<=0.0)
    {
        ja=(long)(3.0+0.1*elev[0]);  /* added (long) to correct */
        jb=np-ja+6;

        for (i=ja-1; i<jb; ++i)
            zsys+=elev[i];

        zsys/=(jb-ja+1);
        q=eno;
    }

    propv.mdvar=12;

    qlrps(frq_mhz,zsys,q,pol,eps_dielect,sgm_conductivity,&prop); /* quick longley rice - setup */

    qlrpfl(elev,propv.klim,propv.mdvar,&prop,&propa,&propv);  /* quick longley-rice, do the calculation */

    fs=32.45+20.0*log10(frq_mhz)+20.0*log10(prop.dist/1000.0);
    q=prop.dist-propa.dla;

    if (trunc(q)<0.0)
        strcpy(strmode,"Line-Of-Sight Mode");
    else
    {
        if (trunc(q)==0.0)
            strcpy(strmode,"Single Horizon");

        else if (trunc(q)>0.0)
            strcpy(strmode,"Double Horizon");

        if (prop.dist<=propa.dlsa || prop.dist <= propa.dx)
            strcat(strmode,", Diffraction Dominant");

        else if (prop.dist>propa.dx)
            strcat(strmode, ", Troposcatter Dominant");
    }

    *dbloss=avar(zr,0.0,zc,&prop,&propv)+fs; /* analysis of variants */
    *errnum=prop.kwx;
}


/******************************************************************************
  point_to_point()

  Note that point_to_point_two has become point_to_point 
  for drop-in interface to splat.cpp.  
  The new variable inputs,
  double enc_ncc_clcref, 
  double clutter_height, 
  double clutter_density, 
  double delta_h_diff, and 
  int mode_var)
  have been given fixed values below. 

pol:
0-Horizontal, 1-Vertical, 2-Circular

radio_climate:
1-Equatorial, 2-Continental Subtropical,
3-Maritime Tropical, 4-Desert, 5-Continental Temperate,
6-Maritime Temperate, Over Land, 7-Maritime Temperate,
Over Sea

conf, rel: .01 to .99

elev[]: [num points - 1], [delta dist(meters)],
[height(meters) point 1], ..., [height(meters) point n]

clutter_height      25.2 meters for compatibility with ITU-R P.1546-2.

clutter_density     1.0 for compatibility with ITU-R P.1546-2.

delta_h_diff        optional delta h for beyond line of sight. 90 m. average.
setting to 0.0 will default to use of original internal
use of delta-h for beyond line-of-sight range.

mode_var        set to 12; or to 1 for FCC ILLR;  see documentation

enc_ncc_clcref         clutter refractivity; 1000 N-units to match ITU-R P.1546-2

eno=eno_ns_surfref    atmospheric refractivity at sea level; 301 N-units nominal
(ranges from 250 for dry, hot day to 450 on hot, humid day]
(stabilizes near 301 in cold, clear weather)

errnum: 0- No Error.
1- Warning: Some parameters are nearly out of range.
Results should be used with caution.
2- Note: Default parameters have been substituted for
impossible ones.
3- Warning: A combination of parameters is out of range.
Results are probably invalid.
Other-  Warning: Some parameters are out of range.
Results are probably invalid.

 *****************************************************************************/
void point_to_point(const elev_t elev[], double tht_m, double rht_m, double eps_dielect, double sgm_conductivity, double eno_ns_surfref, double frq_mhz, int radio_climate, int pol, double conf, double rel, double *dbloss, char *strmode, int *errnum)
{
    prop_type   prop = {0};
    propv_type  propv = {0};
    propa_type  propa = {0};

    double zsys=0;
    double zc, zr;
    double eno, enso, q;
    long ja, jb, i, np;
    /* double dkm, xkm; */
    double tpd, fs;

    prop.hg[0]=tht_m;
    prop.hg[1]=rht_m;
    propv.klim=radio_climate;
    prop.kwx=0;
    prop.mdp=-1;
    prop.ptx=pol;
    prop.thera=0.0;
    prop.thenr=0.0;
    zc=qerfi(conf);
    zr=qerfi(rel);
    np=(long)elev[0];
    /* dkm=(elev[1]*elev[0])/1000.0; */
    /* xkm=elev[1]/1000.0; */
    eno=eno_ns_surfref;
    enso=0.0;
    q=enso;

    /* PRESET VALUES for Basic Version w/o additional inputs active */

    prop.encc = 1000.00;        /*  double enc_ncc_clcref preset  */          
    prop.cch = 22.5;           /* double clutter_height preset to ILLR calibration.;  
                                  use 25.3 for ITU-P1546-2 calibration */
    prop.cd = 1.00;                 /* double clutter_density preset */ 
    int mode_var = 1;         /* int mode_var set to 1 for FCC compatibility;
                                 normally, SPLAT presets this to 12 */
    prop.dhd= 0.0;            /* delta_h_diff preset */

    if (q<=0.0)
    {
        ja=(long)(3.0+0.1*elev[0]);  
        jb=np-ja+6;

        for (i=ja-1; i<jb; ++i)
            zsys+=elev[i];

        zsys/=(jb-ja+1);
        q=eno;
    }

    propv.mdvar=mode_var;
    qlrps(frq_mhz,zsys,q,pol,eps_dielect,sgm_conductivity,&prop);
    qlrpfl2(elev,propv.klim,propv.mdvar,&prop,&propa,&propv);
    tpd=sqrt((prop.he[0]-prop.he[1])*(prop.he[0]-prop.he[1])+(prop.dist)*(prop.dist));
    fs=32.45+20.0*log10(frq_mhz)+20.0*log10(tpd/1000.0);
    q=prop.dist-propa.dla;

    if (trunc(q)<0.0)
        strcpy(strmode,"L-o-S");
    else
    {
        if (trunc(q)==0.0)
            strcpy(strmode,"1_Hrzn");

        else if (trunc(q)>0.0)
            strcpy(strmode,"2_Hrzn");

        if (prop.dist<=propa.dlsa || prop.dist<=propa.dx)

            if(trunc(prop.dl[1])==0.0)
                strcat(strmode,"_Peak");

            else
                strcat(strmode,"_Diff");

        else if (prop.dist>propa.dx)
            strcat(strmode, "_Tropo");
    }

    *dbloss=avar(zr,0.0,zc,&prop,&propv)+fs;
    *errnum=prop.kwx;
}


/*************************************************************************************************
  point_to_pointMDH_two()

pol: 0-Horizontal, 1-Vertical
radio_climate: 1-Equatorial, 2-Continental Subtropical, 3-Maritime Tropical,
4-Desert, 5-Continental Temperate, 6-Maritime Temperate, Over Land,
7-Maritime Temperate, Over Sea
timepct, locpct, confpct: .01 to .99
elev[]: [num points - 1], [delta dist(meters)], [height(meters) point 1], ..., [height(meters) point n]
propmode:  Value   Mode
-1     mode is undefined
0     Line of Sight
5     Single Horizon, Diffraction
6     Single Horizon, Troposcatter
9     Double Horizon, Diffraction
10     Double Horizon, Troposcatter
errnum: 0- No Error.
1- Warning: Some parameters are nearly out of range.
Results should be used with caution.
2- Note: Default parameters have been substituted for impossible ones.
3- Warning: A combination of parameters is out of range.
Results are probably invalid.
Other-  Warning: Some parameters are out of range.
Results are probably invalid.
 *************************************************************************************************/
void point_to_pointMDH_two (const elev_t elev[], double tht_m, double rht_m,
                            double eps_dielect, double sgm_conductivity, double eno_ns_surfref, 
                            double enc_ncc_clcref, double clutter_height, double clutter_density, 
                            double delta_h_diff, double frq_mhz, int radio_climate, int pol, int mode_var, 
                            double timepct, double locpct, double confpct, 
                            double *dbloss, int *propmode, double *deltaH, int *errnum)

{

    prop_type   prop = {0};
    propv_type  propv = {0};
    propa_type  propa = {0};
    double zsys=0;
    double ztime, zloc, zconf;
    double eno, enso, q;
    long ja, jb, i, np;
    /* double dkm, xkm; */
    double fs;

    *propmode = -1;  /* mode is undefined */
    prop.hg[0] = tht_m;   
    prop.hg[1] = rht_m;
    propv.klim = radio_climate;
    prop.encc=enc_ncc_clcref;
    prop.cch=clutter_height;
    prop.cd=clutter_density;
    prop.dhd=delta_h_diff;
    prop.kwx = 0;
    prop.mdp = -1;
    prop.ptx=pol;
    prop.thera=0.0;
    prop.thenr=0.0;
    ztime = qerfi(timepct);
    zloc = qerfi(locpct);
    zconf = qerfi(confpct);
    np = (long)elev[0];
    /* dkm = (elev[1] * elev[0]) / 1000.0; */
    /* xkm = elev[1] / 1000.0; */
    eno = eno_ns_surfref;
    enso = 0.0;
    q = enso;

    /* PRESET VALUES for Basic Version w/o additional inputs active */

    prop.encc = 1000.00;        /*  double enc_ncc_clcref  */          
    prop.cch = 22.5;           /* double clutter_height */
    prop.cd = 1.00;              /* double clutter_density */ 
    mode_var = 1;         /* int mode_var set for FCC ILLR */

    if(q<=0.0)
    {
        ja =(long) (3.0 + 0.1 * elev[0]); /* to match addition of (long) */
        jb = np - ja + 6;
        for(i=ja-1;i<jb;++i)
            zsys+=elev[i];
        zsys/=(jb-ja+1);
        q=eno;
    }
    propv.mdvar=12;
    qlrps(frq_mhz,zsys,q,pol,eps_dielect,sgm_conductivity,&prop);
    qlrpfl2(elev,propv.klim,propv.mdvar,&prop,&propa,&propv);
    fs=32.45+20.0*log10(frq_mhz)+20.0*log10(prop.dist/1000.0);

    *deltaH = prop.dh;
    q = prop.dist - propa.dla;
    if(trunc(q)<0.0)
        *propmode = 0;  /* L-of-S  */
    else
    { if(trunc(q)==0.0)
        *propmode = 4;  /* 1-Hrzn */
        else if(trunc(q)>0.0)
            *propmode = 8;  /* 2-Hrzn */
        if(prop.dist<=propa.dlsa || prop.dist<=propa.dx)
            *propmode += 1; /* Diff */
        else if(prop.dist>propa.dx)
            *propmode += 2; /* Tropo */
    }
    *dbloss = avar(ztime, zloc, zconf, &prop, &propv) + fs;      /*avar(time,location,confidence) */
    *errnum = prop.kwx;
}


/*************************************************************************************************
  point_to_pointDH()

pol: 0-Horizontal, 1-Vertical
radio_climate: 1-Equatorial, 2-Continental Subtropical, 3-Maritime Tropical,
4-Desert, 5-Continental Temperate, 6-Maritime Temperate, Over Land,
7-Maritime Temperate, Over Sea
conf, rel: .01 to .99
elev[]: [num points - 1], [delta dist(meters)], [height(meters) point 1], ..., [height(meters) point n]
errnum: 0- No Error.
1- Warning: Some parameters are nearly out of range.
Results should be used with caution.
2- Note: Default parameters have been substituted for impossible ones.
3- Warning: A combination of parameters is out of range.
Results are probably invalid.
Other-  Warning: Some parameters are out of range.
Results are probably invalid.
 *************************************************************************************************/
void point_to_pointDH (const elev_t elev[], double tht_m, double rht_m, 
                       double eps_dielect, double sgm_conductivity, double eno_ns_surfref,
                       double enc_ncc_clcref, double clutter_height, double clutter_density, 
                       double delta_h_diff, double frq_mhz, int radio_climate, int pol, 
                       double conf, double rel, double loc,double *dbloss, double *deltaH,
                       int *errnum)
{

    char strmode[100];
    prop_type   prop = {0};
    propv_type  propv = {0};
    propa_type  propa = {0};
    double zsys=0;
    double zc, zr;
    double eno, enso, q;
    long ja, jb, i, np;
    /* double dkm, xkm; */
    double fs;

    prop.hg[0]=tht_m;   
    prop.hg[1]=rht_m;
    propv.klim = radio_climate;
    prop.encc=enc_ncc_clcref;
    prop.cch=clutter_height;
    prop.cd=clutter_density;
    prop.dhd=delta_h_diff;
    prop.kwx = 0;
    prop.mdp = -1;
    prop.ptx=pol;
    prop.thera=0.0;
    prop.thenr=0.0;
    zc = qerfi(conf);
    zr = qerfi(rel);
    np = (long)elev[0];
    /* dkm = (elev[1] * elev[0]) / 1000.0; */
    /* xkm = elev[1] / 1000.0; */
    eno = eno_ns_surfref;
    enso = 0.0;
    q = enso;

    /* PRESET VALUES for Basic Version w/o additional inputs active */

    prop.encc = 1000.00;        /*  double enc_ncc_clcref  */          
    prop.cch = 22.5;           /* double clutter_height */
    prop.cd =  1.00;              /* double clutter_density */ 

    if(q<=0.0)
    {
        ja = (long) (3.0 + 0.1 * elev[0]); /* to match KD2BD addition of (long)  */
        jb = np - ja + 6;
        for(i=ja-1; i<jb; ++i)
            zsys+=elev[i];
        zsys/=(jb-ja+1);
        q=eno;
    }
    propv.mdvar=12;
    qlrps(frq_mhz,zsys,q,pol,eps_dielect,sgm_conductivity,&prop);
    qlrpfl2(elev,propv.klim,propv.mdvar,&prop,&propa,&propv);
    fs=32.45+20.0*log10(frq_mhz)+20.0*log10(prop.dist/1000.0);
    *deltaH = prop.dh;
    q = prop.dist - propa.dla;
    if(trunc(q)<0.0)
        strcpy(strmode,"Line-Of-Sight Mode");
    else
    { if(trunc(q)==0.0)
        strcpy(strmode,"Single Horizon");
        else if(trunc(q)>0.0)
            strcpy(strmode,"Double Horizon");
        if(prop.dist<=propa.dlsa || prop.dist<=propa.dx)
            strcat(strmode,", Diffraction Dominant");
        else if(prop.dist>propa.dx)
            strcat(strmode, ", Troposcatter Dominant");
    }
    *dbloss = avar(zr,0.0,zc,&prop,&propv)+fs;      /*avar(time,location,confidence) */
    *errnum = prop.kwx;
}


/********************************************************
 * Area Mode Calculations
 *
 * This is obsolete. It was a less compute-intensive way
 * of generating results for a given area without doing
 * multiple point-to-point calculations, trading accuracy
 * for speed.
 ********************************************************/

/*************************************************************************************************
  area()

pol: 0-Horizontal, 1-Vertical
TSiteCriteria, RSiteCriteria:
0 - random, 1 - careful, 2 - very careful

radio_climate: 1-Equatorial, 2-Continental Subtropical, 3-Maritime Tropical,
4-Desert, 5-Continental Temperate, 6-Maritime Temperate, Over Land,
7-Maritime Temperate, Over Sea
ModVar: 0 - Single: pctConf is "Time/Situation/Location", pctTime, pctLoc not used
1 - Individual: pctTime is "Situation/Location", pctConf is "Confidence", pctLoc not used
2 - Mobile: pctTime is "Time/Locations (Reliability)", pctConf is "Confidence", pctLoc not used
3 - Broadcast: pctTime is "Time", pctLoc is "Location", pctConf is "Confidence"
pctTime, pctLoc, pctConf: .01 to .99
errnum: 0- No Error.
1- Warning: Some parameters are nearly out of range.
Results should be used with caution.
2- Note: Default parameters have been substituted for impossible ones.
3- Warning: A combination of parameters is out of range.
Results are probably invalid.
Other-  Warning: Some parameters are out of range.
Results are probably invalid.
NOTE: strmode is not used at this time.

See ITWOM-SUB-ROUTINES.pdf p460
 *************************************************************************************************/
void area(long ModVar, double deltaH, double tht_m, double rht_m, double dist_km, int TSiteCriteria, int RSiteCriteria, double eps_dielect, double sgm_conductivity, double eno_ns_surfref, double enc_ncc_clcref, double clutter_height, double clutter_density, double delta_h_diff, double frq_mhz, int radio_climate, int pol, int mode_var,
          double pctTime, double pctLoc, double pctConf, double *dbloss, char *strmode, int *errnum)
{

    prop_type prop = {0};
    propv_type propv = {0};
    propa_type propa = {0};
    double zt, zl, zc, xlb;
    double fs;
    long ivar;
    double eps, eno, sgm;
    long ipol;
    int kst[2];

    kst[0]=(int)TSiteCriteria;
    kst[1]=(int)RSiteCriteria;
    zt=qerfi(pctTime/100.0);
    zl=qerfi(pctLoc/100.0);
    zc=qerfi(pctConf/100.0);
    eps=eps_dielect;
    sgm=sgm_conductivity;
    eno=eno_ns_surfref;
    prop.dh=deltaH;
    prop.hg[0]=tht_m;
    prop.hg[1]=rht_m;
    propv.klim=(long)radio_climate;
    prop.encc=enc_ncc_clcref;
    prop.cch=clutter_height;
    prop.cd=clutter_density;
    prop.dhd=delta_h_diff;
    prop.ens=eno;
    prop.kwx=0;
    ivar=(long)ModVar;
    ipol=(long)pol;
    qlrps(frq_mhz, 0.0, eno, ipol, eps, sgm, &prop);
    qlra(kst, propv.klim, ivar, &prop, &propv);

    lrprop2(dist_km*1000.0, &prop, &propa);
    fs=32.45+20.0*log10(frq_mhz)+20.0*log10(prop.dist/1000.0);

    xlb=fs+avar(zt, zl, zc, &prop, &propv);
    *dbloss=xlb;

    if (prop.kwx==0)
        *errnum=0;
    else
        *errnum=prop.kwx;
}


double ITMAreadBLoss(long ModVar, double deltaH, double tht_m, double rht_m,
                     double dist_km, int TSiteCriteria, int RSiteCriteria, 
                     double eps_dielect, double sgm_conductivity, double eno_ns_surfref,
                     double enc_ncc_clcref,double clutter_height,double clutter_density,
                     double delta_h_diff, double frq_mhz, int radio_climate, int pol, 
                     int mode_var, double pctTime, double pctLoc, double pctConf)
{
    char strmode[200];
    int errnum;
    double dbloss;
    area(ModVar,deltaH,tht_m,rht_m,dist_km,TSiteCriteria,RSiteCriteria, 
         eps_dielect,sgm_conductivity,eno_ns_surfref,
         enc_ncc_clcref,clutter_height,clutter_density,
         delta_h_diff,frq_mhz,radio_climate,pol,mode_var,pctTime,
         pctLoc,pctConf,&dbloss,strmode,&errnum);
    return dbloss;
}

double ITWOMVersion()
{
    return 3.0;
}
