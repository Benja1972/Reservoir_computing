
#include<complex>
#include<math.h>
#include<cmath>
#include<cstdlib>
#include<iostream>

/****** uncomment the following line to run generated source code directly ***************/
//#define MANUAL

#define false 0;
#define true 1;

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

#define SQUARE(X) ((X) * (X))

double Heavi(double tin) {
    if(tin>=0)
        return 1.0;
    else
        return 0.0;
}

typedef std::complex<double> Complex;
const Complex ii(0.0, 1.0);
double t;
double *pt_t_ar;
long unsigned int SIM_n;
long unsigned int SIM_size;
#define ABSx(X) (fabs(X))

/* positions of lags in the arrays */
long unsigned int n_c3c07dadf2e76469f0fae41fb666e561;
long unsigned int n_2d41adf0c5db14956faaa058a97dcb04;
/* arrays holding the variabls */
double *pt_x_ar;
double *pt_x_Var;

inline double hermite_x(const double &, const double &, const double &, const double &, 
                                                   const double &, const double &, const double &);
inline double dt_hermite_x(const double &, const double &, const double &, const double &, 
                                                   const double &, const double &, const double &);


            
double interp_x(double t, long unsigned int *n0) 
{
    while(pt_t_ar[*n0+1] < t)
        (*n0)++;
    while(pt_t_ar[*n0] > t) 
        (*n0)--;

    return hermite_x(t, pt_t_ar[*n0],   pt_x_ar[*n0],   pt_x_Var[*n0], 
                              pt_t_ar[*n0+1], pt_x_ar[*n0+1], pt_x_Var[*n0+1]);
}

//interpolation of d(x)/dt
double dt_interp_x(double t, long unsigned int *n0) 
{
    while(pt_t_ar[*n0+1] < t)
        (*n0)++;
    while(pt_t_ar[*n0] > t) 
        (*n0)--;

    return dt_hermite_x(t, pt_t_ar[*n0],   pt_x_ar[*n0],   pt_x_Var[*n0], 
                              pt_t_ar[*n0+1], pt_x_ar[*n0+1], pt_x_Var[*n0+1]);
}

//hermite interpolation
inline double hermite_x(const double &t, const double &tn, const double &Xn, const double &Vn, 
                                     const double &tnp1, const double &Xnp1, const double &Vnp1) 
{
    double h = tnp1 - tn;
    double s = (t - tn) / h;

    return (1.0 + 2 * s) * SQUARE(s-1.0) * Xn + (3.0 - 2 * s) * SQUARE(s) * Xnp1 + h * s * SQUARE(s-1.0) *Vn + h * (s - 1) * SQUARE(s) *Vnp1;
}

inline double dt_hermite_x(const double &t, const double &tn, const double &Xn, const double &Vn, 
                                     const double &tnp1, const double &Xnp1, const double &Vnp1) 
{
    double h = tnp1 - tn;
    double s = (t - tn) / h;

    return (1.0-4.0*s+3.0*SQUARE(s))*Vn + s*(3.0*s-2.0)*Vnp1 + 6.0*(s-1.0)*s*(Xn-Xnp1)/h;
}

        


#ifdef MANUAL
int main()
{
double dtmin = 0.0001;
double dtmax = 0.1;
double tfinal = 100.0;
double RelTol = 1.0E-3;
double AbsTol = 1.0E-6;
int chunk = 10000;
int nstart = 101;
double dt0 = 0.01;
double maxdelay = ...; // Set the maximum delay here !!!
long unsigned int MaxIter = 10000000;
int NumOfDiscont = 4;
double discont[4] = {maxdelay, 2*maxdelay, 3*maxdelay, tfinal};
double PARtau = 0.0122;
double PARbeta = 1.2;
double PARdT = 0.4;
double PART = 24.35;
double PARPhi0 = 0.722566310326;

#endif

t = 0.0; 
long unsigned int i;
long unsigned int NumberOfMinSteps = 0;
int nextdsc = 0;
int hitdsc = false;
double dist;
int TakingMinStep = 0;
SIM_n = nstart; 
SIM_size = SIM_n + chunk;
double RelErr;
double thresh = AbsTol/RelTol;
double dt = dt0;
srand((unsigned)RSEED);
n_c3c07dadf2e76469f0fae41fb666e561 = SIM_n;
n_2d41adf0c5db14956faaa058a97dcb04 = SIM_n;


pt_x_ar  = (double*) malloc((SIM_n+chunk) * sizeof(double));
pt_x_Var = (double*) malloc((SIM_n+chunk) * sizeof(double));

            
pt_t_ar = (double *) malloc((SIM_n+chunk) * sizeof(double));

#ifndef MANUAL
        
for(i = 0; i < nstart+1; i++) {
    pt_x_ar[i]  = histx_ar[i];
    pt_x_Var[i] = Vhistx_ar[i];
}
            
for(i = 0; i < nstart+1; i++) 
    pt_t_ar[i] = Thist_ar[i];
#endif

#ifdef MANUAL
for(i = 0; i < nstart+1; i++) 
    pt_t_ar[i]  = -maxdelay*(nstart-i)/nstart; 

/* set the history here when running generated code directly */
 
for(i = 0; i < nstart+1; i++) {
    pt_x_ar[i]  = 0.2; // history value for the variable
    pt_x_Var[i] = 0.0; // history of the derivatives (0.0 if constant history)
}

#endif

        
double k1x, k2x, k3x, k4x;
double TEMPx;
double ERRORx;
//k1 need to be calculated only for the first step
//due to the FSAL property k1(n+1)=k4(n)
	k1x = (1/ PARtau ) *( PARbeta *pow(cos(interp_x(t - PART, &n_c3c07dadf2e76469f0fae41fb666e561) -interp_x(t - PART - PARdT, &n_2d41adf0c5db14956faaa058a97dcb04) + PARPhi0 ),2) - pt_x_ar[SIM_n] );

while((t <= tfinal || hitdsc) && SIM_n-nstart <= MaxIter) {

	t += dt * 0.5;
	k2x = (1/ PARtau ) *( PARbeta *pow(cos(interp_x(t - PART, &n_c3c07dadf2e76469f0fae41fb666e561) -interp_x(t - PART - PARdT, &n_2d41adf0c5db14956faaa058a97dcb04) + PARPhi0 ),2) - (pt_x_ar[SIM_n] + 0.5 * dt * k1x) );
	t += dt * 0.25;
	k3x = (1/ PARtau ) *( PARbeta *pow(cos(interp_x(t - PART, &n_c3c07dadf2e76469f0fae41fb666e561) -interp_x(t - PART - PARdT, &n_2d41adf0c5db14956faaa058a97dcb04) + PARPhi0 ),2) - (pt_x_ar[SIM_n] + 0.75 * dt * k2x) );
	TEMPx = pt_x_ar[SIM_n] + dt * 1.0/9.0 * (2.0*k1x + 3.0*k2x + 4.0*k3x);
	t += dt * 0.25;
	k4x = (1/ PARtau ) *( PARbeta *pow(cos(interp_x(t - PART, &n_c3c07dadf2e76469f0fae41fb666e561) -interp_x(t - PART - PARdT, &n_2d41adf0c5db14956faaa058a97dcb04) + PARPhi0 ),2) - TEMPx );
	ERRORx = dt/72.0 * (-5.0*k1x + 6.0*k2x + 8.0*k3x - 9.0*k4x);
	RelErr = 0.0;
	ERRORx = ERRORx/MAX( MAX(ABSx(pt_x_ar[SIM_n]), ABSx(TEMPx)), thresh);

	RelErr = MAX(RelErr, ABSx(ERRORx));
	if(RelErr <= RelTol || TakingMinStep ) {
		pt_x_ar[SIM_n+1] = TEMPx;
		pt_x_Var[SIM_n+1] = k4x;
		k1x = k4x; //FSAL

        pt_t_ar[SIM_n+1] = t;
        SIM_n++;
        if(SIM_n - nstart > MaxIter) {
            std::cerr << "Warning: MaxIter reached! EndTime of simulation: " << t << std::endl << std::flush;
        }
        //std::cout << "pass: " << pow(RelTol/RelErr, 0.3333333333333) << std::endl;
        dt = dt * MAX(0.5, 0.8*pow(RelTol/RelErr, 0.3333333333333));

        // hit discontinuities
        hitdsc = false;
        if(nextdsc<NumOfDiscont) {
            dist = discont[nextdsc] - t;
            //std::cout << t << "	" << discont[nextdsc] << "	" << dist << std::endl;
            if(dist <= MIN(1.1*dt, dtmax)) {
                dt = dist;
                nextdsc++;
                hitdsc=true;
            }
            else if(dist <= 2*dt) {
                dt = 0.5 * dist;
            }
        }

    } else {
        //std::cout << "not passed" << std::endl;
        t-= dt;
        dt = 0.5*dt;
    }
    
    if(dt < dtmin && !hitdsc)  {
        //TODO: fix this    
        //if(NumberOfMinSteps == 0)
        //    std::cerr << "Warning: step size very small" << std::endl;
        NumberOfMinSteps++;
        TakingMinStep = true;
        dt = dtmin;
    } else {
        TakingMinStep = false;
    }
    
    if(dt > dtmax) 
        dt = dtmax;

            
    //grow arrays if they are too small
    if(SIM_n+1 == SIM_size)
    {
        SIM_size += chunk;

        double *p;
        p = (double *) realloc(pt_t_ar, (SIM_size) * sizeof(double));
        if (!p) {
            std::cout << "realloc fail, there is probably not enough memory" << std::endl;
            exit(1);
        } 
        pt_t_ar = p;
            
        double *px;
        px = (double *) realloc(pt_x_ar, (SIM_size) * sizeof(double));
        if (!px) {
            std::cout << "realloc fail, there is probably not enough memory" << std::endl;
            exit(1);
        } 
        pt_x_ar = px;

        px = (double *) realloc(pt_x_Var, (SIM_size) * sizeof(double));
        if (!px) {
            std::cout << "realloc fail, there is probably not enough memory" << std::endl;
            exit(1);
        } 
        pt_x_Var = px;
            
    }
}
double *p;
p = (double *) realloc(pt_t_ar, (SIM_n) * sizeof(double));
if (!p) 
    std::cout << "realloc fail when shrinking arrays" << std::endl;
else
    pt_t_ar = p;


double *px;
px = (double *) realloc(pt_x_ar, (SIM_n) * sizeof(double));
if (!px) 
    std::cout << "realloc fail when shrinking arrays" << std::endl;
else
    pt_x_ar = px;

px = (double *) realloc(pt_x_Var, (SIM_n) * sizeof(double));
if (!px) 
    std::cout << "realloc fail when shrinking arrays" << std::endl;
else
    pt_x_Var = px;

            
#ifndef MANUAL
PyObject *erg;
erg = PyDict_New();

PyObject *myarray;
double *array_buf;
npy_intp dim = SIM_n;

myarray = PyArray_SimpleNew(1, &dim, NPY_DOUBLE);

array_buf = (double *) PyArray_DATA(myarray);

for(i = 0; i < SIM_n; i++) 
    *array_buf++ = pt_t_ar[i];
free(pt_t_ar);

PyDict_SetItem(erg, PyString_InternFromString("t"), myarray);
Py_DECREF(myarray);

        
PyObject *myarray_x;
//PyObject *myarray_Vx;
double *array_bufx; 
//double *array_bufVx;

myarray_x = PyArray_SimpleNew(1, &dim, NPY_DOUBLE);
//myarray_Vx = PyArray_SimpleNew(1, &dim, NPY_DOUBLE);

array_bufx = (double *) PyArray_DATA(myarray_x);
//array_bufVx = (double *) PyArray_DATA(myarray_Vx);

for(i = 0; i < SIM_n; i++) {
    *array_bufx++ = pt_x_ar[i];
    //*array_bufVx++ = pt_x_Var[i];
}
free(pt_x_ar);
free(pt_x_Var);

PyDict_SetItem(erg, PyString_InternFromString("x"), myarray_x);
//PyDict_SetItem(erg, PyString_InternFromString("Vx"), myarray_Vx);
Py_DECREF(myarray_x);
//Py_DECREF(myarray_Vx);

            
if(NumberOfMinSteps>1)
    std::cerr << "Number of Minimum steps taken: " << NumberOfMinSteps << std::endl;

return_val =  erg;
Py_DECREF(erg);
#endif
 
#ifdef MANUAL
for(i = 0; i < SIM_n; i++) 
    std:: cout << pt_t_ar[i] << "	" << pt_x_ar[i] << std::endl;
    return 0;
}
#endif

// This hash is included to detect a change in
// the support code and recompile in this case
// a17eca6d2b3fd09be573424776663c09
