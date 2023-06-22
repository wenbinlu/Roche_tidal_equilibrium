
#ifndef MYPARA_H_
#define MYPARA_H_

// defining global variables

extern double npoly;  // polytropic index
extern double pi;
extern double Lmax; // maximum domain size [machine units]
extern int Nresz;  // resolution in the z direction (same as in y, half of x)
extern char fdir[];  // directory to save the data

extern double ***rhoarr;
extern double Mstar;   // total stellar mass [machine units]

#endif  /* MYPARA_H_ */
