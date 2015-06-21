#include <stdio.h>
#include <math.h>
#include <complex.h>

void diff(double _Complex *e_org, double _Complex *e_res, int rh_snum, double rh_max);
void rho_2nd_diff(double _Complex *e_org, double _Complex *e_res, int rh_snum, double rh_max);

int main(int argc, char *argv[]){
  //default parameters
  double a = 0.002; //radius of probe beam
  double b = 3; //radius of thermal lens
  double rh_max = 4*a*a; //maximum radius square of calculated space
  double z_max = 0.2; //maximum z of calculated space
  int rh_snum = 1000; //step number of radius square
  int z_snum = 2000; //step number of z
  double k = 2*M_PI/(8.0*pow(10, -7)); //wave number of probe beam
  double p = b*b*k/2/0.01; //phase change at the radius of thermal lens (*pi)

  int i, j;
  double _Complex e_1[rh_snum];
  double _Complex e_2[rh_snum];
  for(i=0; i<rh_snum; ++i){
    double r = sqrt(i*(rh_max/rh_snum));
    //e_1[i] = cexp(-pow(r/a, 2) + I*p*M_PI*exp(-pow(r/b, 2)));
    e_1[i] = cexp(-(r/a)*(r/a)-I*k*r*r/2/0.2*0);
  }

  for(i=0; i<z_snum; ++i){
    rho_2nd_diff(e_1, e_2, rh_snum, rh_max);
    for(j=0; j<rh_snum; ++j){
      e_1[j] += 2*I*e_2[j]/k*(z_max/z_snum);
      printf("%f\t", cabs(e_1[j]));
    }
    printf("\n");
  }
  return 0;
}

void diff(double _Complex *e_org, double _Complex *e_res, int rh_snum, double rh_max){
  int i;
  double drh = rh_max/rh_snum;
  e_res[0] = (e_org[1]-e_org[0])/drh;
  e_res[rh_snum-1] = (e_org[rh_snum-1]-e_org[rh_snum-2])/drh;
  for(i=1; i<rh_snum-1; ++i){
    e_res[i] = (e_org[i+1] - e_org[i-1])/(2*drh);
  }
}

void rho_2nd_diff(double _Complex *e_org, double _Complex *e_res, int rh_snum, double rh_max){
  int i;
  double drh = rh_max/rh_snum;
  double _Complex e_tmp[rh_snum];
  diff(e_org, e_tmp, rh_snum, rh_max);
  for(i=0; i<rh_snum; ++i){
    e_tmp[i] *= drh*i;
  }
  diff(e_tmp, e_res, rh_snum, rh_max);
}
