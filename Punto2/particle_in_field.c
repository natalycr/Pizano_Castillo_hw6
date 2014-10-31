#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define e 1.602176565e-19 //Elementary charge (Coulomb)
#define mpr 1.672621777e-27 // Proton mass (kg)
#define mel 9.10938291e-31 // Electron mass (kg)
#define c 299792458 // speed of light (m/s)
#define B0 3.07e-5 // Tesla
#define Re 6378137 // meter (Earth radius)
#define tiempo  100.0 // cantidad de segundos que se desea modelar
#define pi 3.14159

double r(double x ,double y, double z);

double v_x(double x2, double y2, double z2);
double v_y(double x2, double y2, double z2);
double v_z(double x2, double y2, double z2);

double deri_vx(double x, double y, double z, double x2, double y2, double z2, double cte);
double deri_vy(double x ,double y, double z, double x2 ,double y2 ,double z2, double cte);
double deri_vz(double x ,double y, double z, double x2 ,double y2 ,double z2, double cte);

void runge_kutta4(double X_old, double Y_old,double Z_old, double X2_old,double Y2_old, double Z2_old, double h, double *x, double *y, double *z,double*x2 ,double *y2 ,double *z2, double cte);

int main (int argc, char **argv){ // (){//
    
    FILE *arch;
    double Ek, pitch, EkJ, vo, gamma, n, h, tmax, tmin, cte, r, x, y, z, x2, y2, z2 , v0; //, num_puntos
    double  *X, *Y, *Z, *X2, *Y2, *Z2, *T;
    int i, j;

    Ek = atof(argv[1]); // Energìa en MeV 1000; //
    pitch= atof(argv[2]); // ángulo 30; //
    EkJ = e*Ek; //Energìa en Joules
    tmin=0.0; //segundos
    h=0.001; // paso
    n=(tiempo-tmin)/n; //numero puntos 
    
    //Abrir archivo para escribir datos
    char nombre[100];

    strcpy(nombre, "trayectoria_");
    strcat(nombre,  argv[1] );
    strcat(nombre,"_");
    strcat(nombre, argv[2]);
    strcat(nombre,".dat");

    // velocidad - gamma - cte:
    vo = sqrt(1-(1-pow(((Ek/mpr)+1),2))); //(EkJ+(mpr*(pow(c,2))))/sqrt(mpr*EkJ*(pow(c,2))*(EkJ+(2*(pow(c,2)))));
    gamma = 1/(sqrt(1-((pow(vo,2))/(pow(c,2)))));
    cte=(e*B0*pow(Re,3))/(gamma*mpr);


     /*
     Condiciones iniciales 
     */


    T=malloc(2*sizeof(float));
    R=malloc(3*sizeof(float)); // Guarda las posiciones  en x,y,z
    V=malloc(3*sizeof(float)); // guarda las velo en vx, vy,vz 

    //matrices  
     R[0]=x;
     R[1]=y;
     R[2]=z;
     V[0]=x2;
     V[1]=y2;
     V[2]=z2;


    x=2*Re;
    y=0.0;
    z=0.0;
    vx=0.0;
    vy = v0*sin(pitch*pi/180.0);
    vz = v0*cos(pitch*pi/180.0);


    for(i=1;i<n;i++){
      x=0.0;
      y=0.0;
      z=0.0;
      x2=0.0;
      y2=0.0;
      z2=0.0;

      //runge_kutta4( X[i-1],  Y[i-1], Z[i-1],  X2[i-1],  Y2[i-1], Z2[i-1], h, &x, &y, &z, &x2 , &y2 , &z2, cte);

      X[i]=x;
      Y[i]=y;
      Z[i]=z;
      X2[i]=x2;
      Y2[i]=y2;
      Z2[i]=z2;
}
    printf("%s  despues hoa  \n", nombre);

    /* 
       genera archivo
    */

    arch =fopen(nombre,"w");
    
    for (j=0;j<n ;j++){

      fprintf( arch, "%f \t %f \t %f \t %f \n", T[j] , X[j], Y[j], Z[j]) ;
    }
    fclose(arch);

    return 0;
}
