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

//cont=(e*B0*pow(Re,3))/(gamma*mpr*pow(r,5));

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
    double Ek, pitch, EkJ, vo, gamma, n, h, tmax, tmin, cte, r, x, y, z, x2, y2,z2 , v0; //, num_puntos
    double  *X, *Y, *Z, *X2, *Y2, *Z2, *T;
    int i, j;

    Ek = atof(argv[1]); // Energìa en MeV 1000; //
    pitch= atof(argv[2]); // ángulo 30; //
    EkJ = e*Ek*10E6; //Energìa en Joules
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

    //file= malloc(50*sizeof(char));

    /*
    if(!arch)
    {
        printf("problems opening the file %s\n", file);
        exit(1);
    }
    */
 
    // velocidad - gamma - cte:
    vo = (EkJ+(mpr*(pow(c,2))))/sqrt(mpr*EkJ*(pow(c,2))*(EkJ+(2*(pow(c,2)))));
    gamma = 1/(sqrt(1-((pow(vo,2))/(pow(c,2)))));
    cte=(e*B0*pow(Re,3))/(gamma*mpr);

     /*
     Condiciones iniciales 
     */

    //num_puntos=100.0;

    T=malloc(n*sizeof(double));
    X=malloc(n*sizeof(double));
    Y=malloc(n*sizeof(double));
    Z=malloc(n*sizeof(double));

    X2=malloc(n*sizeof(double));
    Y2=malloc(n*sizeof(double));
    Z2=malloc(n*sizeof(double));
    X[0]=2*Re;

    printf("%f\n", X[0]);

    for (i=1;i<n;i++){
      T[i]=i*(1.0/n);
      X[i]=0.0;
      Y[i]=0.0;
      Z[i]=0.0;
      X2[i]=0.0;
      Y2[i]=0.0;
      Z2[i]=0.0;
      //printf("%f %f dentro  for \n", X[i], Y[i]);
    }

    printf("%s  despues hoa  \n", nombre);


    /*
    if(pitch==90){
      Z2[0] = 0;  // ahorrar calculos innecesarios
    }
    else
      {
        Z2[0] = v0*cos(pitch*pi/180);
      }
    */

    Y2[0] = v0*sin(pitch*pi/180);
    Z2[0] = v0*cos(pitch*pi/180);

    for(i=1 ; i<n; i++){
      
      x=X[i];
      y=Y[i];
      z=Z[i];
      x2=X2[i];
      y2=Y2[i];
      z2=Z2[i];

      runge_kutta4( X[i-1],  Y[i-1], Z[i-1],  X2[i-1],  Y2[i-1], Z2[i-1], h, &x, &y, &z, &x2 , &y2 , &z2, cte);

      X[i]=x;
      Y[i]=y;
      Z[i]=z;
      X2[i]=x2;
      Y2[i]=y2;
      Z2[i]=z2;
    }

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
void runge_kutta4(double X_old, double Y_old,double Z_old, double X2_old,double Y2_old, double Z2_old, double h, double *x, double *y, double *z,double*x2 ,double *y2 ,double *z2, double cte){

  double k1_x, k1_y, k1_z, k1_vx, k1_vy, k1_vz, x_1, y_1, z_1, x2_1, y2_1, z2_1, k2_x, k2_y, k2_z, k2_vx, k2_vy, k2_vz, x_2, y_2, z_2, x2_2, y2_2, z2_2, k3_x, k3_y, k3_z, k3_vx, k3_vy, k3_vz, x_3, y_3, z_3, x2_3, y2_3, z2_3, k4_x, k4_y, k4_z, k4_vx, k4_vy, k4_vz, x_4, y_4, z_4, x2_4, y2_4, z2_4, promk_x, promk_y, promk_z,  promk_vx,  promk_vy, promk_vz;

  k1_x= v_x(X2_old, Y2_old, Z2_old);
  k1_vx = deri_vx( X_old,  Y_old, Z_old,X2_old, Y2_old, Z2_old, cte);

  k1_y= v_y(X2_old, Y2_old, Z2_old);
  k1_vy = deri_vy( X_old,  Y_old, Z_old,X2_old, Y2_old, Z2_old, cte);
  
  k1_z= v_z(X2_old, Y2_old, Z2_old);
  k1_vz = deri_vz( X_old,  Y_old, Z_old,X2_old, Y2_old, Z2_old, cte);

  //_________________primer paso __________________

  x_1= X_old + (h/2.0)* k1_x;
  x2_1=X2_old + (h/2.0)* k1_vx;

  y_1= Y_old + (h/2.0)* k1_y;
  y2_1= Y2_old + (h/2.0)* k1_vx;

  z_1=Z_old + (h/2.0)* k1_z;
  z2_1=Z2_old + (h/2.0)* k1_vx;

  k2_x=v_x(x_1, y_1, z_1);
  k2_vx=deri_vx(x_1, y_1, z_1, x2_1, y2_1, z2_1, cte);

  k2_y=v_y(x_1, y_1, z_1);
  k2_vy=deri_vy(x_1, y_1, z_1, x2_1, y2_1, z2_1, cte);

  k2_z=v_z(x_1, y_1, z_1);
  k2_vz=deri_vz(x_1, y_1, z_1, x2_1, y2_1, z2_1, cte);
  
  //_________________Segundo  paso __________________

  x_2= X_old + (h/2.0)* k2_x;
  x2_2= X2_old + (h/2.0)* k2_vx;

  y_2= Y_old + (h/2.0)* k2_y;
  y2_2= Y2_old + (h/2.0)* k2_vx;

  z_2=Z_old + (h/2.0)* k2_z;
  z2_2=Z2_old + (h/2.0)* k2_vx;

  k3_x=v_x(x_2, y_2, z_2);
  k3_vx=deri_vx(x_2, y_2, z_2, x2_2, y2_2, z2_2, cte);

  k3_y=v_y(x_2, y_2, z_2);
  k3_vy=deri_vy(x_2, y_2, z_2, x2_2, y2_2, z2_2, cte);

  k3_z=v_z(x_2, y_2, z_2);
  k3_vz=deri_vz(x_2, y_2, z_2, x2_2, y2_2, z2_2, cte);

  //_________________Tercer paso __________________

  x_3= X_old + (h/2.0)* k3_x;
  x2_3= X2_old + (h/2.0)* k3_vx;

  y_3= Y_old + (h/2.0)* k3_y;
  y2_3= Y2_old + (h/2.0)* k3_vx;

  z_3 = Z_old + (h/2.0)* k3_z;
  z2_3 = Z2_old + (h/2.0)* k3_vx;

  k4_x = v_x(x_3, y_3, z_3);
  k4_vx = deri_vx(x_3, y_3, z_3, x2_3, y2_3, z2_3, cte);

  k4_y = v_y(x_3, y_3, z_3);
  k4_vy = deri_vy(x_3, y_3, z_3, x2_3, y2_3, z2_3, cte);

  k4_z = v_z(x_3, y_3, z_3);
  k4_vz = deri_vz(x_3, y_3, z_3, x2_3, y2_3, z2_3, cte);

  //_________________Cuarto paso __________________

  promk_x=(1.0/6.0)*(k1_x + 2.0*k2_x + 2.0*k3_x + k4_x);
  promk_y=(1.0/6.0)*(k1_y + 2.0*k2_y + 2.0*k3_y + k4_y);
  promk_z=(1.0/6.0)*(k1_z + 2.0*k2_z + 2.0*k3_z + k4_z);

  promk_vx=(1.0/6.0)*(k1_vx + 2.0*k2_vx + 2.0*k3_vx + k4_vx);
  promk_vy=(1.0/6.0)*(k1_vy + 2.0*k2_vy + 2.0*k3_vy + k4_vy);
  promk_vz=(1.0/6.0)*(k1_vz + 2.0*k2_vz + 2.0*k3_vz + k4_vz);

  *x= X_old + h* promk_x;
  *y= Y_old + h* promk_y;
  *z= Z_old + h* promk_z;

  *x2= X2_old + h* promk_vx; 
  *y2= Y2_old + h* promk_vy; 
  *z2= Z2_old + h* promk_vz;
}

double r(double x ,double y, double z){
  double denRa;
  denRa=pow(x,2)+pow(y,2)+pow(z,2);
  return pow(denRa, 0.5);
}
/*
  velocidad 
*/

//x1
double v_x(double x2, double y2, double z2){ 
  //dx/dt pero no depende de t por eso no lo coloque 
  return x2;
}
//y1
double v_y(double x2, double y2, double z2){ 
  //dx/dt pero no depende de t por eso no lo coloque 
  return y2;
}
//z1
double v_z(double x2, double y2, double z2){ 
  //dx/dt pero no depende de t por eso no lo coloque 
  return z2;
}

/*
 aceleracion 
*/

//x2
double deri_vx(double x, double y, double z, double x2, double y2, double z2, double cte){ 
  //dx/dt pero no depende de t por eso no lo coloque 
  double rp;
  rp=r(x,y,z);// radio de la particula
  return (cte/(pow(rp,5)))*(x2*(2*pow(z,2)-pow(x,2)-pow(y,2))- (3*z2*y*z));
}


//y2
double deri_vy(double x ,double y, double z, double x2 ,double y2 ,double z2, double cte){ 
  double rp;
  //dy/dt pero no depende de t por eso no lo coloque 
  rp=r(x,y,z);// radio de la particula
  return (cte/(pow(rp,5)))*((z2*3*x*z)-(x2*(2*pow(z,2)-pow(x,2)-pow(y,2))));
}

//z2
double deri_vz(double x ,double y, double z, double x2 ,double y2 ,double z2, double cte){
  //dy/dt pero no depende de t por eso no lo coloque 
  double rp;
  rp=r(x,y,z);// radio de la particula
  return (cte/(pow(rp,5)))*((x2*3*y*z)-(y2*3*x*z));
}
