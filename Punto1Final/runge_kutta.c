# include <stdio.h>
# include  <stdlib.h>
# include <math.h>
#include <string.h>

double x_der1(double x, double y, double A, double B);
double y_der1(double x ,double y, double  C, double  D);
void runfge_kutta4orden(double x_old, double y_old, double A, double B, double C, double D , double h, double *x_new, double *y_new);

double main (int argc, char **argv){
  FILE *archivo;
  double x0, y0, h, x, y;
  double  *X, *Y, *T;
  double  A,B,C,D;
  int i, j, num_puntos;

  /*
    Condiciones  Iniciales
  */

  x0=atof(argv[1]); // 30.0; //
  y0=atof(argv[2]); //20.0; //
  h=1E-4; // tamano del paso 

  //printf("%f,%f \n", x0,y0);

  A=20.0;
  B=1.0;
  C=30.0;
  D=1.0;
  num_puntos=1.0/h;
  /*
    Arrays 
  */

  T=malloc(num_puntos*sizeof(double));
  X=malloc(num_puntos*sizeof(double));
  Y=malloc(num_puntos*sizeof(double));

  for (i=1;i<num_puntos;i++){
    T[i]=i*(1.0/num_puntos);
    X[i]=0.0;
    Y[i]=0.0;
    //printf("%f t \n", T[i]);

}

  X[0]=x0;
  Y[0]=y0;

  for(i=1;i<num_puntos; i++){
    
    x= 0.0;
    y= 0.0;
    //printf("%f %f %f linea ant \n", T[i-1], X[i-1], Y[i-1] );
    
    runfge_kutta4orden(X[i-1], Y[i-1], A , B , C, D, h, &x, &y);
    
    X[i]=x;
    Y[i]=y;
    //printf("%f %f X[i]\n", X[i], Y[i]);
}
  //printf("fordep runge");

  /*

    Archivo que  guarfa los datos de la evolicion
*/

  char nombre[40];

  strcpy(nombre, "poblaciones_");
  strcat(nombre,  argv[1] );
  strcat(nombre,"_");
  strcat(nombre, argv[2]);
  strcat(nombre,".dat");
  
  printf("%s \n",nombre);

  archivo = fopen(nombre, "w");

  for (j=0;j<num_puntos ;j++){

    fprintf( archivo, "%f \t %f \t %f \n", T[j] , X[j], Y[j]);
  }
  fclose(archivo);

  return 0;

}

void runfge_kutta4orden(double x_old, double y_old, double A, double B, double C, double D , double h , double *x_new, double *y_new){

  double  k1_der1x,  k1_der1y, x1, y1, x2, y2, x3, y3, prom_x, k2_x, k3_x, k4_x, prom_y, k2_y, k3_y, k4_y;


  k1_der1x= x_der1(x_old, y_old, A, B);
  k1_der1y= y_der1(x_old, y_old, C, D);

  //printf("%f %f k1_x \n", k1_der1x, k1_der1y);  

  // Primer paso

  x1= x_old + (h/2.0)* k1_der1x; 
  y1= y_old+ (h/2.0) *k1_der1y;

  k2_x =x_der1(x1, y1, A, B);
  k2_y= y_der1(x1, y1, C, D);

  // SEGUNDO PASO 

  x2= x_old + (h/2.0) * k2_x; 
  y2= y_old + (h/2.0) * k2_y;

  k3_x = x_der1(x2, y2, A, B);
  k3_y = y_der1(x2, y2, C, D);

  // Tercer paso

  x3= x_old + (h/2.0)* k3_x;
  y3= y_old + (h/2.0) * k3_y;

  k4_x = x_der1(x3, y3, A, B);
  k4_y = y_der1(x3, y3, C, D);

  // Cuarto paso

  prom_x=(1.0/6.0)*(k1_der1x + 2.0*k2_x + 2.0*k3_x + k4_x);
  prom_y=(1.0/6.0)*(k1_der1y + 2.0*k2_y + 2.0*k3_y + k4_y);

  *x_new=x_old + h* prom_x;
  *y_new=y_old + h* prom_y;
  
  //printf("%f %f aqui \n", *x_new , *y_new);
 
}

double x_der1(double x, double y, double A, double B){
  //dx/dt pero no depende de t por eso no lo coloque 
  return A*x*1.0 - B*x*y*1.0;
}

double y_der1(double x ,double y, double C, double D){
  //dy/dt pero no depende de t por eso no lo coloque 
  return -C*y + D*x*y*1.0;
}




