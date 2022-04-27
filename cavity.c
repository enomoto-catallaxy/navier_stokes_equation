#include<stdio.h>
#include <eigen3/Eigen/Core>

int nx = 64;
int ny = 64;

double u[nx][ny];
double v[nx][ny];
double p[nx][ny];

double us[nx][ny];
double vs[nx][ny];
int L = 1; //length

int Re = 100;
int U = 100;
double dx = L / nx, dy = L / ny;
double dt = 0.01;
double a = dt / dx;
double b = dt / dx / dx ;

int main()
{
  int i, j;
  for (i = 0; i<= nx -1; i++ ){ // set Arrays
    for(j = 0; j <= ny; j++ ){
    u[i][j] = 0;
    v[i][j] = 0;
    p[i][j] = 0;
    us[i][j] = 0;
    vs[i][j] = 0;
    }
  }

  for(j = 1; j < ny; j++){
    for ( i = 0; i < nx; i++)
    {
      v_here = 0.25 * (v[i-1][j] + v[i-1][j+1] + v[i][j] + v[i][j+1]);
      us[i][j] = u[i][j] + b * ((u[i-1][j] - 2 * u[i][j] + u[i+1][j]) + (u[i][j-1] - 2 * u[i][j] + u[i][j+1])) - u[i][j] * (u[i+1][j] - u[i-1][j]) *0.5 /dx - v_here * (u[i][j+1] - u[i][j-1]) * 0.5 /dy;
    }
  }

  
}
