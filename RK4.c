#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double* xs;
double* ys;

double dydx(double x, double y);
void RK4(double x0, double y0, double xn, double h);

int main() {
    double h = 0.001;
    double xn = 10;
    double x0 = 1;
    double y0 = 0;
    RK4(x0, y0, xn, h);
    char * commands[] = {"set title \"f(x)\"", "plot 'data.temp' with line"};
    FILE * temp = fopen("data.temp", "w");
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    
    for (int i=0; i < ceil((xn/h)); i++) {
        fprintf(temp, "%lf %lf \n", xs[i], ys[i]); //Write the data to a temporary file
    }

    for (int i=0; i < 2; i++) {
        fprintf(gnuplotPipe, "%s \n", commands[i]); //Send commands to gnuplot one by one.
    }
    free(xs);
    free(ys);
}


void RK4(double x0, double y0, double xn, double h) {
    double x=x0;
    double y=y0;
    xs = (double*)calloc(ceil((xn/h)), sizeof(double));
    ys = (double*)calloc(ceil((xn/h)), sizeof(double));
    for (int i = 0; x < xn; i++) {
        double m1 = dydx(x0,y0);
        double m2 = dydx((x0+h/2.0),(y0+m1*h/2.0));
        double m3 = dydx((x0+h/2.0),(y0+m2*h/2.0));
        double m4 = dydx((x0+h),(y0+m3*h));
        double m = ((m1+2*m2+2*m3+m4)/6);
        y = y + m * h;
        x += h;

        xs[i] = x;
        ys[i] = y;
        // printf("(%f, %f)\n", xs[i], ys[i]);
    }
}

double dydx(double x, double y) {
    return (3 * x - 2 * y) / x;
}