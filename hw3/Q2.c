#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>


#define D 0.75 // diameter of disc
#define W 1 // distance between parallel lines
#define SAMPLES 1000 // number of discs to drop

#define GRAPH // define GRAPH to display graph of data

int observations[SAMPLES];
double probabilities[SAMPLES];

int x[SAMPLES];

/**
 * Generates a random floating point number 
 * between a and b.
 * Source: https://stackoverflow.com/questions/5289613/generate-random-float-between-two-floats
 */
float RandomFloat(float a, float b) {
    float random = ((float) rand()) / (float) RAND_MAX;
    float diff = b - a;
    float r = random * diff;
    return a + r;
}

/**
 * Simulate a disk toss
 */
int toss_disk(void) {
    double distance = RandomFloat(0, (float)W/2);
    if (distance < (D/2))
        return 1;
    else
        return 0;
}

int main(void) {
    srand(time(NULL));
    probabilities[0] = 0;
    int sum = 0;
    for (int i = 1; i < SAMPLES; i++) {
        if (i % 10000000 == 0) printf("%d\n", i);
        x[i] = i;
        observations[i] =  toss_disk();
        sum += observations[i];
        probabilities[i] = (float)sum / i;
    }

    #ifdef GRAPH
    if (SAMPLES > 10000) exit(1);
    char * commands[] = {"set title \"Disc Toss Probabilities\"", "plot 'data.temp' with points"};
    FILE * temp = fopen("data.temp", "w");
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    
    for (int i=0; i < SAMPLES; i++) {
        fprintf(temp, "%lf %lf \n", (float)x[i], probabilities[i]); //Write the data to a temporary file
    }

    for (int i=0; i < 2; i++) {
        fprintf(gnuplotPipe, "%s \n", commands[i]); //Send commands to gnuplot one by one.
    }
    #endif
}
