#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>


#define D 0.75 // diameter of disc
#define W 1 // distance between parallel lines
#define SAMPLES 1984444444 // number of discs to toss

#if SAMPLES <= 100000
#define GRAPH // define GRAPH to display graph of data
/** NOTE - due to the use of a file buffer for storing the data for graphing,
 *  large sample sizes will not be graphed to prevent filling up of disc space
 */
#endif


// Seeds for pseudo-RNG 
static double seed0;
static double seed1;

/**
 * Initialize seeds for pseudo-RNG
 */
void init_seeds(void) {
    seed0 = 0.15924135;
    seed1 = (double)time(NULL);
    while (seed1 >=1 ) seed1 /= 10.0; // shift right a decimal place until > 1
}

/**
 * Generate a random floating point number between a and b
 * Algorithm sourced from Prof. Deng's lecture 2 notes
 */
double random_float(double a, double b) {
    double temp = seed0;
    seed0 = seed1;
    if (seed1 + temp > 1.0)
        seed1 = fmod((seed1 + temp) - 1.0, 1.0);
    else
        seed1 = fmod((seed1 + temp), 1.0);
    return a + (seed1 * (b-a));
}

/**
 * Simulate a disk toss
 */
int toss_disk(void) {
    double distance = random_float(0, (double)W/2); // Generate random floating point number between 0 and W/2
    if (distance < (D/2))
        return 1;
    else
        return 0;
}


int main(void) {
    init_seeds();
    time_t starting_time = time(NULL);
    #ifdef GRAPH // only save probabilities if they will be graphed at end
    double* probabilities = malloc(sizeof(double) * SAMPLES+1);
    #endif
    int sum = 0;
    for (int i = 0; i < SAMPLES; i++) {
        sum += toss_disk();
        #ifdef GRAPH
        probabilities[i] = (double)sum / (i+1);
        #endif
    }
    printf("Final estimate probability: %f (%lu second(s))\n", (double)sum / SAMPLES, time(NULL)-starting_time);

    #ifdef GRAPH
    if (SAMPLES > 100000) exit(EXIT_SUCCESS);
    char * commands[] = {"set title \"Disc Toss Probabilities\"", "plot 'data.temp' with points"};
    FILE * temp = fopen("data.temp", "w");
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    
    for (int i=0; i < SAMPLES; i++) {
        fprintf(temp, "%lf %lf \n", (double)i, probabilities[i]); //Write the data to a temporary file
    }

    for (int i=0; i < 2; i++) {
        fprintf(gnuplotPipe, "%s \n", commands[i]); //Send commands to gnuplot one by one.
    }
    free(probabilities);
    #endif
}
