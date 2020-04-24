#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>


#define D 0.75 // diameter of disc
#define W 1 // distance between parallel lines
#define SAMPLES 1984444444 // number of discs to drop

#if SAMPLES < 100000
#define GRAPH // define GRAPH to display graph of data
/** NOTE - due to the use of a file buffer for storing the data for graphing,
 *  large sample sizes will not be graphed to prevent filling up of disc space
 */
#endif

/**
 * Generates a random floating point number 
 * between a and b.
 * Source: https://stackoverflow.com/questions/5289613/generate-random-float-between-two-floats
 */
float RandomFloat(float a, float b) {
    float random = ((float) rand()) / (float) RAND_MAX; // generate random number between 0 and 1
    float diff = b - a;
    float r = random * diff;
    return a + r;
}

/**
 * Simulate a disk toss
 */
int toss_disk(void) {
    float distance = RandomFloat(0, (float)W/2); // Generate random floating point number between 0 and W/2
    if (distance < (D/2))
        return 1;
    else
        return 0;
}


int main(void) {
    time_t starting_time = time(NULL);
    #ifdef GRAPH // only save probabilities if they will be graphed at end
    float* probabilities = malloc(sizeof(float) * SAMPLES+1);
    #endif
    int sum = 0;
    for (int i = 0; i < SAMPLES; i++) {
        sum += toss_disk();
        #ifdef GRAPH
        probabilities[i] = (float)sum / (i+1);
        #endif
    }
    printf("Final estimate probability: %f (%lu second(s)\n", (float)sum / SAMPLES, time(NULL)-starting_time);

    #ifdef GRAPH
    if (SAMPLES > 100000) exit(EXIT_SUCCESS);
    char * commands[] = {"set title \"Disc Toss Probabilities\"", "plot 'data.temp' with points"};
    FILE * temp = fopen("data.temp", "w");
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    
    for (int i=0; i < SAMPLES; i++) {
        fprintf(temp, "%lf %lf \n", (float)i, probabilities[i]); //Write the data to a temporary file
    }

    for (int i=0; i < 2; i++) {
        fprintf(gnuplotPipe, "%s \n", commands[i]); //Send commands to gnuplot one by one.
    }
    free(probabilities);
    #endif
}
