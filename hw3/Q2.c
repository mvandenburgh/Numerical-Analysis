#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>


#define D 0.75 // diameter of disc
#define W 1 // distance between parallel lines
#define SAMPLES 100000000 // number of discs to drop

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


double avg(int upto) {
    int sum = 0;
    for (int i = 0; i < upto; i++) {
        sum += observations[i];
    }
    return (float)sum / upto;
}

void printObs(void) {
    for (int i = 0; i < SAMPLES; i++) {
        printf("%d ", observations[i]);
    }
}

void printProbs(void) {
    for (int i = 0; i < SAMPLES; i++) {
        printf("%f ", probabilities[i]);
        if (i % 20 == 0) printf("\n");
    }
}


//rand() % (max_number + 1 - minimum_number) + minimum_number

int main(void) {
    srand(time(NULL));
    probabilities[0] = 0;
    for (int i = 0; i < SAMPLES; i++) {
        x[i] = i;
        observations[i] =  toss_disk();
        probabilities[i] = avg(i);
    }
    // printProbs();
    // printObs();

    char * commands[] = {"set title \"Disc Toss Probabilities\"", "plot 'data.temp' with points"};
    FILE * temp = fopen("data.temp", "w");
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    
    for (int i=0; i < SAMPLES; i++) {
        fprintf(temp, "%lf %lf \n", (float)x[i], probabilities[i]); //Write the data to a temporary file
    }

    for (int i=0; i < 2; i++) {
        fprintf(gnuplotPipe, "%s \n", commands[i]); //Send commands to gnuplot one by one.
    }
}
