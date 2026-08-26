#include <stdio.h>
#include <stdlib.h>

#define C13_CONV  0.251449530
#define N15_CONV  0.101329118

// specfrq_v2.c

// calculates the frequency for zero ppm for 1H, 13C and 15N
// also the chemical shifts for 15N and 13C at the carrier
// position
// Patrik Lundstrom 011126
//
// Non-interactve version

int main(int argc, char *argv[]) {

	double sfrq, dfrq, dfrq2, h1ppm, c13ppm, n15ppm;
	double sfrq0, dfrq0, dfrq20;

	if (argc != 5) {
		printf("usage: %s <sfrq> <h2o_ppm> <dfrq> <dfrq2>\n", argv[0]);
		return 1;
	}
	sfrq = atof(argv[1]);
	h1ppm = atof(argv[2]);
	dfrq = atof(argv[3]);
	dfrq2 = atof(argv[4]);	

	sfrq0  = sfrq / (1.0 + h1ppm*1e-6);
	dfrq0 = sfrq0*C13_CONV;
	dfrq20 = sfrq0*N15_CONV;

	c13ppm = (dfrq-dfrq0)/dfrq0*1e6;
	n15ppm = (dfrq2-dfrq20)/dfrq20*1e6;

	printf("H1  sfrq(0):  %.7f\t shift: %f\n", sfrq0, h1ppm);
	printf("C13 dfrq(0):  %.7f\t shift: %f\n", dfrq0, c13ppm);
	printf("N15 dfrq2(0): %.7f\t shift: %f\n", dfrq20, n15ppm);
return 0;
}
