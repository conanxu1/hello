#include "fftw3.h"

int main()
{ 

    fftw_complex *in, *out;
    fftw_plan p;
    int N= 9;
    int i;
    int j;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	
	
        for( i=0; i < N; i++)
    {
            in[i][0] = 1.0;
            in[i][1] = 0.0;
            printf("%6.2f ",in[i][0]);
    }
	
    printf("\n");
	
    p=fftw_plan_dft_1d(N,in,out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p); /* repeat as needed*/
        
	for(j = 0;j < N;j++)
    {
            printf("%6.2f ",out[j][0]);
    }
    printf("\n");
    
	
	fftw_destroy_plan(p);
    fftw_free(in); 
    fftw_free(out);
    return 0;
}


