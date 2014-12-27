#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <nmmintrin.h>

// http://stackoverflow.com/questions/16619953/parallelizing-the-gaussian-blur-algorithm-with-openmp

float convolve(const float *kernel, const float *ringbuf, const int ksize, const int bufi0) {
    float sum = 0.0f;
    int i;
    for(i=0; i<ksize; i++) {
        sum += kernel[i]*ringbuf[i];
    }
    return sum;
}


void gaussian_blur(float *src, float *dst, int w, int h, float sigma, int nthreads)
{
    int x, y, i;
    int ksize = (int)(sigma * 2.f * 4.f + 1) | 1;
    int halfk = ksize / 2;
    printf("ksize %d\n", ksize);
    float scale = -0.5f/(sigma*sigma);
    float sum = 0.f;
    float *kernel, *ringbuf;
    int xmax = w - halfk;
    int ymax = h - halfk;

    // if sigma too small, just copy src to dst
    if (ksize <= 1)
    {
        for (y = 0; y < h; y++)
            for (x = 0; x < w; x++)
                dst[y*w + x] = src[y*w + x];
        return;
    }

    // create Gaussian kernel
    //kernel = malloc(ksize * sizeof(float));
    kernel =  (float*)_mm_malloc(ksize * sizeof(float),16);
    //ringbuf = malloc(ksize * sizeof(float));
    ringbuf = (float*)_mm_malloc(nthreads*ksize * sizeof(float),16);

    #pragma omp parallel for reduction(+ : sum) if(nthreads>1)
    for (i = 0; i < ksize; i++)
    {
        float x = (float)(i - halfk);
        float t = expf(scale * x * x);
        kernel[i] = t;
        sum += t;
    }

    scale = 1.f / sum;
    #pragma omp parallel for if(nthreads>1)
    for (i = 0; i < ksize; i++)
        kernel[i] *= scale;

    // blur each row
    #pragma omp parallel for if(nthreads>1)
    for (y = 0; y < h; y++)
    {
        int ithread = omp_get_thread_num();
        int x1;
        int bufi0 = ksize-1;
        float tmp = src[y*w + 0];
        for (x1 = 0; x1 < halfk  ; x1++) ringbuf[ksize*ithread + x1] = tmp;
        for (; x1 < ksize-1; x1++) ringbuf[ksize*ithread + x1] = src[y*w + x1-halfk];
        for (x1 = 0; x1 < w; x1++)
        {
            const int ibufi0_fix = (x1+ksize-1)%ksize;

            if(x1 < xmax)
                ringbuf[ksize*ithread + ibufi0_fix] = src[y*w + x1+halfk];
            else
                ringbuf[ksize*ithread + ibufi0_fix] = src[y*w + w-1];
            if (bufi0 == ksize) bufi0 = 0;
            dst[y*w + x1] = convolve(kernel, &ringbuf[ksize*ithread], ksize, bufi0);
        }
    }
    // blur each column
    #pragma omp parallel for
    for (x = 0; x < w; x++)
    {
        int ithread = omp_get_thread_num();
        int y1;
        int bufi0 = ksize-1;
        float tmp = dst[0*w + x];
        for (y1 = 0; y1 < halfk  ; y1++) ringbuf[ksize*ithread + y1] = tmp;
        for (     ; y1 < ksize-1; y1++) ringbuf[ksize*ithread + y1] = dst[(y1-halfk)*w + x];

        for (y1 = 0; y1 < h; y1++)
        {
            const int ibufi0_fix = (y1+ksize-1)%ksize;
            if(y1 < ymax)
                ringbuf[ibufi0_fix] = dst[(y1+halfk)*w + x];
            else
                ringbuf[ibufi0_fix] = dst[(h-1)*w + x];

            if (bufi0 == ksize) bufi0 = 0;
            dst[y1*w + x] = convolve(kernel, &ringbuf[ksize*ithread], ksize, bufi0);
        }
    }

    // clean up
    _mm_free(kernel);
    _mm_free(ringbuf);
}

int compare(float *dst1, float *dst2, const int n) {
    int error = 0;
    int i;
    for(i=0; i<n; i++) {
        if(*dst1 != *dst2) error++;
    }
    return error;
}


int main(int argc, char *argv[]) {
    const int w = 20;
    const int h = 20;
    int i;

    float *src =  (float*)_mm_malloc(w*h*sizeof(float),16);
    float *dst1 =  (float*)_mm_malloc(w*h*sizeof(float),16);
    float *dst2 =  (float*)_mm_malloc(w*h*sizeof(float),16);
    for(i=0; i<w*h; i++) {
        src[i] = i;
    }

    gaussian_blur(src, dst1, w, h, 1.0f, 1);
    gaussian_blur(src, dst2, w, h, 1.0f, 4);
    int error = compare(dst1, dst2, w*h);
    printf("error %d\n", error);
    _mm_free(src);
    _mm_free(dst1);
    _mm_free(dst2);

    return 0;
}
