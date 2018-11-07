bool IsPrintEnabled;
  double *u;
int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    printf("USAGE: GZ_Par1.exe N\n");
    return 1;
  }
  if (argc == 2)
    IsPrintEnabled = 0;
  else
    IsPrintEnabled = atoi(argv[2]);

  int N = atoi(argv[1]);

  GZ_Par(N);  

  return 0;
}

int GZ_Par(int N)
{
  int i, j;
  clock_t begin,end;
  
  double temp;
  double dmax;
  double dm; 
  double d;
  
  int Step;
  omp_lock_t dmax_lock;

  printf("The Gauss - Seidel algorithm. OpenMP Parallel version.\r\n");
  
  if((u = CreateMatrix(N)) == NULL)
    return 1;
  
  Step = 0;
  
  printf("Matrix size: %d\r\n", N);
  begin = clock();
  
  omp_init_lock(&dmax_lock);
  
  do
  {
    dmax = 0;
#pragma omp parallel for shared(u,N,dmax) private(i)
    for (i = 1; i < N - 1; i++)
#pragma omp parallel for shared(u,N,dmax) private(j,temp,d)
      for(j = 1; j < N - 1; j++)
      {
        temp = u[N * i + j];
        u[N * i + j] = 0.25 * (u[N * i + j + 1] + u[N * i + j - 1] + 
          u[N * (i + 1) + j] + u[N * (i - 1) + j]);
        d = fabs(u[N * i + j] - temp);
        
        omp_set_lock(&dmax_lock);
        if ( dmax < d ) 
          dmax = d;
        omp_unset_lock(&dmax_lock);
      }    
      Step++;
  }
  while (dmax > EPS);
  
  end = clock();
  
  long tck = end - begin;
  
  printf("Results:\r\n- current eps = %f\r\n- required eps = %f\r\n- matrix Size = %d\r\n- iterations = %d\r\n- time = %f\r\n",
    dmax, EPS, N, Step, ((double)tck)/CLOCKS_PER_SEC);
  
  if (IsPrintEnabled)
    PrintMatrix(u, N);

  DeleteMatrix(u);
  
  return 0;
}

