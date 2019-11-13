#include <stdio.h>
#include <stdlib.h>
#include "vptree.h"
#include  <time.h>

struct timespec start, finish;
double elapsed;



/* ... */




void printTheTree(vptree * tree,int dims);

int main()
{
    // printf("Hello world!\n");
    // return 0;
    //=NULL;
    srand(time(NULL));
    int n=50000;
    int d=400;
    int i ,j ;
    clock_t t;
    double time_taken;

    //wanna create an n*d array dynamically
    // double **ptr=(double **)malloc(n * sizeof(double *));
    double  * ptr = (double * ) malloc( n*d * sizeof(double) );
    for (int i=0;i<n*d;i++)
      ptr[i]=rand()%10;

    // for (i = 0; i < n; i++) {
    //   ptr[i]=(double *)malloc(d * sizeof(double));
    // }
    // for ( i = 0; i < n; i++) {
    //   for ( j = 0; j < d; j++) {
    //     ptr[i][j]=(double)(rand()%10);//(double)rand()/((double)RAND_MAX/1) ;       (double)(rand()%10);
    //
    //   }
    // }

    // printf("Given dataset:\n");
    // for ( i = 0; i < n; i++) {
    //   printf("[ ");
    //   for ( j = 0; j < d; j++) {
    //     printf("%lf ",ptr[i*d+j]);
    //   }
    //   printf("]\n");
    // }

    //casting the pointer to fit the arguments
    double *pointer = (double *) ptr;
    vptree * VPtree= buildvp(pointer, n,d);


//
// printTheTree(VPtree,d);

//GET the running time of the treebuild
clock_gettime(CLOCK_MONOTONIC, &start);
 buildvp(pointer, n,d);
clock_gettime(CLOCK_MONOTONIC, &finish);
elapsed = (finish.tv_sec - start.tv_sec);
elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

printf("The build of the tree with a dataset %dx%d  took %lf seconds! \n",n,d, elapsed );
}


void printTheTree(vptree * tree,int dims){
  if(tree == NULL){
    return;
  }
  printf("Tree data:\n");
  printf("[");
  printf("Median Distance: %lf ",getMD(tree));
  printf("\n Vantage Point: ");
  for (int i = 0; i < dims; i++) {
    printf("%lf ",(getVP(tree)[i]));
  }
  printf("\n IDX of the Vantage Point: %d ",getIDX(tree));


  printf("]\n \n");

  printf("Inner");
  printTheTree(getInner(tree),dims);

  printf("Outer");
  printTheTree(getOuter(tree),dims);

}
