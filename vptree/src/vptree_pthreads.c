#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <stdbool.h>
#include "vptree.h"
#include "limits.h"
/*FOR THREADING WE USE THE PEER METHOD ,
MAIN THREAD CALLS OTHER THREADS(WORKER THREADS)
BUT IT DOES ALSO PART OF THE WORK */
// Number of worker threads that will help for the computation of the distance.
#define N_OF_DISTANCE_WORKER_THREADS 2
//Number of all Active threads in the programm ( main thread included)
#define MAX_N_OF_ACTIVE_THREADS 4 // Works perfectly for my 8 year old laptop with i5-2410m :')
//Defining a threashold for the distance parallelism
//(bellow that number of spots no parallelism happen in distance calculation)
#define DISTANCE_THRESHOLD 2000
//Defining a threashold for the subtrees build parallelism
//(bellow that number of spots no parallelism happen in subtrees building)
#define PARTS_THRESHOLD 100


//Struct used for the computation of distances in parallel
typedef struct distanceCalcData{
  int numberOfPoints;
  int dims;
  int offset;
  int start_of_computation;
  int end_of_computation;
  double * Xdata;
  double * ed;
}distanceCalcData;

//Struct used for the subtrees building in parallel
typedef struct partOfTree_Data{
  int numberOfPoints;
  int dims;
  int partSize;
  double * parentX;
  vptree *partParent;
  double * ed;
  int * parentIDX;
}partOfTree_Data;

/*Declaring the number of the actives threads as a global variable, in order
to be accesed by all threads , it has a start value of 1 due to the main thread*/
int numberOfActiveThreads=1;
pthread_mutex_t lockforNumberOfThreads;

///////////////////////////
//////////////////////////
//Declaration of functions I created, they are defined lower in this file//
/*The Recursive function that helps build the vantage point tree*/
vptree * vpt_pthreads(double *Xarray, int n, int d, int *idxArray);

/*! Function for computing the Median of the Euclidean Distances of the points
 to the vantage point (ed array), calling quickselect algorithm WITHOUT
 changing the original array's order of items by copying the original array .*/
double computeMedian(double * arr,int n);

/*! Sets the euclidean distances of the points to the vantage array the proper values*/
void computeEuclideanDistances(double *arr ,int numberOfSpots,int dimensions, vptree * tree,double * Xarr);

/*The parallel with pthreads version of computeEuclideanDistances() function*/
void computeEuclideanDistances_pthreads(double *arr ,int numberOfSpots,int dimensions, vptree * tree,double * Xarr);

/*!Sets the inner and the outer subtree of the current tree.*/
void innerAndOuterTree(double *euclideanDistances ,int numberOfSpots,int dimensions, vptree * tree,double * Xarr,int * idxArray);
/*The parallel with pthreads version of  innerAndOuterTree() function*/
void innerAndOuterTree_pthreads(double *euclideanDistances ,int numberOfSpots,int dimensions, vptree * tree,double * Xarr,int * idxArray);

//////////////////////QUICKSELECT//////////////////////////////
//Swap for using in the quickselect algorithm.
void swap(double * xp, double * yp);
// Standard partition process of QuickSort().
// It considers the last element as pivot
// and moves all smaller element to left of it and greater elements to right
int partition(double arr[], int l, int r);
// This function returns k'th smallest
// element in arr[l..r] using QuickSort
// based method.  ASSUMPTION: ALL ELEMENTS
// IN ARR[] ARE DISTINCT
double kthSmallest(double array[], int l, int r, int k);
///////////////////////////////////////////////////////////////
/*Declaration ENDS*/
///////////////////////////
//////////////////////////
//Definitions
vptree * buildvp(double *X, int n, int d){


  //IdxArr is the array that help as get the index of the vantage point of the vptree( IDX)
  //on the input data array vector
  int * IdxArr=(int * )malloc( n*sizeof(int) );
  vptree * root ;
  for (int i = 0; i < n; i++) {
    IdxArr[i]=i;
    }
  root= vpt_pthreads(X, n, d, IdxArr );

  free(IdxArr);
  return root;


}
//The pthreads version of the recursive function that build every node of the tree
vptree * vpt_pthreads(double *Xarray, int n, int d, int *idxArray){
  //If there is no points , (node)-->NULL
  if(n==0){
    return NULL;
  }

  //declaring the vantage point tree
  vptree * T = (vptree *)malloc(sizeof(vptree));

  //Vector ed[n] stores the euclidean distances of each point to the vantage point.
  //ed[i] is the Euclidean distance of X[i-st] point from the vantage point
  double ed[n-1] ;

  /////////////////////////////////////////////
  //0. SETTING THE IDX OF THE ROOT AS THE LAST POINT OF THE IDXARRAY OF THE TREE.
  //We set the vantage point to be the last in order point of a given data points.
  T->idx=idxArray[n-1];

  /////////////////////////////////////////////
  //1.INITIALIZING VANTAGE POINT (VP) OF THE TREE.
  //We set the vantage point to be the last in order point of a given data points.
  T->vp=(double *)malloc(d * sizeof(double));
  for (int i = 0; i < d; i++) {
        T->vp[i]=Xarray[(n-1)*d + i];
  }

  /////////////////////////////////////////////
  //2.
  //Compute and setting the euclidean distances of the points in the ed vector.
  /*If we are above distance threshold and there are active threads(with respect to the max we set)
   left we compute distances  in parallel , otherwise the sequential implementation runs*/
  if((n>=DISTANCE_THRESHOLD)&&(numberOfActiveThreads<=(MAX_N_OF_ACTIVE_THREADS-N_OF_DISTANCE_WORKER_THREADS)) )
    computeEuclideanDistances_pthreads(ed,n,d,T,Xarray);
  else
    computeEuclideanDistances(ed,n,d,T,Xarray);

  //Getting the Median of the Euclidean Distances of the points
  //to the vantage point (median of ed array)
  T->md=computeMedian(ed,n);

  /////////////////////////////////////////////
  //3.GETTING THE INNER AND THE OUTER TREES OF THE VPTREE.
  /*If we are above parts threshold and there are active threads(with respect to the max we set)
   left ,we build subtrees in parallel , otherwise we follow the sequential implementation*/
  if((n>=PARTS_THRESHOLD)&&(numberOfActiveThreads<MAX_N_OF_ACTIVE_THREADS))
    innerAndOuterTree_pthreads(ed,n,d,T,Xarray,idxArray);
  else
    innerAndOuterTree(ed,n,d,T,Xarray,idxArray);
  /////////////////////////////////////////////
  return T;

}


vptree * getInner(vptree * T){
  return T->inner;
}
vptree * getOuter(vptree * T){
  return T->outer;
}
double getMD(vptree * T){
  return T->md;
}
double * getVP(vptree * T){
  return T->vp;

}
int getIDX(vptree * T){
  return T->idx;
}

double computeMedian(double * euclideanDistances,int numberOfSpots){

  //Coppying the euclideanDistances array into a new array in order not to alter
  //the original  euclideanDistances(ed) array.
  double array[numberOfSpots-1];
  for (size_t i = 0; i < numberOfSpots-1; i++) {
      array[i]=euclideanDistances[i];
  }

  //If there is only one point in the given points set,
  //of course there is no median distance.is no median distance.
  if(numberOfSpots==1){
    return 0;
  }

  /*  //If n is an even then n-1(euclideanDistances[] lentgh and number of points
        left (without vp)) is an odd number,
    //Otherwise n-1 (euclideanDistances[] lentgh and number of points
        left (without vp) ) is an even number.
  */
  //If n-1 is odd median=euclideanDistances[n/2].
  else if(numberOfSpots%2==0){
    return kthSmallest(array, 0, numberOfSpots-2, numberOfSpots/2);

  }
  //else if n-1 is even, median=((euclideanDistances[n-1]/2) + (euclideanDistances[n+1]/2))
  else{
    return (kthSmallest(array, 0, numberOfSpots-2, (numberOfSpots-1)  /2)+
    kthSmallest(array, 0, numberOfSpots-2, (numberOfSpots+1)/2) )  /2  ;
  }
}

//Declaration of the function,worker threads use to compute parts of the euclidean distances
void * computePartOfDistance(void * arg);

//Definitions of the serial computeEuclideanDistances() function
void computeEuclideanDistances(double *euclideanDistances ,int numberOfSpots,
  int dimensions, vptree * tree,double * Xarr){
  for (int i = 0; i < (numberOfSpots-1); i++) {
    double temp=0;
    for (int j = 0; j < dimensions; j++) {
      temp+=pow(  (Xarr[i*dimensions+j]-(getVP(tree))[j]), 2 );
    }
    euclideanDistances[i]=sqrt(temp);
  }
}
/*In our parallel distance calcularion implementations we break the euclidean
distances array in N parts ,N depends on the numbber of threads we set to compute
the euclidean distances,and each thread (included the caller thread-peer method )
compute and sets its part of the euclidean distances  */

//Definitions of the pthreads parallel computeEuclideanDistances() function
void computeEuclideanDistances_pthreads(double *euclideanDistances ,int numberOfSpots,
  int dimensions, vptree * tree,double *Xarr){
  int offset;
  //Allocating memory to store the data that will be used by the distance threads
  distanceCalcData * data= (distanceCalcData *)malloc(N_OF_DISTANCE_WORKER_THREADS*sizeof(distanceCalcData));

  pthread_t DistanceThreads[N_OF_DISTANCE_WORKER_THREADS];
  int rc[N_OF_DISTANCE_WORKER_THREADS];
  offset=(numberOfSpots-1)/(N_OF_DISTANCE_WORKER_THREADS);
  if(offset==0){
    offset=numberOfSpots-1;
  }

  //set the properly data for a thread to calcute parts of the euclidean distances.
  for (int i = 0; i < N_OF_DISTANCE_WORKER_THREADS; i++) {
    (data+i)->Xdata=Xarr;
    (data+i)->numberOfPoints=numberOfSpots;
    (data+i)->dims=dimensions;
    (data+i)->ed=euclideanDistances;
    (data+i)->start_of_computation=(1+i)*offset;
    (data+i)->end_of_computation=(2+i)*offset;
    //calling each thread to calculate its part of the distacnes
    rc[i]=pthread_create(DistanceThreads+i,NULL, computePartOfDistance ,(void*)(data+i));
    if(rc[i]){
      printf("ERROR; return code from pthread_create() is %d\n", rc[i]);
      exit(-1);
    }
    //locking the numberOfActiveThreads, for avoiding the race condition
    pthread_mutex_lock(&lockforNumberOfThreads);
    numberOfActiveThreads++;
    pthread_mutex_unlock(&lockforNumberOfThreads);
  }

  //main thread calcutes its part of distances.
  for (int i = 0; i < offset; i++) {
    double temp=0;
    for (int j = 0; j < dimensions; j++) {
      temp+=pow(  (Xarr[i*dimensions+j]-(getVP(tree))[j]), 2 );
    }
    euclideanDistances[i]=sqrt(temp);
  }

  //waiting for every thread to finish its work
  for (int i = 0; i < N_OF_DISTANCE_WORKER_THREADS; i++) {
      rc[i]=pthread_join(DistanceThreads[i],NULL);
      if (rc[i]) {
        printf("ERROR; return code from pthread_join() is %d\n", rc[i]);
        exit(-1);
      }
    }
    //free allocated memory from the heap
    free(data);
  }

//Definition of the function, threads use to compute parts of the euclidean distances
void * computePartOfDistance(void * arg){
  //casting the argument in order to be used by the thread function
  distanceCalcData * data =(distanceCalcData *)arg;


  //getting the values we need for the calculation
  int dimensions=data->dims;
  int start = (data->start_of_computation);
  int end;
  if( (data->end_of_computation)<=(data->numberOfPoints-1) )
    end=(data->end_of_computation);
  else
    end=(data->numberOfPoints-1);
  if( (data->start_of_computation)>(data->numberOfPoints-1) ){
    pthread_mutex_lock(&lockforNumberOfThreads);
    --numberOfActiveThreads;
    pthread_mutex_unlock(&lockforNumberOfThreads);
    pthread_exit(NULL);
  }
    // if((data->numberOfPoints)-1<start){
  //   //locking the numberOfActiveThreads, for avoiding the race condition
  //   pthread_mutex_lock(&lockforNumberOfThreads);
  //   --numberOfActiveThreads;
  //   pthread_mutex_unlock(&lockforNumberOfThreads);
  //   pthread_exit(NULL);
  // }

  //main thread calcutes its part of distances.
  for (int k = start; k < end ; k++) {
    double temp=0;
    for (int j = 0; j < dimensions; j++) {
      temp+=pow(  ((data->Xdata[k*(dimensions)+j])-data->Xdata[((data->numberOfPoints)-1)*dimensions+j]), 2 );
    }
    data->ed[k]=sqrt(temp);
  }
  //locking the numberOfActiveThreads, for avoiding the race condition
  pthread_mutex_lock(&lockforNumberOfThreads);
  --numberOfActiveThreads;
  pthread_mutex_unlock(&lockforNumberOfThreads);
  pthread_exit(NULL);
}

//Definition of the sequential innerAndOuterTree() function
void innerAndOuterTree(double *euclideanDistances ,int numberOfSpots,int dimensions, vptree * tree,double * Xarr, int * idxArray ){
  //innerX and outerX are used as dataset inputs for the inner and outer sub(vp)tree
  double *innerX,*outerX;

  //idxInner and idxOuter are used as IdxArr for the inner and outer subtrees
  //innerSize and outerSize are used for determining the inner and outer subtrees size
  //innerCount and outerCount are used for scaning given iputs and filling inner and outer trees
  int *idxInner=NULL, *idxOuter=NULL, innerSize=0, innerCount=0, outerSize=0, outerCount=0;


    //Getting the size of the inner and the outer subtrees.
    for(int i=0;i<(numberOfSpots-1);i++){
      if(euclideanDistances[i]<=getMD(tree))
        innerSize++;
    }
    outerSize=(numberOfSpots-1)-innerSize;

    //Allocating memory for inner and outer subtrees datasets
    innerX=(double *)malloc( innerSize * dimensions * sizeof(double));
    outerX=(double *)malloc( outerSize * dimensions * sizeof(double));


    //Allocating memory for innerIdx and outerIdx
    idxInner=(int *)malloc( (innerSize)* sizeof(int));
    idxOuter=(int *)malloc( (outerSize)* sizeof(int));

  //Getting the idxInner and idxOuter and filling inner and outer subtree
  //by scanning the euclidean distances of n-1 remaining spots
  for(int i=0;i<(numberOfSpots-1);i++){

    if(euclideanDistances[i]<=getMD(tree)){
        for(int j=0;j<dimensions;j++){
            innerX[innerCount*dimensions+j]=Xarr[i*dimensions+j];
        }
          idxInner[innerCount++]=idxArray[i];
    }

    else if(euclideanDistances[i]>getMD(tree)){
      for(int j=0;j<dimensions;j++){
        outerX[outerCount*dimensions+j]=Xarr[i*dimensions+j];
      }
        idxOuter[outerCount++]=idxArray[i];
    }
  }
  //Creating the current tree's , inner subtree
  tree->inner=vpt_pthreads(innerX,innerCount,dimensions,idxInner);
  //free allocated memory from the heap
  free(innerX);
  free(idxInner);
  //Creating the current tree's , outer subtree
  tree->outer=vpt_pthreads(outerX,outerCount,dimensions,idxOuter);
  //free allocated memory from the heap
  free(outerX);
  free(idxOuter);


}
/*This function is used to build the outer subtree of a vptree but it can easily
 be configured,with minor changes in its code and in
 innerAndOuterTree_pthreads() code to build the inner part */
//Definition of the function, threads use to build a subtree of the vptree
void * createPartOfTree( void * arg){
  //locking the numberOfActiveThreads, for avoiding the race condition
  pthread_mutex_lock(&lockforNumberOfThreads);
  ++numberOfActiveThreads;
  pthread_mutex_unlock(&lockforNumberOfThreads);
  //casting the argument in order to be used by the thread function
  partOfTree_Data *part =(partOfTree_Data *)arg;

  int count=0;
  int n=part->numberOfPoints;
  int d=part->dims;
  double * parentX= part->parentX;
  int size=part->partSize;
  double * ed=part->ed;
  int * parentIDX=part->parentIDX;

  //Allocating memory for outerTree
  int *idx=(int *)malloc( (size)* sizeof(int));
  double * partX=(double *)malloc( size * d * sizeof(double));

  /* ///this works if i want to build the inner subtree with this function
  // if((part->isInner)==true){
  //   for(int i=0;i<(n-1);i++){
  //     if(ed[i]<=getMD(part->partParent)){
  //     for(int j=0;j<d;j++){
  //       partX[count*d+j]=parentX[i*d+j];
  //     }
  //     idx[count++]=parentIDX[i];
  //     }
  //   }
  // }
  // else if((part->isInner)!=true) {*/
    //scanning parentX and getting outer part X and idx array of the outer subtree
    for(int i=0;i<(n-1);i++){
      if(ed[i]>getMD(part->partParent)){
      for(int j=0;j<d;j++){
        partX[count*d+j]=parentX[i*d+j];
      }
      idx[count++]=parentIDX[i];
      }
    }

    //Creating the current tree's , outer subtree
    (part->partParent)->outer=vpt_pthreads(partX,size,d,idx);

    free(idx);
    free(partX);
    //locking the numberOfActiveThreads, for avoiding the race condition
    pthread_mutex_lock(&lockforNumberOfThreads);
    --numberOfActiveThreads;
    pthread_mutex_unlock(&lockforNumberOfThreads);

    pthread_exit(NULL);
}


//Definition of the parallel with pthreads innerAndOuterTree() function
void innerAndOuterTree_pthreads(double *euclideanDistances ,int numberOfSpots,int dimensions, vptree * tree,double * Xarr, int * idxArray ){

  pthread_t outerSubtreeThread;
  partOfTree_Data * part= (partOfTree_Data *)malloc(sizeof(partOfTree_Data));

  //idxInner is used as IdxArr for the inner subtree
  //innerSize is used for determining the inner subtree size
  //innerCount is used for scaning given iputs and filling inner tree
  int innerSize=0,outerSize=0;
    //Getting the size of the inner and the outer subtrees.
    for(int i=0;i<(numberOfSpots-1);i++){
      if(euclideanDistances[i]<=getMD(tree))
        innerSize++;
    }
    outerSize=(numberOfSpots-1)-innerSize;

    //Data of the part(subtree) the thread needs to implement it.
    (part)->numberOfPoints=numberOfSpots;
    (part)->dims=dimensions;
    (part)->partSize=outerSize;
    (part)->parentX=Xarr;
    (part)->ed=euclideanDistances;
    (part)->partParent=tree;
    (part)->parentIDX=idxArray;
    //calling the thread to build the outer tree
    int rc =pthread_create(&outerSubtreeThread,NULL, createPartOfTree , (void *)(part) );
    if(rc){
      printf("ERROR; return code from pthread_create() is %d\n", rc);
      exit(-1);
    }

    //Allocating memory for inner subtree dataset
    double * innerX=(double *)malloc( innerSize * dimensions * sizeof(double));
    //Allocating memory for innerIdx
    int * idxInner=(int *)malloc( (innerSize)* sizeof(int));
    int innerCount=0;


    //Gettin the IdxInner and filling inner subtree data (InnerX)
    //by scanning the euclidean distances of n-1 remaining spots
    for(int i=0;i<(numberOfSpots-1);i++){
      if(euclideanDistances[i]<=getMD(tree)){
          for(int j=0;j<dimensions;j++){
              innerX[innerCount*dimensions+j]=Xarr[i*dimensions+j];
          }
            idxInner[innerCount++]=idxArray[i];
      }
    }
    //Creating the current tree's , inner subtree , by the main (caller) thread
    tree->inner=vpt_pthreads(innerX,innerCount,dimensions,idxInner);
    rc= pthread_join(outerSubtreeThread,NULL);
    if (rc) {
      printf("ERROR; return code from pthread_join() is %d\n", rc);
      exit(-1);
    }

    //free allocated memory from the heap
    free(part);
    free(innerX);
    free(idxInner);

}




//////////////////////QUICKSELECT//////////////////////////////
/* ( source: https://www.geeksforgeeks.org/),  editted it for fitting
the needs of  this program.)
*/
//Swap for using in the quickselect algorithm
void swap(double * xp, double * yp)
{
    double  temp = *xp;
    *xp = *yp;
    *yp = temp;
}
// Standard partition process of QuickSort().
// It considers the last element as pivot
// and moves all smaller element to left of
// it and greater elements to right
int partition(double arr[], int l, int r)
{
    double x = arr[r];
    int i = l;
    for (int j = l; j <= (r - 1); j++) {
        if (arr[j] <= x) {
            swap(&arr[i],&arr[j]);
            i++;
        }
    }
    swap(&arr[i],&arr[r]);
    return i;
}

// This function returns k'th smallest
// element in arr[l..r] using QuickSort
// based method.  ASSUMPTION: ALL ELEMENTS
// IN ARR[] ARE DISTINCT
double kthSmallest(double arr[], int l, int r, int k)
{
    // If k is smaller than number of
    // elements in array
    if ( (k > 0) && k <=( r - l + 1) ) {

        // Partition the array around last
        // element and get position of pivot
        // element in sorted array
        int index = partition(arr, l, r);

        // If position is same as k
        if ((index - l) == (k - 1))
            return arr[index];

        // If position is more, recur
        // for left subarray
        if ((index - l )> (k - 1))
            return kthSmallest(arr, l, index - 1, k);

        // Else recur for right subarray
        return kthSmallest(arr, index + 1, r,k - index + l - 1);
    }

    // If k is more than number of
    // elements in array
    return INT_MAX;
}
