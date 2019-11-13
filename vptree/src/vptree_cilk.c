#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vptree.h"
#include "limits.h"
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
//Defining a threashold for the distance parallelism
//(bellow that number of spots no parallelism happen in distance calculation)
#define DISTANCE_THRESHOLD 0
//there is no real reason to use threshold for the distances with cilk
//actually it makes our code run slightly slower is some tests, but we let this option


///////////////////////////
//////////////////////////
//Declaration of functions I created, they are defined lower in this file//

/*Function that helps build the vantage point tree*/
vptree * vpt(double *Xarray, int n, int d, int *idxArray);

/*! Function for computing the Median of the Euclidean Distances of the points
 to the vantage point (ed array), calling quickselect algorithm WITHOUT
 changing the original array's order of items by copying the original array .*/
double computeMedian(double * arr,int n);

/*! Sets the euclidean distances of the points to the vantage array the proper values*/
void computeEuclideanDistances(double *arr ,int numberOfSpots,int dimensions, vptree * tree,double * Xarr);

/*!Sets the inner and the outer subtree of the current tree.*/
void innerAndOuterTree(double *euclideanDistances ,int numberOfSpots,int dimensions, vptree * tree,double * Xarr,int * idxArray);

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
  root= vpt(X, n, d, IdxArr );

  free(IdxArr);
  return root;


}
//The recursive function that build every node of the tree
vptree * vpt(double *Xarray, int n, int d, int *idxArray){

  //If there is no points , there is no a vp(sub)tree (node).
  if(n==0){
    return NULL;
  }
  //declaring the vantage point tree
  vptree * T = (vptree *)malloc(sizeof(vptree));

  //Vector ed[n] stores the euclidean distances of each point to the vantage point.
  //ed[i] is the Euclidean distance of X[i][] point from the vantage point
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
  computeEuclideanDistances(ed,n,d,T,Xarray);

  //Getting the Median of the Euclidean Distances of the points
  //to the vantage point (median of ed array)
  T->md=computeMedian(ed,n);

  /////////////////////////////////////////////
  //3.GETTING THE INNER AND THE OUTER TREES OF THE VPTREE.
  innerAndOuterTree(ed,n,d,T,Xarray,idxArray);

  /////////////////////////////////////////////
  return T;
  }


vptree * getInner(vptree * T){
  return (T->inner);
}


vptree * getOuter(vptree * T){
  return (T->outer);
}


double getMD(vptree * T){
  return (T->md);
}


double * getVP(vptree * T){
  return (T->vp);

}


int getIDX(vptree * T){
  return (T->idx);

}

////////////////////////////////////////////////////////////////////////////////
//MY FUNCTIONS Definition
double computeMedian(double * euclideanDistances,int numberOfSpots){

  //Coppying the euclideanDistances array into a new array in order not to alter
  //the original  euclideanDistances(ed) array.
  double array[numberOfSpots-1];
  for (size_t i = 0; i < numberOfSpots-1; i++) {
      array[i]=euclideanDistances[i];
  }

  //If there is only one point in the given points set,
  //of course there is no median distance.
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

//the sequential version of computeEuclideanDistances() functions with cilk anotations
//to run in parallel if there are more points than the treshold( threshold =0 , means no thresshold)
void computeEuclideanDistances(double *euclideanDistances ,int numberOfSpots,
  int dimensions, vptree * tree,double * Xarr){

  if(numberOfSpots>=DISTANCE_THRESHOLD){
    //here we use cilk for set parts of the loop in other threads
    cilk_for(int i = 0; i < (numberOfSpots-1); i++) {
      double temp=0;
      for (int j = 0; j < dimensions; j++) {
        temp+=pow(  (Xarr[i*dimensions+j]-(getVP(tree))[j]), 2 );
      }
      euclideanDistances[i]=sqrt(temp);
    }
  }
  else{
    for(int i = 0; i < (numberOfSpots-1); i++) {
      double temp=0;
      for (int j = 0; j < dimensions; j++) {
        temp+=pow(  (Xarr[i*dimensions+j]-(getVP(tree))[j]), 2 );
      }
      euclideanDistances[i]=sqrt(temp);
    }
  }
}

//the sequential version of innerAndOuterTree() functions with cilk anotations
//to run in parallel if there are more points than the treshold
void innerAndOuterTree(double *euclideanDistances ,int numberOfSpots,int dimensions, vptree * tree,double * Xarr, int * idxArray ){
  //innerX and outerX are used as dataset inputs for the inner and outer sub(vp)tree
  double *innerX,*outerX;

  //idxInner and idxOuter are used as IdxArr for the inner and outer subtrees
  //innerSize and outerSize are used for determining the inner and outer subtrees size
  //innerCount and outerCount are used for scaning given iputs and filling inner and outer trees
  int *idxInner, *idxOuter, innerSize=0, innerCount=0, outerSize=0, outerCount=0;


    //Getting the size of the inner and the outer subtrees.
    for(int i=0;i<(numberOfSpots-1);i++){
      if(euclideanDistances[i]<=getMD(tree))
        innerSize++;
    }
    outerSize=(numberOfSpots-1)-innerSize;

    //Allocating memory for inner and outer subtrees datasets
    innerX=(double *)malloc( innerSize * dimensions * sizeof(double ));
    outerX=(double *)malloc( outerSize * dimensions * sizeof(double ));


    //Allocating memory for innerIdx and outerIdx
    idxInner=(int *)malloc( (innerSize)* sizeof(int));
    idxOuter=(int *)malloc( (outerSize)* sizeof(int));

  //Gettin the idxInner and idxOuter andfilling inner and outer subtree
  //by scanning the euclidean distances of n-1 remaining spots
  //(ln. 231 changed after the end of deadline (cilk_for -> for) , we did this minor change to avoid race Condition )
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


  /*Bellow we create a continuation to permit parallelism using cilk spawn.
  We avoid to add a cilk_spawn attribute to the second recursive call to vpt()
   because that will create an empty continuation.*/
  //Creating the current tree's , inner subtree
  tree->inner=cilk_spawn vpt(innerX,innerCount,dimensions,idxInner);

  //Creating the current tree's , outer subtree
  tree->outer=vpt(outerX,outerCount,dimensions,idxOuter);
  cilk_sync;
  //free allocated memory from the heap
  free(innerX);
  free(idxInner);
  free(outerX);
  free(idxOuter);
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
