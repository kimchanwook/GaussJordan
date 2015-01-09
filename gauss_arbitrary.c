#include<stdio.h>
#include<stdlib.h>
#include<math.h>

/*swap the rows*/
void change(double* a, double* b);
void changerows(double* A, int M, int row1, int row2);
/*find the pivot: largest element in the column*/
int pivot(double* A, int M, int r);
/*into Triangular Matrix*/
void nopivoting_triangular(double* A, int M, double*b);//no pivoting, making matrix triangular form
void pivoting_triangular(double* A, int M, double* b);//partial pivoting
/*BACK SUBSTITUTE*/
void backsubstitute(double* A, int M, double* b, double* x);
/*solving the matrix: putting all together*/
void nopartialpivot(double* A, int M, double* b, double* x);
void partialpivot(double* A, int M, double* b, double* x);


/*MAIN FUNCTION*/
int main(int argc, char** argv){

//////////////////////////////////////////////////
  /*Assign arguments to contruct A
    First argument is the dimention of the matrix
    Then, put in matrix elements for A
    Then, put in elements for b[]*/
//////////////////////////////////////////////////

  int N = atoi(argv[1]) * atoi(argv[1]);
  double A[N];
  for(int i=0; i<N; i++) A[i] = atof(argv[i+2]);

  //Construct b
  int M = atoi(argv[1]);//M is the dim of our matrix
  double b[M];
  for(int i=0; i<M; i++) b[i] = atof(argv[N+i+2]);

  //x
  double x[M];

  //print out the matrix A before Gauss
  printf("\n The matrix A \n\n");
  for(int i=0; i<M; i++){
    for(int j=0; j<M; j++)
      printf(" %f ", A[i*M+j]);
      printf("\n");}
      printf("\n");

////////////////////////////////////////////////////////////////////////////////

  //nopartialpivot function to get x
  nopartialpivot(A, M, b, x);

  //print out the matrix A after Gauss
  printf("\n The matrix A after No Pivot  Gauss Elimination \n\n");
  for(int i=0; i<M; i++){
    for(int j=0; j<M; j++)
      printf(" %f ", A[i*M+j]);
    printf("\n");}
  printf("\n");

  //print out result x
  printf("\n solution vector x \n\n");
  for(int i=0; i<M; i++) printf("%f ",x[i]);
  printf("\n");
  printf("\n");

  //compute det A
  double detA = 1;
  for(int i=0; i<M; i++){
    for(int j=0; j<M; j++){
      if(A[i*M+j] != 0){
        detA = detA * A[i*M+j];}}}

  //print out det A
  printf("determinant of the matrix is: %f \n", detA);

///////////////////////////////////////////////////////////////////////////////

  //partialpivot function to get x
  partialpivot(A, M, b, x);

  //print out the matrix A after pivoting + Gauss
  printf("\n The matrix A after Pivoting  Gauss Elimination \n\n");
  for(int i=0; i<M; i++){
    for(int j=0; j<M; j++)
      printf(" %f ", A[i*M+j]);
    printf("\n");}
  printf("\n");

  //print out result x
  printf("\n solution vector x \n\n");
  for(int i=0; i<M; i++) printf("%f ",x[i]);
  printf("\n");
  printf("\n");

  //compute det A
  double detA2 = 1;
  for(int i=0; i<M; i++){
    for(int j=0; j<M; j++){
      if(A[i*M+j] != 0){
        detA2 = detA2 * A[i*M+j];}}}

  //print out det A
  printf("determinant of the matrix is: %f \n", detA2);

return 0;
}


/////////////////////
/*Functions Details*/
/////////////////////

/*swap the rows*/
void change(double* a, double* b){
  double c = *a; //save the address of a into a new variable c
  *a = *b; //store address of b into address of a
  *b = c;} //store address of a into address of b

void changerows(double* A, int M, int row1, int row2){
  for(int i=0; i<M; i++){
    change(A + row1*M + i, A + row2*M + i);}} //use change function M times to change each elements in 
                                          //the row1 to row2 and vice versa
                                          //POINTER = ARRAY!!!!!!!!!!!!!!!! 

/*find the pivot: largest element in the column*/
int pivot(double* A, int M, int r){
  int pivot = r;
  double diagonal = fabs(A[r*M+r]);//magnitude of diagonal element
  for(int nrsc = r+1; nrsc < M; ++nrsc){//nrsc = next row same column
    if(fabs(A[nrsc*M + r]) > diagonal ){
      pivot = nrsc;
      diagonal = fabs(A[nrsc*M + r]);}}//set it to compare it to next nrsc
  return pivot;}

/*into Triangular Matrix*/
void nopivoting_triangular(double* A, int M, double*b){//no pivoting, making matrix triangular form
  for(int r=0; r<M-1; r++){
    double pivot = A[r*M + r];
    if(fabs(pivot) < 1e-10){/*IMPORTANT: testing to see if the matrix is linearly independent*/
      printf("matrix is linearly independent: not solvable");
      abort();}
    /*subtracting the elements with right */
    for(int nrsc = r+1; nrsc < M; ++nrsc){//nrsc=next row same column
      double ratio = -A[nrsc*M + r] / pivot;//need this ratio to subtract the nrsc to be zero
      /*for loop to iterate through each columns*/
      for(int c=r; c<M; c++){
        A[nrsc*M + c] = A[nrsc*M + c] + ratio*A[nrsc*M + c];
        b[nrsc] = b[nrsc] + ratio*b[nrsc]; }}}}

void pivoting_triangular(double* A, int M, double* b){//partial pivoting
  for(int r=0; r<M-1; r++){
    /*pivoting part*/
    int p = pivot(A,M,r);//find the pivot: largest element in the column
    if(p != r){
      changerows(A, M, r, p);//change the rows if the largest element is not our orig. diag. element
      change(b+r, b+p);}//change the elements in b[] also USE b+r, b+p NOTATION!!!

    double pivot = A[r*M + r];
    if(fabs(pivot) < 1e-10){/*IMPORTANT: testing to see if the matrix is linearly independent*/
      printf("matrix is linearly independent: not solvable");
      abort();}
    /*subtracting the elements with right */
    for(int nrsc = r+1; nrsc < M; ++nrsc){//nrsc=next row same column
      double ratio = -A[nrsc*M + r] / pivot;//need this ratio to subtract the nrsc to be zero
      /*for loop to iterate through each columns*/
      for(int c=r; c<M; c++){
        A[nrsc*M + c] = A[nrsc*M + c] + ratio*A[nrsc*M + c];
        b[nrsc] = b[nrsc] + ratio*b[nrsc]; }}}}

/*BACK SUBSTITUTE*/
void backsubstitute(double* A, int M, double* b, double* x){
  for(int r=M-1; r > -1; --r){//going UPWARD: x_n to solve for x_n-1
    double temp = b[r];
    for(int c = r+1; c<M; ++c) temp -= A[r*M+c] * x[c];//just rearranged our equation
    x[r] = temp / A[r*M+r];}}//this give our solution to the system

/*solving the matrix: putting all together*/
void nopartialpivot(double* A, int M, double* b, double* x){
  nopivoting_triangular(A, M, b);
  backsubstitute(A, M, b, x);}

void partialpivot(double* A, int M, double* b, double* x){
  pivoting_triangular(A, M, b);
  backsubstitute(A, M, b, x);}

