/*
 * CX 4220 / CSE 6220 Introduction to High Performance Computing
 *              Programming Assignment 1
 * 
 *  Serial polynomial evaluation algorithm function implementations goes here
 * 
 */

double poly_evaluator(const double x, const int n, const double* constants){
    //Implementation
    constants = constants + n;
    double result = *constants;
    for (int i=0; i<n; i++) {
        constants--;
        result = result*x + *constants;
    }
    return result;
}
