/*
 * CX 4220 / CSE 6220 Introduction to High Performance Computing
 *              Programming Assignment 1
 * 
 *  MPI polynomial evaluation algorithm function implementations go here
 * 
 */

#include "mpi_evaluator.h"
#include "const.h"
#include <cmath>
#include <iostream>

void scatter(const int n, double* scatter_values, int &n_local, double* &local_values, int source_rank, const MPI_Comm comm){
    //Implementation
    int size, rank, n_local_original;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    if(rank == 0) {
        int load = ((n % size) == 0)? 0 : 1;
        n_local = (n / size) + load;
        local_values = (double*)malloc(sizeof(double)*n_local);

        n_local_original = n_local;
        for(int i=0;i < n_local;i++) {
            local_values[i] = scatter_values[i];   //assign to master processor         
        }

        int n_left =n-n_local; // left n number need to scatter.

        for(int i=1;i <size;i++) {
            if(n_left > n_local) {
                MPI_Send(&n_local,1,MPI_INT,i,1,comm);
                n_left = n_left - n_local;
                MPI_Send(&scatter_values[n_local*i], n_local, MPI_DOUBLE, i, 1, comm);
            }
            else {
                n_local = n_left;
                MPI_Send(&n_left,1,MPI_INT,i,1,comm);
                n_left = n_local-n_left  ;
                MPI_Send(&scatter_values[n_local_original*i], n_local, MPI_DOUBLE, i, 1, comm);  
            }  
        }
    }

    else {
        MPI_Status stat;
        MPI_Recv(&n_local,1,MPI_INT,0,1,comm, &stat);
        local_values = (double*)malloc(sizeof(double)*n_local);
        MPI_Recv(&local_values[0], n_local, MPI_DOUBLE, 0, 1, comm, &stat);       
    }

    if(rank == 0) {
        n_local = n_local_original;        
    }
}

double broadcast(double value, int source_rank, const MPI_Comm comm){
    //Implementation
    int size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    if(rank == source_rank) { 
        for(int i=1;i <size;i++) {   
            MPI_Send(&value,1,MPI_DOUBLE,i,11,comm);
        }
    }
    else {
        MPI_Status status;
        MPI_Recv(&value, 1, MPI_DOUBLE,source_rank,11,comm,&status); 
    }
    return value;
}

void parallel_prefix(const int n, const double* values, double* prefix_results, const int OP, const MPI_Comm comm){
    //Implementation
    int size, rank;
    int total_size,power=0;       
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    while(1) {
        int x=pow(2,power); // to get lg(n) size
        if(x > size){
            total_size = x;
            break;
        }   
        else {
            power+=1;            
        }
    }
    double *local_values = (double*)malloc(n *sizeof(double));
    double local_result;
    double global_result;
          
    if (OP == 1) { //PREFIX_OP_SUM
        local_result = 0;
        local_values[0] = values[0];

        for (int i=1; i<n; i++) {
            local_values[i] = local_values[i-1] + values[i];
        }

        global_result = local_values[n-1];
        for (int i = 0; total_size != 1; i++, total_size >>= 1) {
            int des_rank = rank ^ (1 << i);
            double recv;
            if ((des_rank > rank) && (des_rank < size)) {
                MPI_Send(&global_result,1,MPI_DOUBLE, des_rank,111,comm);
                MPI_Status status;
                MPI_Recv(&recv,1,MPI_DOUBLE, des_rank,111,comm,&status);
                global_result = recv;
            }
            else {   
                if(des_rank < size) {
                    MPI_Status status;
                    MPI_Recv(&recv,1,MPI_DOUBLE,des_rank,111,comm,&status);
                    global_result += recv;
                    local_result += recv;
                    MPI_Send(&global_result,1,MPI_DOUBLE,des_rank,111,comm);
                }
            }
        }

        for (int i = 0; i < n; i++) {
            prefix_results[i] = local_values[i] + local_result;
        }
    }
         //PREFIX_OP_PRODUCT
    else if (OP == 2) {
        local_result = 1;
        local_values[0] = values[0];

        for (int i = 1; i < n; i++) {
            local_values[i] = local_values[i-1]*values[i];
        }

        global_result = local_values[n-1];
        for (int i = 0; total_size != 1; i++, total_size >>= 1) {
            int des_rank = rank ^ (1 << i);
            double recv;
            if ((des_rank > rank) && (des_rank < size)) {
                MPI_Send(&global_result,1,MPI_DOUBLE,des_rank,111,comm);
                MPI_Status status;
                MPI_Recv(&recv,1,MPI_DOUBLE, des_rank,111, comm,&status);
                global_result = recv;
            }
            else {   
                if(des_rank < size){
                    MPI_Status status;
                    MPI_Recv(&recv,1,MPI_DOUBLE, des_rank,111,comm,&status);
                    global_result *= recv;
                    local_result *= recv;
                    MPI_Send(&global_result,1,MPI_DOUBLE, des_rank,111,comm);
                }
            }
        }

        for (int i = 0; i < n; i++) {
            prefix_results[i] = local_values[i] * local_result;
        }
      }


}

double mpi_poly_evaluator(const double x, const int n, const double* constants, const MPI_Comm comm){
    //Implementation
    int size, rank;
    double final_value;
    double* values = (double*)malloc(n *sizeof(double));
    double* results = (double*)malloc(n *sizeof(double));
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    values[0] = (rank == 0)? 1 : x;
    for (int i=1;i<n;i++) {
        values[i] = x;
    }

    parallel_prefix(n, values, results, 2, comm);
    for (int i=0;i<n;i++) {
        values[i] = results[i] * constants[i];
    }
    parallel_prefix(n, values, results, 1, comm);

    if (rank == size-1) {
        final_value = results[n -1];
        MPI_Send(&final_value,1,MPI_DOUBLE,0,1111,comm);
    }
    if (rank == 0) {
        MPI_Status status;
        MPI_Recv(&final_value,1,MPI_DOUBLE,size-1,1111,comm,&status);
    }
    final_value= broadcast(final_value,0,comm);
    return final_value;
}
