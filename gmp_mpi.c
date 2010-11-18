#include <gmp.h>
#include <string.h>
#include "aks.h"

#define mpz_buffersize (z) \
    (sizeof(int) * 2 + sizeof(mp_limb_t) * z->_mp_alloc)


int mpi_mpz_send(mpz_t z, int dest, int tag, MPI_Comm comm) {
    void *buf = NULL;
    int buf_size;

    buf = mpz_get_str(NULL, 10, z);
    if (!buf) {
        fprintf(stderr, "mpi_mpz_send:  mpz_get_str() failed.\n");
        return -1;
    }

    buf_size = sizeof(char) * (strlen(buf) + 1);

    MPI_Send(&buf_size, 1, MPI_INT, dest, tag, comm);
    MPI_Send(buf, buf_size, MPI_BYTE, dest, tag, comm);
    snd_cnt++;
    return MPI_SUCCESS;
}



int mpi_mpz_recv(mpz_t z, int src, int tag, MPI_Comm comm, MPI_Status *status) {
    int buf_size;
    int rc = 0;
    void *buf;

    MPI_Recv(&buf_size, 1, MPI_INT, src, tag, comm, status);

    buf = (void *)malloc(sizeof(unsigned char) * buf_size);
    MPI_Recv(buf, buf_size, MPI_BYTE, src, tag, comm, status);

    rc = mpz_set_str(z, buf, 10);
    if (rc) {
        fprintf(stderr, "mpi_mpz_bcast:  mpz_set_str(%s) failed.\n",
            (char*)buf);
        return 0;
    }
 
    recv_cnt++;
    return 1;
}



int mpi_mpz_bcast(mpz_t z, int root, MPI_Comm comm, int id) {
    int buf_size;
    int rc = 0;
    void *buf;

    if (root == id) {
        buf = mpz_get_str(NULL, 10, z);
        if (!buf) {
            fprintf(stderr, "mpi_mpz_bcast:  mpz_get_str() failed.\n");
            return 0;
        }
        buf_size = sizeof(char) * (strlen(buf) + 1);
    }

    MPI_Bcast(&buf_size, 1, MPI_INT, root, comm);

    if (root != id)
        buf = (void*)malloc(sizeof(unsigned char) * buf_size);

    MPI_Bcast(buf, buf_size, MPI_BYTE, root, comm);

    if (root != id) {
        rc = mpz_set_str(z, buf, 10);
        if (rc) {
            fprintf(stderr, "mpi_mpz_bcast:  mpz_set_str(%s) failed.\n",
                (char*)buf);
            return 0;
        }
    }

    return 1;
}

