#ifndef TRANSACTION_H
#define TRANSACTION_H

#include "armadillo"

using namespace arma;
using namespace std;


class transaction
{
public:
    transaction(int N,int n_transactions,int m0,double lambda,double alpha,double gamma);
    vec wealth;
};

#endif // TRANSACTION_H
