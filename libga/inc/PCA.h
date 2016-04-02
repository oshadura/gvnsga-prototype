#ifndef __PCA__
#define __PCA__

#include "mlpack/core.hpp"
#include "mlpack/methods/pca/pca.hpp"
#include "mlpack/methods/kernel_pca/kernel_pca.hpp"

class PCA
{
public:
	PCA();
	~PCA();
	void CorrelationCheckPopulation();
	void InitialisePCA();
	void TransformPCA();
	void InitialiseKPCA();
	void TransformKPCA();
private:

};

#endif
