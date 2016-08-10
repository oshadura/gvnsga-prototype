
/*

// Implementation of the MO-CMA-ES
#include <shark/Algorithms/DirectSearch/MOCMA.h>
// Access to benchmark functions
#include <shark/ObjectiveFunctions/Benchmarks/Benchmarks.h>
		
int main( int argc, char ** argv ) {

	// Adjust the floating-point format to scientific and increase output precision.
	std::cout.setf( std::ios_base::scientific );
	std::cout.precision( 10 );

	// Instantiate both the problem and the optimizer.
	shark::DTLZ2 dtlz2;
	dtlz2.setNumberOfVariables( 3 );

	shark::MOCMA mocma;

	// Initialize the optimizer for the objective function instance.
	dtlz2.init();
	mocma.init( dtlz2 );

	// Iterate the optimizer
	while( dtlz2.evaluationCounter() < 25000 ) {
		mocma.step( dtlz2 );
	}

	// Print the optimal pareto front
	for( std::size_t i = 0; i < mocma.solution().size(); i++ ) {
		for( std::size_t j = 0; j < dtlz2.numberOfObjectives(); j++ ) {
			std::cout<< mocma.solution()[ i ].value[j]<<" ";
		}
		std::cout << std::endl;
	}
}
*/