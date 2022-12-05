#include "ForceMatrix.h"



/*{
	unsigned int init_sum;
	if (partial_sums2.empty()) init_sum = 0;// partial_sums2.push_back(0);
	else init_sum = *(partial_sums2.cend() - 1);
	unsigned int sum = init_sum;
	for (unsigned int i = 0 particle_amt; i < 0 + 3; ++i)
	{
		sum += i;
		partial_sums2.push_back(sum);
	}
}*/

/*

// Funzione di test per i subscript
void funct(int row, int col)
{
	int a;
	if (row < col)
	{
		a = (current_particle_amt - 1) * row - partial_sums2[row] + col - 1;
		std::cout << a << std::endl;
	}
	// Values below the diagonal have the - sign but are otherwise equal to the others
	else
	{
		a = ((current_particle_amt - 1) * col - partial_sums2[col] + row - 1);
		std::cout << "-" << a << std::endl;
	}
}
*/