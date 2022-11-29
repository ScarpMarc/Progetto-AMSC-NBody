#include <iostream>
#include "Particle.h"
#include "Vector.h"

using namespace std;

int main()
{
	Vector<3> vec3({ 1.6, 2.7, 3.8 });
	Vector<2> vec2({ 7.4, 1.8 });
	cout << vec3 << endl << vec2 << endl;
}