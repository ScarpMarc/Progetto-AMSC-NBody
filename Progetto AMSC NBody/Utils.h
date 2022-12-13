#pragma once

inline unsigned int sumUpTo(const unsigned int& start, const unsigned int& end)
{
	unsigned int sum = 0;
	for (unsigned int i = start; i < end; ++i) sum += i;
	return sum;
}
