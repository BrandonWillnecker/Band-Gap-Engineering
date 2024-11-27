#pragma once
#include <cstdint>

class Random
{
public:
	Random(unsigned int seed):m_seed(seed){}

	void reseed(unsigned int seed){m_seed = seed;}

	float get_float()
	{
		m_seed = pcg_hash(m_seed);
		return m_seed * SCALE;
	}

	unsigned int get_uint()
	{
		m_seed = pcg_hash(m_seed);
		return m_seed;
	}

private:
	const float SCALE = 1.0f / UINT32_MAX;
	unsigned int m_seed = 0u;

	unsigned int pcg_hash(unsigned int input)
	{
		m_seed = m_seed * 747796405u + 2891336453u;
		unsigned int state = m_seed;
		unsigned int word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
		return (word >> 22u) ^ word;
	}
};