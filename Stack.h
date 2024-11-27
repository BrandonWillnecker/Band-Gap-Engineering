#pragma once
#include "Materials.h"

class Stack
{
public:
	Stack()
	{
		layer_count = 0;
		for (uint32_t i = 0; i < MAX_LAYER_COUNT; i++) {
			layer_types[i] = 0;
			layer_sizes[i] = 0.0f;
		}
	}

	Stack(uint32_t num_layers)
	{
		layer_count = num_layers;
		for (uint32_t i = 0; i < MAX_LAYER_COUNT; i++) {
			layer_types[i] = 0;
			layer_sizes[i] = 0.0f;
		}
	}

	Stack(const Stack& stack)
	{
		layer_count = stack.layer_count;
		cost = stack.cost;
		energy_transition_matrix = stack.energy_transition_matrix;

		for (uint32_t i = 0; i < layer_count; i++) {
			layer_types[i] = stack.layer_types[i];
			layer_sizes[i] = stack.layer_sizes[i];
		}
		
		for (uint32_t i = layer_count; i < MAX_LAYER_COUNT; i++) {
			layer_types[i] = 0;
			layer_sizes[i] = 0.0f;
		}
	}

	//Stack(Stack&& stack)
	//{
	//	layer_count = stack.layer_count;
	//	cost = stack.cost;

	//	for (uint32_t i = 0; i < layer_count; i++) {
	//		layer_types[i] = stack.layer_types[i];
	//		layer_sizes[i] = stack.layer_sizes[i];
	//	}

	//	for (uint32_t i = layer_count; i < MAX_LAYER_COUNT; i++) {
	//		layer_types[i] = 0;
	//		layer_sizes[i] = 0.0f;
	//	}
	//}



	void save_to_file(const std::string& filename)
	{

	}

	std::string str()
	{
		std::stringstream stack_str;
		stack_str << "(" << layer_count << ") ";
		for (uint32_t layer_index = 0; layer_index < layer_count - 1; layer_index++) 
			stack_str << material_names[layer_types[layer_index]] << ':' << layer_sizes[layer_index] << ", ";
		stack_str << material_names[layer_types[layer_count -1]] << ':' << layer_sizes[layer_count-1];
		return stack_str.str();

	}

	void set_energy_transition_matrix()
	{
		const float mu = 0.03f;

		float total_size = 0.0f;
		std::vector<float> cumulative_thickness(layer_count,0.0f);
		const float dz = 0.1f;

		for (uint32_t i = 0; i < layer_count;i++) {
			total_size += layer_sizes[i];
			cumulative_thickness[i] += total_size;
		}

		uint32_t PSI_SIZE = (uint32_t)(total_size / dz);
		Eigen::MatrixXf H_electron = Eigen::MatrixXf::Zero(PSI_SIZE, PSI_SIZE);
		Eigen::MatrixXf H_hole = Eigen::MatrixXf::Zero(PSI_SIZE, PSI_SIZE);

		float V_electron_max = -FLT_MAX;
		float V_hole_min = FLT_MAX;

		size_t layer_index = 0;
		for (size_t i = 0; i < PSI_SIZE; i++) {
			if (i * dz >= cumulative_thickness[layer_index]) layer_index++;

			float mu_electron = -mu / (dz * dz * material_properties[layer_types[layer_index]][MATERIAL_PROPS::ELECTRON_MASS]);
			float mu_hole = mu / (dz * dz * material_properties[layer_types[layer_index]][MATERIAL_PROPS::HOLE_MASS]);

			float V_electron = -material_properties[layer_types[layer_index]][MATERIAL_PROPS::ELECTRON_AFFINITY];
			float V_hole = V_electron - material_properties[layer_types[layer_index]][MATERIAL_PROPS::BAND_GAP_ENERGY];

			//V_electron += material_properties[6][MATERIAL_PROPS::BAND_GAP_ENERGY] + material_properties[6][MATERIAL_PROPS::ELECTRON_AFFINITY];
			//V_hole += material_properties[6][MATERIAL_PROPS::BAND_GAP_ENERGY] + material_properties[6][MATERIAL_PROPS::ELECTRON_AFFINITY];

			if (V_electron > V_electron_max) V_electron_max = V_electron;
			if (V_hole < V_hole_min) V_hole_min = V_hole;

			if (i >= 1) H_electron(i, i - 1) = mu_electron;
			H_electron(i, i) = -2.0f * mu_electron + V_electron;
			if (i + 1 < PSI_SIZE) H_electron(i, i + 1) = mu_electron;

			if (i >= 1) H_hole(i, i - 1) = mu_hole;
			H_hole(i, i) = -2.0f * mu_hole + V_hole;
			if (i + 1 < PSI_SIZE) H_hole(i, i + 1) = mu_hole;
		}

		Eigen::VectorXf electron_energy_levels = H_electron.eigenvalues().real();
		Eigen::VectorXf hole_energy_levels = H_hole.eigenvalues().real();

		std::vector<float> electron_energy_levels_sorted;
		std::vector<float> hole_energy_levels_sorted;

		//std::cout << "Electron states\n";
		int electron_state_index = 0;
		for (size_t i = 0; i < PSI_SIZE; i++) {
			if (electron_energy_levels(i) < V_electron_max) {
				//std::cout << electron_state_index << ' ' << electron_energy_levels(i) << '\n';
				electron_energy_levels_sorted.push_back(electron_energy_levels(i));
				electron_state_index++;
			}
		}

		//std::cout << "Hole states\n";
		int hole_state_index = 0;
		for (size_t i = 0; i < PSI_SIZE; i++) {
			if (hole_energy_levels(i) > V_hole_min) {
				//std::cout << hole_state_index << ' ' << hole_energy_levels(i) << '\n';
				hole_energy_levels_sorted.push_back(hole_energy_levels(i));
				hole_state_index++;
			}
		}

		std::sort(electron_energy_levels_sorted.begin(), electron_energy_levels_sorted.end());
		std::sort(hole_energy_levels_sorted.begin(), hole_energy_levels_sorted.end());

		size_t n_electron_states = electron_energy_levels_sorted.size();
		size_t n_hole_states = hole_energy_levels_sorted.size();

		energy_transition_matrix = Eigen::MatrixXf::Zero(n_electron_states, n_hole_states);

		for (size_t i = 0; i < n_electron_states; i++)
			for (size_t j = 0; j < n_hole_states; j++)
				energy_transition_matrix(i, j) = abs(hole_energy_levels_sorted[j] - electron_energy_levels_sorted[i]);

	}


public:
	static const uint32_t MAX_LAYER_COUNT = 5;
	static const uint32_t MIN_LAYER_COUNT = 3;
	uint32_t layer_types[MAX_LAYER_COUNT];
	float layer_sizes[MAX_LAYER_COUNT];
	uint32_t layer_count = 0;
	Eigen::MatrixXf energy_transition_matrix;
	float cost = 0.0f;
};


//constexpr float mu = 0.038102817f;