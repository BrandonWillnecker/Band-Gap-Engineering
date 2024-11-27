#pragma once
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <iostream>

std::vector<std::array<float, 4>> material_properties;
std::vector<std::string> material_names;
size_t n_materials = 0;

enum MATERIAL_PROPS
{
	BAND_GAP_ENERGY = 0,
	ELECTRON_AFFINITY,
	ELECTRON_MASS,
	HOLE_MASS
};

void load_material_properties()
{
	std::ifstream material_properties_file("material_properties.txt");
	material_properties_file >> n_materials;
	std::cout << n_materials << " materials in file\n";
	for (uint32_t i = 0; i < n_materials; i++) {
		std::string material_name;
		material_properties_file >> material_name;
		material_names.push_back(material_name);
		std::cout << i + 1 << ' ' << material_name << ' ';
		material_properties.push_back(std::array<float, 4>());
		for (uint32_t j = 0; j < 4; j++) {
			material_properties_file >> material_properties[i][j];
			std::cout << material_properties[i][j] << ' ';
		}
		std::cout << '\n';
	}
}
