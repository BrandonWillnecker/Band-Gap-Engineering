#pragma once
constexpr uint32_t GENERATION_SIZE = 150;
constexpr uint32_t MAX_GENERATION_STEPS = 200;
constexpr uint32_t ELITISM_SIZE = 5;
constexpr uint32_t TOURNEMENT_SIZE = 10;
constexpr float MUTATION_RATE = 0.9f;

std::array<Stack, GENERATION_SIZE> generation;

struct index_pair
{
    uint32_t index1;
    uint32_t index2;
};


void set_cost(Stack& stack, Eigen::MatrixXf& energy_transition_matrix)
{
    for (uint32_t i = 0; i < stack.layer_count-1; i++) {
        if (stack.layer_types[i] == stack.layer_types[i + 1]) {
            stack.cost = FLT_MAX;
            return;
        }
    }

    stack.set_energy_transition_matrix();
    stack.cost = 0.0f;
    const float mismatch_factor = 100.0f;

    if (stack.energy_transition_matrix.rows() == 0 || stack.energy_transition_matrix.cols() == 0) {
        stack.cost = FLT_MAX;
        return;
    }

    uint32_t min_rows = (uint32_t)std::min(stack.energy_transition_matrix.rows(), energy_transition_matrix.rows());
    uint32_t min_cols = (uint32_t)std::min(stack.energy_transition_matrix.cols(), energy_transition_matrix.cols());

    for (size_t row = 0; row < min_rows; row++) {
        for (size_t col = 0; col < min_cols; col++) {
            if (energy_transition_matrix(row, col) > 0) {
                float diff = stack.energy_transition_matrix(row, col) - energy_transition_matrix(row, col);
                stack.cost += diff * diff;
            }
        }
    }

    stack.cost += mismatch_factor * (abs(stack.energy_transition_matrix.rows() - energy_transition_matrix.rows()) + abs(stack.energy_transition_matrix.cols() - energy_transition_matrix.cols()));
}

void init_generation(Random& random, Eigen::MatrixXf& energy_transition_matrix)
{
    for (uint32_t generation_index = 0; generation_index < GENERATION_SIZE; generation_index++)
    {
        uint32_t num_layers = Stack::MIN_LAYER_COUNT + random.get_uint() % (Stack::MAX_LAYER_COUNT-Stack::MIN_LAYER_COUNT + 1);
        generation[generation_index] = Stack(num_layers);
        for (size_t layer_index = 0; layer_index < num_layers; layer_index++)
        {
            generation[generation_index].layer_types[layer_index] = random.get_uint() % n_materials;
            generation[generation_index].layer_sizes[layer_index] = 5.0f + random.get_float() * 10.0f;
        }
        set_cost(generation[generation_index],energy_transition_matrix);
        std::cout << "[" << generation_index << "] " << generation[generation_index].str() << " " << generation[generation_index].cost << "\n";
    }
}


void mutate(Stack& stack, Random& random)
{
    uint32_t rand_layer = random.get_uint() % stack.layer_count;
    if (random.get_float() < 0.2f) {
        //Change the type
        stack.layer_types[rand_layer] = random.get_uint() % n_materials;
    }
    else {
        //Change the thickness
        stack.layer_sizes[rand_layer] += 0.2f * (2.0f * random.get_float() - 1.0f);
        stack.layer_sizes[rand_layer] = std::max(stack.layer_sizes[rand_layer], 5.0f);
    }

    if (random.get_float() < 0.1f && stack.layer_count < Stack::MAX_LAYER_COUNT) {
        //Add a layer
        stack.layer_types[stack.layer_count] = random.get_uint() % n_materials;
        stack.layer_sizes[stack.layer_count] = 5.0f  + random.get_float() * 10.0f;
        stack.layer_count++;
    }
}

index_pair make_parent_selection(Random& random)
{
    //First parent
    uint32_t index1 = 0;
    float cost1 = FLT_MAX;
    for (uint32_t i = 0; i < TOURNEMENT_SIZE; i++){
        uint32_t index = random.get_uint() % GENERATION_SIZE;
        if (generation[index].cost < cost1) {
            cost1 = generation[index].cost;
            index1 = index;
        }
    }

    //second parent
    uint32_t index2 = 0;
    float cost2 = FLT_MAX;
    for (uint32_t i = 0; i < TOURNEMENT_SIZE; i++) {
        uint32_t index = random.get_uint() % GENERATION_SIZE;
        if (generation[index].cost < cost1) {
            cost2 = generation[index].cost;
            index2 = index;
        }
    }

    return { index1,index2 };
}

Stack crossover(index_pair index_pair, Random& random)
{
    Stack child1(generation[index_pair.index1]);
    Stack child2(generation[index_pair.index2]);

    uint32_t min_layer_count = std::min(generation[index_pair.index1].layer_count, generation[index_pair.index2].layer_count);
    uint32_t rand_layer = random.get_uint() % min_layer_count;

    uint32_t tmp_layer_type = generation[index_pair.index1].layer_types[rand_layer];
    float tmp_layer_size = generation[index_pair.index1].layer_sizes[rand_layer];

    child1.layer_sizes[rand_layer] = child2.layer_sizes[rand_layer];
    child1.layer_types[rand_layer] = child2.layer_types[rand_layer];

    child2.layer_sizes[rand_layer] = tmp_layer_size;
    child2.layer_types[rand_layer] = tmp_layer_type;

    if (random.get_float() < 0.5) return child1;
    return child2;
}