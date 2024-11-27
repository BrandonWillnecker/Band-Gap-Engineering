#pragma once
#include <iostream>
#include "Eigen/Dense"
#include "Stack.h"
#include "Random.h"
#include "GA_functions.h"
#include <thread>
#include <functional>

int main()
{
    std::thread::id threadId = std::this_thread::get_id();
    std::hash<std::thread::id> hasher;
    uint32_t threadIdHash = (uint32_t)hasher(threadId);
    uint32_t seed = time(0);
    std::cout << seed << "\n";
    seed += threadIdHash;
    std::cout << seed << "\n";
    Random random(seed);

    Eigen::MatrixXf energy_transition_matrix(3,5);
    //energy_transition_matrix << 0.0f, 0.0f, 1.828683f, 0.0f, 1.652509f,
    //                            0.0f, 2.317309f, 0.0f, 2.062825f, 0.0f,
    //                            0.0f, 0.0f, 2.717855f, 0.0f, 2.541682f;

    energy_transition_matrix << 2.04583665f, 1.90913914f, 1.828683f, 1.7142124f, 1.652509f,
                                2.40676002f, 2.317309f, 2.15735699f, 2.062825f, 2.02533618f,
                                2.87265923f, 2.73596172f, 2.717855f, 2.54103498f, 2.541682f;

    std::cout << "Energy Transition matrix :\n" << energy_transition_matrix << "\n\n";

    std::cout << "Loading matetrials\n";
    load_material_properties();
    std::cout << "\n\n";

    std::cout << "Initializing the generation: " << GENERATION_SIZE << " individuals\n";
    init_generation(random, energy_transition_matrix);
    std::cout << "\n\n";

    std::cout << "Sort the generation\n";
    std::sort(generation.begin(), generation.end(),
        [](Stack& stack1, Stack& stack2) {return stack1.cost < stack2.cost; }
    );

    std::cout << "Initial costs\n";
    for (uint32_t i = 0; i < GENERATION_SIZE; i++) {
        std::cout << generation[i].str() << " " << generation[i].cost << "\n";
    }
    std::cout << "\n\n";

    std::cout << "GA starting\n";
    for (uint32_t generation_step = 0; generation_step < MAX_GENERATION_STEPS; generation_step++)
    {
        std::cout << generation_step << " " << generation[0].str() << " " << generation[0].cost << "\n";

        //Sort the generation
        std::sort(generation.begin(), generation.end(),
            [](Stack& stack1, Stack& stack2) {return stack1.cost < stack2.cost; }
        );

        //Start new generation with the elites
        std::array<Stack,GENERATION_SIZE> new_generation;
        for (uint32_t i = 0; i < ELITISM_SIZE; i++)
            new_generation[i] = generation[i];

        //crossover and mutation
        for (uint32_t i = ELITISM_SIZE; i < GENERATION_SIZE; i++)
        {
            Random random(seed + generation_step + i);
            //parent selection
            index_pair selection_pair = make_parent_selection(random);
            Stack child = crossover(selection_pair, random);

            //child mutation
            if (random.get_float() < MUTATION_RATE) mutate(child,random);

            //Set cost
            set_cost(child, energy_transition_matrix);

            //add to new generation
            new_generation[i] = child;
        }

        std::move(new_generation.begin(), new_generation.end(), generation.begin());
    }

    std::cout << generation[0].str() << "\n";
    std::cout << generation[0].cost << "\n";
    std::cout << generation[0].energy_transition_matrix << "\n";
}
