/******************************************************************************
 * main_brkga.cpp: Main Process.
 *
 * Author: Carlos Eduardo de Andrade <ce.andrade@gmail.com>
 *
 * (c) Copyright 2012-2018, Carlos Eduardo de Andrade. All Rights Reserved.
 *
 * Created on : Jun 06, 2012 by andrade
 * Last update: Apr 26, 2014 by andrade
 *
 * This code is released under LICENSE.md
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *****************************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <cmath>
//#include <iterator>
//#include <cstdlib>
#include <vector>
#include <boost/timer/timer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

using namespace std;
using namespace boost;

#include "k_imdmtsp_instance.hpp"
#include "concorde_wrapper.hpp"
#include "k_imdmtsp_decoder.hpp"
#include "brkga.hpp"
#include "mtrand.hpp"

// Mauro code
#include "./mauro_lib/imdmtsp.h"

//-------------------------[ Some control constants ]-------------------------//

// Controls stop criteria
enum StopRule{GENERATIONS = 'G', TARGET = 'T', IMPROVEMENT = 'I'};

// Scenarios
enum Scenario{ST = 'T', SL = 'L', SQ = 'Q', FREE = 'R'};

ostream& operator<<(ostream &os, const Scenario &scenario) {
    switch(scenario) {
    case ST:
        os << "ST";
        break;
    case SL:
        os << "SL";
        break;
    case SQ:
        os << "SQ";
        break;
    case FREE:
        os << "free";
        break;
    }
    return os;
}


//---------------------------[ Initial chromosome ]---------------------------//

// Build a chromosome based on k-means heuristic
void chromosomeFromKMeans(Population::Chromosome *chromosome,
                          K_IMDMTSP_Instance &k_instance,
                          MTRand &rng, const unsigned max_cycle_size, ConcordeWrapper *cw) {
    // Mauro's instance.
    Instance *instance = new Instance(k_instance.graph, k_instance.coords, k_instance.dist);
    instance->K = k_instance.K;
    instance->has_coords = k_instance.has_coords;
    instance->n = instance->graph.nodeNum();

    // Mauro's solution.
    Solution *solution = kmeans_heuristic(instance);
    solution->concorde = cw;
    solution->improve();

    // Fix the cycles sizes
    for(vector<vector<FullGraph::Node> >::iterator itCycle = solution->subcycles.begin();
        itCycle != --solution->subcycles.end(); ++itCycle) {

        while(itCycle->size() > max_cycle_size) {
            int index;

            do index = rng.randInt(k_instance.K - 1);
            while(!(solution->subcycles[index].size() < max_cycle_size));

            solution->subcycles[index].push_back(itCycle->back());
            solution->cycle_with_node[itCycle->back()] = index;
            itCycle->pop_back();
        }
    }

    // Take some random keys and sort them.
    Population::Chromosome temp(instance->n);
    for(unsigned i = 0; i < temp.size(); ++i)
        temp[i] = rng.randInt(instance->n + 1, numeric_limits< Population::Allele >::max());

    sort(temp.begin(), temp.end(), greater<Population::Allele>());

    // Let's build the chromosome.
    unsigned i = 0;
    unsigned offset = k_instance.K;
    for(vector<vector<FullGraph::Node> >::iterator itCy = solution->subcycles.begin();
        itCy != (--solution->subcycles.end()); ++itCy, ++i) {

        (*chromosome)[i] = itCy->size();

        vector<FullGraph::Node>::iterator itNode = itCy->begin();
        (*chromosome)[k_instance.K + k_instance.graph.id(*itNode)] = temp[i];

        ++itNode;
        for(; itNode != itCy->end(); ++itNode) {
            (*chromosome)[k_instance.K + k_instance.graph.id(*itNode)] = temp[offset++];
        }
    }

    /////////////////////////////
    // Check chromosome sanity //
    /////////////////////////////
    for(unsigned j = 0; j < k_instance.K; ++j) {
        if((*chromosome)[j] == 0)
            (*chromosome)[j] = 1;
    }

    for(unsigned j = k_instance.K; j < unsigned(instance->n + k_instance.K); ++j) {
        if((*chromosome)[j] < unsigned(instance->n + 1))
            (*chromosome)[j] = rng.randInt(instance->n + 1, numeric_limits< Population::Allele >::max());
    }

    delete solution;
    delete instance;
}

//--------------------------------[ Main ]------------------------------------//

int main(int argc, char* argv[]) {
    // First read parameters from command line:
    #ifndef TUNING
    if(argc < 14) {
    #else
    if(argc < 24) {
    #endif
        cerr << "usage: " << argv[0]
             << " <config-file> <seed> <start-run> <runs> <stop-rule> <stop-arg>"
             << " <max-time> <instance-file> <scenario> <K> <max-cycle-size> <alpha> <results-file>"
             #ifdef TUNING
             << " <max_population_size> <elite-percentage> <mutants-percentage> <biasing-rhoe>"
             << " <number-of-populations> <exchange_interval> <num_exchange_indivuduals>"
             << " <reset_interval> <kick_type> <stallcount>"
             #endif
             << "\nwhere: "
             << "\n - <config-file>: parameters of BRKGA algorithm"
             << "\n - <seed>: seed for random generator"
             << "\n - <start-run>: number of the start run"
             << "\n - <runs>: number of runs to be executed"
             << "\n - <stop-rule> <stop-arg>: stop rule and its arguments where"
             << "\n\t+ Generations <number_generations>: the algorithm runs until <number_generations>"
             << "\n\t+ Target <value of expected target>: runs until obtains the target value"
             << "\n\t+ Iterations <max generations without improvement>: runs until the solutions don't"
             << "\n\t  improve by max generations"
             << "\n - <max-time>: max running time (in seconds). If 0, the algorithm stops on chosen stop rule."
             << "\n - <instance-file>: instance file"
             << "\n - <scenario>: indicates which scenario uses:"
             << "\n\t+ ST: set the number of depots to K=0.2*n and maximum external cycle size to C=n/K;"
             << "\n\t+ SL: set the number of depots to K=0.2*n and maximum external cycle size to C=2n/K;"
             << "\n\t+ SQ: set the number of depots and maximum external to K=C=sqrt(n);"
             << "\n\t+ Free: set the sizes according the following parameters;"
             << "\n - <K>: number of depots / size of internal cycle. If \"free\" is set, must be >= 3"
             << "\n - <max-cycle-size>: Maximum external cycle size."
             << "\n - <alpha>: discount factor to internal/depot cycle."
             << "\n - <results-files-name>: results files prefix. "
             << "Save the results in \".txt\", the logs in \".log\","
             << "\n\t solutions plots in \".eps\" and the random number generator state in \".rng\""
             << "\n\n ALL PARAMETERS ARE MANDATORY\n"
             << endl;
        return 64;  // BSD usage error code.
    }

    // Loading parameters from command line
    const char* configFile = argv[1];
    unsigned long seed;
    unsigned start_run;
    unsigned runs;
    char stop_rule;
    double stop_arg;
    double max_time;
    const char* instance_filename = argv[8];
    Scenario scenario;
    unsigned K;
    unsigned max_cycle_size;
    double alpha;
    const char* results_name = argv[13];

    try {
         seed = lexical_cast<unsigned long>(argv[2]);
         start_run = lexical_cast<unsigned>(argv[3]);
         runs = lexical_cast<unsigned>(argv[4]);
         stop_rule = StopRule(toupper(argv[5][0]));
         stop_arg = lexical_cast<double>(argv[6]);
         max_time = lexical_cast<double>(argv[7]);
         scenario = Scenario(toupper(argv[9][1]));
         K = lexical_cast<unsigned>(argv[10]);
         max_cycle_size  = lexical_cast<unsigned>(argv[11]);
         alpha = lexical_cast<double>(argv[12]);

         if(stop_rule != GENERATIONS && stop_rule != TARGET && stop_rule != IMPROVEMENT)
             throw std::runtime_error("Incorrect stop rule.");

         if(scenario != ST && scenario != SL &&
            scenario != SQ && scenario != FREE)
             throw std::runtime_error("Incorrect scenario supplied.");
    }
    catch(std::exception& e) {
        cerr << "*** Bad arguments. Verify them!"
             << "\n*** " << e.what() << endl;
        return 64; // BSD usage error code
    }

    // Take the file names to hold the experiment results
    stringstream result_name_stream;
    result_name_stream
        << results_name
        << "_" << scenario
        << "_K" << K
        << "_MS" << setfill ('0') << setw(2) << max_cycle_size
        << "_alpha" << setiosflags(ios::fixed) << setprecision(1) << alpha;

    #ifndef TUNING
    string results_filename(result_name_stream.str());
    results_filename += ".txt";
    string log_filename(result_name_stream.str());
    log_filename += ".log";
    string rng_filename(result_name_stream.str());
    rng_filename += ".rng";
    #else
    string results_filename("/dev/null");
    string log_filename(results_filename);
    string rng_filename(results_filename);
    #endif

    fstream results_file;
    fstream log_file;
    fstream rng_file;

    // Parameters read from configuration file:
    double pe;          // % of elite items in the population
    double pm;          // % of mutants introduced at each generation
    double rhoe;        // prob. of inheriting each allele from elite parent
    unsigned ind_pop;   // number of independent populations
    unsigned MAX_POP;   // maximum size of each population
    unsigned MAX_THR;   // number of threads in parallel decoding
    unsigned J;         // interval at which elite chromosomes are exchanged (0 means no exchange)
    unsigned M;         // number of elite chromosomes to obtain from each population
    unsigned E;         // interval at which the populations are reset (0 means no reset)

    // Loading algorithm parameters from config file (code from rtoso).
    ifstream fin(configFile, std::ios::in);
    if(!fin) {
        cerr << "Cannot open configuration file: " << configFile << endl;
        return 66; // BSD file not found error code
    }
    else {
        fin.exceptions(ifstream::eofbit | ifstream::failbit | ifstream::badbit);
    }

    try {
        string line;
        fin >> pe;      getline(fin, line); // read the parameter and ignore the rest of the line
        fin >> pm;      getline(fin, line);
        fin >> rhoe;    getline(fin, line);
        fin >> ind_pop; getline(fin, line);
        fin >> MAX_POP; getline(fin, line);
        fin >> MAX_THR; getline(fin, line);
        fin >> J;       getline(fin, line);
        fin >> M;       getline(fin, line);
        fin >> E;
        fin.close();
    }
    catch(ifstream::failure& e) {
        cerr << "Failure when reading configfile: " << e.what() << endl;
        fin.close();
        return 65;  // BSD file read error code
    }

    //-----------------------------------------//
    // Tuning
    //-----------------------------------------//
    #ifdef TUNING
    MAX_POP = boost::lexical_cast<unsigned>(argv[13]);
    pe = boost::lexical_cast<double>(argv[14]);
    pm = boost::lexical_cast<double>(argv[15]);
    rhoe = boost::lexical_cast<double>(argv[16]);
    ind_pop = boost::lexical_cast<unsigned>(argv[17]);
    J = boost::lexical_cast<unsigned>(argv[18]);
    M = boost::lexical_cast<unsigned>(argv[19]);
    E = boost::lexical_cast<unsigned>(argv[20]);
    const int kick_type = boost::lexical_cast<int>(argv[21]);
    const int stallcount = boost::lexical_cast<int>(argv[22]);
    #else
    const int kick_type = 2;
    const int stallcount = 1000;
    #endif

    //-----------------------------------------//

    // Log some information.
    log_file.open(log_filename.c_str(), ios::out|ios::app);
    log_file << "\n------------------------------------------------------\n"
             << "> Experiment started at " << posix_time::second_clock::local_time()
             << "\n> Instance " << instance_filename
             << "\n> Scenario: " << scenario;

    switch(scenario) {
    case ST:
        log_file << "\n> K: 0.2 * n (20% of nodes)"
                 << "\n> Maximum cycle size: C=n/K";
        break;

    case SL:
        log_file << "\n> K: 0.2 * n (20% of nodes)"
                 << "\n> Maximum cycle size: C=2n/K";
        break;

    case SQ:
        log_file << "\n> K: sqrt(n)"
                 << "\n> Maximum cycle size: sqrt(n)";
        break;

    case FREE:
        log_file << "\n> K: " << K
                 << "\n> Maximum cycle size: " << max_cycle_size;
        break;
    }

    log_file << "\n> Parameters"
             << "\n> Config file: " << configFile
             << "\n>    + % of Elite: " << pe
             << "\n>    + % of Mutants: " << pm
             << "\n>    + Prob. inheritance (rhoe): " << rhoe
             << "\n>    + # of independent populations: " << ind_pop
             << "\n>    + maximum size of each population: " << MAX_POP
             << "\n>    + # of threads: " << MAX_THR
             << "\n>    + interval of chromossome exchange: " << J
             << "\n>    + # of elite chromossome of each population: " << M
             << "\n>    + reset interval: " << E
             << "\n> Seed: " << seed
             << "\n> Runs: " << runs
             << "\n> Stop Rule: [" << stop_rule << "] "
             << (stop_rule == 'G' ? "Generations -> " :
                (stop_rule == 'T' ? "Target -> " : "Improvement -> "))
             << stop_arg
             << "\n> Max Time: " << max_time
             << endl;
    log_file.close();

    // Try to load the random number generator status.
    MTRand rng(seed);

    rng_file.open(rng_filename.c_str(), ios::in);
    if(rng_file) {
        rng_file >> rng;
        rng_file.close();
    }

    try {
        // Load instance and prepare the Concorde Wrapper.
        K_IMDMTSP_Instance instance;
        ConcordeWrapper concorde(seed, MAX_THR, NULL, NULL, -1, -1, kick_type, stallcount);

        // Adjust some instance data.
        instance.name = instance_filename;
        string::size_type pos = instance.name.rfind('/');
        if(pos != string::npos)
            instance.name = instance.name.substr(pos + 1);

        pos = instance.name.find('.');
        if(pos != string::npos)
            instance.name = instance.name.substr(0, pos);

        // Load the graph from a TSPLIB file.
        concorde.loadTSP(const_cast<char*>(instance_filename), &instance.graph,
                         &instance.dist, &instance.coords, &instance.has_coords);

        // TODO: modify this to degenerated instances.
        // Verify K.
        if(scenario == SQ) {
            K = ceil(sqrt(instance.graph.nodeNum()));
        }
        else
        if(scenario == ST || scenario == SL) {
            K = ceil(0.2 * instance.graph.nodeNum());
        }
        else if((K < 3) || (int)K > (instance.graph.nodeNum() - 1)) {
            stringstream error_msg;
            error_msg << "\n*** Parameter K isn't compatible with the size of the graph:"
                      << "\n*** K: " << K
                      << "\n*** Number of nodes: " << instance.graph.nodeNum();

            throw runtime_error(error_msg.str());
        }

        // Set the parameter K.
        instance.K = K;

        // TODO: modify this to degenerated instances.
        // Verify the maximum size of external cycles.
        if(scenario == SQ) {
            max_cycle_size = ceil(sqrt(instance.graph.nodeNum()));
        }
        else
        if(scenario == ST) {
            max_cycle_size = ceil(instance.graph.nodeNum()/(double)K);
        }
        else
        if(scenario == SL) {
            max_cycle_size = 2 * ceil(instance.graph.nodeNum()/(double)K);
        }
        else
        if(max_cycle_size < ceil(instance.graph.nodeNum()/(double)K) ||
           max_cycle_size > (instance.graph.nodeNum() - K + 1)) {
            stringstream error_msg;
            error_msg << "\n*** Parameter \"max-cycle-size\" isn't compatible with the size of the graph:"
                      << "\n*** max-cycle-size: " << max_cycle_size
                      << "\n*** Number of nodes: " << instance.graph.nodeNum()
                      << "\n*** Minimum required: " << ceil(instance.graph.nodeNum()/(double)K)
                      << "\n*** Maximum allowed: " << (instance.graph.nodeNum() - K + 1);

            throw runtime_error(error_msg.str());
        }

        // Set the maximum size of external cycles.
        instance.max_cycle_size = max_cycle_size;

        // Set the discount factor.
        instance.alpha = alpha;

        // Log some information.
        log_file.open(log_filename.c_str(), ios::out|ios::app);
        log_file << "\n**> Actual K: " << K
                 << "\n**> Actual MS: " << max_cycle_size
                 << "\n\nResults in \"" << results_name << "\" files."
                 << endl;
        log_file.close();

        // Take the population size.
        unsigned population_size = 10 * instance.graph.nodeNum();
        if(population_size > MAX_POP)
            population_size = MAX_POP;

        // Take parameters to control allele generation.
        const unsigned cut_point = instance.K;
        const unsigned left_lb = 1;
        const unsigned left_ub = max_cycle_size;

        concorde.updateData(&instance.graph, &instance.dist);

        K_IMDMTSP_Decoder decoder(instance, rng, concorde);

        // The BRKGA algorithm object.
        BRKGA< K_IMDMTSP_Decoder, MTRand, less >
            algorithm(instance.K + instance.graph.nodeNum(), population_size,
                      pe, pm, rhoe, decoder, rng, ind_pop, MAX_THR, cut_point,
                      left_lb, left_ub, instance.graph.nodeNum() + 1);

        // The timer.
//        Timer timer;
        boost::timer::cpu_timer timer;

        log_file.open(log_filename.c_str(), ios::out|ios::app);

        Population::Chromosome tmp(instance.K + instance.graph.nodeNum(), 0);
        vector< Population::Chromosome > initial_pop;

        // Try to obtain a solution from k-means heuristic.
        chromosomeFromKMeans(&tmp, instance, rng, max_cycle_size, &concorde);

//        IMDMTSP_Solution solution_kmeans(instance);
//        decoder.getSolutionFromChromosome(tmp, &solution_kmeans);
//        solution_kmeans.draw("k_means_no_bug.eps");
        for(unsigned run_number = start_run; run_number < runs; ++run_number) {
            cout << "Run " << (run_number + 1) << " of " << runs << endl;
            // Hold best results
            IMDMTSP_Solution best_solution(instance);
            Population::Chromosome best_chromosome(instance.K + instance.graph.nodeNum(), 0);

            double best_fitness = numeric_limits< double >::max();

            // Optimization info.
            unsigned last_update_iteration = 0;
            boost::timer::cpu_times last_update_time;
            unsigned update_offset = 0;
            unsigned large_offset = 0;
            boost::timer::cpu_times init_population_time;

            // For control the optimization
            unsigned iteration = 1;
            bool run = true;

            log_file << "\n\n> Run " << (run_number+1) << endl;
            log_file << setiosflags(ios::fixed) << setprecision(6);

            // Initialize properly and take elapsed time to do it.
            timer.start();

            initial_pop.clear();
            initial_pop.push_back(tmp);
            algorithm.setInitialPopulation(initial_pop);

            if(run_number != start_run)
                algorithm.reset(true);
            else
                algorithm.initialize();

            init_population_time = timer.elapsed();

            log_file << "\n\n-----------------------------"
                     << "\n>>>> Optimizing...\n"
                     << "Improvement Iteration | Best Current | Best Overall | "
                     << "Current Time (CPU) | Current Time (Wall)" << endl;

            // Reset the timer and do the iterations
            timer.start();
            while(run) {
                // Run one iteration
                algorithm.evolve();

                // Take the best solution
                double fitness = algorithm.getBestFitness();

                if((best_fitness - fitness) > 0.0000001) {
                    timer.stop();
                    last_update_time = timer.elapsed();

                    log_file << "* "
                             << iteration << " | "
                             << fitness << " | "
                             << ((best_fitness < numeric_limits< double >::max())?
                                  best_fitness : -1.0) << " | "
                             << boost::timer::format(last_update_time, 2, "%u") << " | "
                             << boost::timer::format(last_update_time, 2, "%w")
                             << endl;

                    // Take the best chromosome.
                    const Population::Chromosome &chr = algorithm.getBestChromosome();
                    copy(chr.begin(), chr.end(), best_chromosome.begin());

                    best_fitness = fitness;
                    decoder.getSolutionFromChromosome(best_chromosome, &best_solution);

                    update_offset = iteration - last_update_iteration;
                    last_update_iteration = iteration;

                    if(large_offset < update_offset)
                        large_offset = update_offset;

                    timer.resume();
                }

                // Elite-exchange:
                if(ind_pop > 0 && J > 0 && iteration % J == 0) {
                    algorithm.exchangeElite(M);
                    if(stop_rule != TARGET) {
                        log_file << "Exchanged " << M
                                 << " solutions from each population at iteration "<< iteration
                                 << "; best so far: " << algorithm.getBestFitness()
                                 << endl;
                    }
                }

                // Time to reset?
                unsigned iterWithoutImprovement = iteration - last_update_iteration;
                if(E > 0 && iterWithoutImprovement > 0 && iterWithoutImprovement % E == 0) {
                    algorithm.reset();
                    if(stop_rule != TARGET)
                        log_file << "-- Reset chromosomes with random keys after " << E
                                 << " iterations without improvement." << endl;
                }

                // Time to stop?
                switch(stop_rule) {
                    case GENERATIONS:
                        if(iteration == (unsigned)stop_arg)
                            run = false;
                        break;

                    case IMPROVEMENT:
                        if(iterWithoutImprovement >= (unsigned)stop_arg)
                            run = false;
                        break;

                    case TARGET:
                        if(best_fitness < stop_arg + 0.00001)
                            run = false;
                        break;
                }

                if(((timer.elapsed().wall / 1e9) > max_time) && !(max_time - 1e-6 < 0.0))
                    run = false;

                ++iteration;
            } // end of while

            // Take elapsed time to optimization.
            boost::timer::cpu_times elapsed_time(timer.elapsed());

            // Take the best solution.
            decoder.getSolutionFromChromosome(best_chromosome, &best_solution);

            stringstream result_message;

            result_message
                << instance.name << " & "
                << scenario << " & "
                << K << " & "
                << max_cycle_size << " & "
                << setiosflags(ios::fixed) << setprecision(2)
                << alpha << " & "
                << setiosflags(ios::fixed) << setprecision(2)
                << best_solution.value << " & "
                << setprecision(2)
                << --iteration << " & "            // only to correct the iteration number
                << last_update_iteration << " & "
                << boost::timer::format(last_update_time, 2, "%u") << " & "
                << boost::timer::format(last_update_time, 2, "%w") << " & "
                << update_offset << " & "
                << large_offset << " & "
                << boost::timer::format(init_population_time, 2, "%w") << " & "
                << boost::timer::format(elapsed_time, 2, "%u") << " & "
                << boost::timer::format(elapsed_time, 2, "%w");

            log_file << "\nBest Solution:" << best_solution
                     << "\n\nInstance & Scenario & K & Max cycle size & Discount Factor & Cost & Iterations & Last Update Iteration & LUTCpu & LUTWall & "
                     << "Update Offset & Large Offset & Init. Time (s) & TotalTimeCpu & TotalTimeWall \n"
                     << result_message.str()
                     << endl;
            log_file.flush();

            // The results on the file.
            results_file.open(results_filename.c_str(), ios::out|ios::app);
            if(results_file) {
                results_file << result_message.str() << endl;
                results_file.close();
            }
            else {
                stringstream error_msg;
                error_msg << "\n*** Cannot open results file: "
                          << results_filename;
                throw runtime_error(error_msg.str());
            }

            #ifndef TUNING
            stringstream plot_name;
            plot_name
                << result_name_stream.str()
                << "_run" << setfill ('0') << setw(2) << run_number;

            best_solution.draw(plot_name.str().c_str());
            #else
            cout << best_solution.value;
            #endif

            // Save random number generator status.
            #ifndef FULLDEBUG
            rng_file.open(rng_filename.c_str(), ios::out);
            if(rng_file) {
                rng_file << rng;
            }
            else {
                stringstream error_msg;
                error_msg << "\n*** Cannot open random number generator file: "
                          << rng_filename;
                throw runtime_error(error_msg.str());
            }
            rng_file.close();
            #endif
        } // end of for

        log_file.close();
    }
    // If somethig goes wrong....
    catch(std::exception& e) {
        stringstream message;
        message << "\n***********************************************************"
                << "\n*** Exception Occured: " << e.what()
                << "\n***********************************************************"
                << endl;

        log_file.open(log_filename.c_str(), ios::out|ios::app);
        log_file << message.str();
        log_file.close();

        cerr << message.str();
        return 70; // BSD software internal error code
    }

    // Log final info.
    log_file.open(log_filename.c_str(), ios::out|ios::app);
    log_file << "\n\n> Experiment finished at " << posix_time::second_clock::local_time()
             << "\n------------------------------------------------------" << endl;
    log_file.close();

     return 0;
}

