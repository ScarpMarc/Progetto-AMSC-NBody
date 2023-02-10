#include "Profiling.h"
#include "Constants.h"
#include <omp.h>

#include <fstream>

using namespace std;

void save_profiler_data_text_file(const std::filesystem::path& path_to_file)
{
	ofstream file_out(path_to_file, ios_base::trunc);
	file_out << "NBody profiling" << endl;
#if defined DEBUG || defined _DEBUG
	file_out << "!!! DEBUG MODE !!!";
#endif // DEBUG || defined _DEBUG
	file_out << "--- SIMULATION INFO ---" << endl;
	file_out << "Problem dimension: " << DIM << "D" << endl;
	file_out << "Particle amount: " << std::setw(15) << total_particles << endl;
	file_out << "Number of iterations: " << std::setw(15) << max_ticks << endl;
	file_out << "Simulation time resolution: " << std::setw(15) << ticks_per_second << " ticks per second" << endl;
	file_out << "Simulation elapsed time: " << std::setw(15) << max_ticks * ticks_per_second << " seconds" << endl;
	file_out << "Graphics refresh time: " << std::setw(15) << screen_refresh_millis << " milliseconds" << endl;
	file_out << "Screen resolution: " << std::setw(5) << screenResX << "x" << std::setw(5) << screenResY << endl;
	file_out << endl;
	file_out << "--- TIMING AND THREADS INFO ---" << endl;
	file_out << "Total duration: " << std::setw(15) << total_sim_duration / 1000 << " ms" << endl;
	file_out << "Iteration mean duration: " << std::setw(15) << matrixComp_mean_duration << " us" << endl;

	file_out << "Garbage collecting mean duration: " << std::setw(15) << garbage_collecting_mean_duration << " us" << endl;
	file_out << "Cluster update mean duration: " << std::setw(15) << update_mean_duration << " us" << endl;
	file_out << "Clusters eliminated on average: " << mean_eliminated_clusters << endl;
    
	file_out.close();
}
