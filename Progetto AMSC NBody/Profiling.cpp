#include "Profiling.h"
#include "Globals.h"
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
	file_out << "Total duration: " << std::setw(15) << total_sim_duration << endl;
#ifdef ADVANCED_PROFILING
	file_out << "Force computation mean duration: " << std::setw(15) << forceComp_mean_durations_per_tick << " ns" << endl;
	file_out << "Position computation mean duration: " << std::setw(15) << posComp_mean_durations_per_tick << " ns" << endl;
	file_out << "Matrix computation mean duration: " << std::setw(15) << matrixComp_mean_duration << " ns" << endl;
#endif
	file_out.close();
}
