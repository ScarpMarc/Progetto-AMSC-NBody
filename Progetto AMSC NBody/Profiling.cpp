#include "Profiling.h"
#include "Globals.h"
#include "Constants.h"

#include <fstream>

using namespace std;

void save_profiler_data_text_file(const std::filesystem::path& path_to_file)
{
	ofstream file_out(path_to_file);
	file_out << "NBody profiling" << endl;
	file_out << "Total duration: " << std::setw(15) << total_sim_duration << endl;

#ifdef ADVANCED_PROFILING
	file_out << "Force computation mean duration: " << std::setw(15) << forceComp_mean_durations_per_tick << endl;
	file_out << "Position computation mean duration: " << std::setw(15) << posComp_mean_durations_per_tick << endl;
	file_out << "Matrix computation mean duration: " << std::setw(15) << matrixComp_mean_duration << endl;
#endif
}
