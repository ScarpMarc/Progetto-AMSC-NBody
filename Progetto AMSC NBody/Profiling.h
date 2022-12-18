#pragma once

#include <filesystem>
#include <string>

#include "Globals.h"

extern std::string profiling_folder;

extern std::string profiling_file_name;

void save_profiler_data_text_file(const std::filesystem::path& path_to_file);

