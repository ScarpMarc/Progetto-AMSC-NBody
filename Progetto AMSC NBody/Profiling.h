#pragma once

#include <filesystem>

#include "Globals.h"

extern std::filesystem::path profiling_folder;

extern std::filesystem::path profiling_file_name;

void save_profiler_data_text_file(const std::filesystem::path& path_to_file);

