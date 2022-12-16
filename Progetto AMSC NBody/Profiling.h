#pragma once

#include <filesystem>

#include "Globals.h"

std::filesystem::path profiling_folder = "/";

std::filesystem::path profiling_file_name = "Profiler_";

void save_profiler_data_text_file(const std::filesystem::path& path_to_file);

