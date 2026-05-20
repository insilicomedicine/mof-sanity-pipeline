#pragma once

#include <string>

std::string computeBlake2bHash(const void *data, size_t data_size, size_t digest_size);
