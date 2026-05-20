#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <sstream>

#ifdef USE_SSE
#include "sse/blake2.h"  // Include the BLAKE2 reference header
#else
#include "ref/blake2.h"
#endif

std::string to_hex(const uint8_t *digest, size_t digest_size)
{
    std::ostringstream oss;
    oss << std::hex << std::setfill('0');
    for (size_t i = 0; i < digest_size; i++)
    {
        oss << std::setw(2) << static_cast<unsigned int>(digest[i]);
    }
    return oss.str();
}

std::string computeBlake2bHash(const void *data, size_t data_size, size_t digest_size)
{
  if (digest_size < 1 || digest_size > BLAKE2B_OUTBYTES)
  {
      throw std::runtime_error("Invalid digest size for BLAKE2b");
  }

  std::vector<uint8_t> out(digest_size);

  if (blake2b(out.data(), digest_size, data, data_size, nullptr, 0) != 0)
  {
      throw std::runtime_error("Error computing BLAKE2b hash");
  }

  return to_hex(out.data(), digest_size);
}
