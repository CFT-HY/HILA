#ifndef GPUCUB_H_
#define GPUCUB_H_

// hilapp does not read these, because it does not understand cuda/hip etc.
#if !defined(HILAPP)
#if defined(CUDA)
#include <cub/cub.cuh>
namespace gpucub = cub;
#endif

#if defined(HIP)
#include <hipcub/hipcub.hpp>
namespace gpucub = hipcub;
#endif
#endif // HILAPP

#endif
