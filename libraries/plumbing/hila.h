#ifndef HILA_H_
#define HILA_H_

///////////////////////////////////////////////////////////////////////
/// Catch-(almost)all include to get in most of the hila-system .h -files



#include "plumbing/defs.h"
#include "datatypes/cmplx.h"
#include "datatypes/matrix.h"
#include "datatypes/element_by_element.h"
#include "datatypes/sun_matrix.h"
#include "datatypes/u1.h"
#include "datatypes/su2.h"
#include "datatypes/extended.h"

#include "plumbing/globals.h"
#include "plumbing/coordinates.h"
#include "plumbing/lattice.h"
#include "plumbing/site_index.h"
#include "plumbing/field.h"
#include "plumbing/field_io.h"
#include "plumbing/reduction.h"
#include "plumbing/reductionvector.h"
#include "plumbing/site_select.h"

#if __has_include("hila_signatures.h")
#include "hila_signatures.h"
#endif

//#if defined(CUDA) || defined(HIP)
//#include "plumbing/backend_gpu/gpu_reduction.h"
//#endif

#include "plumbing/gaugefield.h"
#include "plumbing/input.h"
#include "plumbing/cmdline.h"
#include "plumbing/fft.h"

#include "plumbing/spectraldensity.h"

#if defined(OPENMP) && !defined(HILAPP)
#include <omp.h>
#endif

#include "plumbing/random.h"

#include "plumbing/shuffle.h"

#endif
