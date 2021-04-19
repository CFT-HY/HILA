#include "plumbing/memalloc.h"


/// Memory allocator -- gives back aligned memory, if ALIGN defined
/// This is a hacky method to insert automatic check to success of
/// memalloc, so that we get the file name and line number where the
/// the failure happened:


void *memalloc(std::size_t size, const char *filename, const unsigned line) {

#ifndef ALIGNED_MEMALLOC

    void *p;

    if ((p = std::malloc(size)) == nullptr) {
        if (filename != nullptr) {
            hila::output << " *** memalloc failure in file " << filename << " at line "
                         << line << '\n';
        } else {
            hila::output << " *** memalloc failure\n";
        }
        hila::output << "     requested allocation size " << size << " bytes\n"
                     << " *********************************" << std::endl;

        exit(1);
    }
    return p;

#else

    void *p;
    // align to 32 bytes (parameter?)
    int e = posix_memalign(&p, (std::size_t)32, size);
    if (e != 0) {
        if (filename != nullptr) {
            hila::output << " *** memalloc failure in file " << filename << " at line "
                         << line << '\n';
        } else {
            hila::output << " *** memalloc failure\n";
        }
        hila::output << "     requested allocation size " << size << " bytes\n"
                     << "     posix_memalign() error code " << e
                     << " *********************************" << std::endl;
        exit(1);
    }
    return p;

#endif
}

// version without the filename
void *memalloc(std::size_t size) { return memalloc(size, nullptr, 0); }


/// d_malloc allocates memory from "device", either from 

void *d_malloc(std::size_t size) { 

#if !defined(CUDA) && !defined(HIP)

    return memalloc(size);

#else

    void * p;
    gpuMalloc(&p, size);
    check_device_error("d_malloc error");
    return p;

#endif

}

void d_free(void *dptr) {

#if !defined(CUDA) && !defined(HIP)

    if (dptr != nullptr) free(dptr);

#else

    gpuFree(dptr);
    check_device_error("CUDA/HIP d_free");

#endif

}
