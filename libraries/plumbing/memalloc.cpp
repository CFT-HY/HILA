/// Memory allocator -- gives back aligned memory, if ALIGN defined
/// This is a hacky method to insert automatic check to success of
/// memalloc, so that we get the file name and line number where the
/// the failure happened:
///

#include "plumbing/memalloc.h"

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
