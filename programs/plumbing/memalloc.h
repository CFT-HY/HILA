/// Memory allocator -- gives back aligned memory, if ALIGN defined
/// This is a hacky method to insert automatic check to success of
/// memalloc, so that we get the file name and line number where the
/// the failure happened:
/// we use just macro with comma operators and assertion instead of separate function
/// points where the allocation failed.

static void * _memalloc__ptr_;   // this one pollutes namespace, but need to define one

#ifndef ALIGN
// _ptr_ is here a tmp variable
#define memalloc(size) ( _memalloc__ptr_ = std::malloc( size),                    \
    assert( _memalloc__ptr_ != nullptr && "memalloc: memory allocation failure"), \
    _memalloc__ptr_ )

#else
#define memalloc(size) ( void * _ptr_,                            \
    int e = std::posix_memalign( &_ptr_, (std::size_t)32, size),  \
    assert( e == 0 && "memalloc: memory allocation failure"),     \
    _ptr_ )

#endif

