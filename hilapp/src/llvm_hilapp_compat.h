#ifndef LLVM_HILAPP_COMPAT_H
#define LLVM_HILAPP_COMPAT_H

/// Common compatibility macros for different llvm versions

// Define compatibility macro to handle llvm version changing `Expr::isCXX11ConstantExpr` signature
#if LLVM_VERSION_MAJOR > 21 || (LLVM_VERSION_MAJOR == 21 && LLVM_VERSION_MINOR >= 1)
#define HILAPP_ARGS_isCXX11ConstantExpr(a, b, c) (a), (b)
#else
#define HILAPP_ARGS_isCXX11ConstantExpr(a, b, c) (a), (b), (c)
#endif

#endif // LLVM_HILAPP_COMPAT_H