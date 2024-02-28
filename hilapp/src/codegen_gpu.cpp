//------------------------------------------------------------------------------
// Generate transformed
// hardware-dependent "kernels".
//
// Uses Clang RecursiveASTVisitor and Rewriter
// interfaces
//
//------------------------------------------------------------------------------
#include <sstream>
#include <string>

#include "hilapp.h"
#include "toplevelvisitor.h"
#include "stringops.h"
#include "clang/AST/ASTLambda.h"

// define max size of an array passed as a parameter to kernels
#define MAX_PARAM_ARRAY_SIZE 40

// NOTE: If you define ALT_VECTOR_REDUCTION you need to define
// it also in gpu_templated_ops.h !

// #define ALT_VECTOR_REDUCTION

extern std::string looping_var;
extern std::string parity_name;


// write __host__ __device__ to function decl
void GeneralVisitor::gpu_loop_function_marker(FunctionDecl *fd) {

    // SourceLocation sl = fd->getSourceRange().getBegin();
    SourceLocation sl = fd->getInnerLocStart();
    srcBuf *sb = get_file_srcBuf(sl);
    if (sb == nullptr) {
        // it's a system file -- should we do something?
        return;
    }

    if(CXXMethodDecl *MD = dyn_cast<CXXMethodDecl>(fd)){
        if(isLambdaCallOperator(MD)){
            CXXRecordDecl *RD = MD->getParent();
            if(RD->isLambda()){
                // reportDiag(DiagnosticsEngine::Level::Remark,fd->getInnerLocStart(),"fd begin\n");
                // reportDiag(DiagnosticsEngine::Level::Remark,RD->getInnerLocStart(),"RD begin\n");
                // fd->getInnerLocStart points to here:------------------------------------------------v
                //                        auto lambda = [capture, this, and, that]   (int i, int j) -> return_type { ... }
                // RD->getInnerLocStart points to here:-^                          ^
                //                                                                 |
                // __device__ __host__ should be inserted in here -----------------|
                // for now ignore, anyway not standard, nvcc --expt-extended-lambda
                return;
            }
        }
    }

    if (!sb->is_edited(sl))
        sb->insert(sl, "__device__ __host__ ", true, true);
}


// Add the __host__ __device__ keywords to functions called a loop
void GeneralVisitor::handle_loop_function_gpu(call_info_struct &ci) {

    if (ci.is_defaulted)
        return; // cuda can take care of these

    gpu_loop_function_marker(ci.funcdecl);
    FunctionDecl *proto;
    if (find_prototype(ci.funcdecl, proto)) {
        gpu_loop_function_marker(proto);
    }
}

void GeneralVisitor::handle_loop_constructor_gpu(call_info_struct &ci) {

    if (ci.is_defaulted) {
        return;
    }

    SourceLocation sl = ci.ctordecl->getSourceRange().getBegin();
    srcBuf *sb = get_file_srcBuf(sl);
    if (sb == nullptr) {
        // it's a system file -- should we do something?
        return;
    }
    sb->insert(sl, "__device__ __host__ ", true, true);
}

/////////////////////////////////////////////////////////////////////////////////////

/// Help routine to write (part of) a name for a kernel
std::string TopLevelVisitor::make_kernel_name() {
    return kernel_name_prefix +
           clean_name(global.currentFunctionDecl->getNameInfo().getName().getAsString()) + "_" +
           std::to_string(TheRewriter.getSourceMgr().
                          // getSpellingLineNumber(global.location.loop));
                          getFileOffset(global.location.loop));
}

std::string TopLevelVisitor::generate_code_gpu(Stmt *S, bool semicolon_at_end, srcBuf &loopBuf,
                                               bool generate_wait_loops) {

    // "Code" is inserted at the location of the loop statement
    // and the kernel is build in "kernel"
    std::stringstream code, kernel;
    // const std::string t = loopBuf.dump();

    // indexing variable
    extern std::string looping_var;


    // Get kernel name - use line number or file offset (must be deterministic)
    std::string kernel_name = TopLevelVisitor::make_kernel_name();

    // Wait for the communication to finish
    for (field_info &l : field_info_list) {
        // If neighbour references exist, communicate them
        if (!l.is_loop_local_dir) {
            for (dir_ptr &d : l.dir_list)
                if (d.count > 0) {
                    code << l.new_name << ".wait_gather(" << d.direxpr_s << ", "
                         << loop_info.parity_str << ");\n";
                }
        } else {
            code << "for (Direction _HILAdir_ = (Direction)0; _HILAdir_ < NDIRS; "
                    "++_HILAdir_) {\n"
                 << l.new_name << ".wait_gather(_HILAdir_," << loop_info.parity_str << ");\n}\n";
        }
    }

    // Set loop lattice
    // if (field_info_list.size() > 0) {
    //     std::string fieldname = field_info_list.front().old_name;
    //     code << "const lattice_struct & loop_lattice = " << fieldname << ".fs->lattice;\n";
    // } else {
    // now no fields in loop - default lattice
    code << "const lattice_struct & loop_lattice = lattice;\n";


    kernel << "\n\n//------------------ start kernel " << kernel_name << "--------------------\n";


    // Figure out the RNG operation status
    // For RNG loops we optionally create loop with lower number of active threads
    // enabling the use of smaller rng generator footprint

    static bool rng_threads_checked = false;
    static int rng_thread_block_number = 0;

    bool use_thread_blocks = false;
    int thread_block_number = 0;

    if (loop_info.contains_random) {

        if (!rng_threads_checked) {
            // now check if the rng thread number is present
            rng_threads_checked = true;
            std::string rng_str;

            if (is_macro_defined("GPU_RNG_THREAD_BLOCKS", &rng_str)) {
                rng_thread_block_number = std::atoi(rng_str.c_str());
            }
        }

        if (rng_thread_block_number > 0) {
            // loop has random and rng thread blocks in use and number known
            use_thread_blocks = true;
            thread_block_number = rng_thread_block_number;
        }

        if (rng_thread_block_number < 0) {
            reportDiag(DiagnosticsEngine::Level::Warning, loop_info.range.getBegin(),
                       "loop contains random number generator calls but GPU random numbers are "
                       "disabled (GPU_RNG_THREAD_BLOCKS < 0).  This may crash the program.");
        }
    }


    // if we have small arrays, encapsulate them in struct
    // struct has to be defined before the kernel call

    // keep track if we have rv here
    bool loop_has_reductionvector_blocks = false;

    for (array_ref &ar : array_ref_list) {
        // llvm::errs() << "ARRAY " << ar.name << " REF TYPE IS " << ar.type << '\n';

        if (ar.type != array_ref::REPLACE && ar.type != array_ref::REDUCTION) {

            if (ar.size > 0 && ar.size <= MAX_PARAM_ARRAY_SIZE) {

                // copy small arrays directly, create new type
                kernel << "// create encapsulating struct for '" << ar.name << "'\n";
                ar.new_name = var_name_prefix + clean_name(ar.name);
                // use here offset to give unique type
                ar.wrapper_type =
                    "struct " + type_name_prefix + clean_name(ar.name) +
                    std::to_string(TheRewriter.getSourceMgr().getFileOffset(global.location.loop));
                kernel << ar.wrapper_type << " {\n";

                // give the type the correct dimensions
                kernel << ar.element_type << " c";
                for (auto d : ar.dimensions) {
                    kernel << '[' << d << ']';
                }
                kernel << ";\n};\n\n";


            } else {

                // larger array or vector, copy it directly -- allocate

                if (ar.dimensions.size() > 1) {
                    if (ar.refs.size() > 0) {
                        reportDiag(DiagnosticsEngine::Level::Error,
                                   ar.refs[0].E->getSourceRange().getBegin(),
                                   "multidimensional arrays with indices which are not loop "
                                   "constant and size larger than "
                                   "'%0' not implemented in GPU code",
                                   std::to_string(MAX_PARAM_ARRAY_SIZE).c_str());
                    } else {
                        llvm::errs() << "Internal bug - multidimensional array in loop\n";
                        exit(1);
                    }
                }

                ar.new_name = var_name_prefix + clean_name(ar.name);

                code << "// copy array/vector '" << ar.name << "' to device\n";

                code << ar.element_type << " * " << ar.new_name << ";\n";
                code << "gpuMalloc( & " << ar.new_name << ", " << ar.size_expr << " * sizeof("
                     << ar.element_type << ") );\n";

                code << "gpuMemcpy(" << ar.new_name << ", (char *)" << ar.data_ptr << ", "
                     << ar.size_expr << " * sizeof(" << ar.element_type << "), "
                     << "gpuMemcpyHostToDevice);\n\n";
            }

        } else if (ar.type == array_ref::REDUCTION) {

            // Now there is reductionvector - do we have special reduction blocks or use atomic?

            std::string rv_block_str;
            if (is_macro_defined("GPU_VECTOR_REDUCTION_THREAD_BLOCKS", &rv_block_str)) {
                int blocks = std::atoi(rv_block_str.c_str());
                if (blocks > 0) {
                    loop_has_reductionvector_blocks = true;
                    use_thread_blocks = true;

                    // we have to take smaller of block numbers - RNG works with smaller number,
                    // but not larger.
                    if (thread_block_number > 0) {
                        thread_block_number = std::min(thread_block_number, blocks);
                    } else {
                        thread_block_number = blocks;
                    }
                }
            }

            ar.new_name = "r_" + var_name_prefix + clean_name(ar.name);

            code << ar.element_type << " * " << ar.new_name << ";\n";

            std::stringstream array_size; // keep size expression in string

            if (loop_has_reductionvector_blocks) {
                code << "// Create reduction array with " << thread_block_number
                     << " * N_threads parallel reductions\n";

                array_size << ar.size_expr << " * N_threads * " << thread_block_number;

            } else {
                code << "// Create reduction array - using atomicAdd for accumulation\n";

                array_size << ar.size_expr;
            }

            code << "gpuMalloc( & " << ar.new_name << ", " << array_size.str() << " * sizeof("
                 << ar.element_type << "));\n";

            if (ar.reduction_type == reduction::SUM) {
                code << "gpu_set_zero(" << ar.new_name << ", " << array_size.str() << ");\n";
            }

            if (ar.reduction_type == reduction::PRODUCT) {
                code << "gpu_set_value(" << ar.new_name << ", 1, " << array_size.str() << ");\n";
            }

            code << "check_device_error(\"allocate_reduction\");\n";
        }
    }


    // Generate the function definition and call
    // "inline" makes cuda complain, but it is needed to avoid multiple definition error
    // use "static" instead??
    // Switch to static here for now
    // Add __launch_bounds__ directive here
    kernel << "static __global__ void __launch_bounds__(N_threads) " << kernel_name
           << "( backend_lattice_struct d_lattice";
    code << "backend_lattice_struct lattice_info = *(lattice.backend_lattice);\n";
    code << "lattice_info.loop_begin = lattice.loop_begin(" << loop_info.parity_str << ");\n";
    code << "lattice_info.loop_end = lattice.loop_end(" << loop_info.parity_str << ");\n";

    code << "int N_blocks = (lattice_info.loop_end - lattice_info.loop_begin + "
            "N_threads - 1)/N_threads;\n";

    if (use_thread_blocks) {
        code << "if (N_blocks > " << thread_block_number << ") N_blocks = " << thread_block_number
             << ";\n";
    }

    // Check for reductions and allocate device memory for each var
    for (reduction_expr &r : reduction_list) {

        code << r.type << " * dev_" << r.reduction_name << ";\n";
        code << "gpuMalloc( & dev_" << r.reduction_name << ","
             << "sizeof(" << r.type << ") * N_blocks );\n";
        if (r.reduction_type == reduction::SUM) {
            code << "gpu_set_zero(dev_" << r.reduction_name << ", N_blocks);\n";
        }
        if (r.reduction_type == reduction::PRODUCT) {
            code << "gpu_set_value(dev_" << r.reduction_name << ", 1, N_blocks);\n";
        }
    }


    // and for selections
    for (selection_info &s : selection_info_list) {
        s.maskname = s.new_name + "_mask";
        s.sitename = s.new_name + "_sites";
        s.valname = s.new_name + "_values";

        if (s.previous_selection == nullptr) {
            // define mask, coord and possible value holder
            code << "char * " << s.maskname << ";\n";
            code << "gpuMalloc( & " << s.maskname << ", lattice.mynode.volume() );\n";
            code << "SiteIndex * " << s.sitename << ";\n";
            code << "gpuMalloc( & " << s.sitename
                 << ", lattice.mynode.volume()*sizeof(SiteIndex) );\n";
            if (s.assign_expr != nullptr) {
                code << s.val_type << " * " << s.valname << ";\n";
                code << "gpuMalloc( & " << s.valname << ", lattice.mynode.volume()*sizeof("
                     << s.val_type << ") );\n";
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // start kernel building
    ////////////////////////////////////////////////////////////////////////////////////////////////

    if (target.cuda) {
        code << kernel_name << "<<< N_blocks, N_threads >>>( lattice_info";
    } else if (target.hip) {
        code << "hipLaunchKernelGGL(" << kernel_name
             << ", dim3(N_blocks), dim3(N_threads), 0, 0, lattice_info";
    } else {
        llvm::errs() << "Internal bug - unknown kernelized target\n";
        exit(1);
    }

    // print field call list
    int i = 0;
    for (field_info &l : field_info_list) {

        kernel << ", ";
        code << ", ";

        if (!l.is_written)
            kernel << "const ";
        kernel << "field_storage" << l.type_template << " " << l.new_name;
        code << l.new_name + ".fs->payload";
        i++;
    }

    i = 0;
    // and non-field vars
    for (var_info &vi : var_info_list) {

        if (vi.is_raw) {
            // pass the raw ptr as is

            kernel << ", " << vi.type << ' ' << vi.name;
            code << ", " << vi.name;

        } else if (!vi.is_loop_local) {

            vi.new_name = "HILA_var_" + std::to_string(i++) + "_";


            if (vi.reduction_type != reduction::NONE) {
                // // Generate a temporary array for the reduction
                // kernel << ", " << vi.type << " * " << vi.new_name;
                // code << ", dev_" << vi.reduction_name;
            } else {
                kernel << ", const " << vi.type << " " << vi.new_name;
                code << ", " << vi.name;
                // Replace references in the loop body
                for (var_ref &vr : vi.refs) {
                    loopBuf.replace(vr.ref, vi.new_name);
                }
            }
        }
    }

    // write reduction parameters
    i = 0;
    for (reduction_expr &r : reduction_list) {
        r.loop_name = name_prefix + "reduction_" + std::to_string(i++) + "_";
        kernel << ", " << r.type << " * " << r.loop_name;
        code << ", dev_" << r.reduction_name;
    }

    // Then loop constant expressions upgraded
    i = 0;
    for (loop_const_expr_ref lcer : loop_const_expr_ref_list) {
        if (lcer.reduction_type == reduction::NONE) {
            // reductions handled separately
            // replacement was done on the general codegen level
            // loop var is lcer.new_name

            kernel << ", const " << lcer.type << ' ' << lcer.new_name;
            code << ", " << lcer.expression;
        }
    }

    // pass the selection fields to kernel
    for (selection_info &s : selection_info_list) {
        std::string looping("[" + looping_var + "]");

        if (s.previous_selection == nullptr) {
            code << ", " << s.maskname << ", " << s.sitename;
            kernel << ", char * " << s.maskname << ", SiteIndex * " << s.sitename;
        }

        loopBuf.remove_semicolon_after(s.MCE);

        std::string repl_string;
        repl_string = "{\n " + s.maskname + looping + " = 1;\n";
        repl_string +=
            s.sitename + looping + " = SiteIndex(d_lattice.coordinates(" + looping_var + "));\n";

        if (s.assign_expr == nullptr) {

            // replace references
            repl_string += "}\n";
            loopBuf.replace(s.MCE, repl_string);

        } else {
            if (s.previous_selection == nullptr) {
                code << ", " << s.valname;
                kernel << ", " << s.val_type << " * " << s.valname;
            }

            repl_string += s.valname + looping + " = (";

            // get range in a.select(X,value)  between beginning of a and value
            SourceRange r(s.MCE->getSourceRange().getBegin(),
                          s.assign_expr->getSourceRange().getBegin().getLocWithOffset(-1));
            loopBuf.replace(r, repl_string);
            int loc = loopBuf.find_original(r.getEnd(), ')');
            loopBuf.insert(loc + 1, ";\n}\n");
        }
    }


    // In kernelized code we need to handle array expressions as well
    for (array_ref &ar : array_ref_list) {
        if (ar.type == array_ref::REPLACE) {

            // In this case we replace array expression with a new variable
            // Rename the expression
            ar.new_name = var_name_prefix + std::to_string(i) + "_";
            i++;
            kernel << ", const " << ar.element_type << " " << ar.new_name;
            code << ", " << get_stmt_str(ar.refs[0].E);

            loopBuf.replace(ar.refs[0].E, ar.new_name);

        } else if (ar.size > 0 && ar.size <= MAX_PARAM_ARRAY_SIZE) {

            // Now we pass the full array to the kernel
            // ar was set already above
            // Cast the data directly
            code << ", *(" << ar.wrapper_type << "*)(void *)" << ar.data_ptr;

            kernel << ", const " << ar.wrapper_type << ' ' << ar.new_name;

            for (bracket_ref_t &br : ar.refs) {
                loopBuf.replace(br.BASE, ar.new_name + ".c");
            }

        } else if (ar.type != array_ref::REDUCTION) {

            // Now pass the array ptr
            code << ", " << ar.new_name;
            kernel << ", const " << ar.element_type << " * RESTRICT " << ar.new_name;

            for (bracket_ref_t &br : ar.refs) {
                loopBuf.replace(br.BASE, ar.new_name);
            }

        } else {

            // Finally, we have reductionvector here
            // substute here a[i] += b;
            // if not reductionvector blocks, substitute with atomicAdd(&a[i],b);
            // with blocks, each thread adds to its own histogram

            // pass also size here
            code << ", " << ar.new_name;
            kernel << ", " << ar.element_type << " * RESTRICT " << ar.new_name;

            // kernel variable name for size of array
            std::string array_size_varname;
            if (loop_has_reductionvector_blocks) {
                // in this case need to pass the size of the vector (original!) to the kernel
                array_size_varname = ar.new_name + "_size";

                code << ", " << ar.size_expr;
                kernel << ", const int " << array_size_varname;
            }


            for (bracket_ref_t &br : ar.refs) {
                // change the name
                loopBuf.replace(br.BASE, ar.new_name);

                SourceLocation oploc, beginloc, endloc;
                beginloc = br.assign_stmt->getSourceRange().getBegin();
                endloc = getSourceLocationAtEndOfRange(br.assign_stmt->getSourceRange());

                CXXOperatorCallExpr *OP = dyn_cast<CXXOperatorCallExpr>(br.assign_stmt);
                if (OP && OP->isAssignmentOp()) {
                    oploc = OP->getOperatorLoc();

                } else if (CompoundAssignOperator *CAO =
                               dyn_cast<CompoundAssignOperator>(br.assign_stmt)) {
                    oploc = CAO->getOperatorLoc();

                } else {
                    llvm::errs() << "hilapp internal error: vector reduction op in "
                                    "codegen_gpu\n";
                    exit(1);
                }

                if (!loop_has_reductionvector_blocks) {
                    // in this case do atomic ops

                    loopBuf.replace(SourceRange(oploc, oploc.getLocWithOffset(1)), ",");

                    if (ar.reduction_type == reduction::SUM) {

                        // Cuda has a bug where double atomicAdd is not defined for
                        // capability < 6.0, but you nevertheless cannot use the name.
                        // use slightly modified name
                        if (ar.element_type == "double")
                            loopBuf.insert(beginloc, "atomic_Add(&", true);
                        else
                            loopBuf.insert(beginloc, "atomicAdd(&", true);

                    } else {
                        loopBuf.insert(beginloc, "atomicMultiply(&", true);
                    }

                    loopBuf.insert_after(endloc, ")", false);

                } else {
                    // Now each thread writes to its own array section.
                    // give thread offset to the index
                    // Here we change expr of type
                    //    s[ind] += val
                    // to
                    //    new_name[static_cast<int>(ind)*N_threads_total + thread_idx] += val
                    // static cast is added if the type of the index is not of int type

                    // Now have to be careful not to insert on the index variable "area",
                    // go in front of it

                    SourceLocation l =
                        br.Idx.at(0)->getSourceRange().getBegin().getLocWithOffset(-1);
                    // get the char there, most likely '['
                    std::string replstr = loopBuf.get(l, 1);

                    if (br.Idx.at(0)->getType().getTypePtr()->isIntegerType()) {
                        replstr.append("( ");
                    } else {
                        replstr.append("static_cast<int>( ");
                    }
                    int idx = loopBuf.get_index(l);

                    loopBuf.replace(idx, idx, replstr);

                    // _HILA_thread_id is defined in kernel start below
                    l = getSourceLocationAtEndOfRange(br.Idx.at(0)->getSourceRange())
                            .getLocWithOffset(1);

#ifndef ALT_VECTOR_REDUCTION
                    // Normal simpler way to accumulate vectorreductions
                    loopBuf.insert(l, " ) + _HILA_thread_id * " + array_size_varname);

#else
                    // ALT of above is below
                    std::stringstream ss;
                    ss << " )*(N_threads * " << thread_block_number << ") + _HILA_thread_id";
                    loopBuf.insert(l, ss.str());
#endif
                }
            }
        }
    }

    // // change the references to field expressions in the body
    // for ( field_ref & le : field_ref_list ) {
    //   //loopBuf.replace( le.nameExpr, le.info->loop_ref_name );
    //   if (le.dirExpr != nullptr) {
    //     loopBuf.replace(le.fullExpr, le.info->loop_ref_name+"_"+le.dirname);
    //   } else {
    //     loopBuf.replace(le.fullExpr, le.info->loop_ref_name);
    //   }
    // }

    // Handle calls to special in-loop functions
    for (special_function_call &sfc : special_function_call_list) {
        std::string repl = sfc.replace_expression; // comes with ( now
        if (sfc.add_loop_var) {
            repl += looping_var;
            if (sfc.argsExpr != nullptr)
                repl += ',';
            if (sfc.args_string.size() > 0)
                repl += ", " + sfc.args_string;
        }
        loopBuf.replace(sfc.replace_range, repl);
    }

    // Begin the function
    kernel << ")\n{\n";

    // kernel << "backend_lattice_struct *loop_lattice = &d_lattice; \n";

    // Declare the shared reduction variable inside loop
    for (reduction_expr &r : reduction_list) {
        // Generate a temporary array for the reduction
        kernel << "__shared__ " << r.type << " " << r.loop_name << "sh[N_threads];\n";
        kernel << r.type << " " << r.loop_name << "sum;\n";
        // Initialize only the local element
        if (r.reduction_type == reduction::SUM) {
            // OLD method - using operator=, which is not necessarily loop func!
            // kernel << r.loop_name << "sum = 0; \n";
            // kernel << r.loop_name << "sh[threadIdx.x] = 0;\n";
            kernel << "_hila_kernel_set_zero(" << r.loop_name << "sum);\n";
            kernel << "_hila_kernel_set_zero(" << r.loop_name << "sh[threadIdx.x]);\n";
        } else if (r.reduction_type == reduction::PRODUCT) {
            kernel << r.loop_name << "sum = 1; \n";
            kernel << r.loop_name << "sh[threadIdx.x] = 1;\n";
        }

        // Replace references in the loop body
        for (Expr *e : r.refs) {
            loopBuf.replace(e, r.loop_name + "sum");
        }
    }


    /////////////////////////////////////////////////////////////////////////////
    // Standard boilerplate in CUDA kernels: calculate site index

    kernel << "int _HILA_thread_id = threadIdx.x + blockIdx.x * blockDim.x;\n";
    kernel << "unsigned " << looping_var << " = _HILA_thread_id + d_lattice.loop_begin;\n";

    if (!use_thread_blocks) {

        // The last block may exceed the lattice size. Do nothing in that case.
        kernel << "if(" << looping_var << " < d_lattice.loop_end) { \n";

    } else {
        // In this case let all threads to iterate over sites
        // use long long type (64 bits) just in case to avoid wrapping
        // the iteration will stop on the
        kernel << "for (long long _HILA_idx_l_ = " << looping_var
               << "; _HILA_idx_l_ < d_lattice.loop_end; _HILA_idx_l_ += N_threads * "
               << thread_block_number << ") {\n"
               << looping_var << " = _HILA_idx_l_;\n";
    }

    // Create temporary field element variables
    for (field_info &l : field_info_list) {
        if (l.is_read_nb) {
            // this field is nn-read
            if (!l.is_loop_local_dir) {
                // "standard" loop extern neighb direction
                for (dir_ptr &d : l.dir_list) {
                    std::string dirname;
                    if (d.is_constant_direction)
                        dirname = d.direxpr_s; // orig. string
                    else
                        dirname = remove_X(
                            loopBuf.get(d.parityExpr->getSourceRange())); // mapped name was
                                                                          // get_stmt_str(d.e);

                    // Check if the Direction is a variable. These have been renamed.
                    // for ( var_info & vi : var_info_list) for ( var_ref & vr : vi.refs
                    // )
                    //   if( vr.ref == d.e )
                    //     dirname = vi.new_name;

                    // Create the temp variable and call the getter
                    kernel << "const " << l.element_type << " " << d.name_with_dir << " = "
                           << l.new_name << ".get(" << l.new_name << ".neighbours[" << dirname
                           << "][" << looping_var << "], d_lattice.field_alloc_size);\n";

                    // and replace references in loop body
                    for (field_ref *ref : d.ref_list) {
                        loopBuf.replace(ref->fullExpr, d.name_with_dir);
                    }
                }
            } else {

                // now local var dependent neighbour
                // Two methods, either fetch all neighbours to vars and use them in loop,
                // or fetch element every time it is used. Clearly case-by-case either one
                // can be better.  Use the latter for now.
                // TODO: make selection automatic? At least tunable by pragma?

#ifdef BUILD_VARIABLE_NB_DIR_ELEMENTS


                std::string loop_array_name = l.new_name + "_dirs";
                kernel << l.element_type << ' ' << loop_array_name << "[NDIRS];\n";
                kernel << "for (int _HILAdir_ = 0; _HILAdir_ < NDIRS; "
                          "++_HILAdir_) {\n"
                       << loop_array_name << "[_HILAdir_] = " << l.new_name << ".get(" << l.new_name
                       << ".neighbours[_HILAdir_][" << looping_var
                       << "], d_lattice.field_alloc_size);\n}\n";

                // and replace references in loop body
                for (dir_ptr &d : l.dir_list) {
                    std::string dirname;
                    if (d.is_constant_direction)
                        dirname = d.direxpr_s; // orig. string
                    else
                        dirname = remove_X(
                            loopBuf.get(d.parityExpr->getSourceRange())); // mapped name was

                    for (field_ref *ref : d.ref_list) {
                        loopBuf.replace(ref->fullExpr, loop_array_name + "[" + dirname + "]");
                    }
                }

#else

                // and variable direction refs - use accessor directly
                for (dir_ptr &d : l.dir_list) {
                    std::string dirname;
                    if (d.is_constant_direction)
                        dirname = d.direxpr_s; // orig. string
                    else
                        dirname = remove_X(
                            loopBuf.get(d.parityExpr->getSourceRange())); // mapped name was

                    for (field_ref *ref : d.ref_list) {
                        loopBuf.replace(ref->fullExpr, l.new_name + ".get(" + l.new_name +
                                                           ".neighbours[" + dirname + "][" +
                                                           looping_var +
                                                           "], d_lattice.field_alloc_size)");
                    }
                }
#endif
            }
        }

        // TODO:
        if (l.is_read_atX || (loop_info.has_conditional && l.is_written)) {
            // local read
            // if var is not changed mark const
            // NOTE: const here swaps some methods from non-const
            // if (!l.is_written)
            //     kernel << "const ";
            kernel << l.element_type << " " << l.loop_ref_name << " = " << l.new_name << ".get("
                   << looping_var << ", d_lattice.field_alloc_size);\n           ";
            // if (l.is_written)
            //    kernel << "// TODO: READ MAY BE UNNECESSARY!  Do more careful
            //    assignment analysis!\n";

            if (loop_info.has_conditional && !l.is_read_atX) {
                kernel << "// Value of var " << l.loop_ref_name
                       << " read in because loop has conditional\n";
                kernel << "// TODO: MAY BE UNNECESSARY, write more careful analysis\n";
            }

        } else if (l.is_written) {
            // and a var which is not read
            kernel << l.element_type << " " << l.loop_ref_name << ";\n";
            kernel << "// Initial value of var " << l.loop_ref_name << " not needed\n";
        }

        // and finally replace references in body
        for (field_ref *ref : l.ref_list)
            if (!ref->is_direction) {
                loopBuf.replace(ref->fullExpr, l.loop_ref_name);
            }

        // // First create temp variables fields gathered from a Direction
        // for( field_ref *r : l.ref_list ) if( r->dirExpr ) {
        //   std::string dirname = get_stmt_str(r->dirExpr);

        //   // Check if the Direction is a variable. These have been renamed.
        //   for ( var_info & vi : var_info_list) for ( var_ref & vr : vi.refs )
        //     if( vr.ref == r->dirExpr ) {
        //       dirname = vi.new_name;
        //   }
        //   // Create the temp variable and call the getter
        //   kernel << l.element_type << " "  << l.loop_ref_name << "_" <<
        //   r->dirname
        //          << "=" << l.new_name << ".get(loop_lattice->d_neighb["
        //          << dirname << "][" << looping_var
        //          << "], loop_lattice->field_alloc_size);\n";
        // }
        // // Check for references without a Direction. If found, add temp variable
        // for( field_ref *r : l.ref_list ) if(r->dirExpr == nullptr){
        //   kernel << l.element_type << " "  << l.loop_ref_name << "="
        //          << l.new_name << ".get(" << looping_var
        //          << ", loop_lattice->field_alloc_size)" << ";\n";
        //   break;  // Only one needed
        // }

    } // ends Field handling

    // if there are site selections, reset the selection flag to 0 at all sites

    for (selection_info &s : selection_info_list) {
        if (s.previous_selection == nullptr) {
            kernel << s.maskname + '[' + looping_var + "] = 0;\n";
        }
    }

    /////////////////////////////////////////////////////////////////////////////

    // Finally, dump the loop body

    kernel << loopBuf.dump();
    if (semicolon_at_end)
        kernel << ';';
    kernel << '\n';

    // Call setters
    for (field_info &l : field_info_list)
        if (l.is_written) {
            std::string type_name = l.type_template;
            type_name.erase(0, 1).erase(type_name.end() - 1, type_name.end());
            kernel << l.new_name << ".set(" << l.loop_ref_name << ", " << looping_var
                   << ", d_lattice.field_alloc_size );\n";
        }

    // end the if ( looping_var < d_lattice.loop_end) or for() {
    kernel << "}\n";

    ///////////////////////////////////////////////////////////////////////////

    // Assign reductions to shared memory
    for (reduction_expr &r : reduction_list) {
        // OLD
        // kernel << r.loop_name << "sh[threadIdx.x] = " << r.loop_name << "sum;\n";
        kernel << "_hila_kernel_copy_var(" << r.loop_name << "sh[threadIdx.x], " << r.loop_name
               << "sum);\n";
    }

    // Handle reductions: Need to sync threads once, then do reduction
    // locally once per block
    bool sync_done = false;

    for (reduction_expr &r : reduction_list) {
        // Do sync (only if there is a reduction)
        if (!sync_done) {
            kernel << "__syncthreads();\n";
            sync_done = true;
        }

        // Now run the thread level reduction
        // kernel << "if( threadIdx.x == 0 ){\n";
        // if (vi.reduction_type == reduction::SUM) {
        // loopBuf.insert(br.Idx.at(0)->getSourceRange().getBegin(), "(");

        //     kernel << vi.new_name << "[blockIdx.x] = 0;\n";
        // } else if (vi.reduction_type == reduction::PRODUCT) {
        //     kernel << vi.new_name << "[blockIdx.x] = 1;\n";
        // }
        kernel << "for( int _H_i=N_threads/2; _H_i>0; _H_i/=2 ){\n";
        if (r.reduction_type == reduction::SUM) {
            kernel << "if(threadIdx.x < _H_i && _H_i +" << looping_var
                   << " < d_lattice.loop_end) {\n";
            // STD
            // kernel << r.loop_name << "sh[threadIdx.x] += " << r.loop_name
            //        << "sh[threadIdx.x+_H_i];\n";
            kernel << "_hila_kernel_add_var(" << r.loop_name << "sh[threadIdx.x], " << r.loop_name
                   << "sh[threadIdx.x + _H_i]);\n";
            kernel << "}\n";
            kernel << "__syncthreads();\n";
        } else if (r.reduction_type == reduction::PRODUCT) {
            kernel << "if(threadIdx.x < _H_i && _H_i +" << looping_var
                   << " < d_lattice.loop_end) {\n";
            kernel << r.loop_name << "sh[threadIdx.x] *= " << r.loop_name
                   << "sh[threadIdx.x+_H_i];\n";
            kernel << "}\n";
            kernel << "__syncthreads();\n";
        }
        kernel << "}\n";
        // kernel << "}\n";

        kernel << "if(threadIdx.x == 0) {\n"
               // << r.loop_name << "[blockIdx.x] = " << r.loop_name << "sh[0];\n";
               << "_hila_kernel_copy_var(" << r.loop_name << "[blockIdx.x], " << r.loop_name
               << "sh[0]);\n";
        kernel << "}\n";
    }
    // Reduction end handling stops here


    kernel << "}\n//------------------- end kernel " << kernel_name << " ---------------------\n\n";

    // Finally, emit the kernel
    // TheRewriter.InsertText(global.location.function,
    // indent_string(kernel),true,true);

    srcBuf *filebuf = get_file_srcBuf(global.location.kernels);
    filebuf->insert(global.location.kernels, // .getLocWithOffset(-1),
                    indent_string(kernel.str()), false, false);

    ////////////////////////////////////////////////////////////////////////////////////
    // Kernel is finished.  Write the end handling on call code

    code << ");\n\n"; // finishes kernel call
    code << "check_device_error(\"" << kernel_name << "\");\n";


    // If arrays were copied free memory

    for (array_ref &ar : array_ref_list) {
        if (ar.type == array_ref::REDUCTION) {

            if (loop_has_reductionvector_blocks) {
                // there are now N_blocks * loop_blocks -- reduce these
                // (code in hila_gpu.cpp)

                code << "sum_blocked_vectorreduction(" << ar.new_name << ", " << ar.size_expr
                     << ", " << thread_block_number << " * N_threads);\n";

                // after this the data can be collected from the array as in non-blocked reduction!
            }

            code << "{\nstd::vector<" << ar.element_type << "> a_v__tmp(" << ar.size_expr << ");\n";
            code << "gpuMemcpy(a_v__tmp.data(), " << ar.new_name << ", " << ar.size_expr
                 << " * sizeof(" << ar.element_type << "), "
                 << "gpuMemcpyDeviceToHost);\n\n";

            code << "for (int _H_tmp_idx=0; _H_tmp_idx<" << ar.size_expr << "; _H_tmp_idx++) "
                 << ar.name << "[_H_tmp_idx]";
            if (ar.reduction_type == reduction::SUM)
                code << " += a_v__tmp[_H_tmp_idx];\n";
            else
                code << " *= a_v__tmp[_H_tmp_idx];\n";

            code << " }\n";
        }

        if (ar.type != array_ref::REPLACE && (ar.size == 0 || ar.size > MAX_PARAM_ARRAY_SIZE)) {
            code << "gpuFree(" << ar.new_name << ");\n";
        }
    }

    // Check reduction variables
    for (reduction_expr &r : reduction_list) {
        // Run reduction
        if (r.reduction_type == reduction::SUM) {
            code << r.reduction_name << " = gpu_reduce_sum( dev_" << r.reduction_name
                 << ", N_blocks"
                 << ");\n";
        } else if (r.reduction_type == reduction::PRODUCT) {
            code << r.reduction_name << " = gpu_reduce_product( dev_" << r.reduction_name
                 << ", N_blocks"
                 << ");\n";
        }
        // Free memory allocated for the reduction
        code << "gpuFree(dev_" << r.reduction_name << ");\n";
    }


    // and the selection vars
    bool first = true;
    for (selection_info &s : selection_info_list) {
        if (s.previous_selection == nullptr) {
            // "reduce" the flagged arrays using cub functi
            if (s.assign_expr == nullptr) {
                code << s.new_name << ".endloop_action(" << s.maskname << ", " << s.sitename
                     << ");\n";
                code << "gpuFree(" << s.sitename << ");\n";
                code << "gpuFree(" << s.maskname << ");\n";
            } else {
                code << s.new_name << ".endloop_action(" << s.maskname << ", " << s.sitename << ", "
                     << s.valname << ");\n";
                code << "gpuFree(" << s.valname << ");\n";
                code << "gpuFree(" << s.sitename << ");\n";
                code << "gpuFree(" << s.maskname << ");\n";
            }
        }
    }


    return code.str();
}