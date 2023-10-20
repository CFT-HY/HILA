#include <sstream>
#include <iostream>
#include <string>

#include "stringops.h"
#include "hilapp.h"
#include "toplevelvisitor.h"


/////////////////////////////////////////////////////////////////////////////
/// Find the prototype for FunctionDecl fd, if exists
/// Clumsy, depends a lot on whether the func is normal/templated/method with
/// std or templated class

bool GeneralVisitor::find_prototype(FunctionDecl *fd, FunctionDecl *&decl) {


    // first, std function prototypes
    if (FunctionDecl *prevdecl = fd->getPreviousDecl()) {
        decl = prevdecl;
        return true;
    }

    // if the function is class method and fd within class declaration, no proto needed
    if (CXXMethodDecl *cxxmd = dyn_cast<CXXMethodDecl>(fd)) {
        auto *parent = cxxmd->getParent();
        if (parent->getSourceRange().getEnd() > cxxmd->getSourceRange().getBegin()) {
            return false;
        }
    }

    // if func is a template
    if (FunctionTemplateDecl *ftd = fd->getPrimaryTemplate()) {

        FunctionDecl *templ_fd = ftd->getTemplatedDecl();

        if (auto *prev = ftd->getPreviousDecl()) {
            decl = prev->getTemplatedDecl();

            // print_decl_info(decl, "PREVIOUS TEMPLATE");
            return true;
        }

        if (fd->getSourceRange().getBegin() != templ_fd->getSourceRange().getBegin()) {
            decl = templ_fd;
            // print_decl_info(templ_fd, "PRIMARY TEMPLATE");
            return true;
        }

        // print_decl_info(templ_fd, "NO TEMPL PROTOTYPE FOR FUNC ");
        // I think all templated funcs are cleared by this, so we need not go further
        return false;
    }

    // finally, if func is a template class (non-template)method we need to use different way
    if (CXXMethodDecl *cxxmd = dyn_cast<CXXMethodDecl>(fd)) {
        auto *parent = cxxmd->getParent();
        if (parent->getSourceRange().getEnd() < cxxmd->getSourceRange().getBegin()) {
            auto *msi = cxxmd->getMemberSpecializationInfo();
            if (msi) {
                NamedDecl *nd = msi->getInstantiatedFrom();
                // print_decl_info(nd, " NAMED DECL FOUND");
                if (nd) {
                    decl = dyn_cast<FunctionDecl>(nd);
                    return true;
                }
            }
        }
    }

    // print_decl_info(fd, "NO PROTOTYPE FOR FUNC ");

    return false;
}
