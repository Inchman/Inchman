/*
 *  Kernel.h
 *
 *  Created on: 02/02/2012
 *  Created by: aidan
 */

#ifndef __gpgmp_Kernel_h__
#define __gpgmp_Kernel_h__


#include "Globals.h"

#include <cassert>
#include <iostream>
#include <string>

namespace gpgmp {


// todo: It's a bit confusing to have a gpgmp::kernel and a cl::kernel class!
// todo: are the OpenCL C++ wrappers actually outdated (version 1.1)?
class Kernel {

public:
    struct SkipArg {}; // handy for use with Kernel::operator()
    
    Kernel();
    Kernel(const cl::Program &program, const char *name);
    Kernel(const cl::Kernel &kernel, const cl::EnqueueArgs &enqueueArgs);
    Kernel(const cl::Program &program, const char *name, const cl::EnqueueArgs &enqueueArgs);
    virtual ~Kernel();
    
    inline cl::Kernel instance() const;
    
    inline const cl::EnqueueArgs &enqueueArgs() const;
    void setEnqueueArgs(const cl::EnqueueArgs &enqueueArgs);
    void printInternalInformation(cl::Buffer *d_errors);
    void setLocalMemory(const int localMemory);

    inline const cl::CommandQueue &commandQueue() const;
    inline const cl::NDRange &globalOffset() const;
    inline const cl::NDRange &globalSize() const;
    inline const cl::NDRange &localSize() const;
    void setCommandQueue(const cl::CommandQueue &queue);
    void setGlobalOffset(const cl::NDRange &offset);
    void setGlobalSize(const cl::NDRange &global);
    void setLocalSize(const cl::NDRange &local);

    template <class A1> void operator() (A1 &arg1);
    
    template <class A1, class A2> void operator() (A1 &arg1, A2 &arg2);
    
    template <class A1, class A2, class A3> void operator() (A1 &arg1, A2 &arg2, A3 &arg3);
    
    template <class A1, class A2, class A3, class A4> void operator() (A1 &arg1, A2 &arg2, A3 &arg3, A4 &arg4);
    
    template <class A1, class A2, class A3, class A4, class A5>
    void operator() (A1 &arg1, A2 &arg2, A3 &arg3, A4 &arg4, A5 &arg5);
    
    template <class A1, class A2, class A3, class A4, class A5, class A6>
    void operator() (A1 &arg1, A2 &arg2, A3 &arg3, A4 &arg4, A5 &arg5, A6 &arg6);
    
    template <class A1, class A2, class A3, class A4, class A5, class A6, class A7>
    void operator() (A1 &arg1, A2 &arg2, A3 &arg3, A4 &arg4, A5 &arg5, A6 &arg6, A7 &arg7);
    
    template <class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8>
    void operator() (A1 &arg1, A2 &arg2, A3 &arg3, A4 &arg4, A5 &arg5, A6 &arg6, A7 &arg7, A8 &arg8);
    
    template <class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9>
    void operator() (A1 &arg1, A2 &arg2, A3 &arg3, A4 &arg4, A5 &arg5, A6 &arg6, A7 &arg7, A8 &arg8, A9 &arg9);
    
    template <class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9, class A10>
    void operator() (A1 &arg1, A2 &arg2, A3 &arg3, A4 &arg4, A5 &arg5, A6 &arg6, A7 &arg7, A8 &arg8, A9 &arg9, A10 &arg10);
    
    template <class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9, class A10, class A11>
    void operator() (A1 &arg1, A2 &arg2, A3 &arg3, A4 &arg4, A5 &arg5, A6 &arg6, A7 &arg7, A8 &arg8, A9 &arg9, A10 &arg10, A11 &arg11);
    
    template <class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9, class A10, class A11, class A12>
    void operator() (A1 &arg1, A2 &arg2, A3 &arg3, A4 &arg4, A5 &arg5, A6 &arg6, A7 &arg7, A8 &arg8, A9 &arg9, A10 &arg10, A11 &arg11, A12 &arg12);

    template <class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9, class A10, class A11, class A12, class A13>
    void operator() (A1 &arg1, A2 &arg2, A3 &arg3, A4 &arg4, A5 &arg5, A6 &arg6, A7 &arg7, A8 &arg8, A9 &arg9, A10 &arg10, A11 &arg11, A12 &arg12, A13 &arg13);

    template <class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9, class A10, class A11, class A12, class A13, class A14>
    void operator() (A1 &arg1, A2 &arg2, A3 &arg3, A4 &arg4, A5 &arg5, A6 &arg6, A7 &arg7, A8 &arg8, A9 &arg9, A10 &arg10, A11 &arg11, A12 &arg12, A13 &arg13, A14 &arg14);

    template <class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9, class A10, class A11, class A12, class A13, class A14, class A15>
    void operator() (A1 &arg1, A2 &arg2, A3 &arg3, A4 &arg4, A5 &arg5, A6 &arg6, A7 &arg7, A8 &arg8, A9 &arg9, A10 &arg10, A11 &arg11, A12 &arg12, A13 &arg13, A14 &arg14, A15&arg15);

    template <class T> void setArg(size_t index, T &value);
    
    inline virtual void enqueueWithCurrentArgs();
    

private:
    cl::Kernel       m_instance;
    cl::EnqueueArgs  m_enqueueArgs;
};
    
    
inline cl::Kernel Kernel::instance() const { return m_instance; } /* note that cl::Kernel is cheap to copy */
    
inline const cl::EnqueueArgs & Kernel::enqueueArgs() const { return m_enqueueArgs; }
    
inline const cl::CommandQueue & Kernel::commandQueue() const { return m_enqueueArgs.queue_; }
inline const cl::NDRange & Kernel::globalOffset() const { return m_enqueueArgs.offset_; }
inline const cl::NDRange & Kernel::globalSize() const { return m_enqueueArgs.global_; }
inline const cl::NDRange & Kernel::localSize() const { return m_enqueueArgs.local_; }


template <class A1>
void Kernel::operator() (A1 &arg1) {
    //std::cout <<"Setting argument 0 for kernel "<<m_instance.getInfo<CL_KERNEL_FUNCTION_NAME>()<<".\n"<<std::flush;
    setArg(0, arg1);
    enqueueWithCurrentArgs();
}

template <class A1, class A2>
void Kernel::operator() (A1 &arg1, A2 &arg2) {
    //std::cout <<"Setting argument 1 for kernel "<<m_instance.getInfo<CL_KERNEL_FUNCTION_NAME>()<<".\n"<<std::flush;
    setArg(1, arg2);
    operator()(arg1);
}
template <class A1, class A2, class A3>
void Kernel::operator() (A1 &arg1, A2 &arg2, A3 &arg3) {
    //std::cout <<"Setting argument 2 for kernel "<<m_instance.getInfo<CL_KERNEL_FUNCTION_NAME>()<<".\n"<<std::flush;
    setArg(2, arg3);
    operator()(arg1, arg2);
}

template <class A1, class A2, class A3, class A4>
void Kernel::operator() (A1 &arg1, A2 &arg2, A3 &arg3, A4 &arg4) {
    //std::cout <<"Setting argument 3 for kernel "<<m_instance.getInfo<CL_KERNEL_FUNCTION_NAME>()<<".\n"<<std::flush;
    setArg(3, arg4);
    operator()(arg1, arg2, arg3);
}

template <class A1, class A2, class A3, class A4, class A5>
void Kernel::operator() (A1 &arg1, A2 &arg2, A3 &arg3, A4 &arg4, A5 &arg5) {
    //std::cout <<"Setting argument 4 for kernel "<<m_instance.getInfo<CL_KERNEL_FUNCTION_NAME>()<<".\n"<<std::flush;
    setArg(4, arg5);
    operator()(arg1, arg2, arg3, arg4);
}

template <class A1, class A2, class A3, class A4, class A5, class A6>
void Kernel::operator() (A1 &arg1, A2 &arg2, A3 &arg3, A4 &arg4, A5 &arg5, A6 &arg6) {
    //std::cout <<"Setting argument 5 for kernel "<<m_instance.getInfo<CL_KERNEL_FUNCTION_NAME>()<<".\n"<<std::flush;
    setArg(5, arg6);
    operator()(arg1, arg2, arg3, arg4, arg5);
}

template <class A1, class A2, class A3, class A4, class A5, class A6, class A7>
void Kernel::operator() (A1 &arg1, A2 &arg2, A3 &arg3, A4 &arg4, A5 &arg5, A6 &arg6, A7 &arg7) {
    //std::cout <<"Setting argument 6 for kernel "<<m_instance.getInfo<CL_KERNEL_FUNCTION_NAME>()<<".\n"<<std::flush;
    setArg(6, arg7);
    operator()(arg1, arg2, arg3, arg4, arg5, arg6);
}

template <class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8>
void Kernel::operator()
(A1 &arg1, A2 &arg2, A3 &arg3, A4 &arg4, A5 &arg5, A6 &arg6, A7 &arg7, A8 &arg8) {
    //std::cout <<"Setting argument 7 for kernel "<<m_instance.getInfo<CL_KERNEL_FUNCTION_NAME>()<<".\n"<<std::flush;
    setArg(7, arg8);
    operator()(arg1, arg2, arg3, arg4, arg5, arg6, arg7);
}

template <class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9>
void Kernel::operator()
(A1 &arg1, A2 &arg2, A3 &arg3, A4 &arg4, A5 &arg5, A6 &arg6, A7 &arg7, A8 &arg8, A9 &arg9) {
    //std::cout <<"Setting argument 8 for kernel "<<m_instance.getInfo<CL_KERNEL_FUNCTION_NAME>()<<".\n"<<std::flush;
    setArg(8, arg9);
    operator()(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
}

template <class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9, class A10>
void Kernel::operator()
(A1 &arg1, A2 &arg2, A3 &arg3, A4 &arg4, A5 &arg5, A6 &arg6, A7 &arg7, A8 &arg8, A9 &arg9, A10 &arg10) {
    //std::cout <<"Setting argument 9 for kernel "<<m_instance.getInfo<CL_KERNEL_FUNCTION_NAME>()<<".\n"<<std::flush;
    setArg(9, arg10);
    operator()(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
}

template <class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9, class A10, class A11>
void Kernel::operator()
(A1 &arg1, A2 &arg2, A3 &arg3, A4 &arg4, A5 &arg5, A6 &arg6, A7 &arg7, A8 &arg8, A9 &arg9, A10 &arg10, A11 &arg11) {
    //std::cout <<"Setting argument 10 for kernel "<<m_instance.getInfo<CL_KERNEL_FUNCTION_NAME>()<<".\n"<<std::flush;
    setArg(10, arg11);
    operator()(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10);
}

template <class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9, class A10, class A11, class A12>
void Kernel::operator()
(A1 &arg1, A2 &arg2, A3 &arg3, A4 &arg4, A5 &arg5, A6 &arg6, A7 &arg7, A8 &arg8, A9 &arg9, A10 &arg10, A11 &arg11, A12 &arg12) {
    //std::cout <<"Setting argument 11 for kernel "<<m_instance.getInfo<CL_KERNEL_FUNCTION_NAME>()<<".\n"<<std::flush;
    setArg(11, arg12);
    operator()(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11);
}

template <class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9, class A10, class A11, class A12, class A13>
void Kernel::operator()
(A1 &arg1, A2 &arg2, A3 &arg3, A4 &arg4, A5 &arg5, A6 &arg6, A7 &arg7, A8 &arg8, A9 &arg9, A10 &arg10, A11 &arg11, A12 &arg12, A13 &arg13) {
    //std::cout <<"Setting argument 12 for kernel "<<m_instance.getInfo<CL_KERNEL_FUNCTION_NAME>()<<".\n"<<std::flush;
    setArg(12, arg13);
    operator()(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12);
}

template <class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9, class A10, class A11, class A12, class A13, class A14>
void Kernel::operator()
(A1 &arg1, A2 &arg2, A3 &arg3, A4 &arg4, A5 &arg5, A6 &arg6, A7 &arg7, A8 &arg8, A9 &arg9, A10 &arg10, A11 &arg11, A12 &arg12, A13 &arg13, A14 &arg14) {
    //std::cout <<"Setting argument 13 for kernel "<<m_instance.getInfo<CL_KERNEL_FUNCTION_NAME>()<<".\n"<<std::flush;
    setArg(13, arg14);
    operator()(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13);
}

template <class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9, class A10, class A11, class A12, class A13, class A14, class A15>
void Kernel::operator()
(A1 &arg1, A2 &arg2, A3 &arg3, A4 &arg4, A5 &arg5, A6 &arg6, A7 &arg7, A8 &arg8, A9 &arg9, A10 &arg10, A11 &arg11, A12 &arg12, A13 &arg13, A14 &arg14, A15 &arg15) {
    //std::cout <<"Setting argument 13 for kernel "<<m_instance.getInfo<CL_KERNEL_FUNCTION_NAME>()<<".\n"<<std::flush;
    setArg(14, arg15);
    operator()(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13, arg14);
}

template <class T>
void Kernel::setArg(size_t index, T &value) {
    m_instance.setArg(index, value);
}

    /*
template <>
inline void Kernel::setArg<cl::Buffer *>(size_t index, cl::Buffer * &buffer) {
    m_instance.setArg(index, *buffer);
}
    */

    /*
template <>
inline void Kernel::setArg<const cl::Buffer *>(size_t index, const cl::Buffer * &buffer) {
    m_instance.setArg(index, *buffer);
}
    */

// handy for use with Kernel::operator() if you've already set the argument in the past
template <>
inline void Kernel::setArg<Kernel::SkipArg>(size_t /*index*/, Kernel::SkipArg & /*value*/)
{} // should be optimized out

inline void Kernel::enqueueWithCurrentArgs()
{
    std::string kernelName = m_instance.getInfo<CL_KERNEL_FUNCTION_NAME>();
    try {
        m_enqueueArgs.queue_.enqueueNDRangeKernel(m_instance, m_enqueueArgs.offset_, m_enqueueArgs.global_, m_enqueueArgs.local_);
        m_enqueueArgs.queue_.finish();
    }  catch (cl::Error error) {
        std::cout <<std::flush;
        std::cerr << "ERROR in kernel " << kernelName<<" (";
        printRange(std::cerr, m_enqueueArgs.global_);
        std::cerr <<", ";
        printRange(std::cerr, m_enqueueArgs.local_);
        std::cerr << "): "
                << error.what()
                << "("
                << oclErrorString(error.err())
                << ")"
                << std::endl
                << "Exiting..."
                << std::endl
                << std::flush;
        exit(1);
    }
    //std::cout <<"done. \n" <<std::flush;
}
    
    
} // namespace gpgmp

#endif // !__gpgmp_Kernel_h__
