/*
 *  Kernel.cpp
 *
 *  Created on: 02/02/2012
 *  Created by: aidan
 */

#include "Kernel.h"


namespace gpgmp {


Kernel::Kernel()
:   m_enqueueArgs(cl::EnqueueArgs(cl::CommandQueue() /* very important to pass this "null" queue, otherwise a default one will be automatically created on us */, cl::NullRange))
{}
    
Kernel::Kernel(const cl::Kernel &instance, const cl::EnqueueArgs &enqueueArgs)
:   m_instance(instance),
    m_enqueueArgs(enqueueArgs)
{}

Kernel::Kernel(const cl::Program &program, const char *name, const cl::EnqueueArgs &enqueueArgs)
:   m_instance(cl::Kernel(program, name)),
    m_enqueueArgs(enqueueArgs)
{}
    
Kernel::Kernel(const cl::Program &program, const char *name)
    : m_instance(cl::Kernel(program, name)),
      m_enqueueArgs(cl::EnqueueArgs(cl::CommandQueue(), cl::NullRange))
{}

Kernel::~Kernel()
{}

/***
  Reads out the internal information structure that kernels can use to pass information to the host system.
  The information is printed to stdout.

  @todo: Have the d_error memory handled completely by the Kernel instead of DiffusionModel
  ***/
void Kernel::printInternalInformation(cl::Buffer *d_errors) {
    std::vector<Real> h_errors(GMP_NUM_ERRORS, 0);

    // read out errors - todo: make that optional in kernel class
    commandQueue().enqueueReadBuffer(*d_errors, CL_TRUE, 0,
                                     GMP_NUM_ERRORS * sizeof(Real),
                                     (void *) (&(h_errors.front())));
    commandQueue().flush();

    std::cout <<"Error codes are for kernel "<<m_instance.getInfo<CL_KERNEL_FUNCTION_NAME>()<<" are:\n";
    for (int err = 0; err<GMP_NUM_ERRORS; err++)
        std::cout <<"Code "<<err<<" is "<<h_errors[err]<<".\n";
    std::cout << std::flush;
}
    
void Kernel::setEnqueueArgs(const cl::EnqueueArgs &enqueueArgs) { m_enqueueArgs = enqueueArgs; }

void Kernel::setCommandQueue(const cl::CommandQueue &queue) { m_enqueueArgs.queue_ = queue; }
void Kernel::setGlobalOffset(const cl::NDRange &offset) { m_enqueueArgs.offset_ = offset; }
void Kernel::setGlobalSize(const cl::NDRange &global) { m_enqueueArgs.global_ = global; }
void Kernel::setLocalSize(const cl::NDRange &local) { m_enqueueArgs.local_ = local; }

/*!
  Sets the local memory argument of the kernel. By default the local memory is the last parameter
  present in the kernel function.

  \param Size of local memory in bytes.
  */
void Kernel::setLocalMemory(const int localMemory) {
    int localMemArgIndex = m_instance.getInfo<CL_KERNEL_NUM_ARGS>()-1;
    //std::cout <<"Setting local memory argument "<<localMemArgIndex<<".\n"<<std::flush;
    m_instance.setArg(localMemArgIndex, localMemory, 0);
}
} // namespace gpgmp
