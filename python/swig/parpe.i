%module parpe

%include <stl.i>
%include "std_unique_ptr.i"

%ignore hdf5MutexGetLock;
%include "hdf5Misc.h"
%{
#include "hdf5Misc.h"
using namespace parpe;
%}

wrap_unique_ptr(OptimizationOptionsPtr, parpe::OptimizationOptions)

%{
#include "optimizationOptions.h"
using namespace parpe;
%}
%include "optimizationOptions.h"

// TODO: must somehow convince swig to std::move unique ptrs for these to work
%ignore HierachicalOptimizationWrapper;
%ignore HierachicalOptimizationProblemWrapper;
%ignore HierarchicalOptimizationReporter;

wrap_unique_ptr(AmiciSummedGradientFunctionPtr, parpe::AmiciSummedGradientFunction<int>)
wrap_unique_ptr(OptimizationResultWriterPtr, parpe::OptimizationResultWriter)
wrap_unique_ptr(AnalyticalParameterProviderPtr, parpe::AnalyticalParameterProvider)
%{
#include "hierarchicalOptimization.h"
using namespace parpe;
%}

%include "hierarchicalOptimization.h"
