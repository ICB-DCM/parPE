%module parpe

%include <stl.i>
%include "std_unique_ptr.i"

%ignore hdf5MutexGetLock;
%include "hdf5Misc.h"
%{
#include "hdf5Misc.h"
using namespace parpe;
%}

wrap_unique_ptr(GradientFunctionPtr, parpe::GradientFunction)

%{
#include "functions.h"
using namespace parpe;
%}
%include "functions.h"

wrap_unique_ptr(OptimizationOptionsPtr, parpe::OptimizationOptions)

%{
#include "optimizationOptions.h"
using namespace parpe;
%}
%include "optimizationOptions.h"

%feature("action") parpe::OptimizationReporter::OptimizationReporter(GradientFunction *gradFun, std::unique_ptr<OptimizationResultWriter> rw){ result = (parpe::OptimizationReporter *)new parpe::OptimizationReporter(arg1,std::move(arg2));}
%feature("action") parpe::OptimizationReporter::OptimizationReporter(GradientFunction *gradFun, std::unique_ptr<OptimizationResultWriter> rw){ result = (parpe::OptimizationReporter *)new parpe::OptimizationReporter(arg1,std::move(arg2));}
// Auto-generated setter would generate code which calls the copy constructor
// This affects getters as well as setters, not sure how to distinguish%feature("action") parpe::OptimizationProblem::costFun
%ignore parpe::OptimizationProblem::costFun;
%ignore parpe::OptimizationReporter::resultWriter;
%ignore parpe::OptimizationReporterPtr::resultWriter;

wrap_unique_ptr(OptimizationResultWriterPtr, parpe::OptimizationResultWriter);
wrap_unique_ptr(OptimizationReporterPtr, parpe::OptimizationReporter);
%{
#include "optimizationProblem.h"
using namespace parpe;
%}
%include "optimizationProblem.h"


// TODO: must somehow convince swig to std::move unique ptrs for these to work
%ignore HierachicalOptimizationWrapper;
%ignore HierachicalOptimizationProblemWrapper;
%ignore HierarchicalOptimizationReporter;

wrap_unique_ptr(AmiciSummedGradientFunctionPtr, parpe::AmiciSummedGradientFunction<int>)
wrap_unique_ptr(AnalyticalParameterProviderPtr, parpe::AnalyticalParameterProvider)
%{
#include "hierarchicalOptimization.h"
using namespace parpe;
%}
%include "hierarchicalOptimization.h"
