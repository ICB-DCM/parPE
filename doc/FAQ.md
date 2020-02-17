# FAQ

**Q: The program is killed due to memory exhaustion, what should I do?**

**A**: When running with MPI, the master process (rank 0) is consuming more memory
than the others. Consider reserving more memory for this one. For LoadLeveler,
this can be done conveniently via `#@ first_node_tasks`.

---

**Q: Building the example models fails with something like**
`AttributeError: 'SwigPyObject' object has no attribute 'getParameterIds'`
or `swig/python detected a memory leak of type 'std::unique_ptr< amici::Model > *', no destructor found.`

**A**: This is most likely due to mismatching AMICI source code and
AMICI Python module. Delete `build/venv` and `build/examples/parpeamici`
and retry. If it still fails, make sure the correct AMICI source
directory is used.

