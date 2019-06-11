# FAQ

Q: The program is killed due to memory exhaustion, what should I do?

A: When running with MPI, the master process (rank 0) is consuming more memory
than the others. Consider reserving more memory for this one. For LoadLeveler,
this can be done conveniently via `#@ first_node_tasks`.
