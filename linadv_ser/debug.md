## Problems
- [ ] Error is large (66%)
- [ ] Due to periodic boundary conditions, `rank 0` and `last rank` need to communicate boundary conditions, which they are communicating.
- [ ] `update_lw!()` is wrong (only for parallel algorithm) since it will t reat u[25] as u[end] since for a single rank such as rank 0 u[25] is the last element but in general u[end] is u[N]. This maybe because that's how code is written for `update_lw!()`


## Solutions  
- [ ] For third problem there can be two solutions:
    - [ ] Halo exchange method (general method)
    - [ ] Alternatively, allgather the solution vector before calling `update_lw!()` (gather + send to all) but each rank will do a partly updating the solution using `update_lw!()` -- it will result in less MPI communication but more memory usage by each rank.  
- [ ] For second problem:
    - [ ]  `MPI Send/Recv` have been added for communicating boundary conditions.


## Debugging strategies
- [x] Check distribution of domain among processes -- *correct*
- [x] check end points of domain (should be similar to serial version) -- *little difference doesn't matter*
- [x] Remove distribution of domain and generate domain then distribute it to all. -- *no difference in both ways*
- [ ] Is there any other problem than difference in end points? -- yes