## Strategy
- Array on each rank that is `u` and `unew` will be of size `N_local+2`.
- The 2 extra values are for ghost values (including boundary conditions for `rank 0 and rank=size`)


## Bugs
- [x] Error is large (66%) -- *resolved*
- [x] When `size=2` error is large (64%)
- [x] Due to periodic boundary conditions, `rank 0` and `last rank` need to communicate boundary conditions, which they are not communicating.
- [x] `update_lw!()` is wrong (only for parallel algorithm) since it will treat `u[25]` (of rank 0) as `u[end]` since for a single rank such as rank 0 `u[25]` is the last element but in general `u[end]` is last element and `u[25]` is in the middle. This is because that's how code is written for serial `update_lw!()`
- [x] `while` loop condition should be `j < t` not `j <= t` -- caused infinite loop


## Solutions  
- For third bug there can be two solutions:
    - [x] Halo exchange method (general method)
    - [ ] Alternatively, allgather the solution vector before calling `update_lw!()` (gather + send to all) but each rank will do a partly updating the solution using `update_lw!()` -- it will result in less MPI communication but more memory usage by each rank.  -- **Not Effective**
- For second bug:
    - [x] `MPI Send/Recv` have been added for communicating boundary conditions in `get_ghost_values!()` function. 
- For 4th bug:
    - [x] `update_lw!()` will be different for serial and parallel versions. Strategy for my serial and parallel versions is different. 


## Debugging strategies
- [x] Check distribution of domain among processes -- *correct*
- [x] check end points of domain (should be similar to serial version) -- *little difference but matters*
- [x] Remove distribution of domain and generate domain then distribute it to all. -- *no difference in both ways*
- [x] Is there any other problem than difference in end points? -- yes (see bug 3)