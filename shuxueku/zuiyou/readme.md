```shell
mpirun --use-hwthread-cpus   -np  4 octave --eval 'pkg load mpi; op()'
```

```matlab
pkg load parallel




a=parcellfun(4,@op,{zeros(2,2),zeros(3,3),zeros(1,1).zeros(2,2)})

```


