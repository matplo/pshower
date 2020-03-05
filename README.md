# DGLAP parton shower / DN

- Original DGLAP shower code by Duff Neill
- Some modifications by Felix Ringer
- Some code fixes and swig-pythonization by matplo

# dependencies

- CMake
- SWIG
- GLS (GLSBLAS)

# compilation

```
./scripts/build_dglap.sh
```

# example

```
cd test
./run_dglap_dn_test.sh --nev=1000
```

OR run the python version

```
cd test
./run_dglap_dn_test.py --nev 1000
```

try `--help` for more options

```
./run_dglap_dn_test.py --help
```




