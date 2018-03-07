Biased Random-Key Genetic Algorithm for the k-Interconnected Multi-Depot
Multi-Traveling Salesmen Problem
===============================

This project implements a biased random-key genetic algorithm (BRKGA) 
for the k-Interconnected Multi-Depot Multi-Traveling Salesmen Problem
(k-IMDMTSP) proposed in

> C.E. Andrade, F.K. Miyazawa, M.C.G Resende. Evolutionary Algorithm for the
> k-Interconnected Multi-Depot Multi-Traveling Salesmen Problem. Proceedings of
> the Fifteen International Conference on Genetic and Evolutionary Computation
> (GECCO' 13). Pages 463-470, New York, NY, USA, 2013. DOI:
> [10.1145/2463372.2463434](http://dx.doi.org/10.1145/2463372.2463434)

and

> M.C. Lopes, C.E. Andrade, T.A. Queiroz, M.G.C. Resende, F.K. Miyazawa.
> Heuristics for a Hub Location-Routing Problem. Networks, volume 68, number 1,
> pages 54-90, 2016. DOI:
> [10.1002/net.21685](http://dx.doi.org/10.1002/net.21685)

Please, refer to the above papers for detailed algorithms.

Before to use this package, code, and any resource attached to which,
and/or its derivatives, you must accept the **[license](LICENSE.md)**.
And, please, do not forget to refer the above papers.

[Contribution guidelines for this project](CONTRIBUTING.md)


Dependencies
-------------------------------

- [Boost libraries >= 1.58;](http://www.boost.org)

- [Lemon libraries >= 1.3.1](http://lemon.cs.elte.hu)

- Modern C++ compiler supporting C++11;
    - This code was tested using GCC 5.2 and GCC 7.2 successfully. Clang 3.9
      also compiles the code without problems, although changes in the
      parameters are necessary. A caveat is that to the date of this document
      was written, Clang does not use multiple threads of OpenMP, the
      technology behind the paralyzation of the decoding process.

- [Doxygen for documentation](http://www.stack.nl/~dimitri/doxygen)
    (not required to build the code).

This package also includes a modified version of the
[brkgaAPI](https://github.com/rfrancotoso/brkgaAPI). 
The core functionality was kept, but changes in the API make the distributed
version and the original version incompatible.

This package also includes (together brkgaAPI) a Mersenne Twister random number
generator licensed under BSD-3.


Compiling and testing
-------------------------------

First, you must set the proper paths to the "includes" and "libs" in your
system. In the `Makefile`, go to the Section "Lib and include definitions" 
(line 109) and

- Set `INCLUDE_PATH` and `USER_LIBDIRS`, pointing to the location of your boost
  and lemon headers and libraries;

Note that depending of the configuration of your system, other adjustments
may be necessary.

To build the BRKGA, just use:
```bash
$ make all
```

For documentation, type (assuming that you have Doxygen installed):
```bash
$ make doc
```


Debugging and tuning
-------------------------------

In the very begin of the Makefile, we can find the optimization, debug, and
tuning flags. When line 29 is active, 
```makefile
OPT = opt
```

it enables the compilation using optimization flags suitable for tuning and
actual experimentation. Commenting line 29, we activate the compiler debug
flags. Such scenario is ideal for tools like 
[GNU GDB](https://www.gnu.org/software/gdb) and 
[Valgrind](http://valgrind.org).

Lines 33 and 34 enable several additional debugging information, which is printed
in the standard output. Lines 33 and 34 are
```makefile
USER_DEFINES += -DDEBUG
USER_DEFINES += -DFULLDEBUG 
```

Finally, line 37 enables the tuning mode. This mode reads all parameters from
the command line, ignoring the configuration file entirely (although it must be
supplied anyway). This capability enables the usage of tuning tools such as
[irace](https://cran.r-project.org/web/packages/irace){:target="_blank"}.
Line 37 is
```makefile
USER_DEFINES += -DTUNING
```

