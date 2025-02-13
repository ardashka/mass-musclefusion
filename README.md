# README #

### Many Active Systems Simulations ###

drylyotropicF.cpp was developed by Aleksandra Ardaševa and Amin Doostmohammadi 
and was used in Yoann Le Toquin, Sushil Dubey, Aleksandra Ardaševa, Lakshmi 
Balasubramaniam, Emilie Delaune, Valérie Morin,  Amin Doostmohammadi,  
Christophe Marcelle,  Benoît Ladoux. 2024. 'Mechanical stresses govern 
myoblast fusion and myotube growth'. https://doi.org/10.1101/2024.11.22.624831

### Compiling ###

The code relies on the boost::program_options library. Once it is installed, a
simple `make` in the main directory should do it.

Can be used on any operating system (Ubuntu, Windows, Linux).
Compilation time <1 min.

### Running ###

The code is run from the command line and a runcard must always be given as the
first argument:

./mass runcard.dat -fco out 

A runcard is a simple file providing the parameters for the run. Every option can 
also be given to the program using the command line as 
`./mass runcard.dat --option=arg`. A complete list of available options can be
obtained by typing `./mass -h`.

Demo runtime: ~2h.

### Output ###

By default the program writes output files in the current directory. This can be
changed using `--output=dir/` or `-o dir/`, where `dir/` is the target
directory. The program also supports compressed output with the option flag
`--compression` or `-c`. When compression is on, the output name does not
denote the target directory but rather the file name of the compressed archive.
Note that with the current implementation, compression is not recommended for
long runs as the time to add a file to the archive grows with the size of the
archive.

Type `./mass -h` or `./mass -m model-name -h` for a list of available options.

The output is saved to .json file. The output contains:
QQxx, QQyx: nematic (cell) tensor components
fQQxx, fQQyx: ECM tensor components
sigmaXX, sigmaYY: stress tensor components
outS: isotropic stress 
