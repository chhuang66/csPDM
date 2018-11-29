# csPDM
These codes are the implementations of colored-synapse population density method (csPDM), used to simulate large-scale network dynamics of one single neuronal population where neurons receive both the exitatatory and inhibitory inputs only from the outside. The details of the csPDM  are described in Huang, Chih-Hsu and Lin, Chou-Ching K. <i>An efficient population density method for modeling neural networks with synaptic dynamics manifesting finite relaxation time and short-term plasticity</i>, accepted by the <i>eNeuro</i> (prepared for publication, 2018).

### Compilation
-------------
>This implementaion is written by C++ and is compiled by the g++6.3 with cmake3.7. The required libraries are "hdf5" for storing simulation results, "blas" and "lapacke" for numerical solving ODEs through the backward Euler method. The compiled program has been tested in the linux OS environment and runs well, but has not been tested in the MS environment.


1. Create a new sub-directory, e.g., build, in the root directory, and move to this sub-directory.
2. Use the following commands to compile the codes
*  <code>cmake ..</code>
*  <code>make</code>
3. Finally, a program, named by _aEIFONETHREE_semiLocal_ is built.

### Usages
- help: <code>aEIFONETHREE_semiLocal --help</code>
- check external input sources available: <code>aEIFONETHREE_semiLocal -ihelp</code>
- check connections available: <code>aEIFONETHREE_semiLocal -chelp</code>
- check neuronal parameters available: <code>aEIFONETHREE_semiLocal --PYparam=help</code>
- check simulation parameters available: <code>aEIFONETHREE_semiLocal -phelp</code>
- initialize a new simulation: <code>aEIFONETHREE_semiLocal -I [Options] filename </code>
