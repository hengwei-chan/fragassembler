# Fragment Assembler (FragAssembler)

FragAssembler aims on the computer-assisted structure elucidation of (un)knowns compounds by re-assembling 
known fragments of a dedicated library, solely based on given 13C NMR information. <br> 

This approach follows the idea by: <br> 
M. Will et al., <br> 
J. Chem. Inf. Comput. Sci. 1996, 36, 221-227 <br>
https://doi.org/10.1021/ci950092p  

## SSC Library 
 
Each fragment is stored as a substructure-subspectrum-correlation (SSC) object in a SSC library and 
consists of its substructure and the assigned chemical shift (13C NMR) for each carbon atom. <br>
Therefore, the NMRShiftDB<sup>1</sup> was used to build such SSC library. Each 
structure in NMRShiftDB was fragmented at each atom by using the hierarchical order of HOSE codes<sup>2</sup> 
with maximum spheres of 2, 3 and 4.

## Process 

Once a 13C NMR query spectrum, incl. multiplicities, of the compound under investigation is given, each SSC in SSC 
library is compared with the query spectrum and is considered to be a hit if all signals in SSC could match to a 
signal in the query spectrum. All SSC hits are then ranked according their spectral similarity to the query spectrum. <br> 
After that the assembling process starts with the first SSC hit and tries to extend its substructure by the 
next ranked SSCs and their substructures if there is a structural overlap (HOSE code) and 
some validation steps are passed successfully. <br>
The result is a list of assembled structure proposals for a query spectrum in a file containing SMILES.


##References 
<sup>1</sup> C. Steinbeck et al., J. Chem. Inf. Comput. Sci.20034361733-1739, https://doi.org/10.1021/ci0341363 <br>
<sup>2</sup> W. Bremser, 1978, Analytica Chimica Acta, 103(4), 355-365., https://doi.org/10.1016/S0003-2670(01)83100-7