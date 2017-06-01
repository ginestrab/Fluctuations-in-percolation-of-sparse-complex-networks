# Fluctuations-in-percolation-of-sparse-complex-networks
  This code implements the message passing algorithm to evaluate fluctuation in percolation of sparse complex networks
  and compairs the result to brute force average over different realizations of the initial damage.
  
  This code can be redistributed and/or modified
  under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or (at
  your option) any later version.
   
  This program is distributed ny the authors in the hope that it will be 
  useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 
   
  If you use this code please cite the following  paper:
 
  G. Bianconi, arXiv preprint, arXiv:1703.05528 (2017) 
 
  (c) G. Bianconi (email: ginestra.bianconi@gmail.com ) 


 
 Using the notation of the paper our code used node-independent probabilities p^{(1)}=p^{(2)}=p 
 that a node is not initially damaged in each of the two realization of the percolation problem
 and a probability 
 p^{[11]} that the node is not initially damaged in both realizations of the percolation problem
 Here we take
 p^{[11]}=p^a with a>=1 to be an imput parameter of the code
