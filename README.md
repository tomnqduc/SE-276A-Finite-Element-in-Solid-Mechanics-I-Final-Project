<h1>Final Project Code</h1>
This repository contains the code for the SE276AL Finite Element Method in Solid Mechanics at UCSD. The repository contain some utility functions, such as construct shape function, Jacobian matrix in 2D problem, generating integration points and cooresponding weights for Gaussian Quadratures, etc.
</br>
</br>
To see the project outcome, run the <code>FinalProject.m</code>, to adjust the refine parameters, i.e (8 x 8 or 16x16), change the following in the beginning of the <code>FinalProject.m</code>


<p style="text-align: center;"> 

  ```code
  noElemX = % Desired mesh size in x-direction 
  noElemY = % Desired mesh size in y-direction
