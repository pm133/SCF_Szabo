#SCF_Szabo

This code has been written to convert the SCF sample code from Appendix B of the Modern Quantum Chemistry - An Introduction to Electronic Structure (A. Szabo and N. Ostlund) from Fortran IV to C.

Comments are welcome by email - juansanshoo@hotmail.co.uk but note that this code is not designed to be either computationally efficient or code efficient. It is written specifically to make it easier to follow the SCF procedure in the book.
To that effect it is limited to the STO-nG (n=1, 2 or 3) basis set and to the 2 atom molecule HeH+.
The matrix diagonalisation procedure depends on the system having just two atoms.

I have included a makefile, an executable and a log file showing the output.

The code originally written by Szabo and Ostlund is a great launchpad for those interested in developing code of this nature for full use in general computational chemistry systems.
I have commented the code to show which equations in the book are being coded at each step of the SCF procedure and I have also numbered the steps involved to help others see what is going on.


