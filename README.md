# The-expulsion-of-dust-particles-by-stellar-radiation-pressure

This is the project in insterstellar medium.

The code is writen in **Python 3**.
Required libraries: [numpy](https://pypi.org/project/numpy/), [miepython](https://pypi.org/project/miepython/).  
How to install them in Ubuntu:

```bash
pip3 install numpy miepython
```
Using the Mie theory and rayleigh approximation the program calculates:

-planck mean radial pressure factors for different type of dust (Ice, amorphous carbon, silicates),

-dust extraction efficiency for different star types.

Refractive indices for  Ice, amorphous carbon, silicates are given in ICE_WAR.RI, AC_RM.ri, SilicateDraine03.RI


```bash
python3 laba2.py
```
 

You may need to install R on your system and Rstudio.
Required libraries: ggplot2, ggpub, stringi,stringr, scales
In Rstudio terminal:
```bash
install.packages(c("ggplot2", "ggpub", "stringi", "stringr", "scales"))
```
Required libraries: [numpy](https://pypi.org/project/numpy/), [miepython](https://pypi.org/project/miepython/).  
```bash
Rscript plot.r
```
