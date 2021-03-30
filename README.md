# IPM code

The book *Integrated Population Modeling* by Michael Schaub and Marc Kéry, due for publication in July 2021, contains lots of R and BUGS code.

The R package `IPMbook` has all the data sets and the custom functions used in the book. Commented code for the functions is on GitHub [here](https://github.com/mikemeredith/IPMbook).

This repository has all the code in the printed book, plus code referred to as "available on the website" but not printed. The aim is to have code which works with current versions of R, JAGS and contributed R packages. The code is regularly tested and updated code inserted, with the original printed code retained but commented out with `#`. Please open an issue if you find other code which does not work.

Additional code and comments not in the printed book are marked off with twiddly lines like this:
```
#~~~~ oldfunction has been replaced with newfunction ~~~~~~~
# oldfunction(foo)
newfunction(foo)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
```
