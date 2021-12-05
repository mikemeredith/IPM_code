# IPM code

The book *Integrated Population Models* by Michael Schaub and Marc Kéry, published in December 2021, contains lots of R and BUGS code.

The R package `IPMbook` on CRAN has all the data sets and the custom functions used in the book. Commented code for the functions is on GitHub [here](https://github.com/mikemeredith/IPMbook).

This repository has all the code in the printed book, plus code referred to as "available on the website" but not printed.

The aim is to have code which works with current versions of R, JAGS and contributed R packages. The code is regularly tested, and some scripts have "test run" versions with fewer iterations for simulations or MCMC runs. These files have `_tr.R` at the end of the name. Remove/comment out the lines flagged `# ~~~ for testing` for a full run.

Where necessary, updated code is inserted, with the original printed code retained but commented out with `#`. Please open an issue if you find other code which does not work.

Additional code and comments not in the printed book are marked off with twiddly lines like this:
```
#~~~~ oldfunction has been replaced with newfunction ~~~~~~~
# oldfunction(foo)
newfunction(foo)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
```
