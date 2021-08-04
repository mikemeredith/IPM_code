
# Checking code in this repository

When the checking script in this directory is `source`d, it will run every "*.R" script in R's current directory and subdirectories, except those with "#" in the name.

It creates a log file in the target directory with the name "#check_<date>.log".

For each file, it logs the name of the file, the date the file was modified, the time the run started, any error messages, and the script execution time if greater than 20 secs. It checks that graphical parameters are restored to the original values and reports discrepancies.

All graphics produced by the script are saved in a PDF file with the same name, but ending in `#.pdf`. After running the script, the workspace is saved in a .RData file, including an object called `sessionInfo` with the output from `sessionInfo()`. The workspace file name ends with `#.RData`, or `_bad#.RData` if the script threw an error.

If there is an appropriately-named `#.RData` file in the folder from a previous run, the script will load it and compare the current output with the last run using `all.equal`. Objects that don't match are listed in the log file.

It cleans up the workspace and detaches packages between running scripts.

When all are done, it records the overall time taken and the output from `sessionInfo()`.

## Installing and updating packages

Before running the scripts, you may want to check that the necessary packages are installed. I suggest first running `update.packages(ask='graphics',checkBuilt=TRUE)` to ensure everything is up to date, including dependencies. You can then install any extra packages needed with:
```
needed <- c("IPMbook", "jagsUI", "scales", "MCMCglmm", "AHMbook", "wiqid", "RColorBrewer",
    "denstrip", "plotrix")
got <- rownames(installed.packages())

( notgot <- needed[!needed %in% got] )

install.packages(notgot)
```


If you want to install latest devel/patched versions of packages from Github, use or adapt the following code:
```
remotes::install_github("mikemeredith/IPMbook")
packageVersion("IPMbook")
remotes::install_github("mikemeredith/AHMbook")
packageVersion("AHMbook")
remotes::install_github("rbchan/unmarked")
packageVersion("unmarked")
remotes::install_github("kenkellner/jagsUI")
packageVersion("jagsUI")
```
