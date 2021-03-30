

# Install or update packages needed to run the code in this repository
# ====================================================================
# You will need the IPMbook package, which is not yet on CRAN. You can install it with
remotes::install_github("mikemeredith/IPMbook")

# All the rest are available from CRAN

# I suggest first running
update.packages(ask='graphics',checkBuilt=TRUE)
# to ensure everything is up to date, including dependencies.

needed <- c("jagsUI", "scales", "MCMCglm", "AHMbook", "wiqid", "RColorBrewer",
    "denstrip", "plotrix")
got <- rownames(installed.packages())

( notgot <- needed[!needed %in% got] )

install.packages(notgot, dependencies=TRUE)

# 'devel' versions of packages
# ----------------------------
# If you want to try out devel versions of packages from GitHub, install these
#  AFTER 'update.packages' as that may "downdate" to the latest CRAN version.
# For example:
# devtools::install_github("mikemeredith/AHMbook")
# packageVersion("AHMbook")
# news(package="AHMbook")

