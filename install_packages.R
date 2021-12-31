

# Install or update packages needed to run the code in this repository
# ====================================================================

# I suggest first running
update.packages(ask='graphics',checkBuilt=TRUE)
# to ensure everything is up to date, including dependencies.

needed <- c("IPMbook", "jagsUI", "scales", "MCMCglmm", "AHMbook", "wiqid", "RColorBrewer",
    "denstrip", "plotrix", "fields")
got <- rownames(installed.packages())

( notgot <- needed[!needed %in% got] )

install.packages(notgot, dependencies=TRUE)

# 'devel' versions of packages
# ----------------------------
# If you want to try out devel versions of packages from GitHub, install these
#  AFTER 'update.packages' as that may "downdate" to the latest CRAN version.
# For example:
# remotes::install_github("mikemeredith/IPMbook")
# packageVersion("IPMbook")
# news(package="IPMbook")
# remotes::install_github("kenkellner/jagsUI")
# packageVersion("jagsUI")
# news(package="jagsUI")
# remotes::install_github("mikemeredith/AHMbook")
# packageVersion("AHMbook")
# news(package="AHMbook")

