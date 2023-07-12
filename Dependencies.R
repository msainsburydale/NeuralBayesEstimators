CRANMIRROR <- "https://cran.csiro.au"

## Find the packages used throughout this repo using the package renv
if (!("devtools" %in% rownames(installed.packages())))
  install.packages("devtools", repos = CRANMIRROR)

# ---- Helper functions for this script ----

## Facilitates user input regardless of how this script was invoked
user_decision <- function(prompt, allowed_answers = c("y", "n")) {

  if (interactive()) {
    answer <- readline(paste0(prompt, "\n"))
  } else {
    cat(prompt)
    answer <- readLines("stdin", n = 1)
  }

  answer <- tolower(answer)

  if (!(answer %in% allowed_answers)) {
    tmp <- paste(allowed_answers, collapse = " or ")
    cat(paste0("Please enter ", tmp, ".\n"))
    answer <- user_decision(prompt, allowed_answers = allowed_answers)
  }

  return(answer)
}

# Installs the dependencies listed in Dependencies.txt
install_dependencies <- function(install_exact_versions) {

  x <- read.table("Dependencies.txt", header = FALSE)
  pkg_versions <- setNames(as.character(x[, 2]), x[, 1])
  rm(x)

  # ---- Non-CRAN packages ----

  ## These packages are treated individually because they are not available on
  ## CRAN, so we need to specify their repos.

  if("ggOceanMapsData" %in% names(pkg_versions) & !("ggOceanMapsData" %in% rownames(installed.packages())))
    #devtools::install_github("MikkoVihtakari/ggOceanMapsData")
    install.packages("ggOceanMapsData", repos = c("https://mikkovihtakari.github.io/drat", "https://cloud.r-project.org"), dependencies = TRUE)

  ## Remove this from the search list so that the script does not
  ## attempt to re-install them
  pkg_versions <- pkg_versions[!(names(pkg_versions) %in% c("ggOceanMapsData"))]

  # ---- CRAN packages ----

  ## Find the packages not yet installed and add them to the list
  installed_idx <- names(pkg_versions) %in% rownames(installed.packages())
  new_packages  <- names(pkg_versions)[!(installed_idx)]

  if (exists("install_exact_versions") && install_exact_versions) {
    ## Find the packages that are installed, but not the correct version
    installed_pkg_versions <- sapply(names(pkg_versions)[installed_idx], function(pkg) as.character(packageVersion(pkg)))
    idx          <- installed_pkg_versions != pkg_versions[installed_idx]
    already_installed_pkgs_different_versions <- names(installed_pkg_versions)[idx]
  }

  ## Now install the new packages: Here, we always install the correct
  ## package version (no reason not to)
  if(length(new_packages)) {
    cat("Package dependencies are being installed automatically using scripts/Dependencies_install.R\n")
    for (pkg in new_packages) {
      devtools::install_version(pkg, version = pkg_versions[pkg],
                                repos = CRANMIRROR, upgrade = "never",
                                dependencies = TRUE)
    }
  }

  # Change the already installed packages to the correct versions IF we have been told to do so
  if(exists("install_exact_versions") && install_exact_versions && length(already_installed_pkgs_different_versions)) {
    for (pkg in already_installed_pkgs_different_versions) {
      devtools::install_version(pkg, version = pkg_versions[pkg],
                                repos = CRANMIRROR, upgrade = "never",
                                dependencies = TRUE)
    }
  }
}


# ---- Install dependencies ----

install_depends <- user_decision("Do you want to automatically install package dependencies? (y/n)")
if (install_depends == "y") {
  install_exact_versions <- user_decision("Do you want to ensure that all package versions are as given in Dependencies.txt (this option is only recommended for use if there is a problem with the latest version of the packages)? (y/n)")
  install_exact_versions <- install_exact_versions == "y" # Convert to Boolean

  if (install_exact_versions) {
    cat("When changing the packages to the versions specified in Dependencies.txt, please use your discretion when answering the question “Which would you like to update?”.  Updating all packages (i.e., option 3) may cause errors.")
  }

  install_dependencies(install_exact_versions)
}
