.onAttach <- function(...) {
  msg <- create_startup_message()
  packageStartupMessage(msg)
}

# Create package startup message
create_startup_message <- function() {
  v_pkg <- create_desc("pkmstate")
  msg <- paste0(
    "Attached pkmstate", v_pkg,
    ". Type ?pkmstate to get started."
  )
  return(msg)
}

# Create package description
create_desc <- function(pkg_name) {
  Lib <- dirname(system.file(package = pkg_name))
  pkgdesc <- suppressWarnings(
    utils::packageDescription(pkg_name, lib.loc = Lib)
  )
  if (length(pkgdesc) > 1) {
    out <- paste0(" ", pkgdesc$Version)
  } else {
    out <- ""
  }
  return(out)
}
