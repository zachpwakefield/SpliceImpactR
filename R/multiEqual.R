':=' <- function(lhs, rhs) {
  frame <- parent.frame() # Capture the environment of the parent frame to perform assignments there
  lhs <- as.list(substitute(lhs)) # Convert the left-hand side argument to a list of symbols

  # If lhs contains more than one symbol, remove the first element (the function name ':=')
  if (length(lhs) > 1) {
    lhs <- lhs[-1]
  }

  # Single assignment case
  if (length(lhs) == 1) {
    # Perform the assignment in the captured environment
    do.call(`=`, list(lhs[[1]], rhs), envir=frame)
    return(invisible(NULL)) # Return invisibly to avoid console output
  }

  # Handle cases where rhs is a function or formula by wrapping it in a list
  if (is.function(rhs) || is(rhs, 'formula')){
    rhs <- list(rhs)
  }
  # If there are fewer rhs values than lhs symbols, extend rhs with NULLs to match the length
  if (length(lhs) > length(rhs)) {
    rhs <- c(rhs, rep(list(NULL), length(lhs) - length(rhs)))
  }
  # Perform assignments for each lhs-rhs pair in the parent frame environment
  for (i in 1:length(lhs)){
    do.call(`=`, list(lhs[[i]], rhs[[i]]), envir=frame)
    }
  return(invisible(NULL)) # Return invisibly to avoid console output
}
