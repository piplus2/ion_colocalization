.vector <- function(..., names)
{
  v <- vector(...)
  names(v) <- names
  return(v)
}
