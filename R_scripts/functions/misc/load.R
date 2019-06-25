# Custom load function

.load <- function(x)
{
  e <- new.env()
  load(x, envir = e)
  return(e)
}
