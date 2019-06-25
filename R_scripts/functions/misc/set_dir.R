.setDir = function(path, createIfNotExist = FALSE)
{
  if (!dir.exists(path))
  {
    if (createIfNotExist)
    {
      dir.create(path = path, recursive = T)
    } else
    {
      stop('path not found.')
    }
  }
  return(path)
}
