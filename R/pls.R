pls <- function(..., method=c("SIMPLS", "kernelPLS"))
{
  method <- match.arg(method)
  mvr(..., method=method)
}
