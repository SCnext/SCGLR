# weighted scale functions
# @param x vector to scale
# @param w weight
# @return scaled vector
wtScale <-function(x,w) {
  xc=x-sum(w*x)
  v=sum(xc^2*w)
  xcr=xc/sqrt(v)
  return(xcr)
}

# utils functions
# @param x vector to center
# @param w weight
# @return centered vector
wtCenter=function(x,w) {
  xc=x-sum(w*x)
  return(xc)
}

is.categorical <- function(x) return(class(x) %in% c("character","factor","logical"))
which.categorical <- function(x) return(which(unlist(lapply(x,is.categorical))))
has.categorical <- function(x) return(length(which.categorical(x))!=0)

checkLossFunction <- function(type) {
  if(!type %in% c("auc","likelihood","aic","aicc","bic","mspe"))
    stop("Unknown loss function!")
}

# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

# # returns a call in which all of the arguments which were supplied or have presets are specified by their full names and supplied or default values.
# # @param definition a function. See \code{\link[base]{match.call}}.
# # @param call an unevaluated call to the function specified by definition. See \code{\link[base]{match.call}}.
# # @param expand.dots logical. Should arguments matching ... in the call be included or left as a ... argument? See \code{\link[base]{match.call}}.
# # @return An object of class call.
# # @author Fabian Scheipl
# # @export
# # @seealso \code{\link[base]{match.call}}
# expand.call <-
# function(definition=sys.function(sys.parent()), call=sys.call(sys.parent(1)), expand.dots = TRUE,envir=parent.frame(2L))
# {
# 	call <- .Internal(match.call(definition, call, expand.dots,envir))
# 	#given args:
# 	ans <- as.list(call)
# 	# ans1 <- ans[[1]]
# 	# ans <- lapply(ans[-1], eval, envir = sys.frame(sys.parent(2)))
# 	# ans <- c(ans1, ans)
#
# 	#possible args:
# 	frmls <- formals(deparse(ans[[1]]))
# 	#remove formal args with no presets:
# 	frmls <- frmls[!sapply(frmls, is.symbol)]
#
# 	add <- which(!(names(frmls) %in% names(ans)))
# 	return(as.call(c(ans, frmls[add])))
# }

condition <- function(subclass, message, call=sys.call(-1), ...) {
  structure(class = c(subclass, "condition"), list(message=message, call=call), ...)
}

custom_stop <- function(subclass, ..., call=sys.call(-1)) {
  message = .makeMessage(...)
  if(message=="")
    message=subclass
  c <- condition(c(subclass, "error"), message=message, call=call)
  stop(c)
}
