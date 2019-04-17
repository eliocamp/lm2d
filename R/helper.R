smart_svd <- function(X, max_eof) {
    if (check_package("irlba", NULL) &
        max_eof < 0.5*min(nrow(X), ncol(X))) {
      irlba::irlba(X, max_eof, max_eof)
    } else {
      svd(X, max_eof, max_eof)
    }
}



#' @import data.table
tidy2matrix <- function(data, formula, value.var, fill = NULL, ...) {
  row.vars <- all.vars(formula[[2]])
  col.vars <- all.vars(formula[[3]])
  data <- data.table::setDT(data)
  data[, row__ := .GRP, by = c(row.vars)]
  data[, col__ := .GRP, by = c(col.vars)]
  if (is.null(fill)){
    fill <- 0
    rowdims <- data[col__ == 1, (row.vars), with = FALSE]
    coldims <- data[row__ == 1, (col.vars), with = FALSE]
  } else {
    rowdims <- unique(data[, (row.vars), with = FALSE])
    coldims <- unique(data[, (col.vars), with = FALSE])
  }

  data.m <- matrix(fill[1], nrow = max(data[["row__"]]),
                   ncol = max(data[["col__"]]))
  data.m[cbind(data[["row__"]], data[["col__"]])] <- data[[value.var]]

  return(list(matrix = data.m,
              coldims = coldims,
              rowdims = rowdims))
}


as.character.formula <- function(x, ...) {
  form <- paste(deparse(x), collapse = " ")
  form <- gsub("\\s+", " ", form, perl = FALSE)
  return(form)
}


enrich_formula <- function(formula) {
  f <- as.character(formula)
  f <- stringr::str_split(f,"~", n = 2)[[1]]

  value.var <- stringr::str_squish(f[!stringr::str_detect(f, "\\|")])

  matrix.vars <- f[stringr::str_detect(f, "\\|")]
  matrix.vars <- stringr::str_split(matrix.vars,"\\|", n = 2)[[1]]

  row.vars <- stringr::str_squish(stringr::str_split(matrix.vars[1], "\\+")[[1]])
  col.vars <- stringr::str_squish(stringr::str_split(matrix.vars[2], "\\+")[[1]])

  dcast.formula <- stringr::str_squish(f[stringr::str_detect(f, "\\|")])
  dcast.formula <- as.formula(stringr::str_replace(dcast.formula, "\\|", "~"))
  value.var <- stringr::str_squish(f[!stringr::str_detect(f, "\\|")])

  return(list(formula = formula,
              dims = dcast.formula,
              value.var = value.var,
              row.vars = row.vars,
              col.vars = col.vars))
}

data_from_formula <- function(formula, data, extra.vars = NULL) {
  env <- attr(formula$formula, ".Environment")
  formula$formula <- Formula::as.Formula(paste(c(as.character(formula$formula), extra.vars), collapse = " + "))
  attr(formula$formula, ".Environment") <- env

  if (is.null(data)) {
    data <- eval(quote(model.frame(Formula::as.Formula(formula$formula),
                                                 data  = data)))
  } else {
    # Check if columns are indata
    all.cols <- c(formula$value.var, formula$row.vars, formula$col.vars, extra.vars)
    missing.cols <- all.cols[!(all.cols %in% colnames(data))]
    if (length(missing.cols) != 0) {
      stop(paste0("Columns not found in data: ", paste0(missing.cols, collapse = ", ")))
    }
    data <- as.data.frame(data)[all.cols]
  }
  return(data)
}



check_package <- function(package,
                          message = paste0("Package ", package, " not installed.")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    if (is.null(message)) {
      return(FALSE)
    } else {
      stop(message, call. = FALSE)
    }
  }
  invisible(TRUE)
}


if(getRversion() >= "2.15.1")  utils::globalVariables(c("col__", "row__"))
