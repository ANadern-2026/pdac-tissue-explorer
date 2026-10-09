# Small helpers shared by the data-preparation scripts: a provenance log per script run.
start_log <- function(name, dir = file.path("data-raw", "logs")) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(dir, paste0(name, "_", format(Sys.Date()), ".txt"))
  cat(name, " started ", format(Sys.time(), tz = "UTC"), " UTC\n", file = f, sep = "")
  f
}
logc <- function(log_file, ...) { msg <- paste0(...); cat(msg, "\n"); cat(msg, "\n", file = log_file, append = TRUE) }
end_log <- function(log_file) {
  cat("finished ", format(Sys.time(), tz = "UTC"), " UTC\n\n", file = log_file, append = TRUE, sep = "")
  capture.output(sessionInfo(), file = log_file, append = TRUE)
  invisible(log_file)
}
