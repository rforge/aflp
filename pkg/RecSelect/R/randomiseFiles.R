#' Copy a set of files with a randomised order.
#' @param path the path of the directory with the files. Either an absolute path or relative to the working directory. If missing a dialog box will one where must choose a file within the directory.
#' @param output the path of the directory were the new files should be placed. Either an absolute path or relative to the working directory. If missing a dialog box will one where must choose a file within the directory.
#' @param extension which file extension should be used? Defaults to 'WAV'. File extensions are not case sensitive.
#' @param recursive should all subdirectories be processed as well? Defaults to TRUE.
#' @param sub.dir Process the files per subdirectory? Defaults to TRUE
#' @param min.difference the minimal difference between two timestamps in seconds.
#' @export
#' @importFrom plyr ddply
#' @importFrom tcltk tk_choose.dir
randomiseFiles <- function(path, output, recursive = TRUE, sub.dir = TRUE, extension = "WAV", min.difference = 5){
  if(missing(path)){
    path <- tk_choose.dir(getwd(), "Choose the main folder with the input files")
  }
  if(missing(output)){
    output <- tk_choose.dir(getwd(), "Choose the main folder for the output files")
  }
  files <- getTimestamp(path = path, extension = extension, recursive = recursive)
  if(sub.dir){
    files <- ddply(files, "path", function(x){
      x$order <- randomiseTimestamp(x$timestamp)
      x
    })
  } else {
    files$order <- randomiseTimestamp(files$timestamp, min.difference = min.difference)
  }
  files$newpath <- paste(
    output, 
    gsub(path, "", files$path), 
    sep = ""
  )
  np <- sample(unique(files$newpath), 1)
  junk <- sapply(unique(files$newpath), function(np){
    if(!file.exists(np)){
      dir.create(np, recursive = TRUE)
    }
  })
  n <- ceiling(log10(max(files$order)))
  fmt <- sprintf("%%s/%%0%ii_%%s", n)
  files$newfile <- sprintf(fmt, files$newpath, files$order, files$file)
  copied <- file.copy(
    from = paste(files$path, files$file, sep = "/"),
    to = files$newfile
  )
  files <- files[, c("newfile", "timestamp", "order")]
  write.csv2(
    files, 
    file = paste(output, "filelist.csv", sep = "/")
  )
  message("Finished. ", sum(copied), " of ", length(copied), " files copied to ", output)
  invisible(files)
}
