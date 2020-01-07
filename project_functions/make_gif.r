#' make_gif.r
#' This functions depends on the magick and purrr packages.
#'
#' @param directory    directory where images to make gif are stored in order. 
#' @param output.path  where to save gif.
#' @param delay        time between gif frames in seconds.
#'
#' @return             returns a gif file.
#' @export
#'
#' @examples
make_gif <- function(directory, output.path, delay = 1){
  cat('Building gif...\n')
  #grab paths.
  gif.paths <- list.files(directory)
  gif.paths <- paste0(directory, gif.paths)
  
  #stack paths with purrr, join with magick.
  images <- map(gif.paths, image_read)
  images <- image_join(images)
  
  #Save output.
  image_write_gif(images, path = output.path, delay = delay)
  
  #report.
  msg <- paste0('gif save to ',output.path,'.\n')
  cat(msg)
  
} #end function.
