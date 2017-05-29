#####################################################################################
# 2017-05-28 
# Objective: Function renames Micasense RedEdge Camera files for database uploading
# New file naming reflects the original MicaSense Camera folder structure.
# Author: singhdj2@ksu.edu
#####################################################################################

### Note: This function assumes the input argument 'directoryPath' has MicaSense file structure i.e.
### Directory Path >> 000NSET Directories (Actual Flights) >> 000N Folders >> '.tif' image files
### Run this function on server after saving a backup copy of your work elsewhere (local machine/Hard Drive)

rename_RE_files <- function(dirPath){
  if (!is.null(dirPath)) setwd(dirPath)
  main_dir_list <- list.dirs(full.names = F, recursive = F)
  cat("*There are",length(main_dir_list),"flight folders in this dataset...", "\n")
  lapply(main_dir_list, FUN = function(i){ 
    cat("**Entering", i, "flight folder...","\n")
    sub_dir_list <- list.dirs(i, full.names = F, recursive = F)
    dat_list <- list.files(i, pattern = ".dat")
    if (length(dat_list) != 0) {
      lapply(dat_list, FUN= function(x){
        old_dat_path = paste0(i,"/",x)
        dat_path = paste0(dirPath,"/",i,"_",x)
        file.rename(from = old_dat_path, to = dat_path)
      })
    } 
    cat("  ***Folder", i, "has", length(sub_dir_list), "sub-folders..","\n")
    lapply(sub_dir_list, FUN = function(j){
      cat("    ****Entering sub-folder", j, "of", i, "now...", "\n")
      file_path <- paste0(i,"/",j)
      file_list <- list.files(file_path, pattern = ".tif")
      cat("      *****Sub-folder", "'", j,"'", "has", length(file_list), "image files..","\n")
      lapply(file_list, FUN = function(k){
        old_name <- paste0(file_path,"/",k)
        #new_name <- paste0(dir,"/",dir,"_",sub_dir,"_",file)
        #file.rename(from = old_name, to = new_name)
        new_name <- paste0(dirPath,"/",i,"_",j,"_",k)
        file.rename(from = old_name, to = new_name) # careful though; rename deletes the files from their original location. Alternatively use file.copy
        #cat(old_name, "\n")
        #cat(new_name, "\n")
      })
    })
    #delete the copied folders from their original location
    unlink(i, recursive = TRUE)
  })
}


## Use: run function on RedEdge flights folder!
system.time(rename_files("~/Downloads/test_rededge1/"))
