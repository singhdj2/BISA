## Objective: create an sql DB query from 2016-17 GID list
# 2017-05-08
# singhdj2

# load and cleanup the GID data: remove duplicates etc..
# load
library(readr)
gid17 <- read_csv("~/Dropbox/1.Research_and_Education/Kansas_State_University_2014-/Doctoral_Research/Projects_data/UAV_pipeline_work/Data_UAV/image_data/India_2017/raw_data/field_layouts/17SABWGSYT_GID_list.csv")

length(unique(gid17$GID))
gt <- unique(gid17$GID)

## prefix GID for 'sample_name' column
# use sample_name instead of tissue_id as this is more consistent!
query_sample <- paste("SELECT * FROM dna WHERE sample_name IN ('",
               paste("GID",gt,collapse = "','", sep = ""),"')",sep = "")
write.table(query_sample, file= "Downloads/gid17_sample_query.txt")

####### Check GID key file has all GIDs ######
#import
gbs_key <- read_csv("~/Dropbox/1.Research_and_Education/Kansas_State_University_2014-/Doctoral_Research/Projects_data/UAV_pipeline_work/Data_UAV/image_data/India_2017/analysis/sequence_data/gbs_gid1617_key_20170508.txt")

q.gid <- paste("GID",gt, sep = "")
length(unique(gbs_key$sample_name)) == length(unique(q.gid)) #compare lengths
length(unique(gbs_key$sample_name)) == length(unique(gid17$GID)) #compare lengths


################################################
##### Repeat the same for 2015-16 GID list #####
################################################
gid1516 <- read_csv("~/Dropbox/Kevin_DJ_orthos_2016/final_data_all/gid_entry_info_GS_India_2016.csv")
View(gid1516)
length(unique(gid1516$GID))
gt16 <- unique(gid1516$GID)
gt16 <- gt16[!is.na(gt16)]

query_sample <- paste("WHERE sample_name IN ('",
                      paste("GID",gt16,collapse = "','", sep = ""),"')",sep = "")
write.table(query_sample, file= "Downloads/gid16_sample_query.txt")

## run sql query on navicat
#import gbs key file
gbs_key <- read_csv("~/Dropbox/1.Research_and_Education/Kansas_State_University_2014-/Doctoral_Research/Projects_data/UAV_pipeline_work/Data_UAV/image_data/India_2016/sequence_data/gbs_key_1516_20170508.txt")

q.gid <- paste("GID",gt16, sep = "")
length(unique(gbs_key$sample_name)) == length(unique(q.gid)) #compare lengths
length(unique(gbs_key$sample_name)) == length(unique(gt16)) #compare lengths

# why not just combine 2014 to 2017 data together!

#### first check 2014 as well

gid1415 <- read_csv("~/Dropbox/1.Research_and_Education/Kansas_State_University_2014-/Doctoral_Research/Projects_data/UAV_pipeline_work/Data_UAV/image_data/India_2015/Ladhowal/field_info_2015/SABWGPYT01-To-10.csv")
View(gid1415)
length(unique(gid1415$GID))
gt15 <- unique(gid1415$GID)
gt15 <- gt15[!is.na(gt15)]

query_sample <- paste("WHERE sample_name IN ('",
                      paste("GID",gt15,collapse = "','", sep = ""),"')",sep = "")
write.table(query_sample, file= "Downloads/gid15_sample_query.txt")

## run sql query on navicat
#import gbs key file
gbs_key <- read_csv("~/Dropbox/1.Research_and_Education/Kansas_State_University_2014-/Doctoral_Research/Projects_data/UAV_pipeline_work/Data_UAV/image_data/India_2015/sequence_data/gbs_key_1415_20170508.txt")

q.gid <- paste("GID",gt15, sep = "")
length(unique(gbs_key$sample_name)) == length(unique(q.gid)) #compare lengths
length(unique(gbs_key$sample_name)) == length(unique(gt15)) #compare lengths
# Almost same length so combine all years together!


###########################
##### combine all years..##
###########################
gid1417 <- read_csv("~/Dropbox/1.Research_and_Education/Kansas_State_University_2014-/Doctoral_Research/Projects_data/UAV_pipeline_work/Data_UAV/image_data/India_2017/raw_data/field_layouts/gid_list_three_years_2015-17.csv")
length(unique(gid1417$GID))
gt <- unique(gid1417$GID)
gt <- gt[!is.na(gt)]

query_sample <- paste("WHERE sample_name IN ('",
                      paste("GID",gt,collapse = "','", sep = ""),"')",sep = "")
write.table(query_sample, file= "Downloads/gid1417_sample_query.txt")

## run sql query on navicat
#import gbs key file
gbs_key <- read_csv("~/Dropbox/1.Research_and_Education/Kansas_State_University_2014-/Doctoral_Research/Projects_data/UAV_pipeline_work/Data_UAV/image_data/India_2016/sequence_data/gbs_key_1516_20170508.txt")

q.gid <- paste("GID",gt16, sep = "")
length(unique(gbs_key$sample_name)) == length(unique(q.gid)) #compare lengths
length(unique(gbs_key$sample_name)) == length(unique(gt16)) #compare lengths


