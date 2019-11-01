# All migrant winter females from 2018 – one way migration to breeding grounds

library(moveVis)
library(move)
library(raster)
library(ggplot2)
library(magrittr)
library(readxl)

move_all_spring_df <- as.data.frame(read_xls("Data/move_all_spring.xlsx"))

move_all_spring_df <- as.data.frame(move_all_spring)

# remove the duplicate locations from the dataset, then convert to a move object
move_all_spring_df <- move_all_spring_df[!duplicated(move_all_spring_df),]
move_all <- df2move(move_all_spring_df,
                    proj = "+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
                    x = "lon", y = "lat", time = "dt", track_id = "id")

#now look at time steps
unique(timestamps(move_all))
timeLag(move_all, unit = "hours") # results in ~ 2 hour locations(with 8 hours missing overnight)


move_all_data <- align_move(move_all, res = 1, digit = 2 , unit = "hours")
#move_all_data

#frames <- frames_spatial(move_data, path_colours = c("red"),
#                        map_service = "osm", map_type = "watercolor", alpha = 0.5)

#extent of map
#ext <- extent(33.948710, 47.848907, 108.804355, 126.147724)

frames_move_all <- frames_spatial(move_all_data, path_colours = c("red", "green", "blue", "yellow", "orange", "pink", "purple"), path_legend = FALSE, path_size = 2,map_service = "mapbox", map_type = "satellite", map_token = "pk.eyJ1Ijoic3NrYWxvczQiLCJhIjoiY2pzbmZocHFkMDFndzN5cnZxNDBuejB3NCJ9.bZYBIy5C8vuMNbzpe7sVcQ") 


# edit frames
frames_move_all <- add_labels(frames_move_all, x = "Longitude", y = "Latitude") %>% 
  add_progress() %>% 
  add_scalebar(colour = "white", height = 0.015, label_margin = 1, position = "bottomright") %>% 
  add_northarrow(colour = "white", position = "bottomleft") %>% 
  add_timestamps(move_all_data, type = "label")

## Have a look at one of the frames:
#length(frames_move_all) #number of frames
frames_move_all[[500]] #look at one frame

animate_frames(frames_move_all, out_file = "~/Desktop/R_Forever/Dissertation/noha-move-hab/Output/moveVis_all_spring_migration_hours2.mov", overwrite = TRUE)
