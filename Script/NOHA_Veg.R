## This script imports clipped location data to the Suisun 2015 Veg Map with extracted habitat types for each location, summarized using the group_by abd tally functions within dplyr##

getwd()

library(dplyr)
library(readxl)

laureen_veg_df <- as.data.frame(read.csv("~/Desktop/R_Forever/RRF/Data/Laureen_Joined_Veg_Table.csv"))
View(laureen_veg_df)

mama_veg_df <- as.data.frame(read.csv("~/Desktop/R_Forever/RRF/Data/Mama_Joined_Veg_Table.csv"))
View(mama_veg_df)

laureen_veg_df %>% group_by(CWHRType) %>% tally() %>% print(n=40)

mama_veg_df %>% group_by(CWHRType) %>% tally() %>% print(n=40)


