MicroBREW_count <- function(filename, save_file){

  min_trap_size = 8000
  max_trap_size = 13000
  data = read.csv(filename)
  data1 = subset(data, data$intensity_table7 < max_trap_size & data$intensity_table7 > min_trap_size)
  data1 = subset(data1, data1$intensity_table1 <= max(data1$intensity_table1))
  first_frame = subset(data1, data1$intensity_table1 == min(data1$intensity_table1)) # define the first good frame as a reference
  first_frame$trap_ID = round(first_frame$intensity_table9, -2) * 100 + round(first_frame$intensity_table8, -2)
  data1$trap_id = round(data1$intensity_table9, -2) * 100 + round(data1$intensity_table8, -2)
  table_xy = data.frame()
  for(i in unique(first_frame$trap_ID)){
    x0 = as.numeric(substr(i, (nchar(i)-3), nchar(i)))
    y0 = (i - x0)/100;
    x1 = x0 - 100;
    x2 = x0 + 100;
    y1 = y0 - 100;
    y2 = y0 + 100;
    temp_xy = subset(data1, data1$intensity_table8 > x1 & data1$intensity_table8 < x2)
    temp_xy = subset(temp_xy, temp_xy$intensity_table9 > y1 & temp_xy$intensity_table9 < y2)
    rtf = temp_xy
    rtf$trap_new_ID = paste0(i)
    table_xy = rbind(table_xy, rtf)
  }

  table_values <- data.frame()
  peak_values <- data.frame()
  keep_data <- data.frame()
  unique_ids <- unique(table_xy$trap_new_ID)
  for(id in unique_ids){
    temp <- subset(table_xy, trap_new_ID == id)
    temp$avg_center <- movavg(temp$intensity_table3, n = 5, type = "s")
    temp$norm_signal <- temp$avg_center / temp$avg_center[temp$intensity_table1 == min(temp$intensity_table1)]
    peaks_mother_cell <- findpeaks(temp$norm_signal, nups = 4, ndowns = 1, minpeakheight = 3)
    mother_cell_entry <- temp$intensity_table1[peaks_mother_cell[1,2]]
    temp1 <- subset(temp, intensity_table1 >= mother_cell_entry)
    temp1$norm_signal <- NULL

    if(length(mother_cell_entry) > 0 & nrow(temp1) > 6){
      temp1$avg_shape1 <- movavg(temp1$intensity_table4, n = 6, type = "s")
      temp1$avg_shape2 <- movavg(temp1$intensity_table5, n = 6, type = "s")
      temp1$avg_shape3 <- movavg(temp1$intensity_table6, n = 6, type = "s")

      lag_difference <- c(temp1$avg_center[1], diff(temp1$avg_center, lag = 1))
      temp1$diff <- lag_difference
      temp1$sum_value <- temp1$intensity_table3 + temp1$diff

      peaks_shape1 <- findpeaks(temp1$avg_shape1, nups = 4, ndowns = 1, minpeakdistance = 6)
      peaks_shape2 <- findpeaks(temp1$avg_shape2, nups = 4, ndowns = 1, minpeakdistance = 6)
      peaks_shape3 <- findpeaks(temp1$avg_shape3, nups = 4, ndowns = 1, minpeakdistance = 6)

      ptf <- data.frame(trap_id = id, cell_entry = mother_cell_entry,
                        shape1 = ifelse(length(peaks_shape1) == 0, 0, length(peaks_shape1)),
                        shape2 = ifelse(length(peaks_shape2) == 0, 0, length(peaks_shape2)),
                        shape3 = ifelse(length(peaks_shape3) == 0, 0, length(peaks_shape3)))

      if(length(peaks_shape1) == 0){
        peaks_shape1 = matrix(data = 0, nrow = 1, ncol = 4)
      }
      else{peaks_shape1 = peaks_shape1}

      if(length(peaks_shape2) == 0){
        peaks_shape2 = matrix(data = 0, nrow = 1, ncol = 4)
      }
      else{peaks_shape2 = peaks_shape2}
      if(length(peaks_shape3) == 0){
        peaks_shape3 = matrix(data = 0, nrow = 1, ncol = 4)
      }
      else{peaks_shape3 = peaks_shape3}

      peaks_all <- rbind(data.frame(peaks_shape1, type = "shape1"),
                         data.frame(peaks_shape2, type = "shape2"),
                         data.frame(peaks_shape3, type = "shape3"))
      peaks_all$entry_frame <- mother_cell_entry
      colnames(peaks_all) <- c("peak_height", "peak_location", "peak_begin", "peak_end", "which_shape", "mother_cell_entry")
      peaks_all$trap_id <- as.character(id)
      temp1$value_pass <- ifelse(temp1$sum_value < 0, "yes", "no")
      temp2 <- subset(temp1, intensity_table1 > (mother_cell_entry + 50) & value_pass == "yes")

      peaks_all$death_frame <- ifelse(nrow(temp2) == 0, max(temp1$intensity_table1), min(temp2$intensity_table1))
      peaks_all <- subset(peaks_all, peak_location <= death_frame)
      keep_data <- rbind(keep_data, temp1)
      table_values <- rbind(table_values, ptf)
      peak_values <- rbind(peak_values, peaks_all)
    } else {
      table_values <- rbind(table_values, data.frame(trap_id = id, cell_entry = 0, shape1 = 0, shape2 = 0, shape3 = 0))
      peak_values <- rbind(peak_values, data.frame(peak_height = 0, peak_location = 0, peak_begin = 0, peak_end = 0, which_shape = 0, mother_cell_entry = 0, trap_id = id, death_frame = 0))
      keep_data <- rbind(keep_data, data.frame(intensity_table1 = 0, intensity_table2 = 0, intensity_table3 = 0, intensity_table4 = 0, intensity_table5 = 0, intensity_table6 = 0,
                                               intensity_table7 = 0, intensity_table8 = 0, intensity_table9 = 0, trap_id = id, trap_new_ID = 0,
                                               avg_center = 0, avg_shape1 = 0, avg_shape2 = 0, avg_shape3 = 0, diff = 0, sum_value = 0, value_pass = "no"))
    }
  }

  result_values = list(table_values, peak_values, keep_data)
  outputdir = paste0(getwd(), "/RMicroBREW_results")
  if (!dir.exists(outputdir)) {dir.create(outputdir)}
  write.csv(peak_values, file = paste0(outputdir, "/", save_file, "peak_all.csv"))

  return(result_values)
}
