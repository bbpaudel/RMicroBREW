modeRLS <- function(data, mode_Frames, plot_it, save_file){
  data = data
  data = data[data$peak_height != 0, ]
  mode_peaks = Mclust(data$peak_height, mode_Frames)
  df = data.frame(mean_height = mode_peaks$parameters$mean, proportions = mode_peaks$parameters$pro)
  data$classification = factor(mode_peaks$classification)

  outputdir = paste0(getwd(), "/RMicroBREW_results")
  if (!dir.exists(outputdir)) {dir.create(outputdir)}
  write.csv(df, file = paste0(outputdir, "/", save_file,"proportions.csv"))
  write.csv(data, file = paste0(outputdir, "/", save_file,"classification_data.csv"))

  results = list(df, data)
  if(plot_it == "yes"){
    print(ggplot(data, aes(x = classification, y = peak_height)) + geom_boxplot(notch = T) + theme_bw() +
            labs(x = "mode_class", y = "budding_size_pixels"))
  } else{""}
  return(results)
}

